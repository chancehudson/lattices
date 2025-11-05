use std::collections::BTreeMap;
use std::collections::HashMap;
use std::sync::Arc;
use std::sync::LazyLock;
use std::sync::RwLock;

use crate::*;

/// Store a cache of (E::CARDINALITY, sigma) keyed to a CDT
///
/// sigma will be stored as sigma * 10^5 (up to 5 decimals precision for std dev)
/// this is independent of the floating point accuracy inside the CDT
static CDT_CACHE: LazyLock<RwLock<HashMap<(u128, u32), Arc<GaussianCDT>>>> =
    LazyLock::new(|| RwLock::new(HashMap::default()));

/// Refuse to do math with values smaller than this.
static MIN_PRECISION: LazyLock<f64> = LazyLock::new(|| 10f64.powi(5) * f64::EPSILON);

/// An instance of a cumulative distribution table for a finite field, with a specific sigma
/// and tail bounds limited by precision. Panics if the distribution cannot be precisely computed.
///
/// Minimum tail = 7 * sigma
/// Maximum tail = 15 * sigma
///
/// Entries in the finite field are referred to by "displacement". Distance from the 0 element,
/// signed to indicate forward or reverse in the field.
/// forward: `x + E::one()`
/// backward: `x - E::one()`
///
/// In the finite field with 101 elements, the element 0 is at displacement 0. Element 100 at
/// displacement -1, element 1 at displacement 1. Displacement is a measure of distance and
/// direction, and therefore does not exist in the field because fields are not partially ordered.
pub struct GaussianCDT {
    /// Cardinality of the field this CDT operates over
    pub cardinality: u128,
    /// standard deviation of the distribution
    pub sigma: f64,
    /// normalized probability paired with displacement
    pub displacements: HashMap<i32, f64>,
    /// sum of gaussian pdf evaluated over all possible output values
    /// evaluated with min 1e-5 precision
    pub normalized_sum: f64,
    /// Computed bounds based on sigma and decimal precision limits
    pub tail_bounds: (i32, i32),
}

impl GaussianCDT {
    /// generate a distinct identifier for the table
    /// `(E::CARDINALITY, sigma_key)`
    pub fn identifier<E: FieldScalar>(sigma: f64) -> (u128, u32) {
        let scale = 10f64.powi(5);

        let sigma_key = sigma * scale;
        assert!(
            sigma_key - sigma_key.floor() < 1.0,
            "CDT: sigma is too precise"
        );
        assert!(sigma_key <= u32::MAX as f64, "CDT: sigma is too large");
        let sigma_key = sigma_key as u32;

        (E::CARDINALITY, sigma_key)
    }

    /// Given a finite field and a standard deviation precompute the probabilities of selecting
    /// elements in a range up to 15 * sigma. This includes ~ 1 - 2^-167 of the probability mass.
    ///
    /// Compute minimum of 7 * sigma based on precision limitations. Panic if cannot precisely
    /// compute.
    fn new<E: FieldScalar>(sigma: f64) -> Self {
        // 15*sigma = ~2^-167 odds of sampling within the range (according to Claude)
        let std_devs = 15f64;
        let tail = std_devs * sigma;
        assert!(
            tail < 2f64.powi(20),
            "CDT refusing to build table larger than 2^20 tail (i32 limits)"
        );
        // the maximum tail we want given arbitrary precision
        let tail = tail as i32;

        // the tail we can safely sample given precision limits
        let mut actual_tail = tail;
        // compute a table of displacements. This table will extend to or beyond the resulting
        // tail_bounds
        let mut displacements = HashMap::<i32, f64>::default();
        log::debug!("CDT MIN_PRECISION: {}", *MIN_PRECISION);
        for disp in -tail..=tail {
            let prob_exp = (disp as f64).powi(2) / (2.0 * sigma * sigma);
            let is_exp_precise =
                disp == 0 || prob_exp.fract() == 0. || prob_exp.fract() > *MIN_PRECISION;
            // value of the distribution function at this point
            let prob = f64::exp(-prob_exp);
            if prob < *MIN_PRECISION || !is_exp_precise {
                assert_ne!(disp, 0);
                if disp.is_negative() {
                    actual_tail = disp.abs() - 1;
                    continue;
                }
                if disp.is_positive() {
                    assert!(actual_tail < disp);
                    break;
                }
                unreachable!();
            }
            log::debug!("CDT displacement: {} pdf: {}", disp, prob);
            displacements.insert(disp, prob);
        }
        if actual_tail < (7.0 * sigma) as i32 {
            panic!(
                "CDT refusing to build unlikely table with sigma: {} tail: {} stddevs: {}",
                sigma,
                actual_tail,
                ((actual_tail as f64 / sigma) * 10.0).floor() / 10.0
            );
        }
        // both sides inclusive
        let tail_bounds = (-actual_tail, actual_tail);
        let mut raw_sum = 0f64;
        for disp in tail_bounds.0..=tail_bounds.1 {
            let prob = displacements
                .get(&disp)
                .expect("CDT displacement {disp} does not exist in table");
            raw_sum += prob;
            log::debug!("CDT sigma {}, disp: {} prob: {}", sigma, disp, prob);
        }
        // we're guaranteed an additional 5 decimals of precision
        // from MIN_PRECISION, we'll give 5 to the raw sum division operation
        assert!(raw_sum < 10f64.powi(5));
        let mut normalized_sum = 0f64;
        for disp in tail_bounds.0..=tail_bounds.1 {
            let prob = displacements
                .get_mut(&disp)
                .expect("CDT displacement {disp} does not exist in table");
            *prob /= raw_sum;
            normalized_sum += *prob;
            *prob = normalized_sum;
        }
        Self {
            cardinality: E::CARDINALITY,
            sigma,
            tail_bounds,
            displacements,
            normalized_sum,
        }
    }

    /// Retrieve or compute a cumulative distribution table.
    pub fn cache_or_init<E: FieldScalar>(sigma: f64) -> Arc<Self> {
        let identifier = Self::identifier::<E>(sigma);
        if let Some(cdt) = CDT_CACHE.read().unwrap().get(&identifier) {
            return cdt.clone();
        }
        let out = Arc::new(Self::new::<E>(sigma));
        CDT_CACHE.write().unwrap().insert(identifier, out.clone());
        out
    }

    /// Iterate over the displacements in the tail bounds of the probability table.
    fn displacements_iter(&self) -> impl Iterator<Item = (i32, f64)> {
        (self.tail_bounds.0..=self.tail_bounds.1).map(|i| {
            self.displacements
                .get(&i)
                .map(|prob| (i, *prob))
                .expect("CDT did not find entry for displacement {i}")
        })
    }

    /// Sample an element from the distribution.
    pub fn sample<E: FieldScalar, R: Rng>(&self, rng: &mut R) -> E {
        let r: f64 = rng.random_range(0.0..1.0);
        let mut out = None;
        // always iterate over the whole space for timing smoothness
        for (disp, prob) in self.displacements_iter() {
            out = out.or(if r < prob {
                Some(E::at_displacement(disp))
            } else {
                None
            });
        }
        out.unwrap_or(E::at_displacement(self.tail_bounds.1))
    }

    /// Sample a constant size array of elements.
    pub fn sample_arr<const N: usize, E: FieldScalar, R: Rng>(&self, rng: &mut R) -> [E; N] {
        let samples: [f64; N] = std::array::from_fn(|_| rng.random_range(0.0..1.0));

        let mut out = [E::zero(); N];
        let mut matched_out = [false; N];
        let mut matched_sample_count = 0;
        for (disp, prob) in self.displacements_iter() {
            for (i, sample) in samples.iter().enumerate() {
                if *sample < prob && !matched_out[i] {
                    matched_out[i] = true;
                    matched_sample_count += 1;
                    out[i] = E::at_displacement(disp);
                }
            }
        }
        assert_eq!(matched_sample_count, N, "CDT not all samples were matched");
        out
    }

    /// Sample a vector of elements of length `len` from the distribution.
    pub fn sample_vec<E: FieldScalar, R: Rng>(&self, len: usize, rng: &mut R) -> Vector<E> {
        let mut samples = BTreeMap::<usize, f64>::default();
        for i in 0..len {
            samples.insert(i, rng.random_range(0.0..1.0));
        }
        let mut out = BTreeMap::<usize, E>::default();
        for (disp, prob) in self.displacements_iter() {
            samples.retain(|i, sample| {
                if *sample < prob {
                    out.insert(*i, E::at_displacement(disp));
                    return false;
                }
                true
            });
            if samples.is_empty() {
                break;
            }
        }
        assert!(samples.is_empty(), "CDT not all samples were matched");
        assert_eq!(out.len(), len, "CDT outputting invalid sample len");
        out.into_values().collect::<Vector<_>>()
    }

    /// Probability of selecting a certain displacement in this CDT.
    pub(crate) fn prob(&self, disp: &i32) -> f64 {
        let prev_prob = self
            .displacements
            .get(&(*disp - 1))
            .copied()
            .unwrap_or_default();
        *self
            .displacements
            .get(&disp)
            .expect("CDT requested probability of element not in table")
            - prev_prob
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// Repeatedly invoke a test function and provide a CDT + 100,000 samples.
    const SAMPLES_PER_CDT: usize = 100_000;
    fn get_cdt_sample_pairs<E: FieldScalar, R: Rng>(
        test_logic: fn(Arc<GaussianCDT>, HashMap<i32, usize>, rng: &mut R),
        rng: &mut R,
    ) {
        // we'll test std deviations 2.0 to 10.0
        let range = 20..=100;
        let sigma_iter = range.map(|i| (i as f64) / 10.0);

        // individual sampling
        for sigma in sigma_iter.clone() {
            // individual retrieval
            let cdt = GaussianCDT::cache_or_init::<E>(sigma);
            let mut samples = HashMap::<i32, usize>::default();

            for _ in 0..SAMPLES_PER_CDT {
                let disp = cdt.sample::<E, _>(rng).displacement();
                *samples.entry(disp as i32).or_default() += 1;
            }
            test_logic(cdt, samples, rng);
        }

        // vectorized retrieval
        for sigma in sigma_iter.clone() {
            let cdt = GaussianCDT::cache_or_init::<E>(sigma);
            let mut samples = HashMap::<i32, usize>::default();

            let batch_size: usize = rand::random_range(1..30);
            let mut sample_count = 0usize;

            loop {
                if sample_count >= SAMPLES_PER_CDT {
                    break;
                }
                let batch_size = batch_size.min(sample_count.abs_diff(SAMPLES_PER_CDT));
                let disps: Vector<E> = cdt.sample_vec(batch_size, rng);
                for disp in disps {
                    *samples.entry(disp.displacement() as i32).or_default() += 1;
                }
                sample_count += batch_size;
            }
            test_logic(cdt, samples, rng);
        }

        // array retrieval
        for sigma in sigma_iter.clone() {
            let cdt = GaussianCDT::cache_or_init::<E>(sigma);
            let mut samples = HashMap::<i32, usize>::default();

            let mut sample_count = 0usize;
            loop {
                if sample_count >= SAMPLES_PER_CDT {
                    break;
                }
                for disp in cdt.sample_arr::<10, E, _>(rng) {
                    sample_count += 1;
                    *samples.entry(disp.displacement() as i32).or_default() += 1;
                    if sample_count >= SAMPLES_PER_CDT {
                        break;
                    }
                }
            }
            test_logic(cdt, samples, rng);
        }
    }

    // each variable gets locally shadowed as f64 casted type
    macro_rules! as_f64 {
        ($($name: ident),*) => {
            $(
                let $name = ($name).clone() as f64;
            )*
        };
    }

    #[test]
    fn cdt_mean() {
        type Field = OxfoiScalar;
        let rng = &mut rand::rng();

        get_cdt_sample_pairs::<Field, _>(
            |cdt, samples, _rng| {
                let mean = samples
                    .iter()
                    .map(|(disp, count)| {
                        as_f64!(disp, count);
                        disp * count
                    })
                    .sum::<f64>()
                    / SAMPLES_PER_CDT as f64;

                let std_err = cdt.sigma / (SAMPLES_PER_CDT as f64).sqrt();
                let tolerance = 3.5 * std_err;

                // println!("sigma: {}, mean: {mean}, tolerance: {tolerance}", cdt.sigma);
                assert!(mean.abs() < tolerance);
            },
            rng,
        );
    }

    #[test]
    fn cdt_std_dev() {
        type Field = OxfoiScalar;
        let rng = &mut rand::rng();

        get_cdt_sample_pairs::<Field, _>(
            |cdt, samples, _rng| {
                let mut sum = 0f64;
                for (disp, count) in &samples {
                    as_f64!(disp, count);
                    sum += disp * count;
                }
                let mean = sum / SAMPLES_PER_CDT as f64;
                let variance = samples
                    .iter()
                    .map(|(disp, count)| {
                        as_f64!(count, disp, mean);
                        count * (disp - mean).powi(2)
                    })
                    .sum::<f64>()
                    / SAMPLES_PER_CDT as f64;
                let std_dev = variance.sqrt();
                let percent_diff = ((std_dev - cdt.sigma) / cdt.sigma).abs();
                // measured std_dev within 1% of sigma
                assert!(percent_diff < 0.01);
            },
            rng,
        );
    }

    #[test]
    fn cdt_symmetry() {
        type Field = OxfoiScalar;
        let rng = &mut rand::rng();
        get_cdt_sample_pairs::<Field, _>(
            |_cdt, samples, _rng| {
                let mut total_neg = 0f64;
                let mut total_pos = 0f64;
                for (disp, count) in samples {
                    if disp < 0 {
                        total_neg += count as f64;
                    } else if disp > 0 {
                        total_pos += count as f64;
                    }
                }
                // negative and positive counts within 3% diff
                assert!((1.0 - total_neg / total_pos).abs() < 0.03);
            },
            rng,
        );
    }

    #[test]
    fn cdt_chi_squared_fit() {
        type Field = OxfoiScalar;
        let rng = &mut rand::rng();

        get_cdt_sample_pairs::<Field, _>(
            |cdt, mut samples, _rng| {
                let mut chi_sq = 0f64;
                for disp in cdt.tail_bounds.0..cdt.tail_bounds.1 {
                    let count = samples.entry(disp).or_default();
                    let expected = cdt.prob(&disp) * SAMPLES_PER_CDT as f64;
                    if expected < 1.0 {
                        continue;
                    }
                    chi_sq += (*count as f64 - expected).powi(2) / expected;
                }
                let df = samples.len() - 1;
                let expected = chi_sq_95(df);
                assert!(
                    chi_sq < expected,
                    "{chi_sq} outside of bound 95% {expected}"
                );
            },
            rng,
        );
    }
}
