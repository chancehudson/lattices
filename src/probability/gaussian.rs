use std::collections::HashMap;
use std::sync::Arc;
use std::sync::RwLock;

use crate::*;

/// Store a cache of (Element::Cardinality, sigma) keyed to a displacement
/// sigma will be stored as sigma * 10^5 (up to 5 decimals precision for sigma keys)
/// this is independent of the floating point accuracy inside the CDT
static CDT_CACHE: LazyLock<RwLock<HashMap<(u128, u32), Arc<GaussianCDT>>>> =
    LazyLock::new(|| RwLock::new(HashMap::default()));
/// How far from the standard deviation we should precompute.
const TAIL_BOUND_MULTIPLIER: f64 = 8.0;

/// An instance of a cumulative distribution table for a finite field, with a specific sigma
/// and constant tail bounds of TAIL_BOUND_MULTIPLIER * σ.
///
/// Entries in the finite field are referred to by "displacement". Distance from the 0 element,
/// signed to indicate forward or reverse in the field.
///
/// In the finite field with 101 elements, the element 0 is at displacement 0. Element 100 at
/// displacement -1, element 1 at displacement 1. Displacement is a measure of distance and
/// direction, and therefore does not exist in the field because fields are not partially ordered.
pub struct GaussianCDT {
    /// Cardinality of the field this CDT operates over
    pub cardinality: u128,
    /// standard deviation of the distribution
    pub sigma: f64,
    /// normalized probability sum paired with displacement
    pub displacements: Vec<(f64, i32)>,
    // sum of the PDF evaluated over all possible output values
    pub normalized_sum: f64,
}

impl GaussianCDT {
    /// sample an element from the distribution
    pub fn sample<F: Element, R: Rng>(&self, rng: &mut R) -> F {
        let r: f64 = rng.random_range(0.0..1.0);
        for i in 0..self.displacements.len() - 1 {
            let (last_prob, disp) = self.displacements[i];
            let (next_prob, _) = self.displacements[i + 1];
            if r >= last_prob && r < next_prob {
                return F::at_displacement(disp);
            }
        }
        panic!("sampled probability is outside CDT");
    }

    /// Sample a vector of elements of length `len`
    pub fn sample_vec<F: Element, R: Rng>(&self, len: usize, rng: &mut R) -> Vector<F> {
        let mut probs = Vec::with_capacity(len);
        for _ in 0..len {
            probs.push(rng.random_range(0.0..1.0));
        }
        let mut out = Vec::with_capacity(len);
        for i in 0..self.displacements.len() - 1 {
            let (last_prob, disp) = self.displacements[i];
            let (next_prob, _) = self.displacements[i + 1];
            for prob in &probs {
                if *prob >= last_prob && *prob < next_prob {
                    out.push(F::at_displacement(disp));
                }
            }
        }
        assert_eq!(
            out.len(),
            len,
            "CDT vector len mismatch, sampled outside the CDT"
        );
        out.into()
    }

    /// Probability of selecting a certain displacement in this CDT.
    pub(crate) fn prob(&self, disp: i32) -> f64 {
        for i in 1..self.displacements.len() {
            let (last_prob, last_disp) = self.displacements[i - 1];
            let (next_prob, _) = self.displacements[i];
            if last_disp == disp {
                return next_prob - last_prob;
            }
        }
        0.0
    }

    /// Compute or retrieve a cumulative distribution table.
    pub fn new<F: Element>(sigma: f64) -> Arc<Self> {
        let theta_key = sigma * 10f64.powi(5);
        assert!(
            theta_key - theta_key.floor() < 1.0,
            "CDT: sigma is too precise"
        );
        assert!(theta_key <= u32::MAX as f64, "CDT: sigma is too large");
        let theta_key = theta_key as u32;
        if let Some(cdt) = CDT_CACHE.read().unwrap().get(&(F::CARDINALITY, theta_key)) {
            return cdt.clone();
        }
        let dist = (TAIL_BOUND_MULTIPLIER * sigma).ceil() as i32;
        assert!(dist >= 1, "sigma is too small");
        log::info!("Building CDT with max {} elements", dist * 2 + 1);
        if dist > 50 {
            log::warn!("Building CDT with more than 100 elements. Consider adjusting tail bounds.");
        }
        let mut displacements = Vec::default();
        let mut total_prob = 0f64;
        for disp in -dist..=dist {
            let prob_exp = (disp as f64).powi(2) / (2.0 * sigma * sigma);
            // value of the distribution function at this point
            let prob = f64::exp(-prob_exp);
            displacements.push((prob, disp));
            total_prob += prob;
            log::debug!("CDT sigma {}, disp: {} prob: {}", sigma, disp, prob);
        }
        log::debug!("CDT actual size: {}", displacements.len());
        let mut normalized_sum = 0f64;
        for (prob, _disp) in displacements.iter_mut() {
            *prob /= total_prob;
            let prob_floor = normalized_sum;
            normalized_sum += *prob;
            *prob = prob_floor;
        }
        let out = Arc::new(Self {
            cardinality: F::CARDINALITY,
            sigma,
            displacements,
            normalized_sum: total_prob,
        });
        CDT_CACHE
            .write()
            .unwrap()
            .insert((F::CARDINALITY, theta_key), out.clone());
        out
    }
}

#[cfg(test)]
mod test {

    use crate::probability::chi_sq::chi_sq_95;

    use super::*;

    /// Repeatedly invoke a test function and provide a CDT + 100,000 samples.
    const SAMPLES_PER_CDT: usize = 100_000;
    fn get_cdt_sample_pairs<E: Element, R: Rng>(
        test_logic: fn(Arc<GaussianCDT>, HashMap<i32, usize>, rng: &mut R),
        rng: &mut R,
    ) {
        // individual retrieval
        for i in 10..100 {
            let sigma = (i as f64) / 10.;
            let cdt = GaussianCDT::new::<E>(sigma);
            let mut samples = HashMap::<i32, usize>::default();

            for _ in 0..SAMPLES_PER_CDT {
                let disp = cdt.sample::<E, _>(rng).displacement();
                *samples.entry(disp as i32).or_default() += 1;
            }
            test_logic(cdt, samples, rng);
        }
        // vectorized retrieval
        for i in 10..100 {
            let sigma = (i as f64) / 10.;
            let cdt = GaussianCDT::new::<E>(sigma);
            let mut samples = HashMap::<i32, usize>::default();

            let batch_size: usize = rand::random_range(1..10);
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
    }

    #[test]
    fn cdt_mean() {
        type Field = OxfoiScalar;
        let rng = &mut rand::rng();

        get_cdt_sample_pairs::<Field, _>(
            |cdt, samples, _rng| {
                let mut sum = 0f64;
                for (disp, count) in samples.iter() {
                    sum += *disp as f64 * *count as f64;
                }
                // check that mean < 3*sigma/sqrt(N)
                assert!(
                    (sum / SAMPLES_PER_CDT as f64).abs()
                        < (3. * cdt.sigma) / (SAMPLES_PER_CDT as f64).sqrt()
                );
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
                for (disp, count) in samples.iter() {
                    sum += *disp as f64 * *count as f64;
                }
                let mean = sum / SAMPLES_PER_CDT as f64;
                let mut variance = 0f64;
                for (disp, count) in samples {
                    variance += count as f64 * (disp as f64 - mean).powi(2);
                }
                let variance = variance / SAMPLES_PER_CDT as f64;
                let std_dev = variance.sqrt();
                let percent_diff = ((std_dev - cdt.sigma) / cdt.sigma).abs();
                // measured std_dev within 1% ofsigma
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
                for disp in ((-cdt.sigma * 10.) as i32)..((cdt.sigma * 10.) as i32) {
                    let count = samples.entry(disp).or_default();
                    let expected = cdt.prob(disp) * SAMPLES_PER_CDT as f64;
                    if expected < 1.0 {
                        continue;
                    }
                    chi_sq += (*count as f64 - expected).powi(2) / expected;
                }
                let df = samples.len() - 1;
                let expected = chi_sq_95(df);
                println!("{} {}", expected, chi_sq);
                assert!(
                    chi_sq < expected,
                    "{chi_sq} outside of bound 95% {expected}"
                );
            },
            rng,
        );
    }
}
