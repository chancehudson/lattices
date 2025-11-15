use std::collections::HashMap;
use std::sync::Arc;
use std::sync::LazyLock;
use std::sync::RwLock;

use crate::*;

const CARDINALITY: u32 = 455 * 2u32.pow(20) * 9 + 1;
const MONT: Montgomery32 = Montgomery32::new(CARDINALITY);

/// Finite field of cardinality 4293918721. This is ~2^32 - 1mil (close to 32 bit boundary)
/// Q is factored by 2^20 thus supports cyclotomic rings up to degree ~1mil.
///
/// This prime is referred to as "milli"
#[derive(Debug, Copy, Clone, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MilliScalarMont(u32);

impl RingElement for MilliScalarMont {
    const CARDINALITY: u128 = CARDINALITY as u128;

    fn sample_uniform<R: Rng>(rng: &mut R) -> Self {
        Self::from(rng.random::<u32>())
    }
}

impl FieldScalar for MilliScalarMont {
    /// Return the finite field element at a certain displacement.
    fn at_displacement(disp: i32) -> Self {
        if disp.unsigned_abs() > CARDINALITY / 2 {
            log::info!(
                "Attempting to initialize a displacement outside the field: {} {}",
                disp,
                Self::CARDINALITY
            );
            #[cfg(not(debug_assertions))]
            panic!("refusing to use displacement outside of field in production");
        }
        if disp >= 0 {
            Self::from(disp as u32)
        } else {
            Self::from(CARDINALITY - disp.unsigned_abs() as u32)
        }
    }

    fn prime_factorization() -> impl Iterator<Item = (Self, usize)> {
        static PRIME_FACTORIZATION: LazyLock<HashMap<u128, usize>> = LazyLock::new(|| {
            log::info!("Building prime factorization for Z/{}", CARDINALITY);
            prime_factorize((CARDINALITY - 1) as u128)
        });
        PRIME_FACTORIZATION
            .iter()
            .map(|(factor, count)| (Self::from(*factor), *count))
    }

    fn generator() -> Self {
        static GENERATOR: LazyLock<RwLock<Option<u128>>> = LazyLock::new(|| RwLock::new(None));
        if let Some(g) = *GENERATOR.read().unwrap() {
            return Self::from(g);
        }

        log::info!("Building generator for Z/{}", Self::Q);
        let g = find_generator(
            Self::Q,
            &Self::prime_factorization()
                .map(|(factor, count)| (factor.into(), count))
                .collect::<HashMap<u128, usize>>(),
        );
        *GENERATOR.write().unwrap() = Some(g);
        return g.into();
    }

    fn unity_root(len: usize) -> Option<Self> {
        static UNITY_ROOTS: LazyLock<RwLock<HashMap<usize, u32>>> =
            LazyLock::new(|| RwLock::new(HashMap::default()));
        if let Some(root) = UNITY_ROOTS.read().unwrap().get(&len) {
            return Some(Self::from(*root));
        }
        log::info!("lettuce: unity root order {} in Z/{}", len, Self::Q);
        let root_maybe = find_unity_root(len, Self::Q, Self::generator().into());
        if let Some(root) = root_maybe {
            UNITY_ROOTS.write().unwrap().insert(len, root as u32);
            Some(Self::from(root))
        } else {
            None
        }
    }

    fn unity_root_powers(root: Self, len: usize) -> Arc<(Self, Vec<Self>, Vec<Self>)> {
        static UNITY_ROOT_POWERS: LazyLock<
            RwLock<
                HashMap<
                    (u32, usize),
                    Arc<(MilliScalarMont, Vec<MilliScalarMont>, Vec<MilliScalarMont>)>,
                >,
            >,
        > = LazyLock::new(|| RwLock::new(HashMap::default()));
        if let Some(v) = UNITY_ROOT_POWERS.read().unwrap().get(&(root.0, len)) {
            return v.clone();
        }
        let root_inv = root.inverse();
        let mut pow = Self::one();
        let powers = (0..len)
            .map(|_| {
                let out = pow;
                pow *= root;
                out
            })
            .collect::<Vec<_>>();
        let mut pow = Self::one();
        let inv_powers = (0..len)
            .map(|_| {
                let out = pow;
                pow *= root_inv;
                out
            })
            .collect::<Vec<_>>();
        UNITY_ROOT_POWERS
            .write()
            .unwrap()
            .insert((root.0, len), Arc::new((root_inv, powers, inv_powers)));
        Self::unity_root_powers(root, len)
    }
}

impl Display for MilliScalarMont {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!("{}", MONT.from_mont(self.0)))?;
        Ok(())
    }
}

impl Into<u128> for MilliScalarMont {
    fn into(self) -> u128 {
        MONT.from_mont(self.0) as u128
    }
}

impl Product for MilliScalarMont {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut out = Self::one();
        for v in iter {
            out *= v;
        }
        out
    }
}

impl Sum for MilliScalarMont {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut out = Self::zero();
        for v in iter {
            out += v;
        }
        out
    }
}

impl Add<u8> for MilliScalarMont {
    type Output = Self;
    fn add(mut self, rhs: u8) -> Self::Output {
        self += rhs;
        self
    }
}

impl AddAssign<u8> for MilliScalarMont {
    fn add_assign(&mut self, rhs: u8) {
        let m = MONT.to_mont(rhs as u32);
        self.0 = MONT.add(self.0, m);
    }
}

impl Add for MilliScalarMont {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl AddAssign for MilliScalarMont {
    fn add_assign(&mut self, rhs: Self) {
        self.0 = MONT.add(self.0, rhs.0);
    }
}

impl Sub for MilliScalarMont {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl SubAssign for MilliScalarMont {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 = MONT.add(self.0, MONT.neg(rhs.0));
    }
}

impl Mul for MilliScalarMont {
    type Output = Self;
    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl MulAssign for MilliScalarMont {
    fn mul_assign(&mut self, rhs: Self) {
        self.0 = MONT.mul(self.0, rhs.0);
    }
}

impl From<u128> for MilliScalarMont {
    fn from(value: u128) -> Self {
        Self(MONT.to_mont((value % Self::CARDINALITY) as u32))
    }
}

impl From<i32> for MilliScalarMont {
    fn from(value: i32) -> Self {
        Self::at_displacement(value)
    }
}

impl From<u32> for MilliScalarMont {
    fn from(value: u32) -> Self {
        Self(MONT.to_mont((value % CARDINALITY) as u32))
    }
}

impl From<u8> for MilliScalarMont {
    fn from(value: u8) -> Self {
        Self(MONT.to_mont(value as u32))
    }
}

impl From<u16> for MilliScalarMont {
    fn from(value: u16) -> Self {
        Self(MONT.to_mont(value as u32))
    }
}

impl From<usize> for MilliScalarMont {
    fn from(value: usize) -> Self {
        Self(MONT.to_mont((value as u128 % Self::CARDINALITY) as u32))
    }
}

impl From<u64> for MilliScalarMont {
    fn from(value: u64) -> Self {
        Self(MONT.to_mont((value as u128 % Self::CARDINALITY) as u32))
    }
}
