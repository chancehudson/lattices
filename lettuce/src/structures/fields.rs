use std::collections::HashMap;
use std::sync::Arc;
use std::sync::LazyLock;
use std::sync::RwLock;

use crate::*;

macro_rules! export_fields {
    () => {
        integer_prime_field!(7753, Seven753Scalar, u16, u32);
        integer_prime_field!(101, LOLScalar, u8, u16);
        integer_prime_field!(7, SevenScalar, u8, u8);
        integer_prime_field!(2, BinaryScalar, u8, u8);

        // 32 bit prime with roots of unity for power of 2 cycles up to 2^20
        // 455 * 2^20 * 3^2 + 1
        // named milli because it's ~1 million from 2^32 and supports NTT with 1m elements
        integer_prime_field!(455 * 2u128.pow(20) * 9 + 1, MilliScalar, u32, u64);
        // 3 * 2^30 + 1
        integer_prime_field!(3u128 * 2u128.pow(30) + 1, CoolScalar, u32, u64);
        // 2^31 - 1 aka Mersenne31
        integer_prime_field!(2u128.pow(31) - 1, Mersenne31Scalar, u32, u64);
        // 2^64 - 2^32 + 1 aka oxfoi aka goldilocks (ambiguous)
        integer_prime_field!(2u128.pow(64) - 2u128.pow(32) + 1, OxfoiScalar, u64, u128);
    };
}

#[test]
fn generator_test() {
    type Field = Seven753Scalar;
    let g = Field::generator();
    let mut v = g;
    for i in 0..(Field::Q - 2) {
        assert_ne!(v, Field::one(), "{i}");
        v *= g;
    }
    assert_eq!(v, Field::one());
}

/// Generate a scalar prime field implementation backed by a type.
///
/// First argument is cardinality of the field (must be prime). Second argument is the type used to
/// store field elements. Third argument is type used to store multiplication of two field
/// elements.
#[macro_export]
macro_rules! integer_prime_field {
    ($cardinality:expr, $name:ident, $data_ty:ty, $sq_data_ty:ty) => {
        // cardinality must fit in data type
        const _: () = assert!($cardinality < (<$data_ty>::MAX as u128));
        // multiplication must fit inside the next biggest type
        const _: () = assert!($cardinality * $cardinality < (<$sq_data_ty>::MAX) as u128);

        /// Automatically generate `FieldScalar` impl.
        #[derive(Debug, Copy, Clone, Default, PartialEq)]
        #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
        pub struct $name($data_ty);

        impl RingElement for $name {
            const CARDINALITY: u128 = $cardinality;
        }

        impl FieldScalar for $name {
            fn prime_factorization() -> impl Iterator<Item = (Self, usize)> {
                static PRIME_FACTORIZATION: LazyLock<HashMap<u128, usize>> = LazyLock::new(|| {
                    log::info!("Building prime factorization for Z/{}", $cardinality);
                    prime_factorize($cardinality - 1)
                });
                PRIME_FACTORIZATION
                    .iter()
                    .map(|(factor, count)| (Self::from(*factor), *count))
            }

            fn generator() -> Self {
                static GENERATOR: LazyLock<RwLock<Option<u128>>> =
                    LazyLock::new(|| RwLock::new(None));
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
                static UNITY_ROOTS: LazyLock<RwLock<HashMap<usize, $data_ty>>> =
                    LazyLock::new(|| RwLock::new(HashMap::default()));
                if let Some(root) = UNITY_ROOTS.read().unwrap().get(&len) {
                    return Some(Self::from(*root));
                }
                log::info!("lettuce: unity root order {} in Z/{}", len, Self::Q);
                let root_maybe = find_unity_root(len, Self::Q, Self::generator().into());
                if let Some(root) = root_maybe {
                    UNITY_ROOTS.write().unwrap().insert(len, root as $data_ty);
                    Some(Self::from(root))
                } else {
                    None
                }
            }

            fn unity_root_powers(root: Self, len: usize) -> Arc<(Self, Vec<Self>, Vec<Self>)> {
                static UNITY_ROOT_POWERS: LazyLock<
                    RwLock<HashMap<($data_ty, usize), Arc<($name, Vec<$name>, Vec<$name>)>>>,
                > = LazyLock::new(|| RwLock::new(HashMap::default()));
                if let Some(v) = UNITY_ROOT_POWERS.read().unwrap().get(&(root.0, len)) {
                    return v.clone();
                }
                let powers = (0..len).map(|i| root.modpow(i as u128)).collect::<Vec<_>>();
                let root_inv = root.inverse();
                let inv_powers = (0..len)
                    .map(|i| root_inv.modpow(i as u128))
                    .collect::<Vec<_>>();
                UNITY_ROOT_POWERS
                    .write()
                    .unwrap()
                    .insert((root.0, len), Arc::new((root_inv, powers, inv_powers)));
                Self::unity_root_powers(root, len)
            }
        }

        impl Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.write_str(&format!("{}", self.0))?;
                Ok(())
            }
        }

        impl Into<u128> for $name {
            fn into(self) -> u128 {
                self.0.into()
            }
        }

        impl Product for $name {
            fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
                let mut out = Self::one();
                for v in iter {
                    out *= v;
                }
                out
            }
        }

        impl Sum for $name {
            fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
                let mut out = Self::zero();
                for v in iter {
                    out += v;
                }
                out
            }
        }

        impl Add<u8> for $name {
            type Output = Self;
            fn add(mut self, rhs: u8) -> Self::Output {
                self += rhs;
                self
            }
        }

        impl AddAssign<u8> for $name {
            fn add_assign(&mut self, rhs: u8) {
                self.0 = ((self.0 as $sq_data_ty + rhs as $sq_data_ty)
                    % $cardinality as $sq_data_ty) as $data_ty;
            }
        }

        impl Add for $name {
            type Output = Self;
            fn add(mut self, rhs: Self) -> Self::Output {
                self += rhs;
                self
            }
        }

        impl AddAssign for $name {
            fn add_assign(&mut self, rhs: Self) {
                self.0 = ((self.0 as $sq_data_ty + rhs.0 as $sq_data_ty)
                    % $cardinality as $sq_data_ty) as $data_ty;
            }
        }

        impl Sub for $name {
            type Output = Self;
            fn sub(mut self, rhs: Self) -> Self::Output {
                self -= rhs;
                self
            }
        }

        impl SubAssign for $name {
            fn sub_assign(&mut self, rhs: Self) {
                let neg_rhs = $cardinality as $data_ty - rhs.0;
                *self += Self(neg_rhs);
            }
        }

        impl Mul for $name {
            type Output = Self;
            fn mul(mut self, rhs: Self) -> Self::Output {
                self *= rhs;
                self
            }
        }

        impl MulAssign for $name {
            fn mul_assign(&mut self, rhs: Self) {
                self.0 = (((self.0 as $sq_data_ty) * (rhs.0 as $sq_data_ty))
                    % $cardinality as $sq_data_ty) as $data_ty;
            }
        }

        impl From<u128> for $name {
            fn from(value: u128) -> Self {
                Self((value % $cardinality) as $data_ty)
            }
        }

        impl From<i32> for $name {
            fn from(value: i32) -> Self {
                Self::at_displacement(value)
            }
        }

        impl From<u32> for $name {
            fn from(value: u32) -> Self {
                Self((value as u128 % $cardinality) as $data_ty)
            }
        }

        impl From<u8> for $name {
            fn from(value: u8) -> Self {
                Self((value as u128 % $cardinality) as $data_ty)
            }
        }

        impl From<u16> for $name {
            fn from(value: u16) -> Self {
                Self((value as u128 % $cardinality) as $data_ty)
            }
        }

        impl From<usize> for $name {
            fn from(value: usize) -> Self {
                Self((value as u128 % $cardinality) as $data_ty)
            }
        }

        impl From<u64> for $name {
            fn from(value: u64) -> Self {
                Self((value as u128 % $cardinality) as $data_ty)
            }
        }
    };
}

export_fields!();
