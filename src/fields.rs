use crate::*;

macro_rules! export_fields {
    () => {
        integer_prime_field!(7753, Seven753Scalar, u16, u32);
        integer_prime_field!(101, LOLScalar, u8, u16);
        integer_prime_field!(7, SevenScalar, u8, u16);

        // 2^64 - 2^32 + 1 aka oxfoi aka goldilocks (incorrect)
        integer_prime_field!(2u128.pow(64) - 2u128.pow(32) + 1, OxfoiScalar, u64, u128);
    };
}

/// Generate a Scalar implementation backed by a u32.
macro_rules! integer_prime_field {
    ($cardinality:expr, $name:ident, $data_ty:ty, $sq_data_ty:ty) => {
        #[derive(Debug, Copy, Clone, Default, PartialEq)]
        pub struct $name {
            // absolute value of displacement
            norm: $data_ty,
        }

        impl RingElement for $name {
            const CARDINALITY: u128 = $cardinality;
            // const BYTE_LEN: usize =
        }

        impl FieldScalar for $name {}

        impl From<u128> for $name {
            fn from(value: u128) -> Self {
                // allow addition inside the base type
                assert!($cardinality < (<$data_ty>::MAX as u128));
                // allow multiplication inside the next biggest type
                assert!($cardinality * $cardinality < (<$sq_data_ty>::MAX) as u128);
                Self {
                    norm: ((value % $cardinality) as $data_ty),
                }
            }
        }

        impl Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.write_str(&format!("{}", self.norm))?;
                Ok(())
            }
        }

        impl From<$data_ty> for $name {
            fn from(value: $data_ty) -> Self {
                Self {
                    norm: value % $cardinality as $data_ty,
                }
            }
        }

        impl Into<u128> for $name {
            fn into(self) -> u128 {
                self.norm.into()
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
                self.norm = ((self.norm as $sq_data_ty + rhs as $sq_data_ty)
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
                self.norm = ((self.norm as $sq_data_ty + rhs.norm as $sq_data_ty)
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
                let neg_rhs = $cardinality as $data_ty - rhs.norm;
                *self += Self { norm: neg_rhs };
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
                self.norm = (((self.norm as $sq_data_ty) * (rhs.norm as $sq_data_ty))
                    % $cardinality as $sq_data_ty) as $data_ty;
            }
        }
    };
}

export_fields!();
