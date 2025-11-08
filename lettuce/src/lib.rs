mod commitments;
pub(crate) mod log;
mod montgomery;
mod ntt;
mod precomputed;
mod probability;
mod structures;

#[cfg(test)]
mod test;

pub use commitments::*;
pub use montgomery::*;
pub use ntt::*;
pub use precomputed::*;
pub use probability::*;
pub use structures::*;

use std::fmt::Display;
use std::iter::Product;
use std::iter::Sum;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Sub;
use std::ops::SubAssign;
use std::sync::Arc;

use anyhow::Result;
use rand::Rng;

/// An element of a commutative ring.
pub trait RingElement:
    Sized
    + Copy
    + Clone
    + Default
    + Display
    + Add<Output = Self>
    + AddAssign
    + Sub<Output = Self>
    + SubAssign
    + Mul<Output = Self>
    + MulAssign
    + PartialEq
    + From<u128>
    + Sum
    + Product
    + Sync
    + Send
{
    const CARDINALITY: u128;
    const Q: u128 = Self::CARDINALITY;

    /// Is the element the additive identity?
    fn is_zero(&self) -> bool {
        *self == Self::default()
    }

    /// Multiplicative identity.
    fn one() -> Self {
        Self::from(1u128)
    }

    /// Multiplicatively flips a value in centered representation.
    /// Additive inverse of the multiplicative identity.
    /// (Negative one)
    fn negone() -> Self {
        Self::from(Self::Q - 1)
    }

    /// Additive identity.
    fn zero() -> Self {
        Self::default()
    }

    /// Uniform randomly sample an element from the ring provided an RNG source.
    fn sample_uniform<R: Rng>(rng: &mut R) -> Self {
        Self::from(rng.random::<u128>())
    }
}

/// A scalar field element. All `FieldScalar` are also `RingElement`.
pub trait FieldScalar:
    RingElement
    + Into<u128>
    + Add<u8, Output = Self>
    + AddAssign<u8>
    + From<i32>
    + From<u8>
    + From<u16>
    + From<u32>
    + From<usize>
    + From<u64>
{
    /// Find a generator element. This is a primitive root of unity with
    /// cycle length equal to Self::Q - 1. Exponentiating moves through the field
    /// in sequence ending with 1 (multiplicative identity). 0 (additive identity)
    /// is excluded from the sequence. This is the multiplicative group Z_q^*
    ///
    /// in Z_3
    /// g = generator() = 2
    /// g * g = 1
    /// g * g * g = 2
    /// g * g * g * g = 2
    /// g * g * g * g * g = 1
    /// g^6 = 2
    /// g^7 = 1
    /// g^900 = 1 (maybe)
    fn generator() -> Self;

    /// Compute the modular inverse of an element.
    ///
    /// Given x, and inverse x_i, x * x_i = 1
    fn inverse(&self) -> Self {
        self.modpow(Self::Q - 2)
    }

    /// Find a root of unity with a certain cycle length, if it exists.
    ///
    /// A root of unity is a cyclic group over the field. Exponentiating by `len`
    /// yields 1 (multiplicative identity)
    ///
    /// root_3 = unity_root(3)
    ///   1 = root_3 * root_3 * root_3
    fn unity_root(len: usize) -> Option<Self>;

    fn unity_root_powers(root: Self, len: usize) -> Arc<(Self, Vec<Self>, Vec<Self>)>;

    /// Return an iterator over the prime factorization of a field element. Items are factors
    /// paired with the number of times the factor occurs.
    ///
    /// e.g. prime factorization of 100 = (5, 2), (2, 2)
    fn prime_factorization() -> impl Iterator<Item = (Self, usize)>;

    /// Sample from a discrete gaussian distribution with standard deviation sigma.
    ///
    /// Internally uses a statically cached cumulative distribution table (CDT). Table is
    /// lazily computed, so first access will incur many f64 operations. Subsequent accesses are
    /// approximately 10 * sigma f64 comparisons.
    fn sample_gaussian<R: Rng>(sigma: f64, rng: &mut R) -> Self {
        GaussianCDT::cache_or_init::<Self>(sigma).sample(rng)
    }

    /// Return the finite field element at a certain displacement.
    fn at_displacement(disp: i32) -> Self {
        if disp.abs() as u128 > Self::CARDINALITY / 2 {
            log::info!(
                "Attempting to initialize a displacement outside the field: {} {}",
                disp,
                Self::CARDINALITY
            );
            #[cfg(not(debug_assertions))]
            panic!("refusing to use displacement outside of field in production");
        }
        if disp >= 0 {
            Self::from(disp as u128)
        } else {
            Self::from(Self::CARDINALITY - disp.abs() as u128)
        }
    }

    /// Determine the displacement of an element from the zero element. In a Z_q field, if this element
    /// is > q/2 returns a negative value q - self.
    ///
    /// Distance is a measurement, and so not a field element.
    fn displacement(self) -> i128 /* TODO: <- i32 */ {
        // distance from the zero element in the positive dimension only
        let dist: u128 = self.into();
        if dist > (Self::CARDINALITY / 2) {
            -((Self::CARDINALITY - dist) as i128)
        } else {
            dist as i128
        }
    }

    /// Raise an element of the field to a power `exp`. `O(log(exp))` multiplications.
    fn modpow(mut self, mut exp: u128) -> Self {
        let mut out = Self::one();
        loop {
            if exp % 2 == 1 {
                out *= self;
            }
            exp >>= 1;
            if exp == 0 {
                break;
            }
            self *= self;
        }
        out
    }

    /// Number of bits necessary to represent an element.
    const BIT_WIDTH: usize = (Self::CARDINALITY.ilog2() + 1) as usize;
    /// Number of bytes necessary to represent an element.
    const BYTE_WIDTH: u8 =
        (Self::BIT_WIDTH / 8 + if Self::BIT_WIDTH % 8 > 0 { 1 } else { 0 }) as u8;

    /// Determine either number of 2^bits elements in a single element, or upper bound of each
    /// chunked element given `bits` chunks.
    fn bits_vec_len(bits: usize) -> usize {
        Self::BIT_WIDTH.div_ceil(bits)
    }

    /// Break into `bits` field elements. Returns `ceil(log2(F)) / 8` field elements, each
    /// containing a value up to `2^bits`.
    fn as_le_bits_vec(&self, bits: usize) -> Vector<Self> {
        let parts_len = Self::BIT_WIDTH.div_ceil(bits);
        let divisor = 1 << bits;
        let mut v: u128 = (*self).into();
        let mut out = Vector::new(parts_len.try_into().expect("base too large"));
        for i in 0..parts_len {
            if v == 0 {
                break;
            }
            let part = v % divisor;
            out[i] = part.into();
            v >>= bits;
        }
        assert_eq!(v, 0);
        out
    }

    /// Take `parts.len()` field elements each at most `2^parts.len()` and convert them into a
    /// single element.
    fn from_le_bits_vec(parts: Vector<Self>) -> Self {
        let bits_len = Self::bits_vec_len(parts.len());
        let mut mult = 1u128 << bits_len;
        let mut out = Self::default();
        for part in parts {
            out += part * mult.into();
            mult <<= bits_len;
        }
        out
    }

    /// Return an element as a vector of bytes. Individual implementations may want to provide
    /// optimized implementations.
    fn as_le_bytes(&self) -> Vec<u8> {
        let bits = 8;
        let parts_len = Self::BIT_WIDTH.div_ceil(bits);
        let divisor = 1 << bits;
        let mut v: u128 = (*self).into();
        let mut out = Vec::with_capacity(parts_len.into());
        for _ in 0..parts_len {
            out.push(0u8);
        }
        for i in 0..parts_len {
            if v == 0 {
                break;
            }
            let part = v % divisor;
            out[i] = part as u8;
            v >>= bits;
        }
        assert_eq!(v, 0);
        out
    }

    /// Parse a field element from a vector of bytes. Panics if the parsed value is greater than or
    /// equal to the field modulus.
    fn from_le_bytes<'a>(bytes: impl Iterator<Item = &'a u8>) -> Self {
        let mut v = 0u128;
        for (i, byte) in bytes.enumerate() {
            v += (*byte as u128) << i;
        }
        assert!(v < Self::CARDINALITY);
        v.into()
    }
}
