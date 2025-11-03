mod commitments;
mod fields;
mod matrix;
mod polynomial;
mod probability;
mod vector;

#[cfg(test)]
mod test;

use commitments::*;
use fields::*;
use matrix::*;
use polynomial::*;
use probability::*;
use vector::*;

use std::fmt::Display;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Index;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Sub;
use std::ops::SubAssign;
use std::sync::LazyLock;

use anyhow::Result;
use rand::Rng;

/// The default value should be the additive identity.
pub trait Element:
    Sized
    + Copy
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
    + Into<u128>
{
    const BIT_WIDTH: usize;
    const CARDINALITY: u128;

    /// Is the element the additive identity?
    fn is_zero(&self) -> bool;

    /// Multiplicative identity.
    fn one() -> Self {
        Self::from(1)
    }

    /// Multiplicatively flips a value in centered representation.
    /// Additive inverse of the multiplicative identity.
    /// (Negative one)
    fn negone() -> Self {
        Self::zero() - Self::one()
    }

    /// Additive identity.
    fn zero() -> Self {
        Self::from(0)
    }

    /// Return the finite field element at a certain displacement.
    fn at_displacement(disp: i32) -> Self {
        if disp.abs() as u128 > Self::CARDINALITY / 2 {
            log::error!(
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
    /// is > q/2 returns the negated value of the element.
    ///
    /// Distance is a measurement, and so not a field element.
    fn displacement(self) -> i128 {
        // distance from the zero element in the positive dimension only
        let dist: u128 = self.into();
        if dist > (Self::CARDINALITY / 2) {
            -((Self::CARDINALITY - dist) as i128)
        } else {
            dist as i128
        }
    }

    fn sample_rand<R: Rng>(rng: &mut R) -> Self;

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
        let mut out = Vec::with_capacity(parts_len);
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
}

pub struct R1CS<E: Element> {
    a: Matrix<E>,
    b: Matrix<E>,
    c: Matrix<E>,
}

impl<E: Element> R1CS<E> {
    pub fn identity(height: usize, width: usize) -> Self {
        let v = Matrix::zero(height, width);
        Self {
            a: v.clone(),
            b: v.clone(),
            c: v.clone(),
        }
    }

    pub fn eval(&self, witness: &Vector<E>) -> Result<Vector<E>> {
        self.assert_consistency()?;

        let ab = (self.a.clone() * witness) * &(self.b.clone() * witness);
        let c = self.c.clone() * witness;

        Ok(ab - c)
    }

    pub fn dimension(&self) -> (usize, usize) {
        self.a.dimension()
    }

    fn assert_consistency(&self) -> Result<()> {
        let dimension = self.a.dimension();
        if self.b.dimension() != dimension {
            anyhow::bail!(
                "R1CS A and B dimension mismatch, expected {:?}, got {:?}",
                dimension,
                self.b.dimension()
            );
        }
        if self.c.dimension() != dimension {
            anyhow::bail!(
                "R1CS A and C dimension mismatch, expected {:?}, got {:?}",
                dimension,
                self.c.dimension()
            );
        }
        Ok(())
    }
}
