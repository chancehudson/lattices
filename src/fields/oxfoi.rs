use rand::distr::Distribution;
use rand::distr::StandardUniform;

use crate::*;

/// 2^64 - 2^32 + 1
const F: u128 = 18446744069414584321u128;

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct OxfoiScalar {
    val: u128,
}

impl Distribution<OxfoiScalar> for StandardUniform {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> OxfoiScalar {
        OxfoiScalar {
            val: rng.random::<u128>() % F,
        }
    }
}

impl From<u128> for OxfoiScalar {
    fn from(value: u128) -> Self {
        Self { val: value % F }
    }
}

impl Into<u128> for OxfoiScalar {
    fn into(self) -> u128 {
        self.val
    }
}

impl Element for OxfoiScalar {
    const CARDINALITY: u128 = F as u128;
    const BIT_WIDTH: usize = 64;

    fn is_zero(&self) -> bool {
        self.val == 0u128
    }

    fn sample_rand<R: Rng>(rng: &mut R) -> Self {
        rng.random()
    }
}

impl Display for OxfoiScalar {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!("{}", self.val))?;
        Ok(())
    }
}

impl Add for OxfoiScalar {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl AddAssign for OxfoiScalar {
    fn add_assign(&mut self, rhs: Self) {
        self.val = (self.val + rhs.val) % F;
    }
}

impl Sub for OxfoiScalar {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl SubAssign for OxfoiScalar {
    fn sub_assign(&mut self, rhs: Self) {
        self.val = ((self.val + F) - rhs.val) % F;
    }
}

impl Mul for OxfoiScalar {
    type Output = Self;
    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl MulAssign for OxfoiScalar {
    fn mul_assign(&mut self, rhs: Self) {
        self.val = (self.val * rhs.val) % F;
    }
}
