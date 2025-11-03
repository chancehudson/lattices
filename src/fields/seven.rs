use rand::distr::Distribution;
use rand::distr::StandardUniform;

use crate::*;

const F: u8 = 7;

#[derive(Debug, Copy, Clone, Default, PartialEq)]
pub struct SevenScalar {
    val: u8,
}

impl Distribution<SevenScalar> for StandardUniform {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> SevenScalar {
        SevenScalar {
            val: rng.random::<u8>() % F,
        }
    }
}

impl From<u8> for SevenScalar {
    fn from(value: u8) -> Self {
        Self { val: value % F }
    }
}

impl From<u128> for SevenScalar {
    fn from(value: u128) -> Self {
        Self {
            val: ((value % (F as u128)) as u8),
        }
    }
}

impl Into<u128> for SevenScalar {
    fn into(self) -> u128 {
        self.val.into()
    }
}

impl Element for SevenScalar {
    const CARDINALITY: u128 = F as u128;
    const BIT_WIDTH: usize = 8;

    fn is_zero(&self) -> bool {
        self.val == 0
    }

    fn sample_rand<R: Rng>(rng: &mut R) -> Self {
        rng.random()
    }
}

impl Display for SevenScalar {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!("{}", self.val))?;
        Ok(())
    }
}

impl Add for SevenScalar {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl AddAssign for SevenScalar {
    fn add_assign(&mut self, rhs: Self) {
        self.val = (self.val + rhs.val) % F;
    }
}

impl Sub for SevenScalar {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl SubAssign for SevenScalar {
    fn sub_assign(&mut self, rhs: Self) {
        self.val = ((self.val + F) - rhs.val) % F;
    }
}

impl Mul for SevenScalar {
    type Output = Self;
    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl MulAssign for SevenScalar {
    fn mul_assign(&mut self, rhs: Self) {
        self.val = (self.val * rhs.val) % F;
    }
}
