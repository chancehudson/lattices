use crate::*;

/// Univariate polynomial in a ring with modulus X^N + 1.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Polynomial<const N: usize, E: Element> {
    coefs: [E; N],
}

impl<const N: usize, E: Element> Polynomial<N, E> {
    /// Get an iterator over all coefficients. This iterator will always contain N+1 entries.
    pub fn coefs(&self) -> impl Iterator<Item = E> {
        self.coefs.iter().copied()
    }

    pub fn coefs_mut(&mut self) -> impl Iterator<Item = &mut E> {
        self.coefs.iter_mut()
    }
}

impl<const N: usize, E: Element> Default for Polynomial<N, E> {
    fn default() -> Self {
        Self {
            coefs: [E::zero(); N],
        }
    }
}

impl<const N: usize, E: Element> Element for Polynomial<N, E> {
    const BIT_WIDTH: usize = N * E::BIT_WIDTH;
    const CARDINALITY: u128 = E::CARDINALITY.pow((N + 1) as u32);

    fn is_zero(&self) -> bool {
        for coef in self.coefs() {
            if !coef.is_zero() {
                return false;
            }
        }
        true
    }

    fn sample_rand<R: Rng>(rng: &mut R) -> Self {
        Self {
            coefs: std::array::from_fn(|_| E::sample_rand(rng)),
        }
    }
}

impl<const N: usize, E: Element> Display for Polynomial<N, E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("polynomial unimplemented")?;
        Ok(())
    }
}

impl<const N: usize, E: Element> Add for Polynomial<N, E> {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<const N: usize, E: Element> AddAssign for Polynomial<N, E> {
    fn add_assign(&mut self, rhs: Self) {
        for (self_coef, other_coef) in self.coefs_mut().zip(rhs.coefs()) {
            *self_coef += other_coef;
        }
    }
}

impl<const N: usize, E: Element> Sub for Polynomial<N, E> {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<const N: usize, E: Element> SubAssign for Polynomial<N, E> {
    fn sub_assign(&mut self, rhs: Self) {
        for (self_coef, other_coef) in self.coefs_mut().zip(rhs.coefs()) {
            *self_coef -= other_coef;
        }
    }
}

impl<const N: usize, E: Element> Mul for Polynomial<N, E> {
    type Output = Self;
    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<const N: usize, E: Element> MulAssign for Polynomial<N, E> {
    fn mul_assign(&mut self, rhs: Self) {
        let mut out = Self::default();
        for (deg, coef) in self.coefs().enumerate() {
            for (other_deg, other_coef) in rhs.coefs().enumerate() {
                let total_deg = deg + other_deg;
                let rem = total_deg % N;
                let div = total_deg / N;
                let mut out_coef = coef * other_coef;
                if out_coef.is_zero() {
                    continue;
                }
                if div % 2 == 1 {
                    out_coef *= E::negone();
                }
                out.coefs[rem] += out_coef;
            }
        }
        self.coefs = out.coefs;
    }
}

impl<const N: usize, E: Element> From<u128> for Polynomial<N, E> {
    fn from(value: u128) -> Self {
        let mut coefs = [E::zero(); N];
        coefs[0] = value.into();
        Self { coefs }
    }
}

impl<const N: usize, E: Element> Into<u128> for Polynomial<N, E> {
    fn into(self) -> u128 {
        let val = self.coefs[0].into();
        for coef in &self.coefs[1..] {
            assert!(coef.is_zero());
        }
        val
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn polynomial_associative_add() {
        type Field = OxfoiScalar;
        type P = Polynomial<64, Field>;
        let rng = &mut rand::rng();

        // (p + q) + r = p + (q + r)
        let p = P::sample_rand(rng);
        let q = P::sample_rand(rng);
        let r = P::sample_rand(rng);

        assert_eq!((p + q) + r, p + (q + r));
    }

    #[test]
    fn polynomial_associative_mul() {
        type Field = OxfoiScalar;
        type P = Polynomial<64, Field>;
        let rng = &mut rand::rng();
        // (p * q) * r = p * (q * r)
        let p = P::sample_rand(rng);
        let q = P::sample_rand(rng);
        let r = P::sample_rand(rng);

        assert_eq!((p * q) * r, p * (q * r));
        assert_eq!((q * p) * r, p * (r * q));
        assert_eq!(r * (q * p), (r * q) * p);
    }

    #[test]
    fn polynomial_distributive() {
        type Field = OxfoiScalar;
        type P = Polynomial<64, Field>;
        let rng = &mut rand::rng();
        // p * (q + r) = p * q + p * r
        let p = P::sample_rand(rng);
        let q = P::sample_rand(rng);
        let r = P::sample_rand(rng);

        assert_eq!(p * (q + r), p * q + p * r);
        assert_eq!((q + r) * p, q * p + r * p);
    }

    #[test]
    fn polynomial_identity_add() {
        type Field = OxfoiScalar;
        type P = Polynomial<64, Field>;
        let rng = &mut rand::rng();
        // additive identity
        let p = P::sample_rand(rng);
        assert_eq!(p + P::zero(), p);
        assert_eq!(P::zero() + p, p);
    }

    #[test]
    fn polynomial_identity_mul() {
        type Field = OxfoiScalar;
        type P = Polynomial<64, Field>;
        let rng = &mut rand::rng();
        // multiplicative identity
        let p = P::sample_rand(rng);
        assert_eq!(p * P::one(), p);
        assert_eq!(P::one() * p, p);
    }
}
