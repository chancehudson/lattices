use crate::*;

/// Univariate polynomial in a ring with modulus X^N + 1.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Polynomial<const N: usize, E: FieldScalar> {
    coefs: [E; N],
}

impl<const N: usize, E: FieldScalar> RingElement for Polynomial<N, E> {
    const CARDINALITY: u128 = E::CARDINALITY.pow(N as u32);

    /// Uniform randomly sample an element from the ring provided an RNG source.
    fn sample_uniform<R: Rng>(rng: &mut R) -> Self {
        Self {
            coefs: std::array::from_fn(|_| E::sample_uniform(rng)),
        }
    }
}

impl<const N: usize, E: FieldScalar> Polynomial<N, E> {
    /// Retrieve the product of monomial factors split by the provided root of unity.
    pub fn split_root(root: E) -> Result<Self> {
        if root.modpow(N as u128) != E::negone() || root.modpow(2 * N as u128) != E::one() {
            anyhow::bail!(
                "Provided root of unity {root} does not split polynomial ring of degree {N}"
            );
        }
        let mut out = Self::one();
        let mut w = root;
        for _ in 0..N {
            let mut monomial = Self::from(w * E::negone());
            monomial.coefs[1] = E::one();
            out *= monomial;
            w *= root * root;
        }
        for (i, coef) in out.coefs().enumerate() {
            println!("coef {i}: {coef}");
        }
        Ok(out)
    }

    /// Get an iterator over all coefficients. This iterator will always contain N+1 entries.
    pub fn coefs(&self) -> impl Iterator<Item = E> {
        self.coefs.iter().copied()
    }

    pub fn coefs_mut(&mut self) -> impl Iterator<Item = &mut E> {
        self.coefs.iter_mut()
    }

    /// The l2 norm of the polynomial. Represents the magnitude.
    ///
    /// The distance of each coefficient from zero is squared and then
    /// summed. The square root of the sum is returned.
    ///
    /// `sqrt(x[0]^2 + x[1]^2 + x[2]^2 ...)`
    pub fn norm_l2(&self) -> f64 {
        self.coefs()
            .map(|coef| coef.displacement().pow(2) as f64)
            .sum::<f64>()
            .sqrt()
    }

    pub fn sample_gaussian<R: Rng>(sigma: f64, rng: &mut R) -> Self {
        let cdt = GaussianCDT::cache_or_init::<E>(sigma);
        Self {
            coefs: cdt.sample_arr::<N, _, _>(rng),
        }
    }

    /// Encoding polynomial ring elements.
    ///
    /// Self described polynomial of degree N over field with cardinality Q;
    /// - 1 byte: L = byte length of Q
    /// - L bytes: Q-1
    /// - L bytes: x^0 term coefficient
    /// - L bytes: x^1 term coefficient
    /// - L bytes: x^2 term coefficient
    /// ...
    /// - L bytes: x^f s.t. f < N term coefficent (must be non-zero)
    pub fn as_le_bytes(&self) -> Vec<u8> {
        let final_nonzero_term_exp = self
            .coefs
            .iter()
            .enumerate()
            .rev()
            .find_map(|(degree, coef)| if coef.is_zero() { Some(degree) } else { None });
        if let Some(exp) = final_nonzero_term_exp {
            use std::iter::once;
            once(vec![E::BYTE_WIDTH])
                .chain(once(E::negone().as_le_bytes()))
                .chain(self.coefs().take(exp + 1).map(|coef| coef.as_le_bytes()))
                .flatten()
                .collect()
        } else {
            vec![]
        }
    }

    /// Decode a polynomial.
    pub fn from_le_bytes(bytes: &[u8]) -> Result<Self> {
        let mut out = Self::default();
        let q_byte_len = bytes[0];
        let mut e_iter = bytes[1..].chunks(q_byte_len.into());
        let negone_bytes = e_iter.next().ok_or(anyhow::anyhow!(
            "Polynomial data did not have negone value encoded"
        ))?;
        if E::from_le_bytes(negone_bytes.iter()) != E::negone() {
            anyhow::bail!("Polynomial refusing to decode a polynomial over a different field");
        }
        for (degree, bytes) in e_iter.enumerate() {
            if degree >= N {
                anyhow::bail!("Polynomial attempting to parse variable outside of ring");
            }
            out.coefs[degree] = E::from_le_bytes(bytes.iter());
        }
        Ok(out)
    }
}

impl<const N: usize, E: FieldScalar> Default for Polynomial<N, E> {
    fn default() -> Self {
        Self {
            coefs: [E::zero(); N],
        }
    }
}

impl<const N: usize, E: FieldScalar> From<[E; N]> for Polynomial<N, E> {
    fn from(coefs: [E; N]) -> Self {
        Self { coefs }
    }
}

impl<const N: usize, E: FieldScalar> From<E> for Polynomial<N, E> {
    fn from(coef: E) -> Self {
        let mut coefs = [E::zero(); N];
        coefs[0] = coef;
        Self { coefs }
    }
}

impl<const N: usize, E: FieldScalar> Display for Polynomial<N, E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("polynomial unimplemented")?;
        Ok(())
    }
}

impl<const N: usize, E: FieldScalar> Add for Polynomial<N, E> {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<const N: usize, E: FieldScalar> AddAssign for Polynomial<N, E> {
    fn add_assign(&mut self, rhs: Self) {
        for (self_coef, other_coef) in self.coefs_mut().zip(rhs.coefs()) {
            *self_coef += other_coef;
        }
    }
}

impl<const N: usize, E: FieldScalar> Sub for Polynomial<N, E> {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<const N: usize, E: FieldScalar> SubAssign for Polynomial<N, E> {
    fn sub_assign(&mut self, rhs: Self) {
        for (self_coef, other_coef) in self.coefs_mut().zip(rhs.coefs()) {
            *self_coef -= other_coef;
        }
    }
}

impl<const N: usize, E: FieldScalar> Mul for Polynomial<N, E> {
    type Output = Self;
    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<const N: usize, E: FieldScalar> MulAssign for Polynomial<N, E> {
    fn mul_assign(&mut self, rhs: Self) {
        let mut out = Self::default();
        for (deg, coef) in self.coefs().enumerate() {
            if coef.is_zero() {
                continue;
            }
            for (other_deg, other_coef) in rhs.coefs().enumerate() {
                if other_coef.is_zero() {
                    continue;
                }
                let total_deg = deg + other_deg;
                let rem = total_deg % N;
                let div = total_deg / N;
                let mut out_coef = coef * other_coef;
                if div % 2 == 1 {
                    out_coef *= E::negone();
                }
                out.coefs[rem] += out_coef;
            }
        }
        self.coefs = out.coefs;
    }
}

impl<const N: usize, E: FieldScalar> From<u128> for Polynomial<N, E> {
    fn from(value: u128) -> Self {
        let mut coefs = [E::zero(); N];
        coefs[0] = value.into();
        Self { coefs }
    }
}

impl<const N: usize, E: FieldScalar> Into<u128> for Polynomial<N, E> {
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
        let p = P::sample_uniform(rng);
        let q = P::sample_uniform(rng);
        let r = P::sample_uniform(rng);

        assert_eq!((p + q) + r, p + (q + r));
    }

    #[test]
    fn polynomial_associative_mul() {
        type Field = OxfoiScalar;
        type P = Polynomial<64, Field>;
        let rng = &mut rand::rng();
        // (p * q) * r = p * (q * r)
        let p = P::sample_uniform(rng);
        let q = P::sample_uniform(rng);
        let r = P::sample_uniform(rng);

        assert_eq!((p * q) * r, p * (q * r));
        assert_eq!((q * p) * r, p * (r * q));
        assert_eq!(r * (q * p), (r * q) * p);
    }

    #[test]
    fn polynomial_distributive() {
        type Field = OxfoiScalar;
        type P = Polynomial<2, Field>;
        let rng = &mut rand::rng();
        // p * (q + r) = p * q + p * r
        let p = P::sample_uniform(rng);
        let q = P::sample_uniform(rng);
        let r = P::sample_uniform(rng);

        assert_eq!(p * (q + r), p * q + p * r);
        assert_eq!((q + r) * p, q * p + r * p);
    }

    #[test]
    fn polynomial_identity_add() {
        type Field = OxfoiScalar;
        type P = Polynomial<64, Field>;
        let rng = &mut rand::rng();
        // additive identity
        let p = P::sample_uniform(rng);
        assert_eq!(p + P::zero(), p);
        assert_eq!(P::zero() + p, p);
    }

    #[test]
    fn polynomial_identity_mul() {
        type Field = OxfoiScalar;
        type P = Polynomial<64, Field>;
        let rng = &mut rand::rng();
        // multiplicative identity
        let p = P::sample_uniform(rng);
        assert_eq!(p * P::one(), p);
        assert_eq!(P::one() * p, p);
    }
}
