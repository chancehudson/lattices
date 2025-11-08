use crate::*;

/// Univariate polynomial in a ring with modulus X^N + 1.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Polynomial<const N: usize, E: FieldScalar> {
    coefs: [E; N],
    is_coef_form: bool,
}

impl<const N: usize, E: FieldScalar> Polynomial<N, E> {
    pub fn to_eval_form(&mut self) {
        if self.is_coef_form {
            ntt_negacyclic(self.coefs_slice_mut()).unwrap();
            self.is_coef_form = false;
        }
    }

    pub fn to_coef_form(&mut self) {
        if !self.is_coef_form {
            intt_negacyclic(self.coefs_slice_mut()).unwrap();
            self.is_coef_form = true;
        }
    }
    pub fn into_eval_form(mut self) -> Self {
        self.to_eval_form();
        self
    }

    pub fn into_coef_form(mut self) -> Self {
        self.to_coef_form();
        self
    }

    /// Evaluate the polynomial at a point.
    pub fn evaluate(&self, x: E) -> E {
        log::debug!("evaluate Z_{}, x = {x}", E::Q);
        let mut out = E::zero();
        for (i, coef) in self.coefs().enumerate() {
            if coef.is_zero() {
                continue;
            }
            out += coef * x.modpow(i as u128);
        }
        log::debug!("{out} = {self}");
        out
    }

    /// Multiply a single polynomial by a vector of polynomials.
    pub fn batch_mul(mut self, rhs: Vector<Self>) -> Vector<Self> {
        self.to_eval_form();
        rhs.into_iter()
            .map(|mut v| {
                v.to_eval_form();
                for (l, r) in v.coefs_mut().zip(self.coefs()) {
                    *l *= r;
                }
                v
            })
            .collect()
    }

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
        Ok(out)
    }

    pub fn coefs_slice_mut(&mut self) -> &mut [E; N] {
        &mut self.coefs
    }

    /// Get an iterator over all coefficients.
    ///
    /// TODO: coefs_nonzero ?
    pub fn coefs(&self) -> impl Iterator<Item = E> + ExactSizeIterator + Clone {
        self.coefs.iter().copied()
    }

    /// Get a mutable iterator over all coefficients.
    pub fn coefs_mut(&mut self) -> impl Iterator<Item = &mut E> + ExactSizeIterator {
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

    /// Sample a polynomial with coefficients in a gaussian distribution
    /// with standard deviation `sigma`.
    ///
    /// Internally uses a statically cached distribution table for each stddev
    /// initialized on first call.
    pub fn sample_gaussian<R: Rng>(sigma: f64, rng: &mut R) -> Self {
        let cdt = GaussianCDT::cache_or_init::<E>(sigma);
        Self {
            coefs: cdt.sample_arr::<N, _, _>(rng),
            is_coef_form: true,
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

impl<const N: usize, E: FieldScalar> RingElement for Polynomial<N, E> {
    const CARDINALITY: u128 = E::CARDINALITY.pow(N as u32);

    /// Uniform randomly sample an element from the ring provided an RNG source.
    fn sample_uniform<R: Rng>(rng: &mut R) -> Self {
        Self {
            coefs: std::array::from_fn(|_| E::sample_uniform(rng)),
            is_coef_form: true,
        }
    }
}

impl<const N: usize, E: FieldScalar> Default for Polynomial<N, E> {
    fn default() -> Self {
        Self {
            coefs: [E::zero(); N],
            is_coef_form: true,
        }
    }
}

impl<const N: usize, E: FieldScalar> From<&Vector<E>> for Polynomial<N, E> {
    fn from(coefs: &Vector<E>) -> Self {
        assert!(coefs.len() <= N);
        Self {
            coefs: std::array::from_fn(|i| coefs.get(i).copied().unwrap_or_default()),
            is_coef_form: true,
        }
    }
}

impl<const N: usize, E: FieldScalar> From<[E; N]> for Polynomial<N, E> {
    fn from(coefs: [E; N]) -> Self {
        Self {
            coefs,
            is_coef_form: true,
        }
    }
}

impl<const N: usize, E: FieldScalar> From<E> for Polynomial<N, E> {
    fn from(coef: E) -> Self {
        let mut coefs = [E::zero(); N];
        coefs[0] = coef;
        Self {
            coefs,
            is_coef_form: true,
        }
    }
}

impl<const N: usize, E: FieldScalar> Display for Polynomial<N, E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let poly_str = self
            .coefs()
            .enumerate()
            .filter_map(|(i, coef)| {
                if coef.is_zero() {
                    None
                } else if i == 0 {
                    Some(format!("{coef}"))
                } else if i == 1 {
                    Some(format!("{coef}x"))
                } else {
                    Some(format!("{}x^{}", coef, i))
                }
            })
            .collect::<Vec<String>>()
            .join(" + ");
        f.write_str(&poly_str)?;
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
    fn add_assign(&mut self, mut rhs: Self) {
        self.to_coef_form();
        rhs.to_coef_form();
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
    fn sub_assign(&mut self, mut rhs: Self) {
        self.to_coef_form();
        rhs.to_coef_form();
        for (self_coef, other_coef) in self.coefs_mut().zip(rhs.coefs()) {
            *self_coef -= other_coef;
        }
    }
}

impl<const N: usize, E: FieldScalar> Mul<E> for Polynomial<N, E> {
    type Output = Self;
    fn mul(mut self, rhs: E) -> Self::Output {
        self.to_coef_form();
        self *= rhs;
        self
    }
}

impl<const N: usize, E: FieldScalar> MulAssign<E> for Polynomial<N, E> {
    fn mul_assign(&mut self, rhs: E) {
        for coef in self.coefs_mut() {
            *coef *= rhs;
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
    fn mul_assign(&mut self, mut rhs: Self) {
        if self.is_coef_form {
            ntt_negacyclic::<N, _>(self.coefs_slice_mut()).unwrap();
        }
        if rhs.is_coef_form {
            ntt_negacyclic::<N, _>(rhs.coefs_slice_mut()).unwrap();
        }
        for (l, r) in self.coefs_mut().zip(rhs.coefs()) {
            *l *= r;
        }
        if self.is_coef_form {
            intt_negacyclic(self.coefs_slice_mut()).unwrap();
        }

        // let mut out = Self::default();
        // for (deg, coef) in self.coefs().enumerate() {
        //     if coef.is_zero() {
        //         continue;
        //     }
        //     for (other_deg, other_coef) in rhs.coefs().enumerate() {
        //         if other_coef.is_zero() {
        //             continue;
        //         }
        //         let total_deg = deg + other_deg;
        //         let rem = total_deg % N;
        //         let div = total_deg / N;
        //         let mut out_coef = coef * other_coef;
        //         if div % 2 == 1 {
        //             out_coef *= E::negone();
        //         }
        //         out.coefs[rem] += out_coef;
        //     }
        // }
        // self.coefs = out.coefs;
    }
}

impl<const N: usize, E: FieldScalar> Product for Polynomial<N, E> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut out = Self::one();
        for v in iter {
            out *= v;
        }
        out
    }
}

impl<const N: usize, E: FieldScalar> Sum for Polynomial<N, E> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut out = Self::zero();
        for v in iter {
            out += v;
        }
        out
    }
}

impl<const N: usize, E: FieldScalar> From<u128> for Polynomial<N, E> {
    fn from(value: u128) -> Self {
        let mut coefs = [E::zero(); N];
        coefs[0] = value.into();
        Self {
            coefs,
            is_coef_form: true,
        }
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
        type Field = MilliScalarMont;
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
        type Field = MilliScalarMont;
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
        type Field = MilliScalarMont;
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
        type Field = MilliScalarMont;
        type P = Polynomial<64, Field>;
        let rng = &mut rand::rng();
        // additive identity
        let p = P::sample_uniform(rng);
        assert_eq!(p + P::zero(), p);
        assert_eq!(P::zero() + p, p);
    }

    #[test]
    fn polynomial_identity_mul() {
        type Field = MilliScalarMont;
        type P = Polynomial<64, Field>;
        let rng = &mut rand::rng();
        // multiplicative identity
        let p = P::sample_uniform(rng);
        assert_eq!(p * P::one(), p);
        assert_eq!(P::one() * p, p);
    }

    #[test]
    fn polynomial_evaluate() {
        let rng = &mut rand::rng();
        type Field = MilliScalarMont;

        let poly = Polynomial::<64, Field>::sample_uniform(rng);

        let out = poly.evaluate(Field::zero());
        assert_eq!(out, poly.coefs().next().unwrap());

        let out = poly.evaluate(Field::one());
        assert_eq!(out, poly.coefs().sum());

        let mut poly = Polynomial::<1, Field>::zero();
        for _ in 0..1000 {
            let x = Field::sample_uniform(rng);
            let coef_i = rng.random_range(0..1);
            println!("i: {}", coef_i);
            let coef_delta = Field::sample_uniform(rng);
            println!("coef_delta: {}", coef_delta);
            *poly.coefs_mut().skip(coef_i).next().unwrap() += coef_delta;
            let expected: Field = poly
                .coefs()
                .enumerate()
                .map(|(i, coef)| coef * x.modpow(i as u128))
                .sum();
            assert_eq!(expected, poly.evaluate(x));
        }
    }
}
