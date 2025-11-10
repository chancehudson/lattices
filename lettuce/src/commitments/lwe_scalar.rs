use crate::*;

use anyhow::Result;

/// A toy implementation of LWE commitments. Use at your own risk.
#[derive(Clone, Debug)]
pub struct LWEScalar<E: FieldScalar> {
    lattice: Matrix<E>,
    commitment: Vector<E>,
}

impl<E: FieldScalar> LWEScalar<E> {
    pub fn lattice_for<R: Rng>(element_len: usize, rng: &mut R) -> Matrix<E> {
        // m value
        let height: usize = element_len * E::BIT_WIDTH;
        Matrix::<E>::random(height, element_len, rng)
    }

    pub fn commit<R: Rng>(val: Vector<E>, lattice: Matrix<E>, rng: &mut R) -> Self {
        let (height, _width) = lattice.dimension();
        let mut err = Vector::new(height);
        for i in 0..height {
            // generate a value between 0 and 2
            let v = rng.random_range(0..=2);
            // move it to the range -1..1 in the field
            err[i] = E::from(v) - E::one();
        }
        let commitment = &lattice * &val + &err;
        Self {
            lattice,
            commitment,
        }
    }

    /// Attempt to open a commitment to a value, with each error less than `max_err` distance from zero. If successful returns the error vector.
    pub fn try_open(&self, val: &Vector<E>, max_err: u128) -> Result<Vector<E>> {
        let maybe_committed_no_err = &self.lattice * val;
        let err = &self.commitment - maybe_committed_no_err;
        for e in err.iter() {
            let disp = e.displacement();
            if disp.unsigned_abs() > max_err {
                anyhow::bail!(
                    "Error opening LWE commitment, error vector contains element {} beyond displacement bound {}",
                    disp,
                    max_err
                );
            }
        }
        Ok(err)
    }
}

impl<E: FieldScalar> Sub<&Self> for LWEScalar<E> {
    type Output = Self;
    fn sub(mut self, rhs: &Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<E: FieldScalar> SubAssign<&Self> for LWEScalar<E> {
    fn sub_assign(&mut self, rhs: &Self) {
        self.commitment -= &rhs.commitment;
    }
}

impl<E: FieldScalar> Add<&Self> for LWEScalar<E> {
    type Output = Self;
    fn add(mut self, rhs: &Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<E: FieldScalar> AddAssign<&Self> for LWEScalar<E> {
    fn add_assign(&mut self, rhs: &Self) {
        self.commitment += &rhs.commitment;
    }
}

impl<E: FieldScalar> Mul<E> for LWEScalar<E> {
    type Output = Self;
    fn mul(mut self, rhs: E) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<E: FieldScalar> MulAssign<E> for LWEScalar<E> {
    fn mul_assign(&mut self, rhs: E) {
        self.commitment *= rhs;
    }
}

impl<E: FieldScalar> Mul<&Vector<E>> for LWEScalar<E> {
    type Output = Self;
    fn mul(mut self, rhs: &Vector<E>) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<E: FieldScalar> MulAssign<&Vector<E>> for LWEScalar<E> {
    fn mul_assign(&mut self, rhs: &Vector<E>) {
        self.commitment *= rhs;
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn should_be_additively_homomorphic() -> Result<()> {
        type Field = OxfoiScalar;
        let rng = &mut rand::rng();

        let lattice = LWEScalar::lattice_for(1, rng);

        let a = Field::sample_uniform(rng);
        let b = Field::sample_uniform(rng);
        let c = b + a;

        let comm_a = LWEScalar::commit(a.into(), lattice.clone(), rng);
        let comm_b = LWEScalar::commit(b.into(), lattice.clone(), rng);
        let comm_c = LWEScalar::commit(c.into(), lattice, rng);

        let e1 = comm_c.try_open(&c.into(), 1)?;

        let comm_c_homomorphic = comm_a + &comm_b;
        let e2 = comm_c_homomorphic.try_open(&c.into(), 2)?;

        let comm_zero = comm_c_homomorphic - &comm_c;
        // try to open to the zero value
        let e_out = comm_zero.try_open(&Field::zero().into(), 3)?;

        // check that the error vectors match after homomorphic operations
        assert_eq!((e2 - e1), e_out);

        Ok(())
    }
}
