use crate::*;

/// Commitments based on the short integer solution problem over a scalar field. Comitted values
/// should be small/of low norm.
#[derive(Clone)]
pub struct SISScalar<E: FieldScalar> {
    lattice: Matrix<E>,
    pub commitment: Vector<E>,
}

impl<E: FieldScalar> SISScalar<E> {
    pub fn lattice_for<R: Rng>(element_len: usize, rng: &mut R) -> Matrix<E> {
        // m value
        let height: usize = element_len * E::BIT_WIDTH;
        Matrix::<E>::random(height, element_len, rng)
    }

    pub fn commit(val: Vector<E>, lattice: Matrix<E>) -> Self {
        // TODO: warn on big vals
        Self {
            commitment: &lattice * &val,
            lattice,
        }
    }

    pub fn try_open(&self, val: &Vector<E>, max_dist: u128) -> Result<()> {
        for v in val.iter() {
            let disp = v.displacement();
            if disp.unsigned_abs() > max_dist {
                anyhow::bail!(
                    "Error opening SIS commitment, value contains element {} beyond bound {}",
                    disp,
                    max_dist
                );
            }
        }
        let expected_commitment = &self.lattice * val;
        if expected_commitment != self.commitment {
            anyhow::bail!("Error opening SIS commitment, commitment mismatch");
        }
        Ok(())
    }
}

impl<E: FieldScalar> Add<&Self> for SISScalar<E> {
    type Output = Self;
    fn add(mut self, rhs: &Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<E: FieldScalar> AddAssign<&Self> for SISScalar<E> {
    fn add_assign(&mut self, rhs: &Self) {
        self.commitment += &rhs.commitment
    }
}

impl<E: FieldScalar> Mul<E> for SISScalar<E> {
    type Output = Self;
    fn mul(mut self, rhs: E) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<E: FieldScalar> MulAssign<E> for SISScalar<E> {
    fn mul_assign(&mut self, rhs: E) {
        self.commitment *= rhs;
    }
}

impl<E: FieldScalar> Mul<&Vector<E>> for SISScalar<E> {
    type Output = Self;
    fn mul(mut self, rhs: &Vector<E>) -> Self::Output {
        self.commitment *= rhs;
        self
    }
}

impl<E: FieldScalar> MulAssign<&Vector<E>> for SISScalar<E> {
    fn mul_assign(&mut self, rhs: &Vector<E>) {
        self.commitment *= rhs;
    }
}

#[cfg(test)]
mod test {
    use crate::*;

    #[test]
    fn should_be_additively_homomorphic() -> Result<()> {
        type Field = OxfoiScalar;
        let rng = &mut rand::rng();
        const PART_BITS: usize = 8;
        // allow a single addition
        const ARITH_MAX: u128 = 1 << (PART_BITS + 1);

        let a = Field::sample_uniform(rng).as_le_bits_vec(PART_BITS);
        let b = Field::sample_uniform(rng).as_le_bits_vec(PART_BITS);
        let c = a.clone() + &b;

        let lattice = SISScalar::lattice_for(a.len(), rng);

        let comm_a = SISScalar::commit(a, lattice.clone());
        let comm_b = SISScalar::commit(b, lattice.clone());
        let comm_c = SISScalar::commit(c.clone(), lattice.clone());

        let comm_c_homomorphic = comm_a + &comm_b;
        assert_eq!(comm_c.commitment, comm_c_homomorphic.commitment);

        comm_c_homomorphic.try_open(&c, ARITH_MAX)?;
        comm_c.try_open(&c, ARITH_MAX)?;

        Ok(())
    }
}
