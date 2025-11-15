use crate::*;

/// Instance of a rank 1 constraint system over a finite field or ring.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct R1CS<E: RingElement> {
    pub a: Matrix<E>,
    pub b: Matrix<E>,
    pub c: Matrix<E>,
}

impl<E: RingElement> R1CS<E> {
    /// Sample a random R1CS instance and a witness that fulfills the instance.
    ///
    /// This instance will be non-trivial and non-sparse. That is, finding a second
    /// witness vector should be difficult (though not impossible). This function is
    /// designed to create test R1CS instances.
    pub fn sample_uniform<R: Rng>(height: usize, width: usize, rng: &mut R) -> (Self, Vector<E>) {
        let wtns = Vector::<E>::sample_uniform(width, rng);
        let a = Matrix::random(height, width, rng);
        let b = Matrix::random(height, width, rng);

        let aw = &a * &wtns;
        let bw = &b * &wtns;
        let abw = aw.clone() * &bw;
        // now we need a c matrix that yields this abw value
        // we're looking for abw[i] = <c[i], w>
        let c = b
            .iter()
            .enumerate()
            .map(|(i, row)| row.clone() * aw[i])
            .collect::<Matrix<_>>();
        debug_assert!(abw == &c * &wtns);
        (Self { a, b, c }, wtns)
    }

    /// Generate an R1CS instance that is fulfilled by all witness vectors.
    pub fn identity(height: usize, width: usize) -> Self {
        let v = Matrix::zero(height, width);
        Self {
            a: v.clone(),
            b: v.clone(),
            c: v.clone(),
        }
    }

    /// Evaluate the R1CS instance for an input witness. Returns the evaluation of
    /// AwBw - Cw, which may be a non-zero vector.
    ///
    /// For R1CS instances the output should be checked against the 0 value. For relaxed R1CS
    /// instances the output should be checked against the e vector.
    pub fn eval(&self, witness: &Vector<E>) -> Result<Vector<E>> {
        self.assert_consistency()?;

        let ab = (self.a.clone() * witness) * &(self.b.clone() * witness);
        let c = self.c.clone() * witness;

        Ok(ab - c)
    }

    /// Number of constraints in the system.
    pub fn height(&self) -> usize {
        self.a.height()
    }

    /// Number of elements in the witness vector.
    pub fn width(&self) -> usize {
        self.a.width()
    }

    /// Dimension of the constraint matrices.
    pub fn dimension(&self) -> (usize, usize) {
        self.a.dimension()
    }

    /// Assert consistency of matrix dimensions.
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

#[test]
fn random_r1cs() {
    type Field = MilliScalar;
    let rng = &mut rand::rng();
    let (r1cs, wtns) = R1CS::<Field>::sample_uniform(10, 20, rng);
    for v in r1cs.eval(&wtns).iter() {
        assert!(v.is_zero());
    }
}
