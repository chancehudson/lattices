use crate::*;

pub struct R1CS<E: RingElement> {
    a: Matrix<E>,
    b: Matrix<E>,
    c: Matrix<E>,
}

impl<E: RingElement> R1CS<E> {
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
