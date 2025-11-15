use crate::*;

/// A multivariate polynomial of bounded degree N over a scalar field E
/// with V independent variables.
///
/// For convenience we'll denote the first 3 independent variables as
/// x, y, z and make them accessible as methods.
///
/// A multivariate polynomial of structure a*x * b*y = c*z + e
#[derive(Default)]
pub struct MultiPolynomial<E: FieldScalar> {
    pub a: E,
    x: Option<E>,
    pub b: E,
    y: Option<E>,
    pub c: E,
    z: Option<E>,
    pub e: E,
}

impl<E: FieldScalar> MultiPolynomial<E> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn zero() -> Self {
        Self::default()
    }

    /// Get the constant term
    pub fn e(&self) -> E {
        self.e
    }

    /// Get the constant term
    pub fn e_mut(&mut self) -> &mut E {
        &mut self.e
    }

    /// Get the constant x coefficient
    pub fn x_coef_mut(&mut self) -> &mut E {
        &mut self.a
    }

    /// Get the constant y coefficient
    pub fn y_coef_mut(&mut self) -> &mut E {
        &mut self.b
    }

    /// Get the constant z coefficient
    pub fn z_coef_mut(&mut self) -> &mut E {
        &mut self.c
    }

    /// Get the linear x coefficient
    pub fn x_mut(&mut self) -> &mut Option<E> {
        &mut self.x
    }

    /// Get the linear y coefficient
    pub fn y_mut(&mut self) -> &mut Option<E> {
        &mut self.y
    }

    /// Get the linear z coefficient
    pub fn z_mut(&mut self) -> &mut Option<E> {
        &mut self.z
    }

    /// Get the linear x coefficient
    pub fn x(&self) -> Option<E> {
        self.x
    }

    /// Get the linear y coefficient
    pub fn y(&self) -> Option<E> {
        self.y
    }

    /// Get the linear z coefficient
    pub fn z(&self) -> Option<E> {
        self.z
    }

    /// Solve the equation returning the value of each variable.
    pub fn solve(&self) -> Result<[E; 3]> {
        let mut unknown_count = 0;
        let inputs = [self.x(), self.y(), self.z()];
        for v in inputs {
            if v.is_none() {
                unknown_count += 1;
            }
        }
        if unknown_count == 0 {
            return Ok([self.x().unwrap(), self.y.unwrap(), self.z.unwrap()]);
        }
        if unknown_count > 1 {
            anyhow::bail!("MultiPolynomial can't solve equation with more than 1 unknown");
        }

        match (self.x(), self.y(), self.z()) {
            (None, Some(y), Some(z)) => Ok([
                self.a.inverse() * ((self.c * z + self.e) - self.b * y),
                y,
                z,
            ]),
            (Some(x), None, Some(z)) => Ok([
                x,
                self.b.inverse() * ((self.c * z + self.e) - self.a * x),
                z,
            ]),
            (Some(x), Some(y), None) => {
                Ok([x, y, self.c.inverse() * (self.a * x + self.b * y - self.e)])
            }
            _ => {
                unreachable!()
            }
        }
    }
}

#[test]
fn should_solve_multipoly() -> Result<()> {
    type E = MilliScalarMont;
    for _ in 0..100 {
        let mut poly = MultiPolynomial::<E>::new();
        *poly.x_coef_mut() = rand::random::<u32>().into();
        *poly.y_coef_mut() = rand::random::<u32>().into();
        *poly.z_coef_mut() = rand::random::<u32>().into();
        *poly.e_mut() = rand::random::<u32>().into();
        *poly.x_mut() = Some(rand::random::<u32>().into());
        *poly.y_mut() = Some(rand::random::<u32>().into());
        *poly.z_mut() = Some(rand::random::<u32>().into());

        let empty_index = rand::random_range(0..3);
        let expected = match empty_index {
            0 => {
                *poly.x_mut() = None;
                poly.a.inverse()
                    * ((poly.e() + poly.c * poly.z.unwrap()) - poly.b * poly.y.unwrap())
            }
            1 => {
                *poly.y_mut() = None;
                poly.b.inverse()
                    * ((poly.e() + poly.c * poly.z.unwrap()) - poly.a * poly.x.unwrap())
            }
            2 => {
                *poly.z_mut() = None;
                poly.c.inverse()
                    * ((poly.a * poly.x.unwrap()) + poly.b * poly.y.unwrap() - poly.e())
            }
            _ => unreachable!(),
        };
        let solution = poly.solve()?;
        assert_eq!(solution[empty_index], expected);
    }
    Ok(())
}
