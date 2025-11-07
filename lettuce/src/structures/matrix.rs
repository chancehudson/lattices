use std::ops::Index;
use std::ops::IndexMut;

use crate::*;

#[derive(Debug, Clone, PartialEq)]
pub struct Matrix<E: RingElement> {
    width: usize,
    height: usize,
    entries: Vec<Vector<E>>,
}

impl<E: RingElement> Matrix<E> {
    pub fn zero(height: usize, width: usize) -> Self {
        Self {
            width,
            height,
            entries: vec![Vector::new(width); height],
        }
    }

    pub fn random<R: Rng>(height: usize, width: usize, rng: &mut R) -> Self {
        let mut entries = Vec::with_capacity(height);
        for _ in 0..height {
            entries.push(Vector::sample_uniform(width, rng));
        }
        Self {
            width,
            height,
            entries,
        }
    }

    /// Returns the (height, width) dimension of the matrix. Also known as (rows, columns).
    pub fn dimension(&self) -> (usize, usize) {
        (self.height, self.width)
    }

    pub fn height(&self) -> usize {
        self.dimension().0
    }

    pub fn width(&self) -> usize {
        self.dimension().1
    }

    /// Create a square identity matrix of size x size dimension.
    pub fn identity(size: usize) -> Self {
        let mut matrix = Self::zero(size, size);
        for i in 0..size {
            matrix[i][i] = E::one();
        }
        matrix
    }

    /// Take two matrices of equal height and create [self, other] by appending each row of `other`
    /// onto each row of `self`.
    ///
    /// Panics if matrices are not of equal height.
    pub fn compose_horizontal(mut self, other: Self) -> Self {
        let (self_height, self_width) = self.dimension();
        let (other_height, other_width) = other.dimension();
        assert_eq!(
            self_height, other_height,
            "Matrix::compose_horizontal cannot compose matrices of unequal height"
        );
        self.width = self_width + other_width;
        for (self_row, other_row) in self.entries.iter_mut().zip(other.entries.into_iter()) {
            self_row.append(other_row);
        }
        self
    }

    /// Take two matrices of equal width and create [self, other] by appending each column of `other`
    /// onto each column of `self`.
    ///
    /// Panics if matrices are not of equal width.
    pub fn compose_vertical(mut self, mut other: Self) -> Self {
        let (self_height, self_width) = self.dimension();
        let (other_height, other_width) = other.dimension();
        assert_eq!(
            self_width, other_width,
            "Matrix::compose_vertical cannot compose matrices of unequal width"
        );
        self.height = self_height + other_height;
        self.entries.append(&mut other.entries);
        self
    }

    /// Get an iterator over each row of the matrix `self`.
    pub fn iter(&self) -> impl Iterator<Item = &Vector<E>> {
        self.entries.iter()
    }
}

impl<E: RingElement> FromIterator<Vector<E>> for Matrix<E> {
    fn from_iter<T: IntoIterator<Item = Vector<E>>>(iter: T) -> Self {
        let mut width = None;
        let entries = iter
            .into_iter()
            .map(|row| {
                if width.is_none() {
                    width = Some(row.len());
                }
                assert_eq!(
                    width.unwrap(),
                    row.len(),
                    "row length mismatch in Matrix contruction from iterator"
                );
                row
            })
            .collect::<Vec<_>>();
        Self {
            width: width.unwrap_or(0),
            height: entries.len(),
            entries,
        }
    }
}

impl<E: RingElement> Index<usize> for Matrix<E> {
    type Output = Vector<E>;
    fn index(&self, index: usize) -> &Self::Output {
        &self.entries[index]
    }
}

impl<E: RingElement> IndexMut<usize> for Matrix<E> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.entries[index]
    }
}

impl<E: RingElement> AddAssign<&Self> for Matrix<E> {
    fn add_assign(&mut self, rhs: &Self) {
        assert_eq!(
            self.width, rhs.width,
            "cannot add matrices of different width"
        );
        assert_eq!(
            self.height, rhs.height,
            "cannot add matrices of different height"
        );
        for self_row in self.entries.iter_mut() {
            for other_row in rhs.entries.iter() {
                *self_row += other_row;
            }
        }
    }
}

impl<E: RingElement> Add<&Self> for Matrix<E> {
    type Output = Self;
    fn add(mut self, rhs: &Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<E: RingElement> MulAssign<&Self> for Matrix<E> {
    fn mul_assign(&mut self, rhs: &Self) {
        assert_eq!(
            self.width, rhs.width,
            "cannot mul matrices of different width"
        );
        assert_eq!(
            self.height, rhs.height,
            "cannot mul matrices of different height"
        );
        for self_row in self.entries.iter_mut() {
            for other_row in rhs.entries.iter() {
                *self_row *= other_row;
            }
        }
    }
}

impl<E: RingElement> Mul<&Self> for Matrix<E> {
    type Output = Self;
    fn mul(mut self, rhs: &Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<E: RingElement> Mul<&Vector<E>> for Matrix<E> {
    type Output = Vector<E>;
    fn mul(self, rhs: &Vector<E>) -> Self::Output {
        assert_eq!(self.width(), rhs.len());
        self.entries
            .into_iter()
            .map(|row| (row * rhs).into_sum())
            .collect::<Vector<_>>()
    }
}

impl<E: RingElement> Mul<&Vector<E>> for &Matrix<E> {
    type Output = Vector<E>;
    fn mul(self, rhs: &Vector<E>) -> Self::Output {
        assert_eq!(self.width(), rhs.len());
        self.entries
            .iter()
            .map(|row| {
                let mut sum = E::zero();
                for v in row.iter().zip(rhs.iter()) {
                    sum += *v.0 * *v.1;
                }
                sum
            })
            .collect::<Vector<_>>()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn matrix_dimension() {
        type Field = OxfoiScalar;

        let width = 50;
        let height = 100;
        let m = Matrix::<Field>::zero(height, width);
        assert_eq!(m.dimension(), (height, width));
        assert_eq!(m.height(), height);
        assert_eq!(m.width(), width);
    }

    #[test]
    fn matrix_identity() {
        type Field = OxfoiScalar;
        let mut rng = rand::rng();

        for s in 1..100 {
            let identity = Matrix::<Field>::identity(s);
            let vec = Vector::sample_uniform(s, &mut rng);
            assert_eq!(vec, &identity * &vec);
        }
    }

    #[test]
    fn matrix_compose_horizontal() {
        type Field = OxfoiScalar;
        let mut rng = rand::rng();

        let width1 = 100;
        let width2 = 50;
        let m1 = Matrix::<Field>::random(200, width1, &mut rng);
        let m2 = Matrix::<Field>::random(200, width2, &mut rng);

        let m_composed = m1.compose_horizontal(m2.clone());
        assert_eq!(m_composed.width(), width1 + width2);
        for row in m_composed.iter() {
            assert_eq!(row.len(), width1 + width2);
            // TODO: check contents of composed vectors
        }
    }

    #[test]
    fn matrix_compose_vertical() {
        type Field = OxfoiScalar;
        let mut rng = rand::rng();

        let height1 = 100;
        let height2 = 50;
        let m1 = Matrix::<Field>::random(height1, 200, &mut rng);
        let m2 = Matrix::<Field>::random(height2, 200, &mut rng);

        let m_composed = m1.clone().compose_vertical(m2.clone());
        assert_eq!(m_composed.height(), height1 + height2);
        assert_eq!(m_composed.entries.len(), height1 + height2);
        for (row_composed, row) in m_composed.iter().zip(m1.iter().chain(m2.iter())) {
            assert_eq!(row_composed, row);
        }
    }
}
