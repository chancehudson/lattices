use std::ops::Index;
use std::ops::IndexMut;

use crate::*;

// TODO: vectors backed by abstract iterators
// can transfer between threads to evaluate as needed
/// A vector of ring elements. Supports arithmetic with other vectors, matrices, and ring elements.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Vector<E: RingElement> {
    entries: Vec<E>,
}

impl<E: RingElement> Vector<E> {
    pub fn new(len: usize) -> Self {
        Self {
            entries: vec![E::default(); len],
        }
    }

    pub fn resize(&mut self, new_len: usize, val: E) {
        assert!(new_len >= self.entries.len());
        self.entries.resize(new_len, val);
    }

    pub fn chunks(&self, len: usize) -> impl Iterator<Item = Vector<E>> {
        self.entries.chunks(len).map(|v| Self::from(v.to_vec()))
    }

    pub fn as_slice(&self) -> &[E] {
        &self.entries
    }

    pub fn get(&self, i: usize) -> Option<&E> {
        self.entries.get(i)
    }

    pub fn sum(&self) -> E {
        let mut out = E::default();
        for v in &self.entries {
            out += *v;
        }
        out
    }

    pub fn into_sum(mut self) -> E {
        let mut out = E::default();
        for v in std::mem::take(&mut self.entries) {
            out += v;
        }
        out
    }

    pub fn sample_uniform<R: Rng>(len: usize, rng: &mut R) -> Self {
        let mut entries = Vec::with_capacity(len);
        for _ in 0..len {
            entries.push(E::sample_uniform(rng));
        }
        Self { entries }
    }

    pub fn is_zero(&self) -> bool {
        for entry in &self.entries {
            if *entry != E::default() {
                return false;
            }
        }
        true
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn iter(&self) -> impl Iterator<Item = &E> + ExactSizeIterator + Clone {
        self.entries.iter()
    }

    /// Take the entries from `other` and append them to the end of `self`.
    pub fn append(&mut self, mut other: Self) {
        self.entries.append(&mut other.entries);
    }
}

impl<E: RingElement> Display for Vector<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = self
            .entries
            .iter()
            .map(|v| format!("{}", v))
            .collect::<Vec<String>>()
            .join(",");
        f.write_str(&format!("{}", str))?;
        Ok(())
    }
}

impl<E: RingElement> From<Vec<E>> for Vector<E> {
    fn from(value: Vec<E>) -> Self {
        Self { entries: value }
    }
}

impl<E: RingElement> From<E> for Vector<E> {
    fn from(value: E) -> Self {
        Self {
            entries: vec![value],
        }
    }
}

impl<E: RingElement> IntoIterator for Vector<E> {
    type Item = E;
    type IntoIter = <Vec<E> as IntoIterator>::IntoIter;
    fn into_iter(self) -> Self::IntoIter {
        self.entries.into_iter()
    }
}

impl<'a, E: RingElement> IntoIterator for &'a Vector<E> {
    type Item = &'a E;
    type IntoIter = std::slice::Iter<'a, E>;
    fn into_iter(self) -> Self::IntoIter {
        self.entries.iter()
    }
}

impl<'a, E: RingElement> FromIterator<E> for Vector<E> {
    fn from_iter<T: IntoIterator<Item = E>>(iter: T) -> Self {
        Self {
            entries: iter.into_iter().collect(),
        }
    }
}

impl<E: RingElement> Index<usize> for Vector<E> {
    type Output = E;
    fn index(&self, index: usize) -> &Self::Output {
        assert!(index < self.len(), "requested index outside of vector");
        self.entries.get(index).unwrap()
    }
}

impl<E: RingElement> IndexMut<usize> for Vector<E> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        assert!(index < self.len(), "requested index outside of vector");
        self.entries.get_mut(index).unwrap()
    }
}

impl<E: RingElement> Mul<E> for Vector<E> {
    type Output = Self;
    fn mul(mut self, rhs: E) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<E: RingElement> MulAssign<E> for Vector<E> {
    fn mul_assign(&mut self, rhs: E) {
        for entry in self.entries.iter_mut() {
            *entry *= rhs;
        }
    }
}

impl<E: RingElement> Mul<&Vector<E>> for Vector<E> {
    type Output = Self;
    fn mul(mut self, rhs: &Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<E: RingElement> MulAssign<&Vector<E>> for Vector<E> {
    fn mul_assign(&mut self, rhs: &Self) {
        assert_eq!(self.len(), rhs.len());
        for (i, entry) in self.entries.iter_mut().enumerate() {
            *entry *= rhs[i];
        }
    }
}

impl<E: RingElement> AddAssign<&Self> for Vector<E> {
    fn add_assign(&mut self, rhs: &Self) {
        assert_eq!(self.len(), rhs.len());
        for (i, entry) in self.entries.iter_mut().enumerate() {
            *entry += rhs[i];
        }
    }
}

impl<E: RingElement> Add<&Self> for Vector<E> {
    type Output = Self;
    fn add(mut self, rhs: &Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<E: RingElement> Add<E> for Vector<E> {
    type Output = Self;
    fn add(mut self, rhs: E) -> Self::Output {
        self += rhs;
        self
    }
}

impl<E: RingElement> AddAssign<E> for Vector<E> {
    fn add_assign(&mut self, rhs: E) {
        for entry in self.entries.iter_mut() {
            *entry += rhs;
        }
    }
}

impl<E: RingElement> SubAssign for Vector<E> {
    fn sub_assign(&mut self, rhs: Self) {
        *self -= &rhs;
    }
}

impl<E: RingElement> SubAssign<&Self> for Vector<E> {
    fn sub_assign(&mut self, rhs: &Self) {
        assert_eq!(self.len(), rhs.len());
        for (i, entry) in self.entries.iter_mut().enumerate() {
            *entry -= rhs[i];
        }
    }
}

impl<E: RingElement> Sub for Vector<E> {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<E: RingElement> Sub<&Vector<E>> for Vector<E> {
    type Output = Vector<E>;
    fn sub(self, rhs: &Vector<E>) -> Self::Output {
        self.into_iter()
            .zip(rhs.iter())
            .map(|e| e.0 - *e.1)
            .collect::<Vector<_>>()
    }
}

impl<E: RingElement> Sub<Vector<E>> for &Vector<E> {
    type Output = Vector<E>;
    fn sub(self, rhs: Vector<E>) -> Self::Output {
        self.iter()
            .zip(rhs.into_iter())
            .map(|e| *e.0 - e.1)
            .collect::<Vector<_>>()
    }
}
