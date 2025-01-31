use std::ops::{Add, Mul, Sub};

use num::{Integer, One, Zero};
use poly_ring_xnp1::Polynomial;
use rand::{distributions::uniform::SampleUniform, Rng};

use crate::polynomial::rand_polynomial_within;

/// A matrix over polynomial rings Z\[x]/(x^n+1).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Mat<T, const N: usize> {
    pub(crate) polynomials: Vec<Vec<Polynomial<T, N>>>,
}

impl<T, const N: usize> Mat<T, N> {
    /// Create a matrix (m x n) from a single polynomial.
    pub fn from_element(m: usize, n: usize, element: Polynomial<T, N>) -> Self
    where
        T: Clone,
    {
        let polynomials = vec![vec![element.clone(); n]; m];
        Mat { polynomials }
    }

    /// Create a matrix (m x 1) from a vector of polynomials.
    pub fn from_vec(polynomials: Vec<Polynomial<T, N>>) -> Self {
        Mat {
            polynomials: polynomials.into_iter().map(|p| vec![p]).collect(),
        }
    }

    /// Create a random matrix (m x n) with polynomials with coefficients in the
    /// range `[-bound, bound]`.
    pub fn rand(rng: &mut impl Rng, m: usize, n: usize, bound: T) -> Self
    where
        T: Integer + SampleUniform + Clone,
    {
        let polynomials = (0..m)
            .map(|_| {
                (0..n)
                    .map(|_| rand_polynomial_within(rng, bound.clone()))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        Mat { polynomials }
    }

    pub fn dim(&self) -> (usize, usize) {
        let m = self.polynomials.len();
        let n = if self.polynomials.is_empty() {
            0
        } else {
            self.polynomials[0].len()
        };
        (m, n)
    }

    #[allow(clippy::needless_range_loop)]
    pub fn dot(&self, other: &Mat<T, N>) -> Mat<T, N>
    where
        T: Zero + One + Clone,
        for<'a> &'a T: Add<Output = T> + Mul<Output = T> + Sub<Output = T>,
    {
        // mxn * nxp = mxp
        let (m, n) = self.dim();
        let (n2, p) = other.dim();
        assert_eq!(n, n2);

        let mut polynomials = vec![vec![Polynomial::<T, N>::zero(); p]; m];
        for i in 0..m {
            for j in 0..p {
                for k in 0..n {
                    polynomials[i][j] = polynomials[i][j].clone()
                        + self.polynomials[i][k].clone() * other.polynomials[k][j].clone();
                }
            }
        }
        Mat { polynomials }
    }

    #[allow(clippy::needless_range_loop)]
    pub fn add(&self, other: &Mat<T, N>) -> Mat<T, N>
    where
        T: Zero + One + Clone,
        for<'a> &'a T: Add<Output = T> + Mul<Output = T> + Sub<Output = T>,
    {
        let (m, n) = self.dim();
        let (m2, n2) = other.dim();
        assert_eq!(m, m2);
        assert_eq!(n, n2);

        let mut polynomials = vec![vec![Polynomial::<T, N>::zero(); n]; m];
        for i in 0..m {
            for j in 0..n {
                polynomials[i][j] =
                    self.polynomials[i][j].clone() + other.polynomials[i][j].clone();
            }
        }
        Mat { polynomials }
    }

    /// Extend the matrix by adding rows.
    /// Original dimensions: m x n;
    /// New dimensions: (m + m') x n
    pub fn extend_rows(&mut self, other: Mat<T, N>)
    where
        T: Clone,
    {
        let (_, n) = self.dim();
        let (_, n2) = other.dim();
        assert_eq!(n, n2);

        self.polynomials.extend(other.polynomials);
    }

    /// Extend the matrix by adding columns.
    /// Original dimensions: m x n;
    /// New dimensions: m x (n + n')
    pub fn extend_cols(&mut self, other: Mat<T, N>)
    where
        T: Clone,
    {
        let (m, _) = self.dim();
        let (m2, _) = other.dim();
        assert_eq!(m, m2);

        self.polynomials
            .iter_mut()
            .zip(other.polynomials)
            .for_each(|(a, b)| a.extend(b));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const N: usize = 4;

    #[test]
    fn test_dot() {
        let a_0_0 = Polynomial::<i32, N>::new(vec![1, 2, 3]);
        let a_0_1 = Polynomial::<i32, N>::new(vec![4, 5, 6]);

        // 1x2 matrix
        let a = Mat {
            polynomials: vec![vec![a_0_0.clone(), a_0_1.clone()]],
        };

        let b_0_0 = Polynomial::<i32, N>::new(vec![1, 2]);
        let b_1_0 = Polynomial::<i32, N>::new(vec![3, 4]);

        // 2x1 matrix
        let b = Mat {
            polynomials: vec![vec![b_0_0.clone()], vec![b_1_0.clone()]],
        };

        let c = a.dot(&b);

        assert_eq!(
            c.polynomials,
            vec![vec![
                a_0_0.clone() * b_0_0.clone() + a_0_1.clone() * b_1_0.clone()
            ]]
        );
    }

    #[test]
    fn test_add() {
        let a_0_0 = Polynomial::<i32, N>::new(vec![1, 2, 3]);
        let a_0_1 = Polynomial::<i32, N>::new(vec![4, 5, 6]);

        // 1x2 matrix
        let a = Mat {
            polynomials: vec![vec![a_0_0.clone(), a_0_1.clone()]],
        };

        let b_0_0 = Polynomial::<i32, N>::new(vec![1, 2, 3]);
        let b_0_1 = Polynomial::<i32, N>::new(vec![4, 5, 6]);

        // 1x2 matrix
        let b = Mat {
            polynomials: vec![vec![b_0_0.clone(), b_0_1.clone()]],
        };

        let c = a.add(&b);

        assert_eq!(
            c.polynomials,
            vec![vec![
                a_0_0.clone() + b_0_0.clone(),
                a_0_1.clone() + b_0_1.clone()
            ]]
        );
    }

    #[test]
    fn test_extend_rows() {
        let a_0_0 = Polynomial::<i32, N>::new(vec![1, 2, 3]);
        let a_0_1 = Polynomial::<i32, N>::new(vec![4, 5, 6]);

        // 1x2 matrix
        let mut a = Mat {
            polynomials: vec![vec![a_0_0.clone(), a_0_1.clone()]],
        };

        let b_0_0 = Polynomial::<i32, N>::new(vec![1, 2, 3]);
        let b_0_1 = Polynomial::<i32, N>::new(vec![4, 5, 6]);

        // 1x2 matrix
        let b = Mat {
            polynomials: vec![vec![b_0_0.clone(), b_0_1.clone()]],
        };

        a.extend_rows(b);

        assert_eq!(
            a.polynomials,
            vec![
                vec![a_0_0.clone(), a_0_1.clone()],
                vec![b_0_0.clone(), b_0_1.clone()]
            ]
        );
    }

    #[test]
    fn test_extend_cols() {
        let a_0_0 = Polynomial::<i32, N>::new(vec![1, 2, 3]);
        let a_0_1 = Polynomial::<i32, N>::new(vec![4, 5, 6]);

        // 1x2 matrix
        let mut a = Mat {
            polynomials: vec![vec![a_0_0.clone(), a_0_1.clone()]],
        };

        let b_0_0 = Polynomial::<i32, N>::new(vec![1, 2, 3]);
        let b_0_1 = Polynomial::<i32, N>::new(vec![4, 5, 6]);

        // 1x2 matrix
        let b = Mat {
            polynomials: vec![vec![b_0_0.clone(), b_0_1.clone()]],
        };

        a.extend_cols(b);

        assert_eq!(
            a.polynomials,
            vec![vec![
                a_0_0.clone(),
                a_0_1.clone(),
                b_0_0.clone(),
                b_0_1.clone()
            ]]
        );
    }
}
