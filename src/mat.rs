use std::ops::{Add, Mul, Sub};

use num::{One, Zero};
use poly_ring_xnp1::Polynomial;

/// A matrix over polynomial rings Z\[x]/(x^n+1).
pub(crate) struct Mat<T, const N: usize> {
    polynomials: Vec<Vec<Polynomial<T, N>>>,
}

impl<T, const N: usize> Mat<T, N> {
    pub fn dim(&self) -> (usize, usize) {
        let m = self.polynomials.len();
        let n = if self.polynomials.is_empty() {
            0
        } else {
            self.polynomials[0].len()
        };
        (m, n)
    }

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
    pub fn extend_rows(&mut self, other: &Mat<T, N>)
    where
        T: Clone,
    {
        let (_, n) = self.dim();
        let (_, n2) = other.dim();
        assert_eq!(n, n2);

        self.polynomials.extend(other.polynomials.clone());
    }

    /// Extend the matrix by adding columns.
    /// Original dimensions: m x n;
    /// New dimensions: m x (n + n')
    pub fn extend_cols(&mut self, other: &Mat<T, N>)
    where
        T: Clone,
    {
        let (m, _) = self.dim();
        let (m2, _) = other.dim();
        assert_eq!(m, m2);

        for i in 0..m {
            self.polynomials[i].extend(other.polynomials[i].clone());
        }
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

        a.extend_rows(&b);

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

        a.extend_cols(&b);

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
