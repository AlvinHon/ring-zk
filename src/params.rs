use std::{
    iter::Sum,
    ops::{Add, Mul, Sub},
};

use num::{integer::Roots, Integer, NumCast, One, Signed, Zero};
use poly_ring_xnp1::Polynomial;

use crate::{mat::Mat, polynomial::norm_2};

#[derive(Clone)]
pub struct Params<I> {
    /// Prime modulus. q = 2*d + 1 (mod 4d) in Lemma 1. Use d = 2 for this library.
    pub(crate) q: I,
    /// Norm bound for honest prover's randomness.
    pub(crate) b: I,

    // k > n >= l
    /// Height of the commitment matrix (a1).
    pub(crate) n: usize,
    /// Width of the commitment matrices
    pub(crate) k: usize,
    /// Dimension of the message space.
    pub(crate) l: usize,

    /// The maximum norm_1 of any element in Challenge Space C.
    pub(crate) kappa: usize,
}

impl<I> Params<I>
where
    I: Integer + Signed + Sum + Roots + Clone + NumCast,
{
    /// Prepare the value for the commitment. The input is a matrix (of size `l` x 1) of integer vectors.
    /// The generic parameter N indicates the maximum length of the integer vector. It must be a power
    /// of two.
    ///
    /// This method is for wrapping the input value into a polynomial which is used as primitive element
    /// in the library.
    ///
    /// ## Panics
    /// Panics if the `value.len()` is not equal to the message length (`l`),
    /// and the constant `N` is not a power of two.
    #[inline]
    pub fn prepare_value<const N: usize>(&self, value: Vec<Vec<I>>) -> Vec<Polynomial<I, N>>
    where
        I: Zero + One,
        for<'a> &'a I: Mul<Output = I> + Sub<Output = I> + Add<Output = I>,
    {
        assert!(value.len() == self.l);
        value
            .into_iter()
            .map(Polynomial::<I, N>::from_coeffs)
            .collect()
    }

    /// The standard deviation used in the zero-knowledge proof.
    pub(crate) fn standard_deviation(&self, deg_n: usize) -> usize {
        // sigma = 11 * kappa * b * sqrt(k*deg_n)
        self.b.to_usize().unwrap() * (11 * self.kappa) * (self.k * deg_n).sqrt()
    }

    /// Check the commitment constraint. norm_2(r_i) must be less or equal to 4*sigma*sqrt(N).
    pub(crate) fn check_commit_constraint<const N: usize>(&self, r: &Mat<I, N>) -> bool {
        let sigma = self.standard_deviation(N);
        let constraint = 4 * sigma * N.sqrt();
        r.polynomials.iter().all(|r_i| {
            r_i.iter()
                .all(|r_ij| norm_2(r_ij).to_usize().unwrap() <= constraint)
        })
    }

    /// Check the constraint for verification in zk protocol. norm_2(r_i) must be less or equal to 2*sigma*sqrt(N).
    pub(crate) fn check_verify_constraint<const N: usize>(&self, r: &Mat<I, N>) -> bool {
        let sigma = self.standard_deviation(N);
        let constraint = 2 * sigma * N.sqrt();
        r.polynomials.iter().all(|r_i| {
            r_i.iter()
                .all(|r_ij| norm_2(r_ij).to_usize().unwrap() <= constraint)
        })
    }
}

/// Parameters Set 1. The message length (l) is 1.
pub fn params_1() -> Params<i64> {
    let q = 3515337053_i64; // 32 bits
    assert!(q % 8 == 5, "q must be 2*d  + 1 mod 4d. Use d = 2.");
    Params {
        q,
        b: 1,
        n: 1,
        k: 3,
        l: 1,
        kappa: 36,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_deviation() {
        let params = params_1();
        let deg_n = 1024;
        let sigma = params.standard_deviation(deg_n);
        assert_eq!(sigma, 21780);
    }
}
