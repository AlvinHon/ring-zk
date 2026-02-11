//! Defines the public parameters for the protocol.

use std::ops::{Add, Mul, Sub};

use num::{integer::Roots, BigUint, FromPrimitive, One, ToPrimitive, Zero};
use poly_ring_xnp1::{zq::ZqI64, Polynomial};
use rand::{distr::uniform::SampleUniform, RngExt};
use serde::{Deserialize, Serialize};

use crate::{mat::Mat, polynomial::norm_2, CommitmentKey};

/// Public parameters for the protocol.
///
/// ## Safety
/// The struct implements Default for instantiation. If you want to use a custom parameter setting,
/// please carefully check the constraints for the parameters (see the comment-doc for each parameters).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Params<I> {
    /// Prime modulus q' divided by 2, where q' = 2*d + 1 (mod 4d). Use d = 2 in this library.
    pub q: I, // The formula is defined in Lemma 1 of the paper.
    /// Norm bound for honest prover's randomness. It is the `beta` value commonly used
    /// in lattice-based cryptography, that indicates how "random" the commitment is. It
    /// is a trade-off between security and efficiency. The value should be a small constant.
    pub b: I,

    /// Height of the commitment matrix (a1). It should be a value s.t. k > n >= l.
    pub n: usize,
    /// Width of the commitment matrices. It should be a value s.t. k > n >= l.
    pub k: usize,
    /// Dimension of the message space. It should be a value s.t. k > n >= l.
    pub l: usize,

    /// The maximum norm_1 of any element in Challenge Space C. It indicates the
    /// "size" of the challenge space. The value should be a small constant.
    pub kappa: usize,
}

impl<I> Params<I>
where
    I: Clone + PartialOrd + Ord + One + Zero + FromPrimitive + ToPrimitive + SampleUniform,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    /// Generate a new commitment key. The generic parameter N indicates the maximum length of the integer vector.
    /// It must be a power of two.
    ///
    /// ## Panics
    /// Panics if the constant `N` is not a power of two.
    #[inline]
    pub fn generate_commitment_key<const N: usize>(
        &self,
        rng: &mut impl RngExt,
    ) -> CommitmentKey<I, N> {
        CommitmentKey::new(rng, self)
    }

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
    pub fn prepare_value<const N: usize>(
        &self,
        value: Vec<Vec<impl Into<I>>>,
    ) -> Vec<Polynomial<I, N>> {
        assert!(value.len() == self.l);
        value
            .into_iter()
            .map(|v| v.into_iter().map(Into::into).collect())
            .map(Polynomial::<I, N>::from_coeffs)
            .collect()
    }

    /// Prepare the scalar for the commitment. The input is a vector of integers.
    /// The generic parameter N indicates the maximum length of the integer vector. It must be a power
    /// of two.
    ///
    /// This method is for wrapping the input scalar into a polynomial which is used as primitive element
    /// in the library.
    ///
    /// ## Panics
    /// Panics if the constant `N` is not a power of two.
    #[inline]
    pub fn prepare_scalar<const N: usize>(&self, scalar: Vec<impl Into<I>>) -> Polynomial<I, N> {
        Polynomial::from_coeffs(scalar.into_iter().map(Into::into).collect::<Vec<I>>())
    }

    /// The standard deviation used in the zero-knowledge proof.
    pub(crate) fn standard_deviation(&self, deg_n: usize) -> usize {
        // The formula defined in Table 1 of the paper:
        // sigma = 11 * kappa * b * sqrt(k*deg_n)
        self.b.to_usize().unwrap() * (11 * self.kappa) * (self.k * deg_n).sqrt()
    }

    /// Check the commitment constraint. norm_2(r_i) must be less or equal to 4*sigma*sqrt(N).
    /// It is used in the commitment scheme.
    pub(crate) fn check_commit_constraint<const N: usize>(&self, r: &Mat<I, N>) -> bool {
        let sigma = self.standard_deviation(N);
        let constraint = BigUint::from(4 * sigma * N.sqrt());
        r.polynomials
            .iter()
            .all(|r_i| r_i.iter().all(|r_ij| norm_2(r_ij) <= constraint))
    }

    /// Check the constraint for verification in zk protocol. norm_2(r_i) must be less or equal to 2*sigma*sqrt(N).
    /// It is used in the verification step in the zk protocol.
    pub(crate) fn check_verify_constraint<const N: usize>(&self, r: &Mat<I, N>) -> bool {
        let sigma = self.standard_deviation(N);
        let constraint = BigUint::from(2 * sigma * N.sqrt());
        r.polynomials
            .iter()
            .all(|r_i| r_i.iter().all(|r_ij| norm_2(r_ij) <= constraint))
    }
}

impl Default for Params<ZqI64<3515337053_i64>> {
    /// This default parameter setting accepts a message of length 1, and
    /// the integer range in the message (32 bits) is [-3515337053/2, 3515337053/2].
    fn default() -> Self {
        // values (except q) are taken from the paper Table 2 (approximately 32 bits).
        let q = ZqI64::from(3515337053_i64 / 2); // divide by 2 for shifting the range to [-q/2, q/2]
        let b = ZqI64::one();

        Params {
            q,
            b,
            n: 1,
            k: 3,
            l: 1,
            kappa: 36,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_deviation() {
        let params = Params::default();
        let deg_n = 1024;
        let sigma = params.standard_deviation(deg_n);
        assert_eq!(sigma, 21780);
    }

    #[test]
    fn test_prepare_scalar() {
        let params = Params::default();
        let scalar = vec![1, 2, 3, 4];
        let p = params.prepare_scalar::<4>(scalar);
        assert_eq!(p.deg(), 3);
    }

    #[test]
    fn test_prepare_value() {
        let params = Params::default();
        let value = vec![vec![1, 2, 3, 4]];
        let p = params.prepare_value::<4>(value);
        assert_eq!(p.len(), 1);
        assert_eq!(p[0].deg(), 3);
    }
}
