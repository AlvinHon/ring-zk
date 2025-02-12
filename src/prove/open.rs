use std::{
    iter::Sum,
    ops::{Add, Mul, Sub},
};

use num::{integer::Roots, Integer, NumCast, Signed};
use poly_ring_xnp1::Polynomial;
use rand::Rng;
use rand_distr::uniform::SampleUniform;

use crate::{
    challenge_space::random_polynomial_from_challenge_set,
    commit::{Commitment, CommitmentKey, Opening},
    mat::Mat,
    params::Params,
    polynomial::random_polynomial_in_normal_distribution,
};

pub struct OpenProofProver<I, const N: usize>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    pub(crate) params: Params<I>,
    pub(crate) ck: CommitmentKey<I, N>,
}

impl<I, const N: usize> OpenProofProver<I, N>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    pub fn new(ck: CommitmentKey<I, N>, params: Params<I>) -> Self {
        Self { params, ck }
    }
}

impl<I, const N: usize> OpenProofProver<I, N>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    pub fn commit(
        &self,
        rng: &mut impl Rng,
        x: Vec<Polynomial<I, N>>,
    ) -> (OpenProofResponseContext<I, N>, OpenProofCommitment<I, N>) {
        let (opening, c) = self.ck.commit(rng, x, &self.params);

        // y <- N^k_sigma
        let y = Mat::<I, N>::new_with(self.params.k, 1, || {
            random_polynomial_in_normal_distribution::<I, N>(
                rng,
                I::zero().to_f64().unwrap(),
                self.params.standard_deviation(N) as f64,
            )
        });

        // t = A1 * y
        let t = self.ck.a1.dot(&y).one_d_mat_to_vec();

        (
            OpenProofResponseContext { opening, y },
            OpenProofCommitment { c, t },
        )
    }

    pub fn create_response(
        &self,
        context: OpenProofResponseContext<I, N>,
        challenge: OpenProofChallenge<I, N>,
    ) -> OpenProofResponse<I, N> {
        // z = y + d * r
        let z = context
            .y
            .add(&context.opening.r.componentwise_mul(&challenge.d));
        OpenProofResponse { z }
    }
}

pub struct OpenProofVerifier<I, const N: usize>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    params: Params<I>,
    ck: CommitmentKey<I, N>,
}

impl<I, const N: usize> OpenProofVerifier<I, N>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    pub fn new(ck: CommitmentKey<I, N>, params: Params<I>) -> Self {
        OpenProofVerifier { params, ck }
    }
    pub fn generate_challenge(
        &self,
        rng: &mut impl Rng,
        commitment: OpenProofCommitment<I, N>,
    ) -> (OpenProofVerificationContext<I, N>, OpenProofChallenge<I, N>) {
        let d = random_polynomial_from_challenge_set(rng, self.params.kappa);
        let (c1, _) = commitment.c.c1_c2(&self.params);
        (
            OpenProofVerificationContext {
                c1,
                t: commitment.t,
                d: d.clone(),
            },
            OpenProofChallenge { d },
        )
    }

    pub fn verify(
        &self,
        response: OpenProofResponse<I, N>,
        context: OpenProofVerificationContext<I, N>,
    ) -> bool {
        if !self.params.check_verify_constraint(&response.z) {
            return false;
        }
        // A1 * z = t + c1 * d
        let lhs = self.ck.a1.dot(&response.z);
        let rhs = Mat::<I, N>::from_vec(context.t).add(&context.c1.componentwise_mul(&context.d));
        lhs == rhs
    }
}

pub struct OpenProofResponseContext<I, const N: usize> {
    pub opening: Opening<I, N>,
    y: Mat<I, N>, // k x 1 matrix
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpenProofCommitment<I, const N: usize> {
    pub c: Commitment<I, N>,
    t: Vec<Polynomial<I, N>>, // n x 1 matrix
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpenProofChallenge<I, const N: usize> {
    d: Polynomial<I, N>,
}

pub struct OpenProofVerificationContext<I, const N: usize> {
    c1: Mat<I, N>,            // n x 1 matrix
    t: Vec<Polynomial<I, N>>, // n x 1 matrix
    d: Polynomial<I, N>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpenProofResponse<I, const N: usize> {
    z: Mat<I, N>, // k x 1 matrix
}
