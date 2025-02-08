use std::{
    iter::Sum,
    ops::{Add, Mul, Neg, Sub},
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

pub struct LinearProofProver<I, const N: usize>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I> + Neg<Output = I>,
{
    pub(crate) params: Params<I>,
    pub(crate) ck: CommitmentKey<I, N>,
}

impl<I, const N: usize> LinearProofProver<I, N>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I> + Neg<Output = I>,
{
    pub fn new(ck: CommitmentKey<I, N>, params: Params<I>) -> Self {
        Self { params, ck }
    }

    pub fn commit_and_prove(
        &self,
        rng: &mut impl Rng,
        g: Polynomial<I, N>,
        x: Vec<Polynomial<I, N>>,
    ) -> (
        LinearProofResponseContext<I, N>,
        LinearProofCommitment<I, N>,
    ) {
        let gx = x
            .iter()
            .cloned()
            .map(|xi| xi.mul(g.clone()))
            .collect::<Vec<_>>(); // g * x
        let (cp, opening_p) = self.ck.commit(rng, &self.params, gx);
        let (c, opening) = self.ck.commit(rng, &self.params, x);

        // y <- N^k_sigma
        let y = Mat::<I, N>::new_with(self.params.k, 1, || {
            random_polynomial_in_normal_distribution::<I, N>(
                rng,
                I::zero().to_f64().unwrap(),
                self.params.standard_deviation(N) as f64,
            )
        });

        // yp <- N^k_sigma
        let yp = Mat::<I, N>::new_with(self.params.k, 1, || {
            random_polynomial_in_normal_distribution::<I, N>(
                rng,
                I::zero().to_f64().unwrap(),
                self.params.standard_deviation(N) as f64,
            )
        });

        // t = A1 * y
        let t = self.ck.a1.dot(&y).one_d_mat_to_vec();

        // tp = A1 * yp
        let tp = self.ck.a1.dot(&yp).one_d_mat_to_vec();

        // u = g * A2 * y - A2 * yp
        let u = self
            .ck
            .a2
            .dot(&y)
            .componentwise_mul(&g)
            .sub(&self.ck.a2.dot(&yp));

        (
            LinearProofResponseContext {
                opening,
                opening_p,
                y,
                yp,
            },
            LinearProofCommitment { c, cp, g, t, tp, u },
        )
    }

    pub fn create_response(
        &self,
        context: LinearProofResponseContext<I, N>,
        challenge: LinearProofChallenge<I, N>,
    ) -> LinearProofResponse<I, N> {
        // z = y + d * r
        let z = context
            .y
            .add(&context.opening.r.componentwise_mul(&challenge.d));
        // zp = yp + d * rp
        let zp = context
            .yp
            .add(&context.opening_p.r.componentwise_mul(&challenge.d));
        LinearProofResponse { z, zp }
    }
}

pub struct LinearProofVerifier<I, const N: usize>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I> + Neg<Output = I>,
{
    params: Params<I>,
    ck: CommitmentKey<I, N>,
}

impl<I, const N: usize> LinearProofVerifier<I, N>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I> + Neg<Output = I>,
{
    pub fn new(ck: CommitmentKey<I, N>, params: Params<I>) -> Self {
        LinearProofVerifier { params, ck }
    }

    pub fn generate_challenge(
        &self,
        rng: &mut impl Rng,
        commitment: LinearProofCommitment<I, N>,
    ) -> (
        LinearProofVerificationContext<I, N>,
        LinearProofChallenge<I, N>,
    ) {
        let d = random_polynomial_from_challenge_set(rng, self.params.kappa);
        let (c1, c2) = commitment.c.c1_c2(&self.params);
        let (c1p, c2p) = commitment.cp.c1_c2(&self.params);
        (
            LinearProofVerificationContext {
                c1,
                c2,
                c1p,
                c2p,
                g: commitment.g,
                t: commitment.t,
                tp: commitment.tp,
                u: commitment.u,
                d: d.clone(),
            },
            LinearProofChallenge { d },
        )
    }

    pub fn verify(
        &self,
        response: LinearProofResponse<I, N>,
        context: LinearProofVerificationContext<I, N>,
    ) -> bool {
        if !self.params.check_verify_constraint(&response.z) {
            return false;
        }
        if !self.params.check_verify_constraint(&response.zp) {
            return false;
        }
        // A1 * z = t + c1 * d
        let lhs = self.ck.a1.dot(&response.z);
        let rhs = Mat::<I, N>::from_vec(context.t).add(&context.c1.componentwise_mul(&context.d));
        if lhs != rhs {
            return false;
        }
        // A1 * zp = tp + c1p * d
        let lhs = self.ck.a1.dot(&response.zp);
        let rhs = Mat::<I, N>::from_vec(context.tp).add(&context.c1p.componentwise_mul(&context.d));
        if lhs != rhs {
            return false;
        }
        // g * A2 * z - A2 * zp = (g * c2 - c2p) * d + u
        let lhs = self
            .ck
            .a2
            .dot(&response.z)
            .componentwise_mul(&context.g)
            .sub(&self.ck.a2.dot(&response.zp));
        let rhs = context
            .c2
            .componentwise_mul(&context.g)
            .sub(&context.c2p)
            .componentwise_mul(&context.d)
            .add(&context.u);
        lhs == rhs
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LinearProofCommitment<I, const N: usize> {
    c: Commitment<I, N>,
    cp: Commitment<I, N>,
    g: Polynomial<I, N>,
    t: Vec<Polynomial<I, N>>,  // n x 1 matrix
    tp: Vec<Polynomial<I, N>>, // n x 1 matrix
    u: Mat<I, N>,              // l x 1 matrix
}

pub struct LinearProofVerificationContext<I, const N: usize> {
    c1: Mat<I, N>, // n x 1 matrix
    c2: Mat<I, N>, // l x 1 matrix

    c1p: Mat<I, N>, // n x 1 matrix
    c2p: Mat<I, N>, // l x 1 matrix

    g: Polynomial<I, N>,

    t: Vec<Polynomial<I, N>>,  // n x 1 matrix
    tp: Vec<Polynomial<I, N>>, // n x 1 matrix
    u: Mat<I, N>,              // l x 1 matrix
    d: Polynomial<I, N>,
}

pub struct LinearProofResponseContext<I, const N: usize> {
    opening: Opening<I, N>,
    opening_p: Opening<I, N>,
    y: Mat<I, N>,  // k x 1 matrix
    yp: Mat<I, N>, // k x 1 matrix
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LinearProofChallenge<I, const N: usize> {
    d: Polynomial<I, N>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LinearProofResponse<I, const N: usize> {
    z: Mat<I, N>,  // k x 1 matrix
    zp: Mat<I, N>, // k x 1 matrix
}
