//! Implementation of Proof of Linear Relation, defined in section 4.4 of the paper.
//!
//! This modules contains struct [LinearProofProver] and [LinearProofVerifier] for proving and verifying
//! opening of commitments ([LinearProofCommitment]) to `x'` and `x` such that `x' = g * x` for scalar `g`.
//! The prover and verifier will exchange messages [LinearProofChallenge] and [LinearProofResponse] to
//! complete the 3-phase Sigma Protocol.
//! The opening is encapsulated in [LinearProofResponseContext] which is created and used by prover in the
//! protocol. The verifier generates the challenge and verifies the response by using the context
//! [LinearProofVerificationContext].
//!
//!
//! ## Example
//!
//! ```rust
//! use ring_zk::{Params, LinearProofProver, LinearProofVerifier};
//!
//! const N: usize = 512;
//!
//! let rng = &mut rand::rng();
//!
//! let params = Params::default();
//! let ck = params.generate_commitment_key(rng);
//! let x = params.prepare_value::<N>(vec![vec![1, 2, 3, 4]]);
//! let g = params.prepare_scalar::<N>(vec![5, 6]);
//!
//! let prover = LinearProofProver::new(ck.clone(), params.clone());
//! let verifier = LinearProofVerifier::new(ck.clone(), params.clone());
//!
//! // 3-phase Sigma Protocol:
//! // - First create commitment with information for proving the linear relationship of the committed value.
//! let (response_ctx, commitment) = prover.commit(rng, g, x);
//! // - Verifier receives commitment and then create a challenge.
//! let (verification_ctx, challenge) = verifier.generate_challenge(rng, commitment);
//! // - Prover receives the challenge and then create a response.
//! let response = prover.create_response(response_ctx, challenge);
//! // - Verifier verifies the response.
//! assert!(verifier.verify(response, verification_ctx));
//! ```

use std::ops::{Add, Mul, Neg, Sub};

use num::{FromPrimitive, One, ToPrimitive, Zero};
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

/// The prover for the proof of linear relation. It is used to prove that the prover knows the
/// openings of commitments to `x'` and `x` such that `x' = g * x` for scalar `g`.
pub struct LinearProofProver<I, const N: usize> {
    params: Params<I>,
    ck: CommitmentKey<I, N>,
}

impl<I, const N: usize> LinearProofProver<I, N>
where
    I: Clone + PartialOrd + Ord + One + Zero + FromPrimitive + ToPrimitive + SampleUniform,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Neg<Output = I> + Sub<Output = I>,
{
    pub fn new(ck: CommitmentKey<I, N>, params: Params<I>) -> Self {
        Self { params, ck }
    }

    /// Create commitments to `x'` and `x` such that `x' = g * x` for scalar `g`.
    /// It returns the response context and the commitment. The response context is used to create
    /// the response in a later phase of the protocol. Note that the context includes the openings
    /// of commitments to `x'` and `x`.
    ///
    /// ## Panics
    /// Panics if the length of `x` is not equal to the length of `l` defined in the `Params` struct.
    pub fn commit(
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
        let (opening_p, cp) = self.ck.commit(rng, gx, &self.params);
        let (opening, c) = self.ck.commit(rng, x, &self.params);

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

    /// Create the response for the challenge received from the verifier. The response is created
    /// using the context that was created during the commitment phase.
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

/// The verifier for the proof of linear relation. It is used to verify that the prover knows the
/// openings of commitments to `x'` and `x` such that `x' = g * x` for scalar `g`.
pub struct LinearProofVerifier<I, const N: usize> {
    params: Params<I>,
    ck: CommitmentKey<I, N>,
}

impl<I, const N: usize> LinearProofVerifier<I, N>
where
    I: Clone + PartialOrd + Ord + One + Zero + FromPrimitive + ToPrimitive + SampleUniform,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Neg<Output = I> + Sub<Output = I>,
{
    pub fn new(ck: CommitmentKey<I, N>, params: Params<I>) -> Self {
        LinearProofVerifier { params, ck }
    }

    /// Generate the challenge for the prover, given the commitments that says the prover knows its
    /// openings to the commitments to values `x'` and `x` such that `x' = g * x` for scalar `g`.
    /// It returns the verification context and the challenge. The verification context is used to
    /// verify the response in a later phase of the protocol.
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

    /// Verify the response from the prover. It returns `true` if the response is valid, otherwise `false`.
    /// The context was created during the challenge phase in the protocol.
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

/// The response created by the prover upon receiving the challenge from the verifier
/// in the protocol of proof of linear relation. It contains the openings of commitments
/// to `x'` and `x` such that `x' = g * x` for scalar `g`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LinearProofResponseContext<I, const N: usize> {
    /// The opening of the commitment to `x` s.t. `x' = g * x`.
    pub opening: Opening<I, N>,
    /// The opening of the commitment to `x'` s.t. `x' = g * x`.
    pub opening_p: Opening<I, N>,
    y: Mat<I, N>,  // k x 1 matrix
    yp: Mat<I, N>, // k x 1 matrix
}

/// Contains the commitments to the values `x'` and `x` such that `x' = g * x`, used in
/// the proof of linear relation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LinearProofCommitment<I, const N: usize> {
    /// Commitment to value `x` s.t. `x' = g * x`.
    pub c: Commitment<I, N>,
    /// Commitment to value `x'` s.t. `x' = g * x`.
    pub cp: Commitment<I, N>,
    /// The scalar `g` in the relation `x' = g * x`.
    pub g: Polynomial<I, N>,
    t: Vec<Polynomial<I, N>>,  // n x 1 matrix
    tp: Vec<Polynomial<I, N>>, // n x 1 matrix
    u: Mat<I, N>,              // l x 1 matrix
}

/// Contains the context for the verification phase of the proof of linear relation.
/// It is used to verify the response from the prover.
#[derive(Clone, Debug, PartialEq, Eq)]
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

/// The challenge created by the verifier in the protocol of proof of linear relation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LinearProofChallenge<I, const N: usize> {
    d: Polynomial<I, N>,
}

/// The response from the prover to the verifier in the protocol of proof of linear relation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LinearProofResponse<I, const N: usize> {
    z: Mat<I, N>,  // k x 1 matrix
    zp: Mat<I, N>, // k x 1 matrix
}
