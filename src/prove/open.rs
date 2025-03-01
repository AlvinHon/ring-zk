//! Implementation of Proof of Opening a Commitment, defined in section 4.4 of the paper.
//!
//! This modules contains struct [OpenProofProver] and [OpenProofVerifier] for proving and verifying
//! opening of a commitment ([OpenProofCommitment]) to a value.
//! The prover and verifier will exchange messages [OpenProofChallenge] and [OpenProofResponse] to
//! complete the 3-phase Sigma Protocol.
//! The opening is encapsulated in [OpenProofResponseContext] which is created and used by prover in the
//! protocol. The verifier generates the challenge and verifies the response by using the context
//! [OpenProofVerificationContext].
//!
//!
//! ## Example
//!
//! ```rust
//! use ring_zk::{Params, OpenProofProver, OpenProofVerifier};
//!
//! const N: usize = 512;
//!
//! let rng = &mut rand::rng();
//!
//! let params = Params::default();
//! let ck = params.generate_commitment_key(rng);
//! let x = params.prepare_value::<N>(vec![vec![1, 2, 3, 4]]);
//!
//! let prover = OpenProofProver::new(ck.clone(), params.clone());
//! let verifier = OpenProofVerifier::new(ck.clone(), params.clone());
//! // 3-phase Sigma Protocol:
//! // - First create commitment with information for proving the opening.
//! let (response_ctx, commitment) = prover.commit(rng, x);
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
use serde::{Deserialize, Serialize};

use crate::{
    challenge_space::random_polynomial_from_challenge_set,
    commit::{Commitment, CommitmentKey, Opening},
    mat::Mat,
    params::Params,
    polynomial::random_polynomial_in_normal_distribution,
};

/// The prover for the proof of linear relation. It is used to prove that the prover knows the
/// opening of commitment to a value.
pub struct OpenProofProver<I, const N: usize>
where
    I: Zero,
{
    params: Params<I>,
    ck: CommitmentKey<I, N>,
}

impl<I, const N: usize> OpenProofProver<I, N>
where
    I: Clone + PartialOrd + Ord + One + Zero + FromPrimitive + ToPrimitive + SampleUniform,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    pub fn new(ck: CommitmentKey<I, N>, params: Params<I>) -> Self {
        Self { params, ck }
    }

    /// Create commitments to the value `x`.
    /// It returns the response context and the commitment. The response context is used to create
    /// the response in a later phase of the protocol. Note that the context includes the openings
    /// of commitments to `x`.
    ///
    /// ## Panics
    /// Panics if the length of `x` is not equal to the length of `l` defined in the `Params` struct.
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

    /// Create the response for the challenge received from the verifier. The response is created
    /// using the context that was created during the commitment phase.
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

/// The verifier for the proof of opening a commitment. It is used to verify that the prover knows
/// the opening of commitment to a value.
pub struct OpenProofVerifier<I, const N: usize>
where
    I: Zero,
{
    params: Params<I>,
    ck: CommitmentKey<I, N>,
}

impl<I, const N: usize> OpenProofVerifier<I, N>
where
    I: Clone + PartialOrd + Ord + One + Zero + FromPrimitive + ToPrimitive + SampleUniform,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Neg<Output = I> + Sub<Output = I>,
{
    pub fn new(ck: CommitmentKey<I, N>, params: Params<I>) -> Self {
        OpenProofVerifier { params, ck }
    }

    /// Generate the challenge for the prover, given the commitments that says the prover knows its
    /// opening to the commitment to a value.
    /// It returns the verification context and the challenge. The verification context is used to
    /// verify the response in a later phase of the protocol.
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

    /// Verify the response from the prover. It returns `true` if the response is valid, otherwise `false`.
    /// The context was created during the challenge phase in the protocol.
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

/// The response created by the prover upon receiving the challenge from the verifier
/// in the protocol of proof of opening a commitment. It contains the opening of commitment
/// to value `x`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct OpenProofResponseContext<I, const N: usize>
where
    I: Zero,
{
    pub opening: Opening<I, N>,
    y: Mat<I, N>, // k x 1 matrix
}

/// Contains the commitment to the values `x`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct OpenProofCommitment<I, const N: usize>
where
    I: Zero,
{
    /// Commitment to value `x`.
    pub c: Commitment<I, N>,
    t: Vec<Polynomial<I, N>>, // n x 1 matrix
}

/// Contains the context for the verification phase of the proof of opening a commitment.
/// It is used to verify the response from the prover.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct OpenProofVerificationContext<I, const N: usize>
where
    I: Zero,
{
    c1: Mat<I, N>,            // n x 1 matrix
    t: Vec<Polynomial<I, N>>, // n x 1 matrix
    d: Polynomial<I, N>,
}

/// The challenge created by the verifier in the protocol of proof of opening a commitment.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct OpenProofChallenge<I, const N: usize>
where
    I: Zero,
{
    d: Polynomial<I, N>,
}

/// The response from the prover to the verifier in the protocol of proof of opening a commitment.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct OpenProofResponse<I, const N: usize>
where
    I: Zero,
{
    z: Mat<I, N>, // k x 1 matrix
}
