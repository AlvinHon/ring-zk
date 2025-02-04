use std::{
    iter::Sum,
    ops::{Add, Mul, Sub},
};

use num::{integer::Roots, Integer, NumCast, Signed};
use poly_ring_xnp1::Polynomial;
use rand::Rng;
use rand_distr::uniform::SampleUniform;

use crate::{
    commit::{Commitment, CommitmentKey, Opening},
    mat::Mat,
    params::Params,
    polynomial::random_polynomial_in_normal_distribution,
};

pub fn prove_open<I, const N: usize>(
    rng: &mut impl Rng,
    ck: &CommitmentKey<I, N>,
    params: &Params<I>,
    x: Vec<Polynomial<I, N>>,
) -> (OpenProofProver<I, N>, OpenProofCommitment<I, N>)
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    let Params { k, .. } = params.clone();
    let (c, opening) = ck.commit(rng, params, x);

    // y <- N^k_sigma
    let y = Mat::<I, N>::new_with(k, 1, || {
        random_polynomial_in_normal_distribution::<I, N>(
            rng,
            I::zero().to_f64().unwrap(),
            params.standard_deviation(N) as f64,
        )
    });

    // t = A1 * y
    let t = ck.a1.dot(&y).to_vec();

    (OpenProofProver { opening, y }, OpenProofCommitment { c, t })
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpenProofCommitment<I, const N: usize> {
    c: Commitment<I, N>,
    t: Vec<Polynomial<I, N>>, // n x 1 matrix
}

impl<I, const N: usize> OpenProofCommitment<I, N> {
    pub fn create_challenge(
        &self,
        rng: &mut impl Rng,
        params: &Params<I>,
    ) -> (OpenProofVerifier, OpenProofChallenge)
    where
        I: Integer + Clone + SampleUniform,
    {
        // todo!()
        (OpenProofVerifier {}, OpenProofChallenge {})
    }
}

pub struct OpenProofVerifier {}

impl OpenProofVerifier {
    pub fn verify(&self, response: OpenProofResponse) -> bool {
        // todo!()
        true
    }
}

pub struct OpenProofChallenge {}

pub struct OpenProofResponse {}

pub struct OpenProofProver<I, const N: usize> {
    opening: Opening<I, N>,
    y: Mat<I, N>, // k x 1 matrix
}

impl<I, const N: usize> OpenProofProver<I, N> {
    pub fn create_response(&self, challenge: OpenProofChallenge) -> OpenProofResponse {
        // todo!()
        OpenProofResponse {}
    }
}

#[cfg(test)]
mod tests {
    use crate::params::params_1;

    use super::*;

    const N: usize = 4;

    #[test]
    fn test_prove_open() {
        let rng = &mut rand::rng();

        let params = params_1();
        let ck = CommitmentKey::new(rng, &params);
        let x = vec![Polynomial::<i64, N>::from_coeffs(vec![1, 2, 3, 4])];

        // 3-phase Sigma Protocol:
        // - First create commitment with information for proving the opening.
        let (prover, commitment) = prove_open(rng, &ck, &params, x);
        // - Verifier receives commitment and then create a challenge.
        let (verifier, challenge) = commitment.create_challenge(rng, &params);
        // - Prover receives the challenge and then create a response.
        let response = prover.create_response(challenge);
        // - Verifier verifies the response.
        assert!(verifier.verify(response));
    }
}
