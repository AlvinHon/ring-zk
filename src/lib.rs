pub(crate) mod challenge_space;
pub mod commit;
pub use commit::{Commitment, CommitmentKey, Opening};
pub(crate) mod mat;
pub mod params;
pub use params::{params_1, Params};
pub(crate) mod polynomial;
pub mod prove;
pub use prove::open::{
    prove_open, OpenProofChallenge, OpenProofCommitment, OpenProofProver, OpenProofResponse,
    OpenProofVerifier,
};

#[cfg(test)]
mod tests {
    use crate::{params_1, prove_open, CommitmentKey};

    const N: usize = 4;

    #[test]
    fn test_prove_open() {
        let rng = &mut rand::rng();

        let params = params_1();
        let ck = CommitmentKey::new(rng, &params);
        let x = vec![poly_ring_xnp1::Polynomial::<i64, N>::from_coeffs(vec![
            1, 2, 3, 4,
        ])];

        // 3-phase Sigma Protocol:
        // - First create commitment with information for proving the opening.
        let (prover, commitment) = prove_open(rng, x, &ck, &params);
        // - Verifier receives commitment and then create a challenge.
        let (verifier, challenge) = commitment.create_challenge(rng, &params);
        // - Prover receives the challenge and then create a response.
        let response = prover.create_response(challenge);
        // - Verifier verifies the response.
        assert!(verifier.verify(response, &ck, &params));
    }
}
