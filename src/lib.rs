pub(crate) mod challenge_space;
pub mod commit;
pub use commit::{Commitment, CommitmentKey, Opening};
pub(crate) mod mat;
pub mod params;
pub use params::Params;
pub(crate) mod polynomial;
pub mod prove;
pub use prove::{
    linear::{
        prove_linear, LinearProofChallenge, LinearProofCommitment, LinearProofProver,
        LinearProofResponse, LinearProofVerifier,
    },
    open::{
        prove_open, OpenProofChallenge, OpenProofCommitment, OpenProofProver, OpenProofResponse,
        OpenProofVerifier,
    },
};

#[cfg(test)]
mod tests {
    use crate::{prove_linear, prove_open, CommitmentKey, Params};

    const N: usize = 4;

    #[test]
    fn test_prove_open() {
        let rng = &mut rand::rng();

        let params = Params::default();
        let ck = CommitmentKey::new(rng, &params);
        let x = params.prepare_value::<N>(vec![vec![1, 2, 3, 4]]);

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

    #[test]
    fn test_prove_linear() {
        let rng = &mut rand::rng();

        let params = Params::default();
        let ck = CommitmentKey::new(rng, &params);
        let x = params.prepare_value::<N>(vec![vec![1, 2, 3, 4]]);
        let g = params.prepare_scalar::<N>(vec![5, 6]);

        // 3-phase Sigma Protocol:
        // - First create commitment with information for proving the linear relationship of the committed value.
        let (prover, commitment) = prove_linear(rng, g, x, &ck, &params);
        // - Verifier receives commitment and then create a challenge.
        let (verifier, challenge) = commitment.create_challenge(rng, &params);
        // - Prover receives the challenge and then create a response.
        let response = prover.create_response(challenge);
        // - Verifier verifies the response.
        assert!(verifier.verify(response, &ck, &params));
    }
}
