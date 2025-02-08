pub(crate) mod challenge_space;
pub(crate) mod commit;
pub use commit::CommitmentKey;
pub(crate) mod mat;
pub mod params;
pub use params::Params;
pub(crate) mod polynomial;
pub mod prove;
pub use prove::{
    linear::{
        LinearProofChallenge, LinearProofCommitment, LinearProofProver, LinearProofResponse,
        LinearProofResponseContext, LinearProofVerificationContext, LinearProofVerifier,
    },
    open::{
        OpenProofChallenge, OpenProofCommitment, OpenProofProver, OpenProofResponse,
        OpenProofResponseContext, OpenProofVerificationContext, OpenProofVerifier,
    },
};

#[cfg(test)]
mod tests {
    use crate::{
        CommitmentKey, LinearProofProver, LinearProofVerifier, OpenProofProver, OpenProofVerifier,
        Params,
    };

    const N: usize = 4;

    #[test]
    fn test_prove_open() {
        let rng = &mut rand::rng();

        let params = Params::default();
        let ck = params.generate_commitment_key(rng);
        let x = params.prepare_value::<N>(vec![vec![1, 2, 3, 4]]);

        let prover = OpenProofProver::new(ck.clone(), params.clone());
        let verifier = OpenProofVerifier::new(ck, params);

        // 3-phase Sigma Protocol:
        // - First create commitment with information for proving the opening.
        let (response_ctx, commitment) = prover.commit_and_prove(rng, x);
        // - Verifier receives commitment and then create a challenge.
        let (verification_ctx, challenge) = verifier.generate_challenge(rng, commitment);
        // - Prover receives the challenge and then create a response.
        let response = prover.create_response(response_ctx, challenge);
        // - Verifier verifies the response.
        assert!(verifier.verify(response, verification_ctx));
    }

    #[test]
    fn test_prove_linear() {
        let rng = &mut rand::rng();

        let params = Params::default();
        let ck = CommitmentKey::new(rng, &params);
        let x = params.prepare_value::<N>(vec![vec![1, 2, 3, 4]]);
        let g = params.prepare_scalar::<N>(vec![5, 6]);

        let prover = LinearProofProver::new(ck.clone(), params.clone());
        let verifier = LinearProofVerifier::new(ck, params);

        // 3-phase Sigma Protocol:
        // - First create commitment with information for proving the linear relationship of the committed value.
        let (response_ctx, commitment) = prover.commit_and_prove(rng, g, x);
        // - Verifier receives commitment and then create a challenge.
        let (verification_ctx, challenge) = verifier.generate_challenge(rng, commitment);
        // - Prover receives the challenge and then create a response.
        let response = prover.create_response(response_ctx, challenge);
        // - Verifier verifies the response.
        assert!(verifier.verify(response, verification_ctx));
    }
}
