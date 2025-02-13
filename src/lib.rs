#![doc = include_str!("../README.md")]

pub(crate) mod challenge_space;
pub(crate) mod commit;
pub use commit::{Commitment, CommitmentKey, Opening};
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
    sum::{
        SumProofChallenge, SumProofCommitment, SumProofProver, SumProofResponse,
        SumProofResponseContext, SumProofVerificationContext, SumProofVerifier,
    },
};
