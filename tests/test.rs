use poly_ring_xnp1::{rand::CoeffsRangeInclusive, Polynomial};
use rand::Rng;
use ring_zk::{
    LinearProofProver, LinearProofVerifier, OpenProofProver, OpenProofVerifier, Params,
    SumProofProver, SumProofVerifier,
};

const N: usize = 16; // power of two. Should be reasonably long.

/// Test the open proof by generating random inputs over numerous iterations.
#[test]
fn test_open_proof() {
    let rng = &mut rand::rng();

    let params = Params::default();
    let bound = params.q.clone().into();

    for _ in 0..100 {
        let ck = params.generate_commitment_key(rng);
        let x = params.prepare_value::<N>(vec![random_value(rng, bound)]);

        let prover = OpenProofProver::new(ck.clone(), params.clone());
        let verifier = OpenProofVerifier::new(ck.clone(), params.clone());

        let (response_ctx, commitment) = prover.commit(rng, x);
        assert!(commitment.c.verify(&response_ctx.opening, &ck, &params));
        let (verification_ctx, challenge) = verifier.generate_challenge(rng, commitment);
        let response = prover.create_response(response_ctx, challenge);
        assert!(verifier.verify(response, verification_ctx));
    }
}

/// Test the linear proof by generating random inputs over numerous iterations.
#[test]
fn test_linear_proof() {
    let rng = &mut rand::rng();

    let params = Params::default();
    let bound = params.q.clone().into();

    for _ in 0..100 {
        let ck = params.generate_commitment_key(rng);
        let x = params.prepare_value::<N>(vec![random_value(rng, bound)]);
        let g = params.prepare_scalar::<N>(random_value(rng, bound));

        let prover = LinearProofProver::new(ck.clone(), params.clone());
        let verifier = LinearProofVerifier::new(ck.clone(), params.clone());

        let (response_ctx, commitment) = prover.commit(rng, g, x);
        assert!(commitment.c.verify(&response_ctx.opening, &ck, &params));
        assert!(commitment.cp.verify(&response_ctx.opening_p, &ck, &params));
        let (verification_ctx, challenge) = verifier.generate_challenge(rng, commitment);
        let response = prover.create_response(response_ctx, challenge);
        assert!(verifier.verify(response, verification_ctx));
    }
}

/// Test the sum proof by generating random inputs over numerous iterations.
#[test]
fn test_sum_proof() {
    let rng = &mut rand::rng();

    let params = Params::default();
    let bound = params.q.clone().into();
    const VL: usize = 4;

    for _ in 0..100 {
        let ck = params.generate_commitment_key(rng);

        let xs = (0..VL)
            .map(|_| params.prepare_value::<N>(vec![random_value(rng, bound)]))
            .collect::<Vec<_>>();
        let gs = (0..VL)
            .map(|_| params.prepare_scalar::<N>(random_value(rng, bound)))
            .collect::<Vec<_>>();

        let prover = SumProofProver::new(ck.clone(), params.clone());
        let verifier = SumProofVerifier::new(ck.clone(), params.clone());

        let (response_ctx, commitment) = prover.commit(rng, gs, xs);
        commitment.cp.verify(&response_ctx.opening_p, &ck, &params);
        commitment
            .cs
            .iter()
            .zip(response_ctx.openings.iter())
            .for_each(|(c, o)| {
                assert!(c.verify(o, &ck, &params));
            });
        let (verification_ctx, challenge) = verifier.generate_challenge(rng, commitment);
        let response = prover.create_response(response_ctx, challenge);
        assert!(verifier.verify(response, verification_ctx));
    }
}

pub(crate) fn random_value(rng: &mut impl Rng, bound: i64) -> Vec<i64> {
    let range = CoeffsRangeInclusive::from(-bound..=bound);
    let p: Polynomial<i64, N> = rng.random_range(range);
    p.iter().copied().take(rng.random_range(1..=N)).collect()
}
