use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};
use rand::rng;
use ring_zk::{OpenProofProver, OpenProofVerifier, Params};

criterion_group! {
    name = open_proof;
    config = Criterion::default().warm_up_time(Duration::from_secs(1)).sample_size(10).measurement_time(Duration::from_millis(1000));
    targets = bench_open_proof_commit, bench_open_proof_generate_challenge, bench_open_proof_create_response, bench_open_proof_verify,
}

criterion_main!(open_proof);

const N: usize = 512;

fn bench_open_proof_commit(c: &mut Criterion) {
    let rng = &mut rng();

    let params = Params::default();
    let ck = params.generate_commitment_key(rng);
    let x = params.prepare_value::<N>(vec![vec![1, 2, 3, 4]]);

    let prover = OpenProofProver::new(ck, params);

    c.bench_function("open_proof_commit", |b| {
        b.iter_batched(
            || x.clone(),
            |x| {
                _ = prover.commit(rng, x);
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

fn bench_open_proof_generate_challenge(c: &mut Criterion) {
    let rng = &mut rng();

    let params = Params::default();
    let ck = params.generate_commitment_key(rng);
    let x = params.prepare_value::<N>(vec![vec![1, 2, 3, 4]]);

    let prover = OpenProofProver::new(ck.clone(), params.clone());
    let verifier = OpenProofVerifier::new(ck.clone(), params.clone());

    let (_, commitment) = prover.commit(rng, x);

    c.bench_function("open_proof_generate_challenge", |b| {
        b.iter_batched(
            || commitment.clone(),
            |commitment| {
                _ = verifier.generate_challenge(rng, commitment);
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

fn bench_open_proof_create_response(c: &mut Criterion) {
    let rng = &mut rng();

    let params = Params::default();
    let ck = params.generate_commitment_key(rng);
    let x = params.prepare_value::<N>(vec![vec![1, 2, 3, 4]]);

    let prover = OpenProofProver::new(ck.clone(), params.clone());
    let verifier = OpenProofVerifier::new(ck.clone(), params.clone());

    let (response_ctx, commitment) = prover.commit(rng, x);
    let (_, challenge) = verifier.generate_challenge(rng, commitment);

    c.bench_function("open_proof_create_response", |b| {
        b.iter_batched(
            || (response_ctx.clone(), challenge.clone()),
            |(response_ctx, challenge)| {
                _ = prover.create_response(response_ctx, challenge);
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

fn bench_open_proof_verify(c: &mut Criterion) {
    let rng = &mut rng();

    let params = Params::default();
    let ck = params.generate_commitment_key(rng);
    let x = params.prepare_value::<N>(vec![vec![1, 2, 3, 4]]);

    let prover = OpenProofProver::new(ck.clone(), params.clone());
    let verifier = OpenProofVerifier::new(ck.clone(), params.clone());

    let (response_ctx, commitment) = prover.commit(rng, x);
    let (verification_ctx, challenge) = verifier.generate_challenge(rng, commitment);
    let response = prover.create_response(response_ctx, challenge);

    c.bench_function("open_proof_verify", |b| {
        b.iter_batched(
            || (verification_ctx.clone(), response.clone()),
            |(verification_ctx, response)| {
                verifier.verify(response, verification_ctx);
            },
            criterion::BatchSize::SmallInput,
        )
    });
}
