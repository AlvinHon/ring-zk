use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};
use poly_ring_xnp1::{rand::CoeffsRangeInclusive, zq::ZqI64, Polynomial};
use rand::{rng, RngExt};
use ring_zk::{
    LinearProofProver, LinearProofVerifier, OpenProofProver, OpenProofVerifier, Params,
    SumProofProver, SumProofVerifier,
};

criterion_group! {
    name = open_proof;
    config = Criterion::default().warm_up_time(Duration::from_secs(1)).sample_size(10).measurement_time(Duration::from_millis(1000));
    targets = bench_open_proof_commit, bench_open_proof_generate_challenge, bench_open_proof_create_response, bench_open_proof_verify,
}

criterion_group! {
    name = linear_proof;
    config = Criterion::default().warm_up_time(Duration::from_secs(1)).sample_size(10).measurement_time(Duration::from_millis(2000));
    targets = bench_linear_proof_commit, bench_linear_proof_generate_challenge, bench_linear_proof_create_response, bench_linear_proof_verify,
}

criterion_group! {
    name = sum_proof;
    config = Criterion::default().warm_up_time(Duration::from_secs(1)).sample_size(10).measurement_time(Duration::from_millis(4000));
    targets = bench_sum_proof_commit, bench_sum_proof_generate_challenge, bench_sum_proof_create_response, bench_sum_proof_verify,
}

criterion_main!(open_proof, linear_proof, sum_proof);

const N: usize = 512;

// ... bench functions for open_proof ...

fn bench_open_proof_commit(c: &mut Criterion) {
    let rng = &mut rng();
    let (params, prover, _) = setup_open_proof_elements();
    let x = params.prepare_value::<N>(vec![random_value(rng, params.q.clone().into())]);

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

    let (params, prover, verifier) = setup_open_proof_elements();
    let x = params.prepare_value::<N>(vec![random_value(rng, params.q.clone().into())]);

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

    let (params, prover, verifier) = setup_open_proof_elements();
    let x = params.prepare_value::<N>(vec![random_value(rng, params.q.clone().into())]);

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

    let (params, prover, verifier) = setup_open_proof_elements();
    let x = params.prepare_value::<N>(vec![random_value(rng, params.q.clone().into())]);

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

// ... bench functions for linear_proof ...

fn bench_linear_proof_commit(c: &mut Criterion) {
    let rng = &mut rng();

    let (params, prover, _) = setup_linear_proof_elements();
    let bound = params.q.clone().into();
    let x = params.prepare_value::<N>(vec![random_value(rng, bound)]);
    let g = params.prepare_scalar::<N>(random_value(rng, bound));

    c.bench_function("linear_proof_commit", |b| {
        b.iter_batched(
            || (g.clone(), x.clone()),
            |(g, x)| {
                _ = prover.commit(rng, g, x);
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

fn bench_linear_proof_generate_challenge(c: &mut Criterion) {
    let rng = &mut rng();

    let (params, prover, verifier) = setup_linear_proof_elements();
    let bound = params.q.clone().into();
    let x = params.prepare_value::<N>(vec![random_value(rng, bound)]);
    let g = params.prepare_scalar::<N>(random_value(rng, bound));

    let (_, commitment) = prover.commit(rng, g, x);

    c.bench_function("linear_proof_generate_challenge", |b| {
        b.iter_batched(
            || commitment.clone(),
            |commitment| {
                _ = verifier.generate_challenge(rng, commitment);
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

fn bench_linear_proof_create_response(c: &mut Criterion) {
    let rng = &mut rng();

    let (params, prover, verifier) = setup_linear_proof_elements();
    let bound = params.q.clone().into();
    let x = params.prepare_value::<N>(vec![random_value(rng, bound)]);
    let g = params.prepare_scalar::<N>(random_value(rng, bound));

    let (response_ctx, commitment) = prover.commit(rng, g, x);
    let (_, challenge) = verifier.generate_challenge(rng, commitment);

    c.bench_function("linear_proof_create_response", |b| {
        b.iter_batched(
            || (response_ctx.clone(), challenge.clone()),
            |(response_ctx, challenge)| {
                _ = prover.create_response(response_ctx, challenge);
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

fn bench_linear_proof_verify(c: &mut Criterion) {
    let rng = &mut rng();

    let (params, prover, verifier) = setup_linear_proof_elements();
    let bound = params.q.clone().into();
    let x = params.prepare_value::<N>(vec![random_value(rng, bound)]);
    let g = params.prepare_scalar::<N>(random_value(rng, bound));

    let (response_ctx, commitment) = prover.commit(rng, g, x);
    let (verification_ctx, challenge) = verifier.generate_challenge(rng, commitment);
    let response = prover.create_response(response_ctx, challenge);

    c.bench_function("linear_proof_verify", |b| {
        b.iter_batched(
            || (verification_ctx.clone(), response.clone()),
            |(verification_ctx, response)| {
                verifier.verify(response, verification_ctx);
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

// ... bench functions for sum_proof ...

const VL: usize = 4; // number of variables

fn bench_sum_proof_commit(c: &mut Criterion) {
    let rng = &mut rng();

    let (params, prover, _) = setup_sum_proof_elements();
    let bound = params.q.clone().into();

    let xs = (0..VL)
        .map(|_| params.prepare_value::<N>(vec![random_value(rng, bound)]))
        .collect::<Vec<_>>();
    let gs = (0..VL)
        .map(|_| params.prepare_scalar::<N>(random_value(rng, bound)))
        .collect::<Vec<_>>();

    c.bench_function("sum_proof_commit", |b| {
        b.iter_batched(
            || (gs.clone(), xs.clone()),
            |(gs, xs)| {
                _ = prover.commit(rng, gs, xs);
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

fn bench_sum_proof_generate_challenge(c: &mut Criterion) {
    let rng = &mut rng();

    let (params, prover, verifier) = setup_sum_proof_elements();
    let bound = params.q.clone().into();

    let xs = (0..VL)
        .map(|_| params.prepare_value::<N>(vec![random_value(rng, bound)]))
        .collect::<Vec<_>>();
    let gs = (0..VL)
        .map(|_| params.prepare_scalar::<N>(random_value(rng, bound)))
        .collect::<Vec<_>>();

    let (_, commitment) = prover.commit(rng, gs, xs);

    c.bench_function("sum_proof_commit", |b| {
        b.iter_batched(
            || commitment.clone(),
            |commitment| {
                _ = verifier.generate_challenge(rng, commitment);
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

fn bench_sum_proof_create_response(c: &mut Criterion) {
    let rng = &mut rng();

    let (params, prover, verifier) = setup_sum_proof_elements();
    let bound = params.q.clone().into();

    let xs = (0..VL)
        .map(|_| params.prepare_value::<N>(vec![random_value(rng, bound)]))
        .collect::<Vec<_>>();
    let gs = (0..VL)
        .map(|_| params.prepare_scalar::<N>(random_value(rng, bound)))
        .collect::<Vec<_>>();

    let (response_ctx, commitment) = prover.commit(rng, gs, xs);
    let (_, challenge) = verifier.generate_challenge(rng, commitment);

    c.bench_function("sum_proof_create_response", |b| {
        b.iter_batched(
            || (response_ctx.clone(), challenge.clone()),
            |(response_ctx, challenge)| {
                _ = prover.create_response(response_ctx, challenge);
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

fn bench_sum_proof_verify(c: &mut Criterion) {
    let rng = &mut rng();

    let (params, prover, verifier) = setup_sum_proof_elements();
    let bound = params.q.clone().into();

    let xs = (0..VL)
        .map(|_| params.prepare_value::<N>(vec![random_value(rng, bound)]))
        .collect::<Vec<_>>();
    let gs = (0..VL)
        .map(|_| params.prepare_scalar::<N>(random_value(rng, bound)))
        .collect::<Vec<_>>();

    let (response_ctx, commitment) = prover.commit(rng, gs, xs);
    let (verification_ctx, challenge) = verifier.generate_challenge(rng, commitment);
    let response = prover.create_response(response_ctx, challenge);

    c.bench_function("sum_proof_verify", |b| {
        b.iter_batched(
            || (verification_ctx.clone(), response.clone()),
            |(verification_ctx, response)| {
                verifier.verify(response, verification_ctx);
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

// ... utility functions ...

fn setup_open_proof_elements() -> (
    Params<ZqI64<3515337053_i64>>,
    OpenProofProver<ZqI64<3515337053_i64>, N>,
    OpenProofVerifier<ZqI64<3515337053_i64>, N>,
) {
    let rng = &mut rng();

    let params = Params::default();
    let ck = params.generate_commitment_key(rng);
    let prover = OpenProofProver::new(ck.clone(), params.clone());
    let verifier = OpenProofVerifier::new(ck.clone(), params.clone());

    (params, prover, verifier)
}

fn setup_linear_proof_elements() -> (
    Params<ZqI64<3515337053_i64>>,
    LinearProofProver<ZqI64<3515337053_i64>, N>,
    LinearProofVerifier<ZqI64<3515337053_i64>, N>,
) {
    let rng = &mut rng();

    let params = Params::default();
    let ck = params.generate_commitment_key(rng);
    let prover = LinearProofProver::new(ck.clone(), params.clone());
    let verifier = LinearProofVerifier::new(ck.clone(), params.clone());

    (params, prover, verifier)
}

fn setup_sum_proof_elements() -> (
    Params<ZqI64<3515337053_i64>>,
    SumProofProver<ZqI64<3515337053_i64>, N>,
    SumProofVerifier<ZqI64<3515337053_i64>, N>,
) {
    let rng = &mut rng();

    let params = Params::default();
    let ck = params.generate_commitment_key(rng);
    let prover = SumProofProver::new(ck.clone(), params.clone());
    let verifier = SumProofVerifier::new(ck.clone(), params.clone());

    (params, prover, verifier)
}

fn random_value(rng: &mut impl RngExt, bound: i64) -> Vec<i64> {
    let range = CoeffsRangeInclusive::from(-bound..=bound);
    let p: Polynomial<i64, N> = rng.random_range(range);
    p.iter().copied().take(rng.random_range(1..=N)).collect()
}
