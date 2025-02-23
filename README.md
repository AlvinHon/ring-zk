# Efficient Zero-Knowledge Proofs for Commitments from RLWE

Rust implementation of the zero-knowledge proof system for commitments from the Ring Learning With Errors (RLWE) assumption. The proof system is based on the paper [More Efficient Commitments from Structured Lattice Assumptions](https://eprint.iacr.org/2016/997).

> .. this work is a construction of an efficient commitment scheme and accompanying zero-knowledge proofs of knowledge for proving relations among committed values.

The ZK proofs provided in this library have several properties:
- Use lattice-based cryptographic assumptions which are **"post-quantum"** replacements for the discrete logarithm and factoring problem.
- Use an **interactive** protocol, but theoretically can be made non-interactive using the Fiat-Shamir transform.
- Prove the knowledge of valid openings for the committed values and their relations, but **not** the knowledge of the committed values themselves.

## Message Space

A message is the commited value to the commitment scheme used in the proof system.

The message space is a matrix of integers. You can consider the matrix is a vector of polynomials of degree `N-1` with coefficients modulo `q`. The parameters `N` and `q` are defined in the `Params` struct, as well as other parameters for the commitment scheme.

In this implementation, the input message is represented as `Vec<Vec<_>>` where the first dimension is the number of polynomials (limited by the parameter `l`) and the second dimension is the coefficients of the polynomial (limited by the parameter `N`).


## Proof for Opening a Commitment

The prover wants to prove that they know the opening of a commitment to a value `x`.

```rust
use ring_zk::{Params, OpenProofProver, OpenProofVerifier};

const N: usize = 512; // maximum size of message. Must be power of two.

let rng = &mut rand::rng();

let params = Params::default();
let ck = params.generate_commitment_key(rng);
let x = params.prepare_value::<N>(vec![vec![1, 2, 3, 4]]);

let prover = OpenProofProver::new(ck.clone(), params.clone());
let verifier = OpenProofVerifier::new(ck.clone(), params.clone());

// 3-phase Sigma Protocol:
// - First create commitment with information for proving the opening.
let (response_ctx, commitment) = prover.commit(rng, x);
// - Verifier receives commitment and then create a challenge.
let (verification_ctx, challenge) = verifier.generate_challenge(rng, commitment);
// - Prover receives the challenge and then create a response.
let response = prover.create_response(response_ctx, challenge);
// - Verifier verifies the response.
assert!(verifier.verify(response, verification_ctx));
```

## Proof of Relation between Commitments

**Proof of Linear Relation**

The prover wants to prove that they know the openings of commitments `c'` and a commitment to values `x'` and `x` s.t. `x' = g * x`, where `g` is a scalar represented as a vector of integers modulo `q`.

The interface will be similar to the proof of opening a commitment, but with modification: we use the struct `LinearProofProver` and `LinearProofVerifier`.

```rust ignore
// ...
let x = params.prepare_value::<N>(vec![vec![1, 2, 3, 4]]);
let g = params.prepare_scalar::<N>(vec![5, 6]);

let prover = LinearProofProver::new(ck.clone(), params.clone());
let verifier = LinearProofVerifier::new(ck.clone(), params.clone());

// 3-phase Sigma Protocol:
// - First create commitment with information for proving the linear relationship of the committed value.
let (response_ctx, commitment) = prover.commit(rng, g, x);
// ... the same interaction ...
```

**Proof of Sum**

The prover wants to prove that they know the openings of commitments `c'` and a vector of commitments to values `x'` and a vector of values s.t. `x' = g_0 * x_0 + g_1 * x_1 + ...` where `g_i` represents a scalar in a vector of integers modulo `q`.

Similar to the proof of linear relation, we use the struct `SumProofProver` and `SumProofVerifier`.

```rust ignore
// ...
let xs = vec![
    params.prepare_value::<N>(vec![vec![1, 2, 3, 4]]),
    params.prepare_value::<N>(vec![vec![5, 6, 7, 8]]),
];
let gs = vec![
    params.prepare_scalar::<N>(vec![5, 6]),
    params.prepare_scalar::<N>(vec![7, 8]),
];

let prover = SumProofProver::new(ck.clone(), params.clone());
let verifier = SumProofVerifier::new(ck.clone(), params.clone());

// 3-phase Sigma Protocol:
// - First create commitment with information for proving the summation relationship of the committed values.
let (response_ctx, commitment) = prover.commit(rng, gs, xs);
// ... the same interaction ...
```

***In general, Proof of Sum can replace Proof of Linear Relation. Proof of Linear Relation is implemented here for respect to the paper.***

## References

- [More Efficient Commitments from Structured Lattice Assumptions](https://eprint.iacr.org/2016/997)
- [Efficient Zero-Knowledge Proofs for Commitments from Learning With Errors over Rings](https://eprint.iacr.org/2014/889)