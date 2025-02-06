use std::{
    iter::Sum,
    ops::{Add, Mul, Sub},
};

use num::{integer::Roots, Integer, NumCast, One, Signed, Zero};
use poly_ring_xnp1::Polynomial;
use rand::Rng;
use rand_distr::uniform::SampleUniform;

use crate::{
    challenge_space::random_polynomial_from_challenge_set,
    commit::{Commitment, CommitmentKey, Opening},
    mat::Mat,
    params::Params,
    polynomial::random_polynomial_in_normal_distribution,
};

pub fn prove_open<I, const N: usize>(
    rng: &mut impl Rng,
    x: Vec<Polynomial<I, N>>,
    ck: &CommitmentKey<I, N>,
    params: &Params<I>,
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
    let t = ck.a1.dot(&y).one_d_mat_to_vec();

    (OpenProofProver { opening, y }, OpenProofCommitment { c, t })
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpenProofCommitment<I, const N: usize>
where
    I: Integer + Signed + Sum + Roots + Clone + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    c: Commitment<I, N>,
    t: Vec<Polynomial<I, N>>, // n x 1 matrix
}

impl<I, const N: usize> OpenProofCommitment<I, N>
where
    I: Integer + Signed + Sum + Roots + Clone + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    pub fn create_challenge(
        &self,
        rng: &mut impl Rng,
        params: &Params<I>,
    ) -> (OpenProofVerifier<I, N>, OpenProofChallenge<I, N>)
    where
        I: Integer + Clone + SampleUniform,
    {
        let d = random_polynomial_from_challenge_set(rng, params.kappa);
        let (c1, _) = self.c.c1_c2(params);
        (
            OpenProofVerifier {
                c1,
                t: self.t.clone(),
                d: d.clone(),
            },
            OpenProofChallenge { d },
        )
    }
}

pub struct OpenProofVerifier<I, const N: usize>
where
    I: Integer + Signed + Sum + Roots + Clone + NumCast + SampleUniform,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    c1: Mat<I, N>,            // n x 1 matrix
    t: Vec<Polynomial<I, N>>, // n x 1 matrix
    d: Polynomial<I, N>,
}

impl<I, const N: usize> OpenProofVerifier<I, N>
where
    I: Integer + Signed + Sum + Roots + Clone + NumCast + SampleUniform,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    pub fn verify(
        &self,
        response: OpenProofResponse<I, N>,
        ck: &CommitmentKey<I, N>,
        params: &Params<I>,
    ) -> bool {
        if !params.check_verify_constraint(&response.z) {
            return false;
        }
        // A1 * z = t + c1 * d
        let lhs = ck.a1.dot(&response.z);
        let rhs = Mat::<I, N>::from_vec(self.t.clone()).add(&self.c1.componentwise_mul(&self.d));
        lhs == rhs
    }
}

pub struct OpenProofProver<I, const N: usize>
where
    I: Zero + One + Clone,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    opening: Opening<I, N>,
    y: Mat<I, N>, // k x 1 matrix
}

impl<I, const N: usize> OpenProofProver<I, N>
where
    I: Zero + One + Clone,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    pub fn create_response(&self, challenge: OpenProofChallenge<I, N>) -> OpenProofResponse<I, N> {
        // z = y + d * r
        let z = self.y.add(&self.opening.r.componentwise_mul(&challenge.d));
        OpenProofResponse { z }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpenProofChallenge<I, const N: usize> {
    d: Polynomial<I, N>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpenProofResponse<I, const N: usize> {
    z: Mat<I, N>, // k x 1 matrix
}
