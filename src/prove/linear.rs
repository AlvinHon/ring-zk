use std::{
    iter::Sum,
    ops::{Add, Mul, Neg, Sub},
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

pub fn prove_linear<I, const N: usize>(
    rng: &mut impl Rng,
    g: Polynomial<I, N>,
    x: Vec<Polynomial<I, N>>,
    ck: &CommitmentKey<I, N>,
    params: &Params<I>,
) -> (LinearProofProver<I, N>, LinearProofCommitment<I, N>)
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I> + Neg<Output = I>,
{
    let Params { k, .. } = params.clone();

    let gx = x
        .iter()
        .cloned()
        .map(|xi| xi.mul(g.clone()))
        .collect::<Vec<_>>(); // g * x
    let (cp, opening_p) = ck.commit(rng, params, gx);
    let (c, opening) = ck.commit(rng, params, x);

    // y <- N^k_sigma
    let y = Mat::<I, N>::new_with(k, 1, || {
        random_polynomial_in_normal_distribution::<I, N>(
            rng,
            I::zero().to_f64().unwrap(),
            params.standard_deviation(N) as f64,
        )
    });

    // yp <- N^k_sigma
    let yp = Mat::<I, N>::new_with(k, 1, || {
        random_polynomial_in_normal_distribution::<I, N>(
            rng,
            I::zero().to_f64().unwrap(),
            params.standard_deviation(N) as f64,
        )
    });

    // t = A1 * y
    let t = ck.a1.dot(&y).one_d_mat_to_vec();

    // tp = A1 * yp
    let tp = ck.a1.dot(&yp).one_d_mat_to_vec();

    // u = g * A2 * y - A2 * yp
    let u = ck.a2.dot(&y).componentwise_mul(&g).sub(&ck.a2.dot(&yp));

    (
        LinearProofProver {
            opening,
            opening_p,
            y,
            yp,
        },
        LinearProofCommitment { c, cp, g, t, tp, u },
    )
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LinearProofCommitment<I, const N: usize>
where
    I: Integer + Signed + Sum + Roots + Clone + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    c: Commitment<I, N>,
    cp: Commitment<I, N>,
    g: Polynomial<I, N>,
    t: Vec<Polynomial<I, N>>,  // n x 1 matrix
    tp: Vec<Polynomial<I, N>>, // n x 1 matrix
    u: Mat<I, N>,              // l x 1 matrix
}

impl<I, const N: usize> LinearProofCommitment<I, N>
where
    I: Integer + Signed + Sum + Roots + Clone + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    pub fn create_challenge(
        &self,
        rng: &mut impl Rng,
        params: &Params<I>,
    ) -> (LinearProofVerifier<I, N>, LinearProofChallenge<I, N>)
    where
        I: Integer + Clone + SampleUniform,
    {
        let d = random_polynomial_from_challenge_set(rng, params.kappa);
        let (c1, c2) = self.c.c1_c2(params);
        let (c1p, c2p) = self.cp.c1_c2(params);
        (
            LinearProofVerifier {
                c1,
                c2,
                c1p,
                c2p,
                g: self.g.clone(),
                t: self.t.clone(),
                tp: self.tp.clone(),
                u: self.u.clone(),
                d: d.clone(),
            },
            LinearProofChallenge { d },
        )
    }
}

pub struct LinearProofVerifier<I, const N: usize>
where
    I: Integer + Signed + Sum + Roots + Clone + NumCast + SampleUniform,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    c1: Mat<I, N>, // n x 1 matrix
    c2: Mat<I, N>, // l x 1 matrix

    c1p: Mat<I, N>, // n x 1 matrix
    c2p: Mat<I, N>, // l x 1 matrix

    g: Polynomial<I, N>,

    t: Vec<Polynomial<I, N>>,  // n x 1 matrix
    tp: Vec<Polynomial<I, N>>, // n x 1 matrix
    u: Mat<I, N>,              // l x 1 matrix
    d: Polynomial<I, N>,
}

impl<I, const N: usize> LinearProofVerifier<I, N>
where
    I: Integer + Signed + Sum + Roots + Clone + NumCast + SampleUniform,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I> + Neg<Output = I>,
{
    pub fn verify(
        &self,
        response: LinearProofResponse<I, N>,
        ck: &CommitmentKey<I, N>,
        params: &Params<I>,
    ) -> bool {
        if !params.check_verify_constraint(&response.z) {
            return false;
        }
        if !params.check_verify_constraint(&response.zp) {
            return false;
        }
        // A1 * z = t + c1 * d
        let lhs = ck.a1.dot(&response.z);
        let rhs = Mat::<I, N>::from_vec(self.t.clone()).add(&self.c1.componentwise_mul(&self.d));
        if lhs != rhs {
            return false;
        }
        // A1 * zp = tp + c1p * d
        let lhs = ck.a1.dot(&response.zp);
        let rhs = Mat::<I, N>::from_vec(self.tp.clone()).add(&self.c1p.componentwise_mul(&self.d));
        if lhs != rhs {
            return false;
        }
        // g * A2 * z - A2 * zp = (g * c2 - c2p) * d + u
        let lhs = ck
            .a2
            .dot(&response.z)
            .componentwise_mul(&self.g)
            .sub(&ck.a2.dot(&response.zp));
        let rhs = self
            .c2
            .componentwise_mul(&self.g)
            .sub(&self.c2p)
            .componentwise_mul(&self.d)
            .add(&self.u);
        lhs == rhs
    }
}

pub struct LinearProofProver<I, const N: usize>
where
    I: Zero + One + Clone,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    opening: Opening<I, N>,
    opening_p: Opening<I, N>,
    y: Mat<I, N>,  // k x 1 matrix
    yp: Mat<I, N>, // k x 1 matrix
}

impl<I, const N: usize> LinearProofProver<I, N>
where
    I: Zero + One + Clone,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    pub fn create_response(
        &self,
        challenge: LinearProofChallenge<I, N>,
    ) -> LinearProofResponse<I, N> {
        // z = y + d * r
        let z = self.y.add(&self.opening.r.componentwise_mul(&challenge.d));
        // zp = yp + d * rp
        let zp = self
            .yp
            .add(&self.opening_p.r.componentwise_mul(&challenge.d));
        LinearProofResponse { z, zp }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LinearProofChallenge<I, const N: usize> {
    d: Polynomial<I, N>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LinearProofResponse<I, const N: usize> {
    z: Mat<I, N>,  // k x 1 matrix
    zp: Mat<I, N>, // k x 1 matrix
}
