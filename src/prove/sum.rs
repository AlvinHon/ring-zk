use std::{
    iter::Sum,
    ops::{Add, Mul, Neg, Sub},
};

use num::{integer::Roots, Integer, NumCast, Signed};
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

pub struct SumProofProver<I, const N: usize>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I> + Neg<Output = I>,
{
    pub(crate) params: Params<I>,
    pub(crate) ck: CommitmentKey<I, N>,
}

impl<I, const N: usize> SumProofProver<I, N>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I> + Neg<Output = I>,
{
    pub fn new(ck: CommitmentKey<I, N>, params: Params<I>) -> Self {
        Self { params, ck }
    }

    pub fn commit(
        &self,
        rng: &mut impl Rng,
        gs: Vec<Polynomial<I, N>>,
        xs: Vec<Vec<Polynomial<I, N>>>,
    ) -> (SumProofResponseContext<I, N>, SumProofCommitment<I, N>) {
        assert!(!gs.is_empty() && gs.len() == xs.len());
        // xp = g_0 * x_0 + g_1 * x_1 + ...
        let xp = xs
            .iter()
            .cloned()
            .map(Mat::<I, N>::from_vec)
            .zip(gs.iter())
            .map(|(x, g)| x.componentwise_mul(g))
            .reduce(|acc, x| acc.add(&x))
            .unwrap()
            .one_d_mat_to_vec();
        let (opening_p, cp) = self.ck.commit(rng, xp, &self.params);
        let (openings, cs) = xs
            .into_iter()
            .map(|x| self.ck.commit(rng, x, &self.params))
            .unzip::<_, _, Vec<_>, Vec<_>>();

        // y <- N^k_sigma for for y_i
        let ys = (0..gs.len())
            .map(|_| {
                Mat::<I, N>::new_with(self.params.k, 1, || {
                    random_polynomial_in_normal_distribution::<I, N>(
                        rng,
                        I::zero().to_f64().unwrap(),
                        self.params.standard_deviation(N) as f64,
                    )
                })
            })
            .collect::<Vec<_>>();

        // yp <- N^k_sigma
        let yp = Mat::<I, N>::new_with(self.params.k, 1, || {
            random_polynomial_in_normal_distribution::<I, N>(
                rng,
                I::zero().to_f64().unwrap(),
                self.params.standard_deviation(N) as f64,
            )
        });

        // t = A1 * y for each y_i
        let ts = ys
            .iter()
            .map(|y| self.ck.a1.dot(y).one_d_mat_to_vec())
            .collect::<Vec<_>>();

        // tp = A1 * yp
        let tp = self.ck.a1.dot(&yp).one_d_mat_to_vec();

        // u = g_0 * A2 * y_0 +  g_1 * A2 * y_1 + ... - A2 * yp
        let u = gs
            .iter()
            .zip(ys.iter())
            .map(|(g, y)| self.ck.a2.dot(y).componentwise_mul(g))
            .reduce(|acc, x| acc.add(&x))
            .unwrap()
            .sub(&self.ck.a2.dot(&yp));

        (
            SumProofResponseContext {
                openings,
                opening_p,
                yp,
                ys,
            },
            SumProofCommitment {
                cp,
                cs,
                gs,
                tp,
                ts,
                u,
            },
        )
    }

    pub fn create_response(
        &self,
        context: SumProofResponseContext<I, N>,
        challenge: SumProofChallenge<I, N>,
    ) -> SumProofResponse<I, N> {
        // z = y + d * r for each y_i
        let zs = context
            .ys
            .iter()
            .zip(context.openings.iter())
            .map(|(y, opening)| y.add(&opening.r.componentwise_mul(&challenge.d)))
            .collect::<Vec<_>>();
        // zp = yp + d * rp
        let zp = context
            .yp
            .add(&context.opening_p.r.componentwise_mul(&challenge.d));

        SumProofResponse { zs, zp }
    }
}

pub struct SumProofVerifier<I, const N: usize>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I> + Neg<Output = I>,
{
    params: Params<I>,
    ck: CommitmentKey<I, N>,
}

impl<I, const N: usize> SumProofVerifier<I, N>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I> + Neg<Output = I>,
{
    pub fn new(ck: CommitmentKey<I, N>, params: Params<I>) -> Self {
        SumProofVerifier { params, ck }
    }

    pub fn generate_challenge(
        &self,
        rng: &mut impl Rng,
        commitment: SumProofCommitment<I, N>,
    ) -> (SumProofVerificationContext<I, N>, SumProofChallenge<I, N>) {
        let d = random_polynomial_from_challenge_set(rng, self.params.kappa);
        let cs = commitment
            .cs
            .iter()
            .map(|c| c.c1_c2(&self.params))
            .collect();
        let (c1p, c2p) = commitment.cp.c1_c2(&self.params);
        (
            SumProofVerificationContext {
                c1p,
                c2p,
                cs,
                gs: commitment.gs.clone(),
                ts: commitment.ts.clone(),
                tp: commitment.tp.clone(),
                u: commitment.u.clone(),
                d: d.clone(),
            },
            SumProofChallenge { d },
        )
    }

    pub fn verify(
        &self,
        response: SumProofResponse<I, N>,
        context: SumProofVerificationContext<I, N>,
    ) -> bool {
        if !response
            .zs
            .iter()
            .all(|z| self.params.check_verify_constraint(z))
        {
            return false;
        }
        if !self.params.check_verify_constraint(&response.zp) {
            return false;
        }
        // check lengths
        if response.zs.len() != context.ts.len() && response.zs.len() != context.cs.len() {
            return false;
        }

        // A1 * z = t + c1 * d for each z_i
        let lhs = response
            .zs
            .iter()
            .map(|z| self.ck.a1.dot(z))
            .collect::<Vec<_>>();
        let rhs = context
            .cs
            .iter()
            .zip(context.ts)
            .map(|((c1, _), t)| Mat::<I, N>::from_vec(t).add(&c1.componentwise_mul(&context.d)))
            .collect::<Vec<_>>();
        if lhs != rhs {
            return false;
        }

        // A1 * zp = tp + c1p * d
        let lhs = self.ck.a1.dot(&response.zp);
        let rhs = Mat::<I, N>::from_vec(context.tp).add(&context.c1p.componentwise_mul(&context.d));
        if lhs != rhs {
            return false;
        }

        // g_0 * A2 * z_0 + g_1 * A2 * z_1 + ... - A2 * zp = (g_0 * c2_0 + g_1 * c2_1 + ... - c2p) * d + u
        let lhs = response
            .zs
            .iter()
            .zip(context.gs.iter())
            .map(|(z, g)| self.ck.a2.dot(z).componentwise_mul(g))
            .reduce(|acc, x| acc.add(&x))
            .unwrap()
            .sub(&self.ck.a2.dot(&response.zp));
        let rhs = context
            .cs
            .iter()
            .zip(context.gs.iter())
            .map(|((_, c2), g)| c2.componentwise_mul(g))
            .reduce(|acc, x| acc.add(&x))
            .unwrap()
            .sub(&context.c2p)
            .componentwise_mul(&context.d)
            .add(&context.u);
        lhs == rhs
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SumProofCommitment<I, const N: usize> {
    /// commitment to x'
    pub cp: Commitment<I, N>,
    /// commitments to x_i where x' = g_0 x_0 + g_1 x_1 + ..
    pub cs: Vec<Commitment<I, N>>,
    gs: Vec<Polynomial<I, N>>,      // vector of scalar g_i
    tp: Vec<Polynomial<I, N>>,      // n x 1 matrix
    ts: Vec<Vec<Polynomial<I, N>>>, // vector of n x 1 matrices
    u: Mat<I, N>,                   // l x 1 matrix
}

pub struct SumProofVerificationContext<I, const N: usize> {
    c1p: Mat<I, N>,                  // n x 1 matrix
    c2p: Mat<I, N>,                  // l x 1 matrix
    cs: Vec<(Mat<I, N>, Mat<I, N>)>, // vector of (n x 1, l x 1) matrices
    gs: Vec<Polynomial<I, N>>,       // vector of scalar g_i
    ts: Vec<Vec<Polynomial<I, N>>>,  // vector of n x 1 matrices
    tp: Vec<Polynomial<I, N>>,       // n x 1 matrix
    u: Mat<I, N>,                    // l x 1 matrix
    d: Polynomial<I, N>,
}

pub struct SumProofResponseContext<I, const N: usize> {
    /// vector of openings of x_i where x' = g_0 x_0 + g_1 x_1 + ..
    pub openings: Vec<Opening<I, N>>,
    /// opening of x'
    pub opening_p: Opening<I, N>,
    yp: Mat<I, N>,      // k x 1 matrix
    ys: Vec<Mat<I, N>>, // vector of k x 1 matrices
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SumProofChallenge<I, const N: usize> {
    d: Polynomial<I, N>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SumProofResponse<I, const N: usize> {
    zp: Mat<I, N>,      // k x 1 matrix
    zs: Vec<Mat<I, N>>, // vector of k x 1 matrices
}
