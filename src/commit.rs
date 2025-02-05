use std::{
    iter::Sum,
    ops::{Add, Mul, Sub},
};

use num::{integer::Roots, Integer, NumCast, One, Signed, Zero};
use poly_ring_xnp1::Polynomial;
use rand::{distr::uniform::SampleUniform, Rng};

use crate::{mat::Mat, params::Params, polynomial::random_polynomial_within};

pub struct CommitmentKey<I, const N: usize>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    pub(crate) a1: Mat<I, N>, // n x k matrix
    pub(crate) a2: Mat<I, N>, // l x k matrix
}

impl<I, const N: usize> CommitmentKey<I, N>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    pub fn new(rng: &mut impl Rng, params: &Params<I>) -> Self {
        let Params { q, n, k, l, .. } = params.clone();

        // a1 is a polynomial matrix of size n x k
        // a1 = [I_n a1'], where a1 is a polynomial matrix of size n x (k-n)
        let a1 = {
            let mut tmp = Mat::<I, N>::diag(n, n, Polynomial::<I, N>::one());
            let a1_prime =
                Mat::<I, N>::new_with(n, k - n, || random_polynomial_within(rng, q.clone()));
            tmp.extend_cols(a1_prime);
            tmp
        };

        // a2 is a polynomial matrix of size l x k
        // a2 = [0_lxn I_l a2], where a2 is a polynomial matrix of size l x (k-n-l)
        let a2 = {
            let mut tmp = Mat::<I, N>::from_element(l, n, Polynomial::<I, N>::zero());
            let i_l = Mat::<I, N>::diag(l, l, Polynomial::<I, N>::one());
            let a2_prime =
                Mat::<I, N>::new_with(l, k - n - l, || random_polynomial_within(rng, q.clone()));
            tmp.extend_cols(i_l);
            tmp.extend_cols(a2_prime);
            tmp
        };

        CommitmentKey { a1, a2 }
    }

    pub fn commit(
        &self,
        rng: &mut impl Rng,
        params: &Params<I>,
        x: Vec<Polynomial<I, N>>,
    ) -> (Commitment<I, N>, Opening<I, N>) {
        let Params { b, n, k, l, .. } = params.clone();
        assert_eq!(l, x.len());

        let x_mat = Mat::<I, N>::from_vec(x.clone());
        let r = {
            let mut tmp;
            loop {
                tmp = Mat::<I, N>::new_with(k, 1, || random_polynomial_within(rng, b.clone()));
                if params.check_commit_constraint(&tmp) {
                    break;
                }
            }
            tmp
        };

        let a = {
            // [a1 a2]
            let mut a1 = self.a1.clone();
            a1.extend_rows(self.a2.clone());
            a1
        };

        let z = {
            // [0_n x]
            let mut tmp = Mat::<I, N>::from_element(n, 1, Polynomial::<I, N>::zero());
            tmp.extend_rows(x_mat);
            tmp
        };

        // [c1 c2] = [a1 a2] * r + [0_n x]
        let c = a.dot(&r).add(&z);

        (Commitment { c }, Opening { x, r, f: None })
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Commitment<I, const N: usize> {
    pub(crate) c: Mat<I, N>,
}

impl<I, const N: usize> Commitment<I, N> {
    pub fn verify(
        &self,
        params: &Params<I>,
        ck: &CommitmentKey<I, N>,
        opening: &Opening<I, N>,
    ) -> bool
    where
        I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
        for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
    {
        let Params { n, .. } = params.clone();
        let Opening { x, r, f } = opening;

        if !params.check_commit_constraint(r) {
            return false;
        }

        let a = {
            // [a1 a2]
            let mut a1 = ck.a1.clone();
            a1.extend_rows(ck.a2.clone());
            a1
        };

        let z = {
            // [0_n x]
            let mut tmp = Mat::<I, N>::from_element(n, 1, Polynomial::<I, N>::zero());
            tmp.extend_rows(Mat::<I, N>::from_vec(x.clone()));
            tmp
        };

        // f * [c1 c2] = [a1 a2] * r + f * [0_n x]
        match f {
            Some(f) => {
                let lhs = self.c.componentwise_mul(f);
                let rhs = a.dot(r).add(&z.componentwise_mul(f));
                lhs == rhs
            }
            None => a.dot(r).add(&z) == self.c,
        }
    }

    /// Split the commitment into two parts: c1 (dim: n x 1) and c2 (dim: l x 1).
    pub(crate) fn c1_c2(&self, params: &Params<I>) -> (Mat<I, N>, Mat<I, N>)
    where
        I: Clone,
    {
        self.c.clone().split_rows(params.n)
    }
}

pub struct Opening<I, const N: usize> {
    /// The committed value.
    pub(crate) x: Vec<Polynomial<I, N>>,
    /// The randomness used in the commit method.
    pub(crate) r: Mat<I, N>,
    /// Additional randomness for randomizing the opening.
    /// The polynomial must be **non-zero** in Challenge Space.
    /// None means the `f` is the `identity`` for verification.
    pub(crate) f: Option<Polynomial<I, N>>,
}

#[cfg(test)]
mod tests {
    use crate::params::params_1;

    use super::*;

    const N: usize = 4;

    #[test]
    fn test_commitment() {
        let mut rng = rand::rng();
        let params = params_1();
        let ck = CommitmentKey::new(&mut rng, &params);

        let x = vec![Polynomial::<i64, N>::from_coeffs(vec![1, 2, 3, 4])];

        let (com, open) = ck.commit(&mut rng, &params, x.clone());
        assert!(com.verify(&params, &ck, &open));
    }
}
