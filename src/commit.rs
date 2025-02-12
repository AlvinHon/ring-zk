//! Definition of the commitment scheme, defined in section 4.1 of the paper.

use std::{
    iter::Sum,
    ops::{Add, Mul, Sub},
};

use num::{integer::Roots, Integer, NumCast, One, Signed, Zero};
use poly_ring_xnp1::Polynomial;
use rand::{distr::uniform::SampleUniform, Rng};

use crate::{mat::Mat, params::Params, polynomial::random_polynomial_within};

/// The commitment key for the commitment scheme. It is used by both the prover and the verifier.
/// The prover uses it to commit to the message while the verifier uses it to verify the commitment.
/// They should use the same commitment key in the protocol.
///
/// The size of the commitment key contains (n + l) x k polynomials, where n, k, and l are the parameters
/// defined in the `Params` struct.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CommitmentKey<I, const N: usize> {
    pub(crate) a1: Mat<I, N>, // n x k matrix
    pub(crate) a2: Mat<I, N>, // l x k matrix
}

impl<I, const N: usize> CommitmentKey<I, N>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    /// Generate a random new commitment key given the parameters.
    pub(crate) fn new(rng: &mut impl Rng, params: &Params<I>) -> Self {
        let Params { q, n, k, l, .. } = params.clone();
        // Defined in equation (5) of the paper:
        // a1 = [I_n a1'], where a1 is a polynomial matrix of size n x (k-n)
        // a1 is a polynomial matrix of size n x k
        let a1 = {
            let mut tmp = Mat::<I, N>::diag(n, n, Polynomial::<I, N>::one());
            let a1_prime =
                Mat::<I, N>::new_with(n, k - n, || random_polynomial_within(rng, q.clone()));
            tmp.extend_cols(a1_prime);
            tmp
        };

        // Defined in equation (6) of the paper:
        // a2 = [0_lxn I_l a2], where a2 is a polynomial matrix of size l x (k-n-l)
        // a2 is a polynomial matrix of size l x k
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

    /// Commit to the message `x` using the commitment key. It returns the opening and the commitment.
    ///
    /// ## Example
    ///
    /// ```rust
    /// use ring_zk::Params;
    ///
    /// const N: usize = 4; // Must be a power of two
    ///
    /// let rng = &mut rand::rng();
    /// let params = Params::default();
    /// let ck = params.generate_commitment_key::<N>(rng);
    ///
    /// let x = params.prepare_value(vec![vec![1, 2, 3, 4]]);
    /// let (open, com) = ck.commit(rng, x, &params);
    /// assert!(com.verify(&open, &ck, &params));
    /// ```
    ///
    /// ## Safety
    /// This method contains a loop that generates a random polynomial `r` until the commitment constraint
    /// defined in the `Params` struct is satisfied. This check is to ensure the comitment can be verified
    /// correctly.
    /// So it is important to ensure that the parameters are set carefully to avoid an infinite loop.
    ///
    /// ## Panics
    /// Panics if the length of `x` is not equal to the length of `l` defined in the `Params` struct.
    pub fn commit(
        &self,
        rng: &mut impl Rng,
        x: Vec<Polynomial<I, N>>,
        params: &Params<I>,
    ) -> (Opening<I, N>, Commitment<I, N>) {
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

        // Defined in equation (7) of the paper:
        // [c1 c2] = [a1 a2] * r + [0_n x]
        let c = a.dot(&r).add(&z);

        (Opening { x, r, f: None }, Commitment { c })
    }
}

/// The commitment in the commitment scheme.
///
/// The size of the commitment contains (n + l) x 1 polynomials, where n and l are the parameters
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Commitment<I, const N: usize> {
    /// The commitment value [c1 c2]. A combined matrix is used here for compactness.
    pub(crate) c: Mat<I, N>,
}

impl<I, const N: usize> Commitment<I, N>
where
    I: Integer + Signed + Sum + Roots + Clone + SampleUniform + NumCast,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    /// Verify the validity of opening r.s.t the commitment.
    ///
    /// ## Example
    ///
    /// ```rust
    /// use ring_zk::Params;
    ///
    /// const N: usize = 4; // Must be a power of two
    ///
    /// let rng = &mut rand::rng();
    /// let params = Params::default();
    /// let ck = params.generate_commitment_key::<N>(rng);
    ///
    /// let x = params.prepare_value(vec![vec![1, 2, 3, 4]]);
    /// let (open, com) = ck.commit(rng, x, &params);
    /// assert!(com.verify(&open, &ck, &params));
    ///
    /// let x2 = params.prepare_value(vec![vec![4, 5, 6, 7]]);
    /// let (open2, com2) = ck.commit(rng, x2, &params);
    /// assert!(com2.verify(&open2, &ck, &params));
    ///
    /// assert!(!com2.verify(&open, &ck, &params));
    /// assert!(!com.verify(&open2, &ck, &params));
    /// ```
    ///
    pub fn verify(
        &self,
        opening: &Opening<I, N>,
        ck: &CommitmentKey<I, N>,
        params: &Params<I>,
    ) -> bool {
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

        // Defined in the method `Open` in section 4.1 of the paper:
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

/// The opening in the commitment scheme.
#[derive(Clone, Debug, PartialEq, Eq)]
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
