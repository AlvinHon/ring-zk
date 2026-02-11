//! Defines the Challenge Space C and its difference C-bar, defined in section 2.2 of the paper.

use std::ops::{Add, Mul, Neg, Sub};

use num::{One, Zero};
use poly_ring_xnp1::Polynomial;
use rand::{distr::uniform::SampleUniform, seq::SliceRandom, RngExt};

/// Create a random polynomial in Challenge Space C.
/// The Challenge Space C is defined as  `{c in R_q | norm_infinity(c) = 1, norm_1(c) = kappa}`.
/// In other words, there exists exactly `kappa` amount of coefficients = 1 or -1, and the rest are 0.
pub(crate) fn random_polynomial_from_challenge_set<I, const N: usize>(
    rng: &mut impl RngExt,
    kappa: usize,
) -> Polynomial<I, N>
where
    I: Clone + One + Zero + SampleUniform,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Sub<Output = I>,
{
    let zero = I::zero();
    let bound = I::one(); // since norm_infinity(c) = 1, the coefficients are in the range [-1, 1]
    let mut coeffs = (0..N).map(|_| zero.clone()).collect::<Vec<_>>();
    coeffs.iter_mut().take(kappa).for_each(|c| {
        // either 1 or -1
        *c = if rng.random_bool(0.5) {
            bound.clone()
        } else {
            &zero - &bound
        };
    });
    coeffs.shuffle(rng);
    Polynomial::new(coeffs)
}

/// Create a random polynomial in Set difference (C-bar) in Challenge Space C.
/// Defines C-bar as `{c - c', where c, c' in C}`. This difference `c - c'` has
/// a special property that the returned polynomial is invertible in `R_q`.
#[allow(unused)]
pub(crate) fn random_polynomial_from_challenge_set_difference<I, const N: usize>(
    rng: &mut impl RngExt,
    kappa: usize,
) -> Polynomial<I, N>
where
    I: Clone + One + Zero + SampleUniform + PartialEq,
    for<'a> &'a I: Add<Output = I> + Mul<Output = I> + Neg<Output = I> + Sub<Output = I>,
{
    let c1 = random_polynomial_from_challenge_set(rng, kappa);
    loop {
        let c2 = random_polynomial_from_challenge_set(rng, kappa);
        if c1 != c2 {
            return c1 - c2;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomial::{norm_1, norm_infinity};
    use num::ToPrimitive;

    const N: usize = 256;

    #[test]
    fn test_random_polynomial_from_challenge_set() {
        let mut rng = rand::rng();
        let kappa = 60;
        let c = random_polynomial_from_challenge_set::<i32, N>(&mut rng, kappa);
        assert_eq!(norm_1(&c).to_usize().unwrap(), kappa);
        assert_eq!(norm_infinity(&c).to_usize().unwrap(), 1);
    }

    #[test]
    fn test_random_polynomial_from_challenge_set_difference() {
        let mut rng = rand::rng();
        let kappa = 60;
        let c = random_polynomial_from_challenge_set_difference::<i32, N>(&mut rng, kappa);
        c.iter().for_each(|c| {
            assert!(c >= &-2i32); // the coefficients are in the range [-1, 1] so the possible max is 2
        });
    }
}
