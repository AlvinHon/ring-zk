use num::Integer;
use poly_ring_xnp1::Polynomial;
use rand::{distributions::uniform::SampleUniform, seq::SliceRandom, Rng};
use std::ops::{Add, Mul, Neg, Sub};

/// Create a random polynomial in Challenge Space C.
/// The Challenge Space C is defined as  `{c in R_q | norm_infinity(c) = 1, norm_1(c) = kappa}`.
/// In other words, there exists exactly `kappa` amount of coefficients = 1 or -1, and the rest are 0.
pub(crate) fn rand_polynomial_from_challenge_set<I, const N: usize>(
    rng: &mut impl Rng,
    kappa: usize,
) -> Polynomial<I, N>
where
    I: Integer + Clone + SampleUniform,
{
    let zero = I::zero();
    let bound = I::one(); // since norm_infinity(c) = 1, the coefficients are in the range [-1, 1]
    let mut coeffs = (0..N).map(|_| zero.clone()).collect::<Vec<_>>();
    coeffs.iter_mut().take(kappa).for_each(|c| {
        // either 1 or -1
        *c = if rng.gen_bool(0.5) {
            bound.clone()
        } else {
            zero.clone() - bound.clone()
        };
    });
    coeffs.shuffle(rng);
    Polynomial::new(coeffs)
}

/// Create a random polynomial in Set difference (C-bar) in Challenge Space C.
/// Defines C-bar as `{c - c', where c, c' in C}`. This difference `c - c'` has
/// a special property that the returned polynomial is invertible in `R_q`.
pub(crate) fn rand_polynomial_from_challenge_set_difference<I, const N: usize>(
    rng: &mut impl Rng,
    kappa: usize,
) -> Polynomial<I, N>
where
    I: Integer + Clone + SampleUniform,
    for<'a> &'a I: Neg<Output = I> + Mul<Output = I> + Sub<Output = I> + Add<Output = I>,
{
    let c1 = rand_polynomial_from_challenge_set(rng, kappa);
    loop {
        let c2 = rand_polynomial_from_challenge_set(rng, kappa);
        if c1 != c2 {
            return c1 - c2;
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::polynomial::{norm_1, norm_infinity};

    use super::*;

    const N: usize = 256;

    #[test]
    fn test_rand_polynomial_from_challenge_set() {
        let mut rng = rand::thread_rng();
        let kappa = 60;
        let c = rand_polynomial_from_challenge_set::<i32, N>(&mut rng, kappa);
        assert_eq!(norm_1(&c), kappa as i32);
        assert_eq!(norm_infinity(&c), 1i32);
    }

    #[test]
    fn test_rand_polynomial_from_challenge_set_difference() {
        let mut rng = rand::thread_rng();
        let kappa = 60;
        let c = rand_polynomial_from_challenge_set_difference::<i32, N>(&mut rng, kappa);
        c.iter().for_each(|c| {
            assert!(c >= &-2i32); // the coefficients are in the range [-1, 1] so the possible max is 2
        });
    }
}
