//! An utilities module for polynomial operations. It provides functions for generating random polynomials and calculating norms.

use std::iter::Sum;

use num::{integer::Roots, Integer, NumCast, Signed};
use poly_ring_xnp1::Polynomial;
use rand::{distr::uniform::SampleUniform, Rng};
use rand_distr::{Distribution, Normal};

/// Returns a random polynomial with coefficients in the range `[-bound, bound]`.
///
/// ## Safety
/// **bound** must be positive.
pub(crate) fn random_polynomial_within<I, const N: usize>(
    rng: &mut impl Rng,
    bound: I,
) -> Polynomial<I, N>
where
    I: Integer + Clone + SampleUniform,
{
    let lower = I::zero() - bound.clone();

    let mut upper = bound;
    upper.inc(); // inclusive bound

    let range = lower.clone()..upper.clone();
    let coeffs = (0..N).map(|_| rng.random_range(range.clone())).collect();

    Polynomial::new(coeffs)
}

/// Returns a random polynomial with coefficients in the normal distribution.
pub(crate) fn random_polynomial_in_normal_distribution<I, const N: usize>(
    rng: &mut impl Rng,
    mean: f64,
    std_dev: f64,
) -> Polynomial<I, N>
where
    I: Integer + Clone + NumCast,
{
    let normal = Normal::new(mean, std_dev).unwrap();
    let coeffs = normal
        .sample_iter(rng)
        .take(N)
        .map(|x| NumCast::from(x).unwrap())
        .collect();

    Polynomial::new(coeffs)
}

/// Returns the 1-norm of the polynomial. It is the sum of the absolute values of the coefficients.
#[allow(unused)]
#[inline]
pub(crate) fn norm_1<I, const N: usize>(p: &Polynomial<I, N>) -> I
where
    I: Integer + Signed + Sum + Clone,
{
    p.iter().map(|c| c.abs()).sum()
}

/// Returns the 2-norm of the polynomial. It is the square root of the sum of the squares of the coefficients.
#[inline]
pub(crate) fn norm_2<I, const N: usize>(p: &Polynomial<I, N>) -> I
where
    I: Integer + Signed + Sum + Clone + Roots,
{
    p.iter().map(|c| c.clone() * c.clone()).sum::<I>().sqrt()
}

/// Returns the infinity-norm of the polynomial. It is the maximum absolute value of the coefficients.
#[allow(unused)]
#[inline]
pub(crate) fn norm_infinity<I, const N: usize>(p: &Polynomial<I, N>) -> I
where
    I: Integer + Signed + Sum + Clone,
{
    p.iter().map(|c| c.abs()).max().unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    const N: usize = 4;

    #[test]
    fn test_random_polynomial_within() {
        let mut rng = rand::rng();
        let bound = 10;
        let p = random_polynomial_within::<_, N>(&mut rng, bound);
        for c in p.iter() {
            assert!(-bound <= *c && *c <= bound);
        }
    }

    #[test]
    fn test_norm_1() {
        let p = Polynomial::<i32, N>::new(vec![1, -2, 3, -4]);
        assert_eq!(norm_1(&p), 10);
    }

    #[test]
    fn test_norm_2() {
        let p = Polynomial::<i32, N>::new(vec![1, -2, 3, -4]);
        assert_eq!(norm_2(&p), 5);
    }

    #[test]
    fn test_norm_infinity() {
        let p = Polynomial::<i32, N>::new(vec![1, -2, 3, -4]);
        assert_eq!(norm_infinity(&p), 4);
    }

    #[test]
    fn test_random_polynomial_in_normal_distribution() {
        let mut rng = rand::rng();
        let mean = 0.0;
        let std_dev = 10.0;
        let p = random_polynomial_in_normal_distribution::<i64, N>(&mut rng, mean, std_dev);
        p.iter().for_each(|c| {
            assert!(c.abs() <= 3 * std_dev as i64); // 99.7% of the data
        });
    }
}
