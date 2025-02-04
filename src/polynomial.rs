use num::{integer::Roots, Integer, Signed};
use poly_ring_xnp1::Polynomial;
use rand::{distr::uniform::SampleUniform, Rng};
use std::iter::Sum;

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

pub(crate) fn norm_1<I, const N: usize>(p: &Polynomial<I, N>) -> I
where
    I: Integer + Signed + Sum + Clone,
{
    p.iter().map(|c| c.abs()).sum()
}

pub(crate) fn norm_2<I, const N: usize>(p: &Polynomial<I, N>) -> I
where
    I: Integer + Signed + Sum + Clone + Roots,
{
    p.iter().map(|c| c.clone() * c.clone()).sum::<I>().sqrt()
}

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
        let p = Polynomial::<_, N>::new(vec![1, -2, 3, -4]);
        assert_eq!(norm_1(&p), 10);
    }

    #[test]
    fn test_norm_2() {
        let p = Polynomial::<_, N>::new(vec![1, -2, 3, -4]);
        assert_eq!(norm_2(&p), 5);
    }

    #[test]
    fn test_norm_infinity() {
        let p = Polynomial::<_, N>::new(vec![1, -2, 3, -4]);
        assert_eq!(norm_infinity(&p), 4);
    }
}
