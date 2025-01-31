use num::{cast::AsPrimitive, integer::Roots, Integer};

#[derive(Clone)]
pub struct Params<I> {
    /// Prime modulus. q = 2*d + 1 (mod 4d) in Lemma 1. Use d = 2 for this library.
    pub(crate) q: I,
    /// Norm bound for honest prover's randomness.
    pub(crate) b: I,

    // k > n >= l
    /// Height of the commitment matrix (a1).
    pub(crate) n: usize,
    /// Width of the commitment matrices
    pub(crate) k: usize,
    /// Dimension of the message space.
    pub(crate) l: usize,

    /// The maximum norm_1 of any element in Challenge Space C.
    pub(crate) kappa: usize,
}

impl<I> Params<I>
where
    I: Integer + Clone + AsPrimitive<usize>,
{
    /// The standard deviation used in the zero-knowledge proof.
    pub fn standard_deviation(&self, deg_n: usize) -> usize {
        // sigma = 11 * kappa * b * sqrt(k*deg_n)
        self.b.as_() * (11 * self.kappa) * (self.k * deg_n).sqrt()
    }
}

/// Parameters Set 1. The message length (l) is 1.
pub fn params_1() -> Params<i64> {
    Params {
        q: 2147483647i64,
        b: 1,
        n: 1,
        k: 3,
        l: 1,
        kappa: 36,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_deviation() {
        let params = params_1();
        let deg_n = 1024;
        let sigma = params.standard_deviation(deg_n);
        assert_eq!(sigma, 21780);
    }
}
