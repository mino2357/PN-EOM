// simple vector (Vector with time).
//
// Implicit assumption: dimensions do not change.
// => No exception handling for operations with different dimensions.
//
use std::ops;

// Simple Vector
// (time, vector)
#[derive(Clone, Debug)]
pub struct SmpVector {
    pub time: f64,
    pub vec: Vec<f64>,
}

impl SmpVector {
    // initialize
    pub fn initialize() -> SmpVector {
        todo!()
    }

    // 0 fill
    pub fn set_zero_vector(dim: usize) -> SmpVector {
        SmpVector {
            time: 0.0,
            vec: vec![0.0; dim],
        }
    }

    // 1 fill
    pub fn set_one_fill(dim: usize) -> SmpVector {
        SmpVector {
            time: 0.0,
            vec: vec![1.0; dim],
        }
    }

    // Harmonic sequence vector (for testing).
    pub fn set_harmonic_vector(dim: usize) -> SmpVector {
        let mut ret = SmpVector::set_zero_vector(dim);
        for i in 0..ret.vec.len() {
            ret.vec[i] = 1.0 / (i as f64 + 1.0);
        }
        ret.time = 0.0;
        ret
    }

    // 2-norm
    #[inline(always)]
    pub fn norm(self: &SmpVector) -> f64 {
        let mut ret = 0.0;
        for i in 0..self.vec.len() {
            ret += self.vec[i] * self.vec[i];
        }
        ret.sqrt()
    }
}

// inner product
impl ops::Mul<&SmpVector> for &SmpVector {
    type Output = f64;
    #[inline(always)]
    fn mul(self, rhs: &SmpVector) -> Self::Output {
        let mut ret = 0.0;
        for i in 0..self.vec.len() {
            ret += self.vec[i] * rhs.vec[i];
        }
        ret
    }
}

// sum
impl ops::Add for &SmpVector {
    type Output = SmpVector;
    #[inline(always)]
    fn add(self, right: &SmpVector) -> Self::Output {
        let mut ret = SmpVector::set_zero_vector(self.vec.len());
        for i in 0..self.vec.len() {
            ret.vec[i] = self.vec[i] + right.vec[i]
        }
        ret.time = self.time;
        ret
    }
}

// difference
impl ops::Sub for &SmpVector {
    type Output = SmpVector;
    #[inline(always)]
    fn sub(self, right: &SmpVector) -> Self::Output {
        let mut ret = SmpVector::set_zero_vector(self.vec.len());
        for i in 0..self.vec.len() {
            ret.vec[i] = self.vec[i] - right.vec[i]
        }
        ret.time = self.time;
        ret
    }
}

// scalar product
impl ops::Mul<f64> for SmpVector {
    type Output = SmpVector;
    #[inline(always)]
    fn mul(self, rhs: f64) -> Self::Output {
        let mut ret = SmpVector::set_zero_vector(self.vec.len());
        for i in 0..ret.vec.len() {
            ret.vec[i] = self.vec[i] * rhs;
        }
        ret.time = self.time;
        ret
    }
}

// scalar product
impl ops::Mul<&SmpVector> for f64 {
    type Output = SmpVector;
    #[inline(always)]
    fn mul(self, lhs: &SmpVector) -> Self::Output {
        let mut ret = SmpVector::set_zero_vector(lhs.vec.len());
        for i in 0..ret.vec.len() {
            ret.vec[i] = self * lhs.vec[i];
        }
        ret.time = lhs.time;
        ret
    }
}

//
// Unit test
//
#[cfg(test)]
pub mod tests {
    extern crate approx;
    use crate::num_core::smp_vector::SmpVector;

    #[test]
    fn vector_add() {
        let a = SmpVector::set_one_fill(10);
        let b = SmpVector::set_one_fill(10);
        let c = 2.0 * &SmpVector::set_one_fill(10);
        approx::assert_abs_diff_eq!((&a + &b).norm(), c.norm());
    }

    #[test]
    fn vector_sub() {
        let a = 2.0 * &SmpVector::set_one_fill(10);
        let b = 3.0 * &SmpVector::set_one_fill(10);
        let c = -1.0 * &SmpVector::set_one_fill(10);
        approx::assert_abs_diff_eq!((&a - &b).norm(), c.norm(), epsilon = 1.0e-14);
    }

    #[test]
    fn vector_zero() {
        let a = 2.0 * &SmpVector::set_zero_vector(10);
        approx::assert_abs_diff_eq!(a.norm(), 0.0, epsilon = 1.0e-14);
    }

    #[test]
    fn inner_product() {
        let a = SmpVector::set_harmonic_vector(1000);
        let b = SmpVector::set_harmonic_vector(1000);
        approx::assert_abs_diff_eq!((6.0 * (&a * &b)).sqrt(), 3.14, epsilon = 1.0e-2);
    }
}
