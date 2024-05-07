use crate::operations::{add_points, scalar_multiply};
use ark_mnt4_298::G1Projective;
use ark_ff::Zero;

// Naive approach to Multi-Scalar Multiplication
pub fn naive_msm(points: &[G1Projective], scalars: &[u32]) -> G1Projective {
    // Ensure points and scalars have the same length
    assert_eq!(points.len(), scalars.len(), "Points and scalars must have the same length");

    let mut result = G1Projective::zero();
    
    for (&scalar, point) in scalars.iter().zip(points) {
        result = add_points(result, scalar_multiply(*point, scalar.into()));
    }

    // Return a single point as MSM result
    result
}
