use msm::naive::naive_msm;
use msm::operations::{add_points,scalar_multiply};
use ark_mnt4_298::G1Projective;
use ark_ff::Zero;
use ark_std::{test_rng, UniformRand};

#[test]
fn test_naive_msm_with_zero_scalars() {
    let mut rng = test_rng();
    let points = vec![G1Projective::rand(&mut rng), G1Projective::rand(&mut rng)];
    let scalars = vec![0, 0];  // Zero scalars
    assert_eq!(naive_msm(&points, &scalars), G1Projective::zero(), "MSM with zero scalars should be the zero point");
}

#[test]
fn test_naive_msm_with_mixed_scalars() {
    let mut rng = test_rng();
    let points = vec![G1Projective::rand(&mut rng), G1Projective::rand(&mut rng)];
    let scalars = vec![1, 2];  // Mixed scalars
    // Compare against result from point addition and scalar multiplication
    let expected_result = add_points(scalar_multiply(points[0], scalars[0].into()), scalar_multiply(points[1], scalars[1].into()));
    assert_eq!(naive_msm(&points, &scalars), expected_result, "MSM with mixed scalars failed");
}

#[test]
fn test_naive_msm_with_all_ones_scalars() {
    let mut rng = test_rng();
    let points = vec![G1Projective::rand(&mut rng), G1Projective::rand(&mut rng)];
    let scalars = vec![1, 1];  // All ones
    // Compare against result from point addition and scalar multiplication
    let expected_result = add_points(points[0], points[1]);
    assert_eq!(naive_msm(&points, &scalars), expected_result, "MSM with all ones scalars failed");
}

#[test]
fn test_naive_msm_with_large_scalars() {
    let mut rng = test_rng();
    let points = vec![G1Projective::rand(&mut rng)];
    let large_scalar = 1u32 << 30;  // Large scalar
    // Compare against result from point addition and scalar multiplication
    let expected_result = scalar_multiply(points[0], large_scalar.into());
    assert_eq!(naive_msm(&points, &[large_scalar]), expected_result, "MSM with large scalar failed");
}

#[test]
fn test_naive_msm_with_empty_lists() {
    let points: Vec<G1Projective> = Vec::new(); // Empty list
    let scalars: Vec<u32> = Vec::new(); // Empty list
    assert_eq!(naive_msm(&points, &scalars), G1Projective::zero(), "MSM with large scalar failed");
}

#[test]
#[should_panic(expected = "Points and scalars must have the same length")]
fn test_naive_msm_with_different_lengths() {
    let mut rng = test_rng();
    let points = vec![G1Projective::rand(&mut rng)];
    let scalars = vec![1, 2];
    let panic_result = naive_msm(&points, &scalars); // This should panic
    assert_eq!(panic_result,G1Projective::zero()) 
}



