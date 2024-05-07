use msm::operations::{add_points, scalar_multiply};
use ark_ff::{Zero, One};
use ark_mnt4_298::{G1Projective, Fr};
use ark_std::{UniformRand, test_rng};

#[test]
// Point added with zero point should return same point
fn test_point_addition_with_zero() {
    let mut rng = test_rng();
    let point1 = G1Projective::rand(&mut rng);
    assert_eq!(point1, add_points(point1, G1Projective::zero()), "Point addition with zero point failed");
}

#[test]
// Point addition should be commutative
fn test_point_addition_commutativity() {
    let mut rng = test_rng();
    let point1 = G1Projective::rand(&mut rng);
    let point2 = G1Projective::rand(&mut rng);
    assert_eq!(add_points(point1, point2), add_points(point2, point1), "Point addition is not commutative");
}

#[test]
// Point addition should be associative
fn test_point_addition_associativity() {
    let mut rng = test_rng();
    let point1 = G1Projective::rand(&mut rng);
    let point2 = G1Projective::rand(&mut rng);
    let point3 = G1Projective::rand(&mut rng);
    assert_eq!(add_points(add_points(point1, point2), point3), add_points(point1, add_points(point2, point3)), "Point addition is not associative");
}

#[test]
// Point added with its negative should return the zero point
fn test_point_addition_invertibility() {
    let mut rng = test_rng();
    let point = G1Projective::rand(&mut rng);
    let neg_point = -point;
    assert_eq!(G1Projective::zero(), add_points(point, neg_point), "Point addition does not satisfy invertability");
}

#[test]
// Point multiplied by one should return same point
fn test_scalar_multiplication_identity() {
    let mut rng = test_rng();
    let point1 = G1Projective::rand(&mut rng);
    assert_eq!(point1, scalar_multiply(point1, Fr::one()), "Scalar multiplication by one failed");
}

#[test]
// Scalar multiplication should be commutative
fn test_scalar_multiplication_commutativity() {
    let mut rng = test_rng();
    let point = G1Projective::rand(&mut rng);
    let scalar1 = Fr::rand(&mut rng);
    let scalar2 = Fr::rand(&mut rng);
    assert_eq!(scalar_multiply(scalar_multiply(point, scalar1), scalar2), scalar_multiply(scalar_multiply(point, scalar2), scalar1), "Scalar multiplication is not commutative");
}

#[test]
// Scalar multiplication should be associative
fn test_scalar_multiplication_associativity() {
    let mut rng = test_rng();
    let point = G1Projective::rand(&mut rng);
    let scalar1 = Fr::rand(&mut rng);
    let scalar2 = Fr::rand(&mut rng);
    assert_eq!(scalar_multiply(scalar_multiply(point, scalar1), scalar2), scalar_multiply(point, scalar1 * scalar2), "Scalar multiplication is not associative");
}

#[test]
// Point multiplied by 0 should return the zero point
fn test_scalar_multiplication_zero() {
    let mut rng = test_rng();
    let point = G1Projective::rand(&mut rng);
    assert_eq!(G1Projective::zero(), scalar_multiply(point, Fr::zero()), "Multiplication by zero failed");
}

#[test]
// Scalar multiplication should be distributive
fn test_scalar_multiplication_distributivity() {
    let mut rng = test_rng();
    let point1 = G1Projective::rand(&mut rng);
    let point2 = G1Projective::rand(&mut rng);
    let scalar = Fr::rand(&mut rng);
    assert_eq!(scalar_multiply(add_points(point1, point2), scalar), add_points(scalar_multiply(point1, scalar), scalar_multiply(point2, scalar)), "Scalar multiplication is not distributive over point addition");
}

#[test]
fn test_scalar_multiplication_large_scalar() {
    let mut rng = test_rng();
    let point = G1Projective::rand(&mut rng);
    let large_scalar = Fr::from(1 << 30); // Example large scalar
    // Simply testing that the operation completes without error
    let _ = scalar_multiply(point, large_scalar);
}
