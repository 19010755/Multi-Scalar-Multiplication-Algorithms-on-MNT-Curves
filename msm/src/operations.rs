extern crate ark_ec;
extern crate ark_ff;
extern crate ark_mnt4_298;

use ark_mnt4_298::{G1Projective, Fr};

// Point Addition - Adds two points on an MNT curve.
pub fn add_points (point1: G1Projective, point2: G1Projective) -> G1Projective {
    point1 + point2
}

// Scalar Multiplication - Multiplies a point on an MNT curve by a scalar.
pub fn scalar_multiply (point: G1Projective, scalar: Fr) -> G1Projective  {
    point * scalar
}