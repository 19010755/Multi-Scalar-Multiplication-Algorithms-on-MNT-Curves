use ark_ff::Zero;
use ark_ec::Group;
use ark_mnt4_298::G1Projective;
use crate::operations::add_points;

// Trivial approach to Multi-Scalar Multiplication using doubling and addition
pub fn trivial_msm(points: &[G1Projective], scalars: &[u32]) -> G1Projective {
    // Ensure points and scalars have the same length
    assert_eq!(points.len(), scalars.len(), "Points and scalars must have the same length");

    let mut result = G1Projective::zero();

    // Iterate over each point and scalar pair
    for (i, point) in points.iter().enumerate() {
        let mut point_contribution = G1Projective::zero();
        let mut scalar = scalars[i];

        // Start with the current point to add if the LSB of the scalar is 1
        let mut current_point = *point;

        // Process each bit from LSB to MSB
        while scalar != 0 {
            if scalar & 1 == 1 {
                // Add the current point to the point_contribution if the current bit is 1
                point_contribution = add_points(point_contribution, current_point);
            }

            // Double current point regardless of bit
            current_point = current_point.double();
            scalar >>= 1; // Shift the scalar to the right for the next bit
        }

        // Add the contribution from this point-scalar pair to the total result
        result = add_points(result, point_contribution);
    }

    result
}