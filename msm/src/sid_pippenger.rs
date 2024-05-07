use crate::operations::add_points;
use ark_ec::Group;
use ark_ff::Zero;
use ark_mnt4_298::G1Projective;
use std::collections::HashMap;
use std::ops::Neg;
// use std::time::Instant;

// Main function for Pippenger with Signed Integer Decomposition Decomposition
pub fn sid_pippenger(points: &[G1Projective], scalars: &[u32], window_size: usize) -> G1Projective {
    assert_eq!(points.len(), scalars.len(), "Points and scalars must have the same length");
    
    let partitions = sid_partition_msm(scalars, window_size);
    let decomposed_partitions = sid_decompose_partitions(&partitions, window_size);
    sid_combine_partitioned_msm(&decomposed_partitions, points, window_size)
}

pub struct SidMsmPartition {
    pub bit_index: usize,
    pub window_values: Vec<u32>,
}

// New struct adjusted for signed integer decomposition
pub struct SidMsmPartitionDecomposed {
    pub bit_index: usize,
    pub window_values: Vec<i64>, 
}


pub fn sid_partition_msm(scalars: &[u32], window_size: usize) -> Vec<SidMsmPartition> {
    let num_partitions = (32 + window_size - 1) / window_size;
    let mut partitions = Vec::new();

    for partition_index in 0..num_partitions {
        let bit_index = partition_index * window_size;
        let window_values: Vec<u32> = scalars.iter().map(|&scalar| {
            (scalar >> bit_index) & ((1 << window_size) - 1)
        }).collect();
        partitions.push(SidMsmPartition { bit_index, window_values });
    }

    partitions
}

// Signed Integer Decomposition Step
pub fn sid_decompose_partitions(partitions: &[SidMsmPartition], window_size: usize) -> Vec<SidMsmPartitionDecomposed> {
    // let start_decomposition = Instant::now();
    let base = 2u32.pow(window_size as u32);
    let threshold = base / 2;

    // Initialise decomposed partitions with the same structure but empty window values
    let mut decomposed_partitions: Vec<SidMsmPartitionDecomposed> = partitions.iter()
        .map(|p| SidMsmPartitionDecomposed {
            bit_index: p.bit_index,
            window_values: vec![0i64; p.window_values.len()], // Initialize with zeros
        })
        .collect();

    // Append an extra partition for overflow handling
    let last_bit_index = partitions.last().unwrap().bit_index + window_size;
    decomposed_partitions.push(SidMsmPartitionDecomposed {
        bit_index: last_bit_index,
        window_values: vec![0i64; partitions[0].window_values.len()], // Initialise with zeros for overflow handling
    });

    // Iterate through each position of window values across all partitions
    for i in 0..partitions[0].window_values.len() {
        let mut carry = 0i64;

        for j in 0..partitions.len() {
            let window_value = partitions[j].window_values[i] as i64;
            let adjusted_value = window_value + carry;
            carry = 0; // Reset carry for the next window value

            if adjusted_value >= threshold as i64 {
                decomposed_partitions[j].window_values[i] = adjusted_value - base as i64;
                // Ensure carry is forwarded to the next partition's same position
                carry = 1;
            } else {
                decomposed_partitions[j].window_values[i] = adjusted_value;
            }
        }

        if carry > 0 {
            decomposed_partitions[partitions.len()].window_values[i] += carry;
        }
    }
    // let duration_decomposition = start_decomposition.elapsed();
    // println!("Decomposition took: {:?}", duration_decomposition);

    decomposed_partitions
}

pub fn sid_compute_msm_for_partition(partition: &SidMsmPartitionDecomposed, points: &[G1Projective], window_size: usize) -> G1Projective {
    let mut buckets: HashMap<u32, Vec<(usize, i64)>> = HashMap::new(); // Use absolute value for keys and keep sign with index for values

    // Assign points to buckets based on the absolute value while keeping track of the original value's sign
    for (index, &value) in partition.window_values.iter().enumerate() {
        if value != 0 {
            let abs_value = value.abs() as u32; // Use absolute value for bucket key
            buckets.entry(abs_value).or_insert_with(Vec::new).push((index, value));
        }
    }

    // Calculate the maximum scalar value based on the absolute values
    let max_scalar_value = 1 << (window_size - 1);

    let mut msm_result = G1Projective::zero();
    let mut temp = G1Projective::zero();

    // Iterating over scalar values in decreasing order
    for scalar_value in (1..=max_scalar_value).rev() {
        if let Some(index_sign_pairs) = buckets.get(&scalar_value) {
            let sum_of_points: G1Projective = index_sign_pairs.iter()
                .map(|&(i, sign)| {
                    let mut point = points[i];
                    if sign < 0 {
                        point = point.neg(); // Negate the point if the original value was negative
                    }
                    point
                })
                .fold(G1Projective::zero(), |acc, p| add_points(acc, p));

            temp = add_points(temp, sum_of_points);
        }

        // Add temp to msm_result after each scalar value iteration
        msm_result = add_points(msm_result, temp);
    }

    msm_result
}


pub fn sid_combine_partitioned_msm(partitions: &[SidMsmPartitionDecomposed], points: &[G1Projective], window_size: usize) -> G1Projective {
    
    let mut final_result = G1Projective::zero();
   
    // Iterating over each partition in reverse to ensure doubling mimics scaling accurately
    for partition in partitions.iter().rev() {
        let partition_msm = sid_compute_msm_for_partition(partition, points, window_size);

        // Double the final result window_size times to mimic scaling by 2^bit_index
        for _ in 0..window_size {
            final_result = final_result.double();
        }

        // Add the iteratively doubled result to the accumulated result
        final_result = add_points(final_result, partition_msm);
    }
    final_result
}