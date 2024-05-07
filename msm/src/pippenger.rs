use crate::operations::add_points;
use ark_ec::Group;
use ark_ff::Zero;
use ark_mnt4_298::G1Projective;
use std::collections::HashMap;
// use std::time::Instant;

// Main pippenger function
pub fn pippenger(points: &[G1Projective], scalars: &[u32], window_size: usize) -> G1Projective {

    // Ensure points and scalars have the same length
    assert_eq!(points.len(), scalars.len(), "Points and scalars must have the same length");
    
    let partitions = partition_msm(scalars, window_size);
    combine_partitioned_msm(&partitions, points, window_size)
}

pub struct MsmPartition {
    pub bit_index: usize,
    pub window_values: Vec<u32>,
}

/// Step 1: Split MSM with b-bit scalars into b/c MSMs with c-bit scalars. c == window_size
pub fn partition_msm(scalars: &[u32], window_size: usize) -> Vec<MsmPartition> {
    // let start_partitioning = Instant::now();
    
    // Calculate the total number of partitions based on window size
    // (32 + window_size - 1) / window_size is used so that if 32 divides window_size, num_partitions will return the divisor
    // But if 32 does not divide window_size, then num_partitions will round up the float instead of round down, the default in Rust
    let num_partitions = (32 + window_size - 1) / window_size;

    // Vector to hold information on partitions
    let mut partitions = Vec::new();

    for partition_index in 0..num_partitions {
        // Calculate the starting bit index for the current partition
        let bit_index = partition_index * window_size;

        // Collect the corresponding window values for each scalar
        let window_values: Vec<u32> = scalars.iter().map(|&scalar| {
            // Shift the scalar to right align the desired bits and mask off the rest
            (scalar >> bit_index) & ((1 << window_size) - 1)
        }).collect();

        // Push the partition information to the list of partitions
        partitions.push(MsmPartition { bit_index, window_values });
        // let duration_partitioning = start_partitioning.elapsed();
        // println!("Partitioning took: {:?}", duration_partitioning);
    }

    // Return the list of partitions, each containing its bit index and window values
    partitions
}

pub fn compute_msm_for_partition(partition: &MsmPartition, points: &[G1Projective], window_size: usize) -> G1Projective {
    // let start_bucketing = Instant::now();
    let mut buckets: HashMap<u32, Vec<usize>> = HashMap::new();
    for (index, &value) in partition.window_values.iter().enumerate() {
        if value != 0 {
            buckets.entry(value).or_insert_with(Vec::new).push(index);
        }
    }
    // let duration_bucketing = start_bucketing.elapsed();
    // println!("Bucketing took: {:?}", duration_bucketing);

    let max_scalar_value = (1 << window_size) - 1;
    let mut msm_result = G1Projective::zero();
    let mut temp = G1Projective::zero();

    // let start_subsum = Instant::now();
    for scalar_value in (1..=max_scalar_value).rev() {
        if let Some(indexes) = buckets.get(&scalar_value) {
            let sum_of_points: G1Projective = indexes.iter()
                .map(|&i| points[i])
                .fold(G1Projective::zero(), |acc, p| add_points(acc, p));
            temp = add_points(temp, sum_of_points);
        }
        msm_result = add_points(msm_result, temp);
    }
    // let duration_subsum = start_iterating.elapsed();
    // println!("Subsum accumulation took: {:?}", duration_subsum);

    msm_result
}



// Step 3: Compute the final MSM result by combining all partitions
pub fn combine_partitioned_msm(partitions: &[MsmPartition], points: &[G1Projective], window_size: usize) -> G1Projective {
    // Variable to store the final MSM result
    let mut final_result = G1Projective::zero();

    // Iterating over each partition in reverse to ensure doubling mimics scaling accurately
    for partition in partitions.iter().rev() {
        // Computing MSM for the current partition
        let partition_msm = compute_msm_for_partition(partition, points, window_size);

        // Double the final result window_size times to mimic scaling by 2^bit_index
        for _ in 0..window_size {
            final_result = final_result.double();
        }

        // Adding the partition MSM to the final result
        final_result = add_points(final_result, partition_msm);
    }

    // Returning the final combined MSM result
    final_result
}