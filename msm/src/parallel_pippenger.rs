use crate::operations::add_points;
use ark_ec::Group;
use ark_ff::Zero;
use ark_mnt4_298::G1Projective;
use std::collections::HashMap;
use std::thread;

// Main function for Pippenger with parallelism
pub fn parallel_pippenger(points: &[G1Projective], scalars: &[u32], window_size: usize) -> G1Projective {
    assert_eq!(points.len(), scalars.len(), "Points and scalars must have the same length");
    
    let partitions = parallel_partition_msm(scalars, window_size);
    parallel_combine_partitioned_msm(&partitions, points, window_size)
}

pub struct ParallelMsmPartition {
    pub bit_index: usize,
    pub window_values: Vec<u32>,
}

// Cloned struct for parallelism
impl Clone for ParallelMsmPartition {
    fn clone(&self) -> ParallelMsmPartition {
        ParallelMsmPartition {
            bit_index: self.bit_index,
            window_values: self.window_values.clone(),
        }
    }
}

pub fn parallel_partition_msm(scalars: &[u32], window_size: usize) -> Vec<ParallelMsmPartition> {
    let num_partitions = (32 + window_size - 1) / window_size;
    let mut partitions = Vec::new();

    for partition_index in 0..num_partitions {
        let bit_index = partition_index * window_size;
        let window_values: Vec<u32> = scalars.iter().map(|&scalar| {
            (scalar >> bit_index) & ((1 << window_size) - 1)
        }).collect();
        partitions.push(ParallelMsmPartition { bit_index, window_values });
    }

    partitions
}

pub fn parallel_compute_msm_for_partition(partition: &ParallelMsmPartition, points: &[G1Projective], window_size: usize) -> G1Projective {
    let mut buckets: HashMap<u32, Vec<usize>> = HashMap::new();
    for (index, &value) in partition.window_values.iter().enumerate() {
        if value != 0 {
        buckets.entry(value).or_insert_with(Vec::new).push(index);
        }
    }

    // Variables to store the computed MSM for this partition
    let mut msm_result = G1Projective::zero();
    let mut temp = G1Projective::zero();

    // Get the maximum scalar value (which is the number of buckets minus 1)
    let max_scalar_value = (1 << window_size) - 1;

    // Iterating over scalar values in decreasing order
    for scalar_value in (1..=max_scalar_value).rev() {
        if let Some(indexes) = buckets.get(&scalar_value) {
            // Summing up the points corresponding to the indexes in the bucket
            let sum_of_points: G1Projective = indexes.iter()
                .map(|&i| points[i])
                .fold(G1Projective::zero(), |acc, p| add_points(acc, p));

            // Add to temp regardless of whether sum_of_points is zero
            temp = add_points(temp, sum_of_points);
        }

        // Add temp to msm_result
        msm_result = add_points(msm_result, temp);
    }

    msm_result
}


pub fn parallel_combine_partitioned_msm(partitions: &[ParallelMsmPartition], points: &[G1Projective], window_size: usize) -> G1Projective {
    let mut handles = Vec::new();

    // Spawn a thread for each partition, iterate through them in reverse to ensure doubling mimics the bit scaling process accurately
    for partition in partitions.iter().rev() {
        let partition_clone = partition.clone();
        let points_clone = points.to_vec();

        let handle = thread::spawn(move || {
            let msm_result = parallel_compute_msm_for_partition(&partition_clone, &points_clone, window_size);
            msm_result
        });

        handles.push(handle);
    }

    // Collect results from each thread and combine
    let mut final_result = G1Projective::zero();
    for handle in handles {
        let partition_result = handle.join().unwrap();
        
        // Double the final result window_size times to mimic scaling by 2^bit_index
        for _ in 0..window_size {
            final_result = final_result.double();
        }

        // Adding the partition MSM to the final result
        final_result = add_points(final_result, partition_result);
    }

    final_result
}