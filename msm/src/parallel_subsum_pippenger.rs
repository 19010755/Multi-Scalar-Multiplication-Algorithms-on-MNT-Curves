use crate::operations::add_points;
use ark_ec::Group;
use ark_ff::Zero;
use ark_mnt4_298::G1Projective;
use std::collections::BTreeMap;
use std::thread;

// Main pippenger function
pub fn parallel_subsum_pippenger(points: &[G1Projective], scalars: &[u32], window_size: usize) -> G1Projective {

    // Ensure points and scalars have the same length
    assert_eq!(points.len(), scalars.len(), "Points and scalars must have the same length");
    
    let partitions = parallel_subsum_partition_msm(scalars, window_size);
    parallel_subsum_combine_partitioned_msm(&partitions, points, window_size)
}

pub struct ParallelSubsumMsmPartition {
    pub bit_index: usize,
    pub window_values: Vec<u32>,
}

// Cloned struct for parallelism
impl Clone for ParallelSubsumMsmPartition {
    fn clone(&self) -> ParallelSubsumMsmPartition {
        ParallelSubsumMsmPartition {
            bit_index: self.bit_index,
            window_values: self.window_values.clone(),
        }
    }
}

// Step 1: Split MSM with 32-bit scalars into 32/c MSMs with c-bit scalars. c == window_size
pub fn parallel_subsum_partition_msm(scalars: &[u32], window_size: usize) -> Vec<ParallelSubsumMsmPartition> {
    
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
            (scalar >> bit_index) & ((1 << window_size) - 1)
        }).collect();

        // Push the partition information to the list of partitions
        partitions.push(ParallelSubsumMsmPartition { bit_index, window_values });
    }

    // Return the list of partitions, each containing its bit index and window values
    partitions
}

// Step 2: Compute MSM for each partition using parallel_subsum accumulation
pub fn parallel_subsum_compute_msm_for_partition(partition: &ParallelSubsumMsmPartition, points: &[G1Projective]) -> G1Projective {
    let mut buckets: BTreeMap<u32, Vec<usize>> = BTreeMap::new();

    // Add an empty bucket with index 0 as requirement for new parallel_subsum accumulation algorithm
    buckets.insert(0, Vec::new());

    // Populate buckets with indexes grouped by their scalar value
    for (index, &value) in partition.window_values.iter().enumerate() {
        if value != 0 {
            buckets.entry(value).or_insert_with(Vec::new).push(index);
        }
    }

    // Collect all the scalar values (keys of the buckets) into a vector
    let scalars: Vec<u32> = buckets.keys().copied().collect();

    // Find the maximum difference between consecutive scalars
    let max_diff = scalars.windows(2)
        .map(|window| window[1] - window[0])
        .max()
        .unwrap_or(1) as usize;
    
    // Initialise tmp array of length max_diff + 1
    let mut tmp = vec![G1Projective::zero(); max_diff + 1];

    // Use a peekable iterator to keep track of the next_scalar logic
    let mut iter = buckets.iter().rev().peekable();

    // Iterate through the sorted buckets in reverse order
    while let Some((&scalar, indexes)) = iter.next() {
        let sum_of_points = indexes.iter()
            .map(|&index| points[index])
            .fold(G1Projective::zero(), |acc, point| add_points(acc, point));

        tmp[0] = add_points(tmp[0], sum_of_points);

        // Peek to see if there is a next scalar and calculate the difference if so
        if let Some((&next_scalar_val, _)) = iter.peek() {
            let k = (scalar - next_scalar_val) as usize;
            if k >= 1 && k < tmp.len() {
                // Add the current sum to tmp[k] based on the gap to the next scalar
                tmp[k] = add_points(tmp[k], tmp[0]);
            }
        }
    }

    let mut temp = G1Projective::zero();
    let mut msm_result = G1Projective::zero();

    // ParallelSubsum accumulation on tmp array
    for i in (1..=max_diff).rev() {
        temp = add_points(temp, tmp[i]);
        msm_result = add_points(msm_result, temp);
    }

    msm_result
}

// Step 3: Compute the final MSM result by combining all partitions
pub fn parallel_subsum_combine_partitioned_msm(partitions: &[ParallelSubsumMsmPartition], points: &[G1Projective], window_size: usize) -> G1Projective {
    let mut handles = Vec::new();

    // Spawn a thread for each partition, iterate through them in reverse to ensure doubling mimics the bit scaling process accurately
    for partition in partitions.iter().rev() {
        let partition_clone = partition.clone();
        let points_clone = points.to_vec();

        let handle = thread::spawn(move || {
            let msm_result = parallel_subsum_compute_msm_for_partition(&partition_clone, &points_clone);
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