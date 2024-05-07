use crate::operations::add_points;
use ark_ec::Group;
use ark_ff::Zero;
use ark_mnt4_298::G1Projective;
use std::collections::BTreeMap;
use std::ops::Neg;
use std::thread;

// Main function for Pippenger with parallelism and Signed Integer Decomposition
pub fn parallel_sid_subsum_pippenger(points: &[G1Projective], scalars: &[u32], window_size: usize) -> G1Projective {
    assert_eq!(points.len(), scalars.len(), "Points and scalars must have the same length");
    
    let partitions = parallel_sid_subsum_partition_msm(scalars, window_size);
    let decomposed_partitions = parallel_sid_subsum_decompose_partitions(&partitions, window_size);
    parallel_sid_subsum_combine_partitioned_msm(&decomposed_partitions, points, window_size)
}

pub struct ParallelSidSubsumMsmPartition {
    pub bit_index: usize,
    pub window_values: Vec<u32>,
}

// Cloned struct for parallelism
impl Clone for ParallelSidSubsumMsmPartition {
    fn clone(&self) -> ParallelSidSubsumMsmPartition {
        ParallelSidSubsumMsmPartition {
            bit_index: self.bit_index,
            window_values: self.window_values.clone(),
        }
    }
}

// New struct adjusted for signed integer decomposition
pub struct ParallelSidSubsumMsmPartitionDecomposed {
    pub bit_index: usize,
    pub window_values: Vec<i64>,
}

// Cloned struct for parallelism
impl Clone for ParallelSidSubsumMsmPartitionDecomposed {
    fn clone(&self) -> ParallelSidSubsumMsmPartitionDecomposed {
        ParallelSidSubsumMsmPartitionDecomposed {
            bit_index: self.bit_index,
            window_values: self.window_values.clone(),
        }
    }
}

pub fn parallel_sid_subsum_partition_msm(scalars: &[u32], window_size: usize) -> Vec<ParallelSidSubsumMsmPartition> {
    let num_partitions = (32 + window_size - 1) / window_size;
    let mut partitions = Vec::new();

    for partition_index in 0..num_partitions {
        let bit_index = partition_index * window_size;
        let window_values: Vec<u32> = scalars.iter().map(|&scalar| {
            (scalar >> bit_index) & ((1 << window_size) - 1)
        }).collect();
        partitions.push(ParallelSidSubsumMsmPartition { bit_index, window_values });
    }

    partitions
}

pub fn parallel_sid_subsum_decompose_partitions(partitions: &[ParallelSidSubsumMsmPartition], window_size: usize) -> Vec<ParallelSidSubsumMsmPartitionDecomposed> {
    let base = 2u32.pow(window_size as u32);
    let threshold = base / 2;

    // Initialise decomposed partitions with the same structure but empty window values
    let mut decomposed_partitions: Vec<ParallelSidSubsumMsmPartitionDecomposed> = partitions.iter()
        .map(|p| ParallelSidSubsumMsmPartitionDecomposed {
            bit_index: p.bit_index,
            window_values: vec![0i64; p.window_values.len()], // Initialise with zeros
        })
        .collect();

    // Append an extra partition for overflow handling
    let last_bit_index = partitions.last().unwrap().bit_index + window_size;
    decomposed_partitions.push(ParallelSidSubsumMsmPartitionDecomposed {
        bit_index: last_bit_index,
        window_values: vec![0i64; partitions[0].window_values.len()], // Initialise with zeros for overflow handling
    });

    // Iterate through each position of window values across all partitions, not window_values in the same bit_index
    for i in 0..partitions[0].window_values.len() {
        let mut carry = 0i64;

        for j in 0..partitions.len() {
            let window_value = partitions[j].window_values[i] as i64;
            let adjusted_value = window_value + carry;
            carry = 0; // Reset carry for the next window value

            if adjusted_value >= threshold as i64 {
                decomposed_partitions[j].window_values[i] = adjusted_value - base as i64;
                // Forward carry to the next partition's same position
                carry = 1;
            } else {
                decomposed_partitions[j].window_values[i] = adjusted_value;
            }
        }

        if carry > 0 {
            decomposed_partitions[partitions.len()].window_values[i] += carry;
        }
    }

    decomposed_partitions
}

pub fn parallel_sid_subsum_compute_msm_for_partition(partition: &ParallelSidSubsumMsmPartitionDecomposed, points: &[G1Projective]) -> G1Projective {
    let mut buckets: BTreeMap<u32, Vec<(usize, i64)>> = BTreeMap::new();

    // Add an empty bucket with index 0 as requirement for new parallel_subsum accumulation algorithm
    buckets.insert(0, Vec::new());

    // Assign points to buckets based on the absolute value while keeping track of the original value's sign
    for (index, &value) in partition.window_values.iter().enumerate() {
        if value != 0 {
            let abs_value = value.abs() as u32;
            buckets.entry(abs_value).or_insert_with(Vec::new).push((index, value));
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
        let sum_of_points: G1Projective = indexes.iter()
                .map(|&(i, sign)| {
                    let mut point = points[i];
                    if sign < 0 {
                        point = point.neg(); // Negate the point if the original value was negative
                    }
                    point
                })
                .fold(G1Projective::zero(), |acc, p| add_points(acc, p));

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

    // Subsum accumulation on tmp array
    for i in (1..=max_diff).rev() {
        temp = add_points(temp, tmp[i]);
        msm_result = add_points(msm_result, temp);
    }

    msm_result
}

pub fn parallel_sid_subsum_combine_partitioned_msm(partitions: &[ParallelSidSubsumMsmPartitionDecomposed], points: &[G1Projective], window_size: usize) -> G1Projective {
    let mut handles = Vec::new();

    // Spawn a thread for each partition, iterate through them in reverse to ensure doubling mimics the bit scaling process accurately
    for partition in partitions.iter().rev() {
        let partition_clone = partition.clone();
        let points_clone = points.to_vec();

        let handle = thread::spawn(move || {
            let msm_result = parallel_sid_subsum_compute_msm_for_partition(&partition_clone, &points_clone);
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