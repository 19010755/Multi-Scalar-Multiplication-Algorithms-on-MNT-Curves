use crate::operations::add_points;
use ark_ec::Group;
use ark_ff::Zero;
use ark_mnt4_298::G1Projective;
use std::collections::BTreeMap;
use std::ops::Neg;

// Main function for Pippenger with Signed Integer Decomposition and New Subsum Accumulation
pub fn sid_subsum_pippenger(points: &[G1Projective], scalars: &[u32], window_size: usize) -> G1Projective {
    assert_eq!(points.len(), scalars.len(), "Points and scalars must have the same length");

    let partitions = sid_subsum_partition_msm(scalars, window_size);
    let decomposed_partitions = sid_subsum_decompose_partitions(&partitions, window_size);
    sid_subsum_combine_partitioned_msm(&decomposed_partitions, points, window_size)
}

pub struct SidSubsumMsmPartition {
    pub bit_index: usize,
    pub window_values: Vec<u32>,
}

// New struct adjusted for signed integer decomposition
pub struct SidSubsumMsmPartitionDecomposed {
    pub bit_index: usize,
    pub window_values: Vec<i64>, 
}

// Step 1: Split MSM with 32-bit scalars into 32/c MSMs with c-bit scalars. c == window_size
pub fn sid_subsum_partition_msm(scalars: &[u32], window_size: usize) -> Vec<SidSubsumMsmPartition> {
    let num_partitions = (32 + window_size - 1) / window_size;
    let mut partitions = Vec::new();

    for partition_index in 0..num_partitions {
        let bit_index = partition_index * window_size;
        let window_values: Vec<u32> = scalars.iter().map(|&scalar| {
            (scalar >> bit_index) & ((1 << window_size) - 1)
        }).collect();
        partitions.push(SidSubsumMsmPartition { bit_index, window_values });
    }

    partitions
}

// Step 1.5: Decompose scalars using Signed Integer Decomposition
pub fn sid_subsum_decompose_partitions(partitions: &[SidSubsumMsmPartition], window_size: usize) -> Vec<SidSubsumMsmPartitionDecomposed> {
    let base = 2u32.pow(window_size as u32); // 2^(window_size) -> can't use this!
    let threshold = base / 2; // 2^(window_size-1) -> can't use this!

    // Initialise decomposed partitions with the same structure but empty window values
    let mut decomposed_partitions: Vec<SidSubsumMsmPartitionDecomposed> = partitions.iter()
        .map(|p| SidSubsumMsmPartitionDecomposed {
            bit_index: p.bit_index,
            window_values: vec![0i64; p.window_values.len()], // Initialise with zeros
        })
        .collect();

    // Append an extra partition for overflow handling
    let last_bit_index = partitions.last().unwrap().bit_index + window_size;
    decomposed_partitions.push(SidSubsumMsmPartitionDecomposed {
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

    decomposed_partitions
}

// Step 2: Compute MSM for each partition using subsum accumulation
pub fn sid_subsum_compute_msm_for_partition(partition: &SidSubsumMsmPartitionDecomposed, points: &[G1Projective]) -> G1Projective {
    let mut buckets: BTreeMap<u32, Vec<(usize, i64)>> = BTreeMap::new(); // Use absolute value for keys and keep sign with index for values

    // Add an empty bucket with index 0 as requirement for new subsum accumulation algorithm
    buckets.insert(0, Vec::new());

    // Assign points to buckets based on the absolute value while keeping track of the original value's sign
    for (index, &value) in partition.window_values.iter().enumerate() {
        if value != 0 {
            let abs_value = value.abs() as u32; // Use absolute value for bucket key
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

// Step 3: Compute the final MSM result by combining all partitions
pub fn sid_subsum_combine_partitioned_msm(partitions: &[SidSubsumMsmPartitionDecomposed], points: &[G1Projective], window_size: usize) -> G1Projective {
    
    let mut final_result = G1Projective::zero();
   
    // Iterating over each partition in reverse to ensure doubling mimics scaling accurately
    for partition in partitions.iter().rev() {
        let partition_msm = sid_subsum_compute_msm_for_partition(partition, points);

        // Double the final result window_size times to mimic scaling by 2^bit_index
        for _ in 0..window_size {
            final_result = final_result.double();
        }

        // Add the iteratively doubled result to the accumulated result
        final_result = add_points(final_result, partition_msm);
    }
    final_result
}