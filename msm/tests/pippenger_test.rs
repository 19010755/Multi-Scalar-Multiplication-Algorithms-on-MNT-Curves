use msm::pippenger::{pippenger, MsmPartition, partition_msm, compute_msm_for_partition, combine_partitioned_msm};
use msm::naive::naive_msm;
use msm::operations::add_points;
use ark_mnt4_298::G1Projective;
use ark_ff::Zero;
use ark_std::{test_rng, UniformRand};
use rand::{Rng, thread_rng};

#[test]
fn test_pippenger_with_zero_scalars() {
    let mut rng = test_rng();
    let points = vec![G1Projective::rand(&mut rng), G1Projective::rand(&mut rng)];
    let scalars = vec![0, 0];
    let window_size = 2; // Example window size
    assert_eq!(pippenger(&points, &scalars, window_size), G1Projective::zero(), "Pippenger with zero scalars should return the zero point");
}

#[test]
fn test_pippenger_with_all_ones_scalars() {
    let mut rng = test_rng();
    let points = vec![G1Projective::rand(&mut rng), G1Projective::rand(&mut rng)];
    let scalars = vec![1, 1];
    let window_size = 2;
    // Compare against result from naive msm
    let expected_result = naive_msm(&points, &scalars);
    assert_eq!(pippenger(&points, &scalars, window_size), expected_result, "Pippenger with all ones scalars failed");
}

#[test]
fn test_pippenger_with_large_scalars() {
    let mut rng = test_rng();
    let points = vec![G1Projective::rand(&mut rng)];
    let large_scalar = 1u32 << 30;
    let window_size = 2;
    // Compare against result from naive msm
    let expected_result = naive_msm(&points, &[large_scalar]);
    assert_eq!(pippenger(&points, &[large_scalar], window_size), expected_result, "Pippenger with large scalar failed");
}

#[test]
fn test_pippenger_with_empty_lists() {
    let points: Vec<G1Projective> = Vec::new();
    let scalars: Vec<u32> = Vec::new();
    let window_size = 2;
    assert_eq!(pippenger(&points, &scalars, window_size), G1Projective::zero(), "Pippenger with empty lists should return the zero point");
}

#[test]
#[should_panic(expected = "Points and scalars must have the same length")]
fn test_pippenger_with_different_lengths() {
    let mut rng = test_rng();
    let points = vec![G1Projective::rand(&mut rng)];
    let scalars = vec![1, 2];
    let window_size = 2;
    let panic_result = pippenger(&points, &scalars, window_size); // This should panic
    assert_eq!(panic_result, G1Projective::zero())
}

// Helper function to generate n points
fn generate_points(num_points: usize) -> Vec<G1Projective> {
    let mut rng = test_rng();
    (0..num_points).map(|_| G1Projective::rand(&mut rng)).collect()
}

// Helper function to generate n random scalars of type u32
fn generate_scalars(num_scalars: usize) -> Vec<u32> {
    let mut rng = thread_rng();
    (0..num_scalars).map(|_| rng.gen_range(0..65536)).collect()
}

#[test]
// Test for Step 1 Part 1: Number of partitions should be 32/c. c == window_size
fn test_partition_msm() {
    let scalars = vec![182, 255, 129];
    let window_size = 2;
    let partitions = partition_msm(&scalars, window_size);

    assert_eq!(partitions.len(), 32/window_size, "Incorrect number of partitions");
}

#[test]
// Test for Step 1 Part 1: Testing behaviour of partition_msm when partitions does not divide window_size
fn test_partition_msm_1() {
    let scalars = vec![182, 255, 129];
    let window_size = 3;
    let partitions = partition_msm(&scalars, window_size);

    // The number of partitions should be (32/window_size).ceil() i.e. 11 when scalars are u32 bit integeres and window_size = 3
    assert_eq!(partitions.len(), 11, "Incorrect number of partitions");
}

#[test]
// Test for Step 1 Part 2: Bit_index should be (i * window_size)
fn test_partition_msm_2() {
    let scalars = vec![182, 255, 129]; 
    let window_size = 2;
    let partitions = partition_msm(&scalars, window_size);

    for (i, partition) in partitions.iter().enumerate() {
        assert_eq!(partition.bit_index, i * window_size, "Incorrect bit_index");
    }
}

#[test]
// Test for Step 1 Part 3: Testing correctness of window_values
fn test_partition_msm_3() {
    let scalars = vec![182, 255, 129];
    let window_size = 2;
    let partitions = partition_msm(&scalars, window_size);

    // Manually calculated expected window values for each scalar
    let expected_values: Vec<Vec<u32>> = vec![
        vec![2, 1, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  // For 182 (10110110)
        vec![3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  // For 255 (11111111)
        vec![1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  // For 129 (10000001)
    ];

    for (i, ..) in partitions.iter().enumerate() {
        for (j, &scalar) in scalars.iter().enumerate() {
            let expected = expected_values[j][i];
            println!("Exp: {:?}", expected);
            let actual = (scalar >> (i * window_size)) & ((1 << window_size) - 1);
            println!("Actual: {:?}", actual);
            assert_eq!(actual, expected, "Incorrect window value for scalar {} in partition {}", scalar, i);
        }
    }
}

#[test]
// Test for Step 2: Compute MSM for each partition
fn test_compute_msm_for_partition() {
    let points = generate_points(10);
    let partition = MsmPartition { bit_index: 0, window_values: vec![1, 0, 1, 0, 1, 0, 1, 0, 1, 0] };
    let window_size = 2;

    let msm_result = compute_msm_for_partition(&partition, &points, window_size);
    // Compare against result by adding points
    let expected_result = points.iter().step_by(2).fold(G1Projective::zero(), |acc, &p| add_points(acc, p));
    assert_eq!(msm_result, expected_result, "MSM computation for partition failed");
}

#[test]
// Test for Step 3: Compute the final MSM result by combining all partitions
fn test_combine_msm() {
    let points = generate_points(10);
    let scalars: Vec<u32> = vec![0b01; 10]; // Scalars with only the second least significant bit set
    let window_size = 2;

    let partitions = partition_msm(&scalars, window_size);
    let combined_result = combine_partitioned_msm(&partitions, &points, window_size);
    // Compare against result from naive msm
    let expected_result = naive_msm(&points, &scalars);
    assert_eq!(combined_result, expected_result, "Combined MSM result is incorrect");
}

#[test]
// "Comprehensive test with 10 points"
fn test_pippenger_algorithm() {
    let points = generate_points(10);
    let scalars: Vec<u32> = generate_scalars(10);
    let window_size = 2;

    let msm_result = pippenger(&points, &scalars, window_size);
    // Compare against result from naive msm
    let expected_result = naive_msm(&points, &scalars);
    assert_eq!(msm_result, expected_result, "Pippenger algorithm did not match expected result");
}

#[test]
// "Comprehensive test with 100 points"
fn test_pippenger_algorithm_1() {
    let points = generate_points(100);
    let scalars: Vec<u32> = generate_scalars(100);
    let window_size = 2;

    let msm_result = pippenger(&points, &scalars, window_size);
    // Compare against result from naive msm
    let expected_result = naive_msm(&points, &scalars);
    assert_eq!(msm_result, expected_result, "Pippenger algorithm did not match expected result");
}

#[test]
// "Comprehensive test with 1000 points"
fn test_pippenger_algorithm_2() {
    let points = generate_points(1000);
    let scalars: Vec<u32> = generate_scalars(1000);
    let window_size = 2;

    let msm_result = pippenger(&points, &scalars, window_size);
    // Compare against result from naive msm
    let expected_result = naive_msm(&points, &scalars);
    assert_eq!(msm_result, expected_result, "Pippenger algorithm did not match expected result");
}