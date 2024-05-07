extern crate ark_ff;
extern crate ark_ec;
extern crate ark_mnt4_298;
extern crate ark_std;
extern crate std;

mod operations;
mod pippenger;
mod parallel_pippenger;
mod naive;
mod sid_pippenger;
mod subsum_pippenger;
mod parallel_sid_pippenger;
mod parallel_subsum_pippenger;
mod sid_subsum_pippenger;
mod parallel_sid_subsum_pippenger;
mod trivial;

use ark_mnt4_298::G1Projective;
use ark_std::{UniformRand, test_rng};
use rand::{Rng, thread_rng};
use std::time::Instant;
use parallel_pippenger::parallel_pippenger;
use pippenger::pippenger;
use naive::naive_msm;
use trivial::trivial_msm;
use sid_pippenger::sid_pippenger;
use subsum_pippenger::subsum_pippenger;
use parallel_sid_pippenger::parallel_sid_pippenger;
use parallel_subsum_pippenger::parallel_subsum_pippenger;
use sid_subsum_pippenger::sid_subsum_pippenger;
use parallel_sid_subsum_pippenger::parallel_sid_subsum_pippenger;

fn main() {

    fn generate_points(num_points: usize) -> Vec<G1Projective> {
        let mut rng = test_rng();
        (0..num_points).map(|_| G1Projective::rand(&mut rng)).collect()
    }
    
    fn generate_scalars(num_points: usize) -> Vec<u32> {
        let mut rng = thread_rng();
        (0..num_points).map(|_| rng.gen_range(0..4294967295)).collect() // 32-bit integers
    }

    let num_points = 100;

    let points = generate_points(num_points);
    let scalars = generate_scalars(num_points);
    let window_size = 2;

    let start = Instant::now();
    let result_naive = naive_msm(&points, &scalars);
    let duration_naive = start.elapsed();
        
    let start = Instant::now();
    let result_trivial = trivial_msm(&points, &scalars);
    let duration_trivial = start.elapsed();

    let start = Instant::now();
    let result_pippenger = pippenger(&points, &scalars, window_size);
    let duration_pippenger = start.elapsed();

    let start = Instant::now();
    let result_parallel = parallel_pippenger(&points, &scalars, window_size);
    let duration_parallel = start.elapsed();

    let start = Instant::now();
    let result_sid = sid_pippenger(&points, &scalars, window_size);
    let duration_sid = start.elapsed();

    let start = Instant::now();
    let result_subsum = subsum_pippenger(&points, &scalars, window_size);
    let duration_subsum = start.elapsed();

    let start = Instant::now();
    let result_parallel_sid = parallel_sid_pippenger(&points, &scalars, window_size);
    let duration_parallel_sid = start.elapsed();

    let start = Instant::now();
    let result_parallel_subsum = parallel_subsum_pippenger(&points, &scalars, window_size);
    let duration_parallel_subsum = start.elapsed();

    let start = Instant::now();
    let result_sid_subsum = sid_subsum_pippenger(&points, &scalars, window_size);
    let duration_sid_subsum = start.elapsed();

    let start = Instant::now();
    let result_parallel_sid_subsum = parallel_sid_subsum_pippenger(&points, &scalars, window_size);
    let duration_parallel_sid_subsum = start.elapsed();

    assert_eq!(result_naive, result_pippenger, "Results of pippenger and naive MSM should match");
    assert_eq!(result_naive, result_trivial, "Results of trivial and naive MSM should match");
    assert_eq!(result_naive, result_parallel, "Results of pippenger with parallelism and naive MSM should match");
    assert_eq!(result_naive, result_sid, "Results of pippenger with Signed Integer Decomposition should match with naive MSM");
    assert_eq!(result_naive, result_subsum, "Results of pippenger with more efficient subsum accumulation should match with naive MSM");
    assert_eq!(result_naive, result_parallel_sid, "Results of pippenger with Signed Integer Decomposition and parallelism should match with naive MSM");
    assert_eq!(result_naive, result_parallel_subsum, "Results of pippenger with more efficient subsum accumulation and parallelism should match with naive MSM");
    assert_eq!(result_naive, result_sid_subsum, "Results of pippenger with Signed Integer Decomposition and more efficient subsum accumulation should match with naive MSM");
    assert_eq!(result_naive, result_parallel_sid_subsum, "Results of pippenger with Signed Integer Decomposition, parallelism and more efficient subsum accumulation should match with naive MSM");

    println!("Naive: {:?}", duration_naive);
    println!("Trivial: {:?}", duration_trivial);
    println!("Pippenger: {:?}", duration_pippenger);
    println!("Parallel: {:?}", duration_parallel);
    println!("SID: {:?}", duration_sid);
    println!("Subsum: {:?}", duration_subsum);
    println!("Parallel SID: {:?}", duration_parallel_sid);
    println!("Parallel Subsum: {:?}", duration_parallel_subsum);
    println!("SID Subsum: {:?}", duration_sid_subsum);
    println!("Parallel SID Subsum: {:?}", duration_parallel_sid_subsum);    
}