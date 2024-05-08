# Multi-Scalar-Multiplication-Algorithms-on-MNT-Curves

This repository is dedicated to the implementation of Multi-Scalar Multiplication (MSM) algorithms on the MNT4-298 curve using Rust. Our suite contains 10 algorithm implementations and an implementation of basic operations on elliptic curves. Our test suite contains one test file for each implementation in the experimental suite.

## System Requirements

To run these implementations, you will need the Rust programming environment set up on your machine. The code has been tested and confirmed to work with Rust edition "2021".

## Dependencies

- ark-mnt4-298: 0.4.0
- ark-ec: 0.4.0
- ark-ff: 0.4.0
- ark-std: 0.4.0
- rand: 0.8.5

Ensure that you have Cargo installed, as it will handle these dependencies automatically

## Installation instructions

1. Install Rust and Cargo using rustup.
2. Clone the repository to your local machine: "git clone https://github.com/19010755/Multi-Scalar-Multiplication-Algorithms-on-MNT-Curves"
3. Navigate to the cloned repository's root directory.
4. Run "cargo build --release" to compile the project.

## Usage

To run the main application and view the runtime outputs of the ten algorithm implementations:

1. Use the command "cargo run". This will execute the main.rs file, where the execution times of the algorithms are calculated and displayed.
2. Within main.rs, you can modify the window_size variable to adjust the window size, and num_points to change the number of point-scalar pairs used in the calculations.
3. The generate_scalar function includes a default maximum scalar value set to 4294967295, which is the maximum for a 32-bit unsigned integer. Feel free to adjust this value as needed to fit your testing requirements.
4. To obtain runtimes for specific stages of the Pippenger bucket method or to assess the additional cost of signed integer decomposition, uncomment the relevant timing lines in pippenger.rs and sid_pippenger.rs. By extension, you may also add these timing lines to any algorithm file to get hold of how long a specific step takes. With variability in variable names, look for lines similar to:

```rust
        use std::time::Instant; // This library will need to be used
        let start_partitioning = Instant::now(); // This will be at the start of the step you'd like to measure
        let duration_partitioning = start_partitioning.elapsed(); // This will be at the end of the step you'd like to measure
        println!("Partitioning took: {:?}", duration_partitioning); // Output when the command "cargo run" is executed
```

## Testing

For testing of the algorithms:

1. Run cargo test to execute the test suites for all implemented algorithms. This will verify the correctness of each algorithm and ensure they are functioning as expected.
# Multi-Scalar-Multiplication-Algorithms-on-MNT-Curves
