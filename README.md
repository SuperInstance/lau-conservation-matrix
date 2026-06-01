# lau-conservation-matrix

> Conservation matrix — doubly stochastic transition matrices, Birkhoff decomposition, Sinkhorn balancing, exchange graphs, conservation flows, and stochastic processes

## What This Does

Conservation matrix — doubly stochastic transition matrices, Birkhoff decomposition, Sinkhorn balancing, exchange graphs, conservation flows, and stochastic processes. Part of the PLATO/LAU ecosystem — a mathematically rigorous framework for building educational agents that learn, teach, and evolve.

## The Key Idea

This crate implements the core abstractions needed for its domain, with a focus on correctness, composability, and conservation guarantees. Every public type is serializable (serde), every algorithm is tested, and every invariant is verified.

## Install

```bash
cargo add lau-conservation-matrix
```

## Quick Start

See the API Reference below for complete usage. Key entry points:

```rust
use lau_conservation_matrix::*;
// See types and methods below for complete usage
```

## API Reference

```rust
pub struct ConservationMatrix 
    pub fn new(matrix: Vec<Vec<f64>>) -> Result<Self, MatrixError> 
    pub fn is_row_stochastic(&self) -> bool 
    pub fn is_column_stochastic(&self) -> bool 
    pub fn is_doubly_stochastic(&self) -> bool 
    pub fn row_sum(&self, i: usize) -> f64 
    pub fn col_sum(&self, j: usize) -> f64 
    pub fn total_sum(&self) -> f64 
    pub fn apply(&self, vector: &[f64]) -> Vec<f64> 
    pub fn trace(&self) -> f64 
    pub fn is_permutation(&self) -> bool 
    pub fn is_identity(&self) -> bool 
    pub fn entropy(&self) -> f64 
pub struct ExchangeGraph 
    pub fn new(agents: usize, resources: usize, exchange_matrix: ConservationMatrix) -> Self 
    pub fn execute(&self, holdings: &[Vec<f64>]) -> Vec<Vec<f64>> 
    pub fn total_conserved(&self, before: &[Vec<f64>], after: &[Vec<f64>]) -> bool 
    pub fn fairness(&self) -> f64 
    pub fn try_new(
pub struct SinkhornBalancer 
    pub fn new(max_iterations: usize, tolerance: f64) -> Self 
    pub fn balance(&mut self, matrix: &[Vec<f64>]) -> Result<ConservationMatrix, MatrixError> 
    pub fn convergence_history(&self) -> Vec<f64> 
pub struct BirkhoffDecomposition 
    pub fn decompose(matrix: &ConservationMatrix) -> Self 
    pub fn reconstruct(&self) -> Result<ConservationMatrix, MatrixError> 
    pub fn num_permutations(&self) -> usize 
pub struct ConservationFlow 
    pub fn new(nodes: usize, edges: Vec<(usize, usize, f64)>) -> Self 
    pub fn flow_matrix(&self) -> Vec<Vec<f64>> 
    pub fn is_conserved_at_node(&self, node: usize) -> bool 
    pub fn is_conserved(&self) -> bool 
    pub fn total_flow(&self) -> f64 
pub struct StochasticProcess 
    pub fn new(transition: ConservationMatrix, initial_state: Vec<f64>) -> Self 
    pub fn step(&mut self) 
    pub fn run(&mut self, steps: usize) -> Vec<Vec<f64>> 
    pub fn stationary_distribution(&self) -> Vec<f64> 
    pub fn mixing_time(&self, epsilon: f64) -> usize 
    pub fn is_reversible(&self) -> bool 
pub enum MatrixError 
```

## How It Works

Read the source in `src/` for full implementation details. All algorithms are documented with inline comments explaining the mathematical foundations.

## The Math

This crate implements formal mathematical constructs. See the source documentation for theorem statements and proofs of correctness.

## Testing

**66 tests** covering construction, serialization, correctness properties, edge cases, and composability with other lau-* crates.

## License

MIT
