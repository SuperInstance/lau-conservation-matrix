# lau-conservation-matrix

**Doubly stochastic transition matrices, Birkhoff decomposition, Sinkhorn balancing, exchange graphs, conservation flows, and stochastic processes.**

A crate for building, analyzing, and decomposing matrices that conserve mass — where rows sum to 1, columns sum to 1, and total energy is preserved through transitions. Covers the full pipeline: construct a matrix → verify conservation → balance it → decompose into permutations → simulate stochastic processes.

## What This Does

This library handles matrices where **what goes in must come out**:

- **ConservationMatrix** — An n×n matrix with non-negative entries. Validates row-stochastic, column-stochastic, and doubly-stochastic properties. Computes entropy, trace, and applies matrix-vector products.
- **SinkhornBalancer** — Iteratively normalizes rows then columns until the matrix converges to doubly-stochastic. Tracks convergence history.
- **BirkhoffDecomposition** — Decomposes any doubly-stochastic matrix into a convex combination of permutation matrices (the Birkhoff–von Neumann theorem). Uses a greedy algorithm with recursive assignment.
- **ExchangeGraph** — Models resource exchange between agents via a conservation matrix. Verifies conservation of total resources and measures fairness of exchange.
- **ConservationFlow** — A directed graph where edge weights represent flow. Checks conservation at each node (inflow = outflow).
- **StochasticProcess** — A Markov chain with power iteration for stationary distributions, mixing time estimation, and detailed balance (reversibility) checking.

## Key Idea

A **doubly stochastic matrix** (also called a bistochastic matrix) is the mathematical formalization of "fair redistribution." Every row sums to 1 (nothing created) and every column sums to 1 (nothing destroyed). The Birkhoff–von Neumann theorem guarantees that every such matrix is a **convex combination of permutation matrices** — meaning any fair redistribution can be decomposed into a lottery over deterministic reassignments.

This library makes that theorem operational: given any non-negative matrix, Sinkhorn balancing projects it to the nearest doubly-stochastic matrix, and Birkhoff decomposition breaks it into permutations with weights.

## Install

```toml
[dependencies]
lau-conservation-matrix = "0.1.0"
```

## Quick Start

### Build and verify a doubly-stochastic matrix

```rust
use lau_conservation_matrix::*;

let matrix = vec![
    vec![0.5, 0.5, 0.0],
    vec![0.0, 0.5, 0.5],
    vec![0.5, 0.0, 0.5],
];
let cm = ConservationMatrix::new(matrix)?;
assert!(cm.is_doubly_stochastic());
assert!(cm.is_row_stochastic());
assert!(cm.is_column_stochastic());
```

### Sinkhorn balancing

```rust
// Start with an arbitrary non-negative matrix
let matrix = vec![
    vec![1.0, 2.0, 3.0],
    vec![4.0, 5.0, 6.0],
    vec![7.0, 8.0, 9.0],
];

let mut balancer = SinkhornBalancer::new(100, 1e-10);
let balanced = balancer.balance(&matrix)?;
assert!(balanced.is_doubly_stochastic());

// Check convergence
let history = balancer.convergence_history();
println!("Converged in {} iterations", history.len());
```

### Birkhoff decomposition

```rust
let cm = ConservationMatrix::new(vec![
    vec![0.5, 0.5, 0.0],
    vec![0.0, 0.5, 0.5],
    vec![0.5, 0.0, 0.5],
])?;

let decomp = BirkhoffDecomposition::decompose(&cm);
println!("Decomposed into {} permutations", decomp.num_permutations());
for (i, w) in decomp.weights.iter().enumerate() {
    println!("  Permutation {} with weight {:.3}", i, w);
}

// Reconstruct (should match original)
let reconstructed = decomp.reconstruct()?;
```

### Stochastic process simulation

```rust
let transition = ConservationMatrix::new(vec![
    vec![0.7, 0.3, 0.0],
    vec![0.0, 0.8, 0.2],
    vec![0.1, 0.0, 0.9],
])?;

let mut process = StochasticProcess::new(transition, vec![1.0, 0.0, 0.0]);
let history = process.run(20);

// Find stationary distribution
let pi = process.stationary_distribution();
println!("Stationary: {:?}", pi);

// Check reversibility (detailed balance)
println!("Reversible: {}", process.is_reversible());
```

### Exchange graph

```rust
let exchange = ConservationMatrix::new(vec![
    vec![0.5, 0.3, 0.2],
    vec![0.2, 0.6, 0.2],
    vec![0.3, 0.1, 0.6],
])?;

let graph = ExchangeGraph::new(3, 1, exchange);
let holdings = vec![vec![10.0, 20.0, 30.0]]; // 1 resource, 3 agents
let after = graph.execute(&holdings);
assert!(graph.total_conserved(&holdings, &after));
println!("Fairness: {:.2}", graph.fairness());
```

## API Reference

### `ConservationMatrix`

| Method | Description |
|---|---|
| `new(matrix)` | Construct with validation (square, non-negative) |
| `is_row_stochastic()` | All rows sum to 1 |
| `is_column_stochastic()` | All columns sum to 1 |
| `is_doubly_stochastic()` | Both row and column stochastic |
| `is_permutation()` | Doubly stochastic with entries ∈ {0, 1} |
| `is_identity()` | Diagonal matrix with 1s |
| `row_sum(i)` / `col_sum(j)` | Sum of row/column |
| `total_sum()` | Sum of all entries |
| `apply(&vector)` | Matrix-vector multiplication |
| `trace()` | Sum of diagonal |
| `entropy()` | Shannon entropy of all entries |

### `SinkhornBalancer`

| Method | Description |
|---|---|
| `new(max_iterations, tolerance)` | Create with convergence parameters |
| `balance(&matrix)` | Project to doubly-stochastic |
| `convergence_history()` | Per-iteration max row error |

### `BirkhoffDecomposition`

| Method | Description |
|---|---|
| `decompose(&matrix)` | Greedy Birkhoff–von Neumann decomposition |
| `reconstruct()` | Rebuild matrix from permutations + weights |
| `num_permutations()` | Count of permutation matrices |

### `ExchangeGraph`

| Method | Description |
|---|---|
| `new(agents, resources, exchange_matrix)` | Create graph |
| `execute(&holdings)` | Apply exchange to resource holdings |
| `total_conserved(&before, &after)` | Verify resource conservation |
| `fairness()` | 1.0 = uniform, lower = biased |

### `ConservationFlow`

| Method | Description |
|---|---|
| `new(nodes, edges)` | Create flow graph with (from, to, capacity) edges |
| `flow_matrix()` | Adjacency matrix of flows |
| `is_conserved_at_node(node)` | Inflow = outflow at one node |
| `is_conserved()` | Conservation at all nodes |
| `total_flow()` | Sum of all edge weights |

### `StochasticProcess`

| Method | Description |
|---|---|
| `new(transition, initial_state)` | Create Markov chain |
| `step()` | One transition step |
| `run(steps)` | Multiple steps with history |
| `stationary_distribution()` | Power iteration to convergence |
| `mixing_time(epsilon)` | Steps to reach stationarity |
| `is_reversible()` | Detailed balance check (πᵢPᵢⱼ = πⱼPⱼᵢ) |

### `MatrixError`

```rust
enum MatrixError {
    NotSquare,
    InvalidEntry(usize, usize, f64),
    NotRowStochastic,
    NotColumnStochastic,
}
```

## How It Works

### Conservation Matrix

A `ConservationMatrix` wraps a `Vec<Vec<f64>>` with validation: must be square and non-negative. It provides fast checks for stochasticity — all sums compared against 1.0 within tolerance 1e-9. The `entropy()` method computes the Shannon entropy of the matrix entries, measuring how "spread out" the transitions are.

### Sinkhorn Balancing

The **Sinkhorn–Knopp algorithm** iteratively normalizes rows and columns:

1. Divide each row by its sum (row normalization)
2. Divide each column by its sum (column normalization)
3. Repeat until convergence

Each iteration reduces the maximum deviation from doubly-stochastic. The balancer tracks this error per iteration in a history vector. Convergence is guaranteed for positive matrices (all entries > 0); for matrices with zeros, convergence depends on the support structure.

### Birkhoff Decomposition

The **Birkhoff–von Neumann theorem** states that every n×n doubly-stochastic matrix can be written as:

```
M = w₁P₁ + w₂P₂ + ... + wₖPₖ
```

where each Pᵢ is a permutation matrix and wᵢ ≥ 0 with Σwᵢ = 1.

The greedy algorithm:
1. Find a permutation (perfect matching) in the remaining matrix using recursive backtracking
2. Extract the minimum weight at permutation positions as w
3. Subtract w·P from the remaining matrix
4. Repeat until the remainder is negligible

The number of permutations needed is at most n² − n + 1 (Marcus–Ree theorem), though in practice far fewer.

### Exchange Graph

An `ExchangeGraph` models resource redistribution among agents. The exchange matrix describes what fraction of each resource each agent sends to others. Executing an exchange applies the matrix to the current holdings vector. Conservation is verified by checking that total resource amounts don't change.

**Fairness** measures how close the exchange matrix is to uniform (all entries = 1/n). Fairness = 1.0 means perfectly uniform redistribution; lower values indicate bias.

### Conservation Flow

A `ConservationFlow` is a directed graph with weighted edges. Conservation at a node means total inflow equals total outflow. This models physical flow networks (fluid, current, traffic) where the law of conservation applies.

### Stochastic Process

A `StochasticProcess` wraps a transition matrix and a state vector. Each `step()` multiplies: **s' = P · s**. The **stationary distribution** π is found by power iteration: start from uniform, repeatedly apply P and renormalize, until convergence.

**Mixing time** is the number of steps from an arbitrary starting state until the state distribution is within ε of π (measured by max absolute deviation).

**Reversibility** checks the detailed balance condition: πᵢPᵢⱼ = πⱼPⱼᵢ for all (i, j). A reversible Markov chain has a symmetric flow structure and is easier to analyze.

## The Math

### Doubly Stochastic Matrices

A matrix M ∈ ℝⁿˣⁿ is doubly stochastic if:

```
Mᵢⱼ ≥ 0           ∀ i, j
Σⱼ Mᵢⱼ = 1         ∀ i     (row sums)
Σᵢ Mᵢⱼ = 1         ∀ j     (column sums)
```

The set of all n×n doubly-stochastic matrices forms the **Birkhoff polytope** Bₙ, a convex polytope in ℝⁿ². Its vertices are the n! permutation matrices.

### Birkhoff–von Neumann Theorem

Every doubly-stochastic matrix M is a convex combination of at most (n−1)² + 1 permutation matrices:

```
M = Σₖ wₖ Pₖ,  wₖ ≥ 0,  Σwₖ = 1
```

This is a consequence of the fact that Bₙ is the convex hull of the permutation matrices, combined with Carathéodory's theorem.

### Sinkhorn's Theorem

For any matrix A with positive entries, the iterative row-column normalization converges to a unique doubly-stochastic matrix D = D₁AD₂, where D₁ and D₂ are diagonal scaling matrices. Convergence is linear.

### Shannon Entropy of a Matrix

```
H(M) = −Σᵢⱼ Mᵢⱼ ln(Mᵢⱼ)
```

For a doubly-stochastic matrix, this ranges from 0 (permutation matrix) to 2ln(n) (uniform matrix with all entries 1/n). Higher entropy = more "mixed" transitions.

### Stationary Distribution

For an irreducible, aperiodic Markov chain with transition matrix P, the stationary distribution π satisfies:

```
πP = π,  Σπᵢ = 1
```

Power iteration converges because P has spectral radius 1 and the eigenvalue 1 has multiplicity 1 (by Perron–Frobenius).

### Detailed Balance

A chain satisfies **detailed balance** if there exists π such that:

```
πᵢPᵢⱼ = πⱼPⱼᵢ  ∀ i, j
```

This means the flow from i to j equals the flow from j to i in stationarity. Such chains are **time-reversible**.

## Testing

66 integration tests covering:

- ConservationMatrix: construction, validation, stochasticity checks, trace, entropy, identity, permutation detection
- SinkhornBalancer: convergence, doubly-stochastic output, empty input, tolerance
- BirkhoffDecomposition: decompose and reconstruct, weight sum = 1, identity decomposes to single permutation
- ExchangeGraph: conservation verification, fairness measurement, resource holding updates
- ConservationFlow: node conservation, total flow, flow matrix
- StochasticProcess: stepping, stationary distribution, mixing time, reversibility
- Serde round-trips for all types
- Edge cases: empty matrices, single-element matrices

Run with `cargo test`.

## License

MIT
