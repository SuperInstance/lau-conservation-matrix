# lau-conservation-matrix

Conservation matrix — doubly stochastic transition matrices and the mathematical structure of fair exchange.

## Modules

- **ConservationMatrix** — doubly stochastic matrix validation, entropy, trace, apply
- **BirkhoffDecomposition** — decompose doubly stochastic into convex combination of permutations (Birkhoff-von Neumann)
- **SinkhornBalancer** — project arbitrary positive matrices to doubly stochastic via Sinkhorn iteration
- **ExchangeGraph** — agents exchanging resources with conservation guarantees
- **ConservationFlow** — flow through networks with conservation at each node
- **StochasticProcess** — Markov chains with doubly stochastic transitions

## Theorems Verified

1. Identity matrix is doubly stochastic
2. Permutation matrices are doubly stochastic
3. Convex combination of doubly stochastic matrices is doubly stochastic
4. Birkhoff-von Neumann decomposition
5. Sinkhorn balancing convergence for positive matrices
6. Doubly stochastic preserves total sum
7. Uniform stationary distribution for doubly stochastic Markov chains
8. Finite mixing time for irreducible, aperiodic chains
9. Exchange graph conserves total resources
10. Conservation flow: inflow = outflow at interior nodes
11. Entropy of permutation matrix = 0
12. Entropy of uniform matrix is maximal
