use lau_conservation_matrix::*;

// ─── ConservationMatrix basics ───

#[test]
fn test_identity_is_doubly_stochastic() {
    let m = ConservationMatrix::new(vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ])
    .unwrap();
    assert!(m.is_doubly_stochastic());
    assert!(m.is_identity());
}

#[test]
fn test_permutation_is_doubly_stochastic() {
    let m = ConservationMatrix::new(vec![
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
        vec![1.0, 0.0, 0.0],
    ])
    .unwrap();
    assert!(m.is_doubly_stochastic());
    assert!(m.is_permutation());
}

#[test]
fn test_uniform_matrix_is_doubly_stochastic() {
    let v = 1.0 / 3.0;
    let m = ConservationMatrix::new(vec![
        vec![v, v, v],
        vec![v, v, v],
        vec![v, v, v],
    ])
    .unwrap();
    assert!(m.is_doubly_stochastic());
    assert!(!m.is_permutation());
}

#[test]
fn test_not_square_error() {
    let err = ConservationMatrix::new(vec![vec![1.0, 0.0], vec![0.0, 1.0], vec![0.0, 0.0]]);
    assert!(matches!(err, Err(MatrixError::NotSquare)));
}

#[test]
fn test_negative_entry_error() {
    let err = ConservationMatrix::new(vec![vec![-0.5, 1.5], vec![1.0, 0.0]]);
    assert!(matches!(err, Err(MatrixError::InvalidEntry(0, 0, -0.5))));
}

#[test]
fn test_row_sum() {
    let m = ConservationMatrix::new(vec![vec![0.3, 0.7], vec![0.6, 0.4]]).unwrap();
    assert!((m.row_sum(0) - 1.0).abs() < 1e-9);
    assert!((m.row_sum(1) - 1.0).abs() < 1e-9);
}

#[test]
fn test_col_sum() {
    let m = ConservationMatrix::new(vec![vec![0.5, 0.5], vec![0.5, 0.5]]).unwrap();
    assert!((m.col_sum(0) - 1.0).abs() < 1e-9);
    assert!((m.col_sum(1) - 1.0).abs() < 1e-9);
}

#[test]
fn test_total_sum() {
    let m = ConservationMatrix::new(vec![
        vec![0.25, 0.75],
        vec![0.75, 0.25],
    ])
    .unwrap();
    assert!((m.total_sum() - 2.0).abs() < 1e-9);
}

#[test]
fn test_apply() {
    let m = ConservationMatrix::new(vec![vec![0.0, 1.0], vec![1.0, 0.0]]).unwrap();
    let result = m.apply(&[3.0, 5.0]);
    assert!((result[0] - 5.0).abs() < 1e-9);
    assert!((result[1] - 3.0).abs() < 1e-9);
}

#[test]
fn test_trace() {
    let m = ConservationMatrix::new(vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ])
    .unwrap();
    assert!((m.trace() - 3.0).abs() < 1e-9);
}

#[test]
fn test_entropy_permutation_is_zero() {
    let m = ConservationMatrix::new(vec![vec![1.0, 0.0], vec![0.0, 1.0]]).unwrap();
    assert!(m.entropy().abs() < 1e-9);
}

#[test]
fn test_entropy_uniform_is_maximal() {
    let v = 0.5;
    let m = ConservationMatrix::new(vec![vec![v, v], vec![v, v]]).unwrap();
    let expected = -4.0 * v * v.ln(); // max entropy for 2x2
    assert!((m.entropy() - expected).abs() < 1e-9);
}

#[test]
fn test_empty_matrix() {
    let m = ConservationMatrix::new(vec![]).unwrap();
    assert!(m.is_doubly_stochastic());
    assert!(m.is_identity());
    assert!(m.is_permutation());
}

#[test]
fn test_is_row_stochastic_true() {
    let m = ConservationMatrix::new(vec![vec![0.2, 0.8], vec![1.0, 0.0]]).unwrap();
    assert!(m.is_row_stochastic());
}

#[test]
fn test_is_column_stochastic_false() {
    let m = ConservationMatrix::new(vec![vec![0.2, 0.8], vec![1.0, 0.0]]).unwrap();
    assert!(!m.is_column_stochastic());
}

// ─── Theorem: Convex combination of doubly stochastic is doubly stochastic ───

#[test]
fn test_convex_combination_doubly_stochastic() {
    let a = ConservationMatrix::new(vec![
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
        vec![1.0, 0.0, 0.0],
    ])
    .unwrap();
    let b = ConservationMatrix::new(vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ])
    .unwrap();
    let alpha = 0.4;
    let n = a.n;
    let mut combined = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            combined[i][j] = alpha * a.matrix[i][j] + (1.0 - alpha) * b.matrix[i][j];
        }
    }
    let c = ConservationMatrix::new(combined).unwrap();
    assert!(c.is_doubly_stochastic());
}

// ─── Doubly stochastic preserves total sum ───

#[test]
fn test_doubly_stochastic_preserves_sum() {
    let m = ConservationMatrix::new(vec![
        vec![0.5, 0.5, 0.0],
        vec![0.0, 0.5, 0.5],
        vec![0.5, 0.0, 0.5],
    ])
    .unwrap();
    let state = vec![2.0, 3.0, 5.0];
    let after = m.apply(&state);
    let before_sum: f64 = state.iter().sum();
    let after_sum: f64 = after.iter().sum();
    assert!((before_sum - after_sum).abs() < 1e-9);
}

// ─── Birkhoff Decomposition ───

#[test]
fn test_birkhoff_identity() {
    let m = ConservationMatrix::new(vec![
        vec![1.0, 0.0],
        vec![0.0, 1.0],
    ])
    .unwrap();
    let decomp = BirkhoffDecomposition::decompose(&m);
    assert!(decomp.num_permutations() >= 1);
    let reconstructed = decomp.reconstruct().unwrap();
    for i in 0..m.n {
        for j in 0..m.n {
            assert!(
                (reconstructed.matrix[i][j] - m.matrix[i][j]).abs() < 1e-9,
                "Mismatch at ({}, {})",
                i, j
            );
        }
    }
}

#[test]
fn test_birkhoff_uniform() {
    let v = 0.5;
    let m = ConservationMatrix::new(vec![vec![v, v], vec![v, v]]).unwrap();
    let decomp = BirkhoffDecomposition::decompose(&m);
    let reconstructed = decomp.reconstruct().unwrap();
    for i in 0..m.n {
        for j in 0..m.n {
            assert!(
                (reconstructed.matrix[i][j] - m.matrix[i][j]).abs() < 1e-9,
                "Mismatch at ({}, {})",
                i, j
            );
        }
    }
}

#[test]
fn test_birkhoff_weights_sum_to_one() {
    let v = 1.0 / 3.0;
    let m = ConservationMatrix::new(vec![
        vec![v, v, v],
        vec![v, v, v],
        vec![v, v, v],
    ])
    .unwrap();
    let decomp = BirkhoffDecomposition::decompose(&m);
    let weight_sum: f64 = decomp.weights.iter().sum();
    assert!((weight_sum - 1.0).abs() < 1e-9);
}

#[test]
fn test_birkhoff_permutations_are_permutations() {
    let v = 1.0 / 3.0;
    let m = ConservationMatrix::new(vec![
        vec![v, v, v],
        vec![v, v, v],
        vec![v, v, v],
    ])
    .unwrap();
    let decomp = BirkhoffDecomposition::decompose(&m);
    for perm in &decomp.permutations {
        let pm = ConservationMatrix::new(perm.clone()).unwrap();
        assert!(pm.is_permutation());
    }
}

// ─── Sinkhorn ───

#[test]
fn test_sinkhorn_positive_matrix() {
    let mut balancer = SinkhornBalancer::new(1000, 1e-12);
    let input = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
    let result = balancer.balance(&input).unwrap();
    assert!(result.is_doubly_stochastic());
}

#[test]
fn test_sinkhorn_convergence() {
    let mut balancer = SinkhornBalancer::new(1000, 1e-12);
    let input = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0], vec![7.0, 8.0, 9.0]];
    let result = balancer.balance(&input).unwrap();
    assert!(result.is_doubly_stochastic());
    let history = balancer.convergence_history();
    assert!(!history.is_empty());
    assert!(*history.last().unwrap() < 1e-6);
}

#[test]
fn test_sinkhorn_already_doubly_stochastic() {
    let mut balancer = SinkhornBalancer::new(100, 1e-10);
    let input = vec![vec![0.5, 0.5], vec![0.5, 0.5]];
    let result = balancer.balance(&input).unwrap();
    assert!(result.is_doubly_stochastic());
}

#[test]
fn test_sinkhorn_empty() {
    let mut balancer = SinkhornBalancer::new(10, 1e-10);
    let result = balancer.balance(&vec![]).unwrap();
    assert!(result.is_doubly_stochastic());
}

// ─── ExchangeGraph ───

#[test]
fn test_exchange_conserves_resources() {
    let eg = ExchangeGraph::try_new(
        3,
        2,
        vec![
            vec![0.5, 0.25, 0.25],
            vec![0.25, 0.5, 0.25],
            vec![0.25, 0.25, 0.5],
        ],
    )
    .unwrap();
    let holdings = vec![vec![10.0, 20.0, 30.0], vec![5.0, 15.0, 10.0]];
    let after = eg.execute(&holdings);
    assert!(eg.total_conserved(&holdings, &after));
}

#[test]
fn test_exchange_fairness_uniform() {
    let v = 1.0 / 3.0;
    let eg = ExchangeGraph::try_new(
        3,
        1,
        vec![vec![v, v, v], vec![v, v, v], vec![v, v, v]],
    )
    .unwrap();
    assert!((eg.fairness() - 1.0).abs() < 1e-9);
}

#[test]
fn test_exchange_fairness_permutation() {
    let eg = ExchangeGraph::try_new(
        3,
        1,
        vec![
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
            vec![1.0, 0.0, 0.0],
        ],
    )
    .unwrap();
    assert!(eg.fairness() < 1.0);
}

#[test]
fn test_exchange_doubly_stochastic_preserves_totals() {
    let eg = ExchangeGraph::try_new(
        2,
        1,
        vec![vec![0.3, 0.7], vec![0.7, 0.3]],
    )
    .unwrap();
    let before = vec![vec![100.0, 50.0]];
    let after = eg.execute(&before);
    let sum_before: f64 = before[0].iter().sum();
    let sum_after: f64 = after[0].iter().sum();
    assert!((sum_before - sum_after).abs() < 1e-9);
}

// ─── ConservationFlow ───

#[test]
fn test_flow_conserved_cycle() {
    let flow = ConservationFlow::new(3, vec![(0, 1, 1.0), (1, 2, 1.0), (2, 0, 1.0)]);
    assert!(flow.is_conserved());
    assert!(flow.is_conserved_at_node(0));
    assert!(flow.is_conserved_at_node(1));
    assert!(flow.is_conserved_at_node(2));
}

#[test]
fn test_flow_not_conserved() {
    let flow = ConservationFlow::new(3, vec![(0, 1, 2.0), (1, 2, 1.0)]);
    assert!(!flow.is_conserved());
}

#[test]
fn test_flow_total_flow() {
    let flow = ConservationFlow::new(3, vec![(0, 1, 2.0), (1, 2, 3.0), (2, 0, 1.5)]);
    assert!((flow.total_flow() - 6.5).abs() < 1e-9);
}

#[test]
fn test_flow_matrix() {
    let flow = ConservationFlow::new(2, vec![(0, 1, 3.0), (1, 0, 2.0)]);
    let fm = flow.flow_matrix();
    assert!((fm[0][1] - 3.0).abs() < 1e-9);
    assert!((fm[1][0] - 2.0).abs() < 1e-9);
}

#[test]
fn test_flow_empty() {
    let flow = ConservationFlow::new(3, vec![]);
    assert!(flow.is_conserved());
    assert!((flow.total_flow() - 0.0).abs() < 1e-9);
}

// ─── StochasticProcess ───

#[test]
fn test_markov_step() {
    let transition = ConservationMatrix::new(vec![
        vec![0.5, 0.5],
        vec![0.5, 0.5],
    ])
    .unwrap();
    let mut sp = StochasticProcess::new(transition, vec![1.0, 0.0]);
    sp.step();
    assert!((sp.state[0] - 0.5).abs() < 1e-9);
    assert!((sp.state[1] - 0.5).abs() < 1e-9);
}

#[test]
fn test_markov_run() {
    let transition = ConservationMatrix::new(vec![
        vec![0.5, 0.5],
        vec![0.5, 0.5],
    ])
    .unwrap();
    let mut sp = StochasticProcess::new(transition, vec![1.0, 0.0]);
    let history = sp.run(3);
    assert_eq!(history.len(), 4); // initial + 3 steps
}

#[test]
fn test_uniform_stationary_for_doubly_stochastic() {
    let transition = ConservationMatrix::new(vec![
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
        vec![1.0, 0.0, 0.0],
    ])
    .unwrap();
    let sp = StochasticProcess::new(transition, vec![1.0, 0.0, 0.0]);
    let stat = sp.stationary_distribution();
    for v in &stat {
        assert!((v - 1.0 / 3.0).abs() < 1e-6);
    }
}

#[test]
fn test_mixing_time_finite() {
    let transition = ConservationMatrix::new(vec![
        vec![0.5, 0.5],
        vec![0.5, 0.5],
    ])
    .unwrap();
    let sp = StochasticProcess::new(transition, vec![1.0, 0.0]);
    let mt = sp.mixing_time(0.01);
    assert!(mt < 100);
}

#[test]
fn test_reversible() {
    let transition = ConservationMatrix::new(vec![
        vec![0.5, 0.5],
        vec![0.5, 0.5],
    ])
    .unwrap();
    let sp = StochasticProcess::new(transition, vec![0.5, 0.5]);
    assert!(sp.is_reversible());
}

#[test]
fn test_not_reversible() {
    // Non-doubly-stochastic transition that is not reversible
    let transition = ConservationMatrix::new(vec![
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
        vec![1.0, 0.0, 0.0],
    ])
    .unwrap();
    let sp = StochasticProcess::new(transition, vec![1.0 / 3.0; 3]);
    // This is actually reversible since uniform pi with cyclic perm...
    // A permutation matrix IS reversible with uniform pi (trivially πᵢPᵢⱼ = πⱼPⱼᵢ for cyclic)
    // Let's check: πᵢ=1/3 for all, P₀₁=1, P₁₀=0 → 1/3*1 ≠ 1/3*0 → not reversible!
    assert!(!sp.is_reversible());
}

// ─── Serde round-trip ───

#[test]
fn test_serde_conservation_matrix() {
    let m = ConservationMatrix::new(vec![vec![0.5, 0.5], vec![0.5, 0.5]]).unwrap();
    let json = serde_json::to_string(&m).unwrap();
    let m2: ConservationMatrix = serde_json::from_str(&json).unwrap();
    assert_eq!(m, m2);
}

#[test]
fn test_serde_matrix_error() {
    let err = MatrixError::InvalidEntry(1, 2, -0.5);
    let json = serde_json::to_string(&err).unwrap();
    let err2: MatrixError = serde_json::from_str(&json).unwrap();
    assert_eq!(err, err2);
}

#[test]
fn test_serde_sinkhorn() {
    let s = SinkhornBalancer::new(50, 0.001);
    let json = serde_json::to_string(&s).unwrap();
    let s2: SinkhornBalancer = serde_json::from_str(&json).unwrap();
    assert_eq!(s.max_iterations, s2.max_iterations);
    assert_eq!(s.tolerance, s2.tolerance);
}

#[test]
fn test_serde_exchange_graph() {
    let eg = ExchangeGraph::try_new(2, 1, vec![vec![0.5, 0.5], vec![0.5, 0.5]]).unwrap();
    let json = serde_json::to_string(&eg).unwrap();
    let eg2: ExchangeGraph = serde_json::from_str(&json).unwrap();
    assert_eq!(eg.agents, eg2.agents);
}

#[test]
fn test_serde_conservation_flow() {
    let cf = ConservationFlow::new(3, vec![(0, 1, 1.0), (1, 2, 2.0)]);
    let json = serde_json::to_string(&cf).unwrap();
    let cf2: ConservationFlow = serde_json::from_str(&json).unwrap();
    assert_eq!(cf.nodes, cf2.nodes);
    assert_eq!(cf.edges, cf2.edges);
}

#[test]
fn test_serde_stochastic_process() {
    let t = ConservationMatrix::new(vec![vec![1.0, 0.0], vec![0.0, 1.0]]).unwrap();
    let sp = StochasticProcess::new(t, vec![0.5, 0.5]);
    let json = serde_json::to_string(&sp).unwrap();
    let sp2: StochasticProcess = serde_json::from_str(&json).unwrap();
    assert_eq!(sp.state, sp2.state);
}

#[test]
fn test_serde_birkhoff() {
    let m = ConservationMatrix::new(vec![vec![1.0, 0.0], vec![0.0, 1.0]]).unwrap();
    let bd = BirkhoffDecomposition::decompose(&m);
    let json = serde_json::to_string(&bd).unwrap();
    let bd2: BirkhoffDecomposition = serde_json::from_str(&json).unwrap();
    assert_eq!(bd.weights, bd2.weights);
}

// ─── Additional theorem verifications ───

#[test]
fn test_theorem_identity_doubly_stochastic() {
    let n = 5;
    let mut m = vec![vec![0.0; n]; n];
    for i in 0..n {
        m[i][i] = 1.0;
    }
    let cm = ConservationMatrix::new(m).unwrap();
    assert!(cm.is_doubly_stochastic());
    assert!(cm.is_identity());
    assert!(cm.is_permutation());
}

#[test]
fn test_theorem_permutation_doubly_stochastic() {
    // Reverse permutation
    let n = 4;
    let mut m = vec![vec![0.0; n]; n];
    for i in 0..n {
        m[i][n - 1 - i] = 1.0;
    }
    let cm = ConservationMatrix::new(m).unwrap();
    assert!(cm.is_doubly_stochastic());
    assert!(cm.is_permutation());
    assert!(!cm.is_identity());
}

#[test]
fn test_theorem_preserves_sum_3x3() {
    let m = ConservationMatrix::new(vec![
        vec![0.2, 0.3, 0.5],
        vec![0.4, 0.1, 0.5],
        vec![0.4, 0.6, 0.0],
    ])
    .unwrap();
    assert!(m.is_doubly_stochastic());
    let state = vec![10.0, 20.0, 30.0];
    let after = m.apply(&state);
    assert!((state.iter().sum::<f64>() - after.iter().sum::<f64>()).abs() < 1e-9);
}

#[test]
fn test_markov_uniform_stationary_general() {
    // Any doubly stochastic transition → uniform stationary
    let m = ConservationMatrix::new(vec![
        vec![0.5, 0.25, 0.25],
        vec![0.25, 0.5, 0.25],
        vec![0.25, 0.25, 0.5],
    ])
    .unwrap();
    assert!(m.is_doubly_stochastic());
    let sp = StochasticProcess::new(m, vec![1.0, 0.0, 0.0]);
    let stat = sp.stationary_distribution();
    for v in &stat {
        assert!((v - 1.0 / 3.0).abs() < 1e-6);
    }
}

#[test]
fn test_entropy_permutation_zero() {
    let perms: Vec<Vec<Vec<f64>>> = vec![
        vec![vec![1.0, 0.0], vec![0.0, 1.0]],
        vec![vec![0.0, 1.0], vec![1.0, 0.0]],
    ];
    for p in perms {
        let m = ConservationMatrix::new(p).unwrap();
        assert!(m.is_permutation());
        assert!(m.entropy().abs() < 1e-9);
    }
}

#[test]
fn test_entropy_uniform_maximal_3x3() {
    let v = 1.0 / 3.0;
    let m = ConservationMatrix::new(vec![vec![v, v, v]; 3]).unwrap();
    let entropy = m.entropy();
    let expected = -9.0 * v * v.ln();
    assert!((entropy - expected).abs() < 1e-9);
    // Check it's larger than a non-uniform doubly stochastic entropy
    let m2 = ConservationMatrix::new(vec![
        vec![0.5, 0.25, 0.25],
        vec![0.25, 0.5, 0.25],
        vec![0.25, 0.25, 0.5],
    ])
    .unwrap();
    assert!(entropy > m2.entropy());
}

#[test]
fn test_sinkhorn_3x3_positive() {
    let mut balancer = SinkhornBalancer::new(1000, 1e-12);
    let input = vec![
        vec![1.0, 2.0, 3.0],
        vec![4.0, 5.0, 6.0],
        vec![7.0, 8.0, 9.0],
    ];
    let result = balancer.balance(&input).unwrap();
    assert!(result.is_row_stochastic());
    assert!(result.is_column_stochastic());
}

#[test]
fn test_birkhoff_3x3_mixed() {
    let m = ConservationMatrix::new(vec![
        vec![0.5, 0.25, 0.25],
        vec![0.25, 0.5, 0.25],
        vec![0.25, 0.25, 0.5],
    ])
    .unwrap();
    let decomp = BirkhoffDecomposition::decompose(&m);
    let reconstructed = decomp.reconstruct().unwrap();
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (reconstructed.matrix[i][j] - m.matrix[i][j]).abs() < 1e-9,
                "at ({}, {})",
                i, j
            );
        }
    }
    let ws: f64 = decomp.weights.iter().sum();
    assert!((ws - 1.0).abs() < 1e-9);
}

#[test]
fn test_flow_conservation_interior_nodes() {
    // Line: 0 → 1 → 2, with 1.0 units each
    let flow = ConservationFlow::new(3, vec![(0, 1, 1.0), (1, 2, 1.0)]);
    // Node 1: inflow=1, outflow=1 → conserved
    assert!(flow.is_conserved_at_node(1));
    // Node 0: inflow=0, outflow=1 → not conserved (source)
    assert!(!flow.is_conserved_at_node(0));
    // Overall not conserved (source and sink)
    assert!(!flow.is_conserved());
}

#[test]
fn test_process_stationary_is_fixed_point() {
    let m = ConservationMatrix::new(vec![
        vec![0.5, 0.25, 0.25],
        vec![0.25, 0.5, 0.25],
        vec![0.25, 0.25, 0.5],
    ])
    .unwrap();
    let stat = {
        let sp = StochasticProcess::new(m.clone(), vec![1.0, 0.0, 0.0]);
        sp.stationary_distribution()
    };
    let after = m.apply(&stat);
    for i in 0..3 {
        assert!((stat[i] - after[i]).abs() < 1e-6);
    }
}

#[test]
fn test_exchange_graph_multiple_resources() {
    let eg = ExchangeGraph::try_new(
        3,
        2,
        vec![
            vec![0.6, 0.2, 0.2],
            vec![0.2, 0.6, 0.2],
            vec![0.2, 0.2, 0.6],
        ],
    )
    .unwrap();
    let holdings = vec![vec![10.0, 20.0, 30.0], vec![5.0, 10.0, 15.0]];
    let after = eg.execute(&holdings);
    assert!(eg.total_conserved(&holdings, &after));
    // Check each resource individually
    for r in 0..2 {
        let before_sum: f64 = holdings[r].iter().sum();
        let after_sum: f64 = after[r].iter().sum();
        assert!((before_sum - after_sum).abs() < 1e-9);
    }
}

#[test]
fn test_mixing_time_decreases() {
    let transition = ConservationMatrix::new(vec![
        vec![0.5, 0.5, 0.0],
        vec![0.0, 0.5, 0.5],
        vec![0.5, 0.0, 0.5],
    ])
    .unwrap();
    let sp = StochasticProcess::new(transition, vec![1.0, 0.0, 0.0]);
    let mt_loose = sp.mixing_time(0.1);
    let mt_tight = sp.mixing_time(0.001);
    assert!(mt_tight >= mt_loose);
}

#[test]
fn test_flow_multi_edge() {
    let flow = ConservationFlow::new(
        2,
        vec![(0, 1, 1.5), (0, 1, 0.5), (1, 0, 2.0)],
    );
    // Node 0: out=2.0, in=2.0
    // Node 1: in=2.0, out=2.0
    assert!(flow.is_conserved());
    assert!((flow.total_flow() - 4.0).abs() < 1e-9);
}

// ─── Additional tests to reach 65+ ───

#[test]
fn test_birkhoff_permutation_matrix() {
    let m = ConservationMatrix::new(vec![
        vec![0.0, 0.0, 1.0],
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
    ]).unwrap();
    let decomp = BirkhoffDecomposition::decompose(&m);
    assert_eq!(decomp.num_permutations(), 1);
    assert!((decomp.weights[0] - 1.0).abs() < 1e-9);
}

#[test]
fn test_sinkhorn_4x4() {
    let mut balancer = SinkhornBalancer::new(1000, 1e-12);
    let input = vec![
        vec![1.0, 2.0, 3.0, 4.0],
        vec![5.0, 6.0, 7.0, 8.0],
        vec![2.0, 3.0, 4.0, 5.0],
        vec![6.0, 7.0, 8.0, 9.0],
    ];
    let result = balancer.balance(&input).unwrap();
    assert!(result.is_doubly_stochastic());
}

#[test]
fn test_apply_preserves_dimension() {
    let m = ConservationMatrix::new(vec![
        vec![0.5, 0.5],
        vec![0.5, 0.5],
    ]).unwrap();
    let result = m.apply(&[1.0, 2.0]);
    assert_eq!(result.len(), 2);
}

#[test]
fn test_stochastic_process_multiple_steps() {
    let transition = ConservationMatrix::new(vec![
        vec![0.8, 0.2],
        vec![0.3, 0.7],
    ]).unwrap();
    let mut sp = StochasticProcess::new(transition, vec![1.0, 0.0]);
    let history = sp.run(200);
    // Should converge toward stationary distribution
    let final_state = history.last().unwrap();
    // Stationary: π₀*0.8 + π₁*0.3 = π₀, π₀+π₁=1 → π₀=0.6
    assert!((final_state[0] - 0.6).abs() < 0.1, "got {}", final_state[0]);
}

#[test]
fn test_exchange_conserves_each_resource() {
    let eg = ExchangeGraph::try_new(
        3,
        3,
        vec![
            vec![0.5, 0.3, 0.2],
            vec![0.2, 0.5, 0.3],
            vec![0.3, 0.2, 0.5],
        ],
    ).unwrap();
    let holdings = vec![
        vec![100.0, 200.0, 300.0],
        vec![10.0, 20.0, 30.0],
        vec![1.0, 2.0, 3.0],
    ];
    let after = eg.execute(&holdings);
    for r in 0..3 {
        let before_sum: f64 = holdings[r].iter().sum();
        let after_sum: f64 = after[r].iter().sum();
        assert!((before_sum - after_sum).abs() < 1e-9, "resource {} not conserved", r);
    }
}

#[test]
fn test_flow_conserved_star_topology() {
    // Center node 0 connects to 1,2,3 with equal flow
    let flow = ConservationFlow::new(
        4,
        vec![(0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0), (1, 0, 1.0), (2, 0, 1.0), (3, 0, 1.0)],
    );
    assert!(flow.is_conserved());
    assert!(flow.is_conserved_at_node(0));
    assert!((flow.total_flow() - 6.0).abs() < 1e-9);
}
