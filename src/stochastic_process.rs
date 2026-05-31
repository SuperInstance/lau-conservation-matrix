use crate::ConservationMatrix;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StochasticProcess {
    pub transition: ConservationMatrix,
    pub state: Vec<f64>,
}

impl StochasticProcess {
    pub fn new(transition: ConservationMatrix, initial_state: Vec<f64>) -> Self {
        Self {
            transition,
            state: initial_state,
        }
    }

    pub fn step(&mut self) {
        self.state = self.transition.apply(&self.state);
    }

    pub fn run(&mut self, steps: usize) -> Vec<Vec<f64>> {
        let mut history = vec![self.state.clone()];
        for _ in 0..steps {
            self.step();
            history.push(self.state.clone());
        }
        history
    }

    /// Find stationary distribution via power iteration.
    pub fn stationary_distribution(&self) -> Vec<f64> {
        let n = self.transition.n;
        if n == 0 {
            return vec![];
        }
        // Start from uniform
        let mut v = vec![1.0 / n as f64; n];
        for _ in 0..1000 {
            let new_v = self.transition.apply(&v);
            let sum: f64 = new_v.iter().sum();
            if sum > 1e-15 {
                v = new_v.iter().map(|x| x / sum).collect();
            }
        }
        v
    }

    /// Steps until state is within epsilon of stationary distribution
    pub fn mixing_time(&self, epsilon: f64) -> usize {
        let n = self.transition.n;
        if n == 0 {
            return 0;
        }
        let uniform = vec![1.0 / n as f64; n];
        let mut state = vec![0.5; n]; // arbitrary start
        let sum: f64 = state.iter().sum();
        state = state.iter().map(|x| x / sum).collect();

        for step in 0..10000 {
            let max_diff = state
                .iter()
                .zip(uniform.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            if max_diff < epsilon {
                return step;
            }
            state = self.transition.apply(&state);
            let s: f64 = state.iter().sum();
            if s > 1e-15 {
                state = state.iter().map(|x| x / s).collect();
            }
        }
        10000
    }

    /// Detailed balance: πᵢ Pᵢⱼ = πⱼ Pⱼᵢ
    pub fn is_reversible(&self) -> bool {
        let pi = self.stationary_distribution();
        let n = self.transition.n;
        for i in 0..n {
            for j in 0..n {
                let lhs = pi[i] * self.transition.matrix[i][j];
                let rhs = pi[j] * self.transition.matrix[j][i];
                if (lhs - rhs).abs() > 1e-9 {
                    return false;
                }
            }
        }
        true
    }
}
