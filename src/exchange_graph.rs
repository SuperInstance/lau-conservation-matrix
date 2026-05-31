use crate::{ConservationMatrix, MatrixError};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExchangeGraph {
    pub agents: usize,
    pub resources: usize,
    pub exchange_matrix: ConservationMatrix,
}

impl ExchangeGraph {
    pub fn new(agents: usize, resources: usize, exchange_matrix: ConservationMatrix) -> Self {
        Self {
            agents,
            resources,
            exchange_matrix,
        }
    }

    /// holdings[r][a] = amount of resource r held by agent a
    /// result[r][a] = new amount of resource r for agent a
    pub fn execute(&self, holdings: &[Vec<f64>]) -> Vec<Vec<f64>> {
        holdings
            .iter()
            .map(|resource_holdings| self.exchange_matrix.apply(resource_holdings))
            .collect()
    }

    pub fn total_conserved(&self, before: &[Vec<f64>], after: &[Vec<f64>]) -> bool {
        for r in 0..self.resources {
            let sum_before: f64 = before[r].iter().sum();
            let sum_after: f64 = after[r].iter().sum();
            if (sum_before - sum_after).abs() > 1e-9 {
                return false;
            }
        }
        true
    }

    /// Fairness: 1.0 = uniform exchange, lower = more biased
    pub fn fairness(&self) -> f64 {
        if self.agents == 0 {
            return 1.0;
        }
        let uniform_val = 1.0 / self.agents as f64;
        let total_deviation: f64 = self
            .exchange_matrix
            .matrix
            .iter()
            .flat_map(|r| r.iter())
            .map(|&v| (v - uniform_val).abs())
            .sum();
        let max_deviation = 2.0 * (self.agents - 1) as f64 * uniform_val;
        if max_deviation < 1e-15 {
            1.0
        } else {
            1.0 - total_deviation / max_deviation
        }
    }
}

impl ExchangeGraph {
    pub fn try_new(
        agents: usize,
        resources: usize,
        exchange_matrix: Vec<Vec<f64>>,
    ) -> Result<Self, MatrixError> {
        let cm = ConservationMatrix::new(exchange_matrix)?;
        Ok(Self::new(agents, resources, cm))
    }
}
