use crate::{ConservationMatrix, MatrixError};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SinkhornBalancer {
    pub max_iterations: usize,
    pub tolerance: f64,
    #[serde(skip)]
    history: Vec<f64>,
}

impl SinkhornBalancer {
    pub fn new(max_iterations: usize, tolerance: f64) -> Self {
        Self {
            max_iterations,
            tolerance,
            history: Vec::new(),
        }
    }

    /// Alternate row and column normalization to project to doubly stochastic.
    pub fn balance(&mut self, matrix: &[Vec<f64>]) -> Result<ConservationMatrix, MatrixError> {
        self.history.clear();
        let n = matrix.len();
        if n == 0 {
            return ConservationMatrix::new(vec![]);
        }
        for row in matrix.iter() {
            if row.len() != n {
                return Err(MatrixError::NotSquare);
            }
        }

        let mut m = matrix.to_owned();

        for _ in 0..self.max_iterations {
            // Row normalize
            for row in m.iter_mut().take(n) {
                let rs: f64 = row.iter().sum();
                if rs > 1e-15 {
                    for val in row.iter_mut() {
                        *val /= rs;
                    }
                }
            }

            // Column normalize
            for j in 0..n {
                let cs: f64 = m.iter().map(|r| r[j]).sum();
                if cs > 1e-15 {
                    for row in m.iter_mut().take(n) {
                        row[j] /= cs;
                    }
                }
            }

            // Compute max row-stochastic error
            let max_row_err = m
                .iter()
                .map(|row| (row.iter().sum::<f64>() - 1.0).abs())
                .fold(0.0_f64, f64::max);

            self.history.push(max_row_err);

            if max_row_err < self.tolerance {
                break;
            }
        }

        ConservationMatrix::new(m)
    }

    pub fn convergence_history(&self) -> Vec<f64> {
        self.history.clone()
    }
}
