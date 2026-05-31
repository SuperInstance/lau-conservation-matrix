use crate::error::MatrixError;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ConservationMatrix {
    pub matrix: Vec<Vec<f64>>,
    pub n: usize,
}

impl ConservationMatrix {
    pub fn new(matrix: Vec<Vec<f64>>) -> Result<Self, MatrixError> {
        let n = matrix.len();
        if n == 0 {
            return Ok(Self { matrix, n: 0 });
        }
        for (i, row) in matrix.iter().enumerate() {
            if row.len() != n {
                return Err(MatrixError::NotSquare);
            }
            for (j, &val) in row.iter().enumerate() {
                if val < 0.0 {
                    return Err(MatrixError::InvalidEntry(i, j, val));
                }
            }
        }
        Ok(Self { matrix, n })
    }

    pub fn is_row_stochastic(&self) -> bool {
        if self.n == 0 {
            return true;
        }
        (0..self.n).all(|i| (self.row_sum(i) - 1.0).abs() <= 1e-9)
    }

    pub fn is_column_stochastic(&self) -> bool {
        if self.n == 0 {
            return true;
        }
        (0..self.n).all(|j| (self.col_sum(j) - 1.0).abs() <= 1e-9)
    }

    pub fn is_doubly_stochastic(&self) -> bool {
        self.is_row_stochastic() && self.is_column_stochastic()
    }

    pub fn row_sum(&self, i: usize) -> f64 {
        self.matrix[i].iter().sum()
    }

    pub fn col_sum(&self, j: usize) -> f64 {
        self.matrix.iter().map(|row| row[j]).sum()
    }

    pub fn total_sum(&self) -> f64 {
        self.matrix.iter().flat_map(|r| r.iter()).sum()
    }

    pub fn apply(&self, vector: &[f64]) -> Vec<f64> {
        let mut result = vec![0.0; self.n];
        for (i, row) in self.matrix.iter().enumerate() {
            for (j, v) in vector.iter().enumerate() {
                result[i] += row[j] * v;
            }
        }
        result
    }

    pub fn trace(&self) -> f64 {
        (0..self.n).map(|i| self.matrix[i][i]).sum()
    }

    pub fn is_permutation(&self) -> bool {
        if !self.is_doubly_stochastic() {
            return false;
        }
        self.matrix
            .iter()
            .flat_map(|r| r.iter())
            .all(|&v| v == 0.0 || (v - 1.0).abs() < 1e-9)
    }

    pub fn is_identity(&self) -> bool {
        if self.n == 0 {
            return true;
        }
        for (i, row) in self.matrix.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                if (val - expected).abs() > 1e-9 {
                    return false;
                }
            }
        }
        true
    }

    pub fn entropy(&self) -> f64 {
        self.matrix
            .iter()
            .flat_map(|r| r.iter())
            .filter(|&&v| v > 0.0)
            .map(|&v| -v * v.ln())
            .sum()
    }
}
