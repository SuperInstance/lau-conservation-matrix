use crate::{ConservationMatrix, MatrixError};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BirkhoffDecomposition {
    pub permutations: Vec<Vec<Vec<f64>>>,
    pub weights: Vec<f64>,
}

impl BirkhoffDecomposition {
    /// Birkhoff-von Neumann decomposition: greedy algorithm.
    pub fn decompose(matrix: &ConservationMatrix) -> Self {
        let n = matrix.n;
        if n == 0 {
            return Self {
                permutations: vec![],
                weights: vec![],
            };
        }

        let mut remaining = matrix.matrix.clone();
        let mut permutations = Vec::new();
        let mut weights = Vec::new();

        let max_iters = n * n;
        for _ in 0..max_iters {
            // Check if anything remains
            let has_positive = remaining
                .iter()
                .flat_map(|r| r.iter())
                .any(|&v| v > 1e-12);
            if !has_positive {
                break;
            }

            // Find a permutation via greedy matching
            let perm = Self::find_permutation(&remaining, n);
            let Some(perm) = perm else {
                break;
            };

            // Weight = minimum positive entry at permutation positions
            let w = (0..n)
                .map(|i| remaining[i][perm[i]])
                .fold(f64::MAX, f64::min);

            if w < 1e-14 {
                break;
            }

            // Build permutation matrix
            let mut perm_matrix = vec![vec![0.0; n]; n];
            for i in 0..n {
                perm_matrix[i][perm[i]] = 1.0;
            }

            // Subtract
            for i in 0..n {
                remaining[i][perm[i]] = (remaining[i][perm[i]] - w).max(0.0);
            }

            permutations.push(perm_matrix);
            weights.push(w);
        }

        Self {
            permutations,
            weights,
        }
    }

    fn find_permutation(remaining: &[Vec<f64>], n: usize) -> Option<Vec<usize>> {
        let mut assignment = vec![None; n];
        let mut col_used = vec![false; n];
        Self::assign_row(remaining, 0, &mut assignment, &mut col_used)
    }

    fn assign_row(
        remaining: &[Vec<f64>],
        row: usize,
        assignment: &mut Vec<Option<usize>>,
        col_used: &mut Vec<bool>,
    ) -> Option<Vec<usize>> {
        let n = assignment.len();
        if row == n {
            return Some(assignment.iter().map(|a| a.unwrap()).collect());
        }
        for col in 0..n {
            if !col_used[col] && remaining[row][col] > 1e-12 {
                col_used[col] = true;
                assignment[row] = Some(col);
                if let Some(result) = Self::assign_row(remaining, row + 1, assignment, col_used) {
                    return Some(result);
                }
                col_used[col] = false;
                assignment[row] = None;
            }
        }
        None
    }

    pub fn reconstruct(&self) -> Result<ConservationMatrix, MatrixError> {
        if self.permutations.is_empty() {
            return ConservationMatrix::new(vec![]);
        }
        let n = self.permutations[0].len();
        let mut result = vec![vec![0.0; n]; n];
        for (k, perm) in self.permutations.iter().enumerate() {
            let w = self.weights[k];
            for (i, row) in perm.iter().enumerate() {
                for (j, val) in row.iter().enumerate() {
                    result[i][j] += w * val;
                }
            }
        }
        ConservationMatrix::new(result)
    }

    pub fn num_permutations(&self) -> usize {
        self.permutations.len()
    }
}
