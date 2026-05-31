use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum MatrixError {
    NotSquare,
    InvalidEntry(usize, usize, f64),
    NotRowStochastic,
    NotColumnStochastic,
}
