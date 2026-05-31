mod birkhoff;
mod conservation_flow;
mod conservation_matrix;
mod error;
mod exchange_graph;
mod sinkhorn;
mod stochastic_process;

pub use birkhoff::BirkhoffDecomposition;
pub use conservation_flow::ConservationFlow;
pub use conservation_matrix::ConservationMatrix;
pub use error::MatrixError;
pub use exchange_graph::ExchangeGraph;
pub use sinkhorn::SinkhornBalancer;
pub use stochastic_process::StochasticProcess;
