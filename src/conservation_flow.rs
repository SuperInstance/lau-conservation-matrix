use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConservationFlow {
    pub nodes: usize,
    pub edges: Vec<(usize, usize, f64)>,
}

impl ConservationFlow {
    pub fn new(nodes: usize, edges: Vec<(usize, usize, f64)>) -> Self {
        Self { nodes, edges }
    }

    pub fn flow_matrix(&self) -> Vec<Vec<f64>> {
        let mut m = vec![vec![0.0; self.nodes]; self.nodes];
        for &(from, to, cap) in &self.edges {
            m[from][to] += cap;
        }
        m
    }

    pub fn is_conserved_at_node(&self, node: usize) -> bool {
        let inflow: f64 = self
            .edges
            .iter()
            .filter(|(_, to, _)| *to == node)
            .map(|(_, _, c)| c)
            .sum();
        let outflow: f64 = self
            .edges
            .iter()
            .filter(|(from, _, _)| *from == node)
            .map(|(_, _, c)| c)
            .sum();
        (inflow - outflow).abs() < 1e-9
    }

    pub fn is_conserved(&self) -> bool {
        // Conservation at all nodes: total inflow = total outflow at every node
        let mut inflow = vec![0.0; self.nodes];
        let mut outflow = vec![0.0; self.nodes];
        for &(from, to, cap) in &self.edges {
            outflow[from] += cap;
            inflow[to] += cap;
        }
        for i in 0..self.nodes {
            if (inflow[i] - outflow[i]).abs() > 1e-9 {
                return false;
            }
        }
        true
    }

    pub fn total_flow(&self) -> f64 {
        self.edges.iter().map(|(_, _, c)| c).sum()
    }
}
