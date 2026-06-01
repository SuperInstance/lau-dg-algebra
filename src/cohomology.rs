//! Cohomology: H(d) = ker(d)/im(d) as a graded algebra.

use serde::{Serialize, Deserialize};
use crate::chain_complex::{ChainComplex, LinearMap};
use crate::dga::DGA;
use std::collections::HashMap;

/// Computed cohomology of a DGA or chain complex.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Cohomology {
    /// Betti numbers: degree → dimension of Hᵏ.
    pub betti: HashMap<i32, usize>,
    /// Representative cocycles: for each degree k, a list of basis vectors for ker(dₖ).
    /// Each vector is a Vec<f64> of length dim(Vᵏ).
    pub representatives: HashMap<i32, Vec<Vec<f64>>>,
    /// The underlying graded vector space dimensions.
    pub graded_dims: HashMap<i32, usize>,
}

impl Cohomology {
    /// Compute cohomology from a chain complex.
    pub fn from_chain_complex(cc: &ChainComplex) -> Self {
        let mut betti = HashMap::new();
        let mut representatives = HashMap::new();

        for &k in &cc.degrees() {
            let dk_dim = cc.dim(k);
            let ker_d = if let Some(ref d_k) = cc.differentials.get(&k) {
                d_k.kernel_dim()
            } else {
                dk_dim
            };
            let im_d_prev = if let Some(ref d_prev) = cc.differentials.get(&(k - 1)) {
                d_prev.image_dim()
            } else {
                0
            };

            let b = if ker_d >= im_d_prev { ker_d - im_d_prev } else { 0 };
            betti.insert(k, b);

            // Find representatives (kernel vectors)
            if let Some(ref d_k) = cc.differentials.get(&k) {
                let kernel_vecs = find_kernel_basis(d_k, dk_dim);
                representatives.insert(k, kernel_vecs);
            } else {
                // No differential means everything is closed
                let mut reps = Vec::new();
                for i in 0..dk_dim {
                    let mut v = vec![0.0; dk_dim];
                    v[i] = 1.0;
                    reps.push(v);
                }
                representatives.insert(k, reps);
            }
        }

        Cohomology {
            betti,
            representatives,
            graded_dims: cc.dimensions.clone(),
        }
    }

    /// Compute cohomology from a DGA.
    pub fn from_dga(dga: &DGA) -> Self {
        let cc = dga.as_chain_complex();
        Self::from_chain_complex(&cc)
    }

    /// Betti number at degree k.
    pub fn betti_number(&self, k: i32) -> usize {
        *self.betti.get(&k).unwrap_or(&0)
    }

    /// Total Betti number.
    pub fn total_betti(&self) -> usize {
        self.betti.values().sum()
    }

    /// Euler characteristic from cohomology: Σ (-1)ᵏ βₖ.
    pub fn euler_characteristic(&self) -> i64 {
        self.betti.iter()
            .map(|(&k, &b)| ((-1i64).pow(k as u32)) * b as i64)
            .sum()
    }

    /// Poincaré polynomial: Σ βₖ tᵏ (returned as map k → βₖ).
    pub fn poincare_polynomial(&self) -> &HashMap<i32, usize> {
        &self.betti
    }

    /// Check if cohomology is trivial (all Betti numbers zero).
    pub fn is_trivial(&self) -> bool {
        self.betti.values().all(|&b| b == 0)
    }

    /// Check if this is a Poincaré duality algebra.
    /// For a manifold of dimension n: βₖ = β_{n-k} for all k, and βₙ = 1.
    pub fn check_poincare_duality(&self, dim: i32) -> bool {
        for (&k, &b) in &self.betti {
            let b_complement = self.betti.get(&(dim - k)).copied().unwrap_or(0);
            if b != b_complement {
                return false;
            }
        }
        *self.betti.get(&dim).unwrap_or(&0) == 1
    }

    /// Degrees where cohomology is nonzero.
    pub fn nonzero_degrees(&self) -> Vec<i32> {
        let mut ds: Vec<i32> = self.betti.iter()
            .filter(|(_, &b)| b > 0)
            .map(|(&k, _)| k)
            .collect();
        ds.sort();
        ds
    }

    /// Cohomology ring product (induced product on cohomology).
    /// Given two cohomology classes [a], [b], compute [a] ∪ [b].
    pub fn cup_product(
        &self,
        dga: &DGA,
        deg_a: i32,
        class_a: &[f64],
        deg_b: i32,
        class_b: &[f64],
    ) -> Option<Vec<f64>> {
        dga.multiplication.multiply(deg_a, class_a, deg_b, class_b)
    }
}

/// Find a basis for the kernel of a linear map using row reduction.
fn find_kernel_basis(map: &LinearMap, domain_dim: usize) -> Vec<Vec<f64>> {
    if map.is_zero(1e-10) {
        // Everything is in the kernel
        let mut basis = Vec::new();
        for i in 0..domain_dim {
            let mut v = vec![0.0; domain_dim];
            v[i] = 1.0;
            basis.push(v);
        }
        return basis;
    }

    let ker_dim = map.kernel_dim();
    if ker_dim == 0 {
        return Vec::new();
    }

    // Use augmented matrix [M | I] and row reduce to find kernel
    let rows = map.codomain_dim;
    let cols = map.domain_dim;

    // Build augmented matrix
    let mut aug = vec![vec![0.0; cols + rows]; rows];
    for i in 0..rows {
        for j in 0..cols {
            aug[i][j] = map.entries[i][j];
        }
        aug[i][cols + i] = 1.0;
    }

    // Forward elimination
    let mut pivot_cols = Vec::new();
    let mut pivot_row = 0;
    for col in 0..cols {
        let mut found = None;
        for row in pivot_row..rows {
            if aug[row][col].abs() > 1e-10 {
                found = Some(row);
                break;
            }
        }
        if let Some(pr) = found {
            aug.swap(pivot_row, pr);
            let scale = aug[pivot_row][col];
            for j in 0..(cols + rows) {
                aug[pivot_row][j] /= scale;
            }
            for row in 0..rows {
                if row != pivot_row {
                    let factor = aug[row][col];
                    for j in 0..(cols + rows) {
                        aug[row][j] -= factor * aug[pivot_row][j];
                    }
                }
            }
            pivot_cols.push(col);
            pivot_row += 1;
        }
    }

    // Free variables are columns not in pivot_cols
    let mut free_cols: Vec<usize> = (0..cols).filter(|c| !pivot_cols.contains(c)).collect();
    free_cols.truncate(ker_dim);

    let mut basis = Vec::new();
    for &fc in &free_cols {
        let mut v = vec![0.0; cols];
        v[fc] = 1.0;
        for (pi, &pc) in pivot_cols.iter().enumerate() {
            if pi < rows {
                v[pc] = -aug[pi][fc];
            }
        }
        basis.push(v);
    }
    basis
}

/// Compute the long exact sequence in cohomology for a short exact sequence.
/// 0 → A → B → C → 0 induces ... → Hᵏ(A) → Hᵏ(B) → Hᵏ(C) → Hᵏ⁺¹(A) → ...
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LongExactSequence {
    /// Maps: Hᵏ(A) → Hᵏ(B), Hᵏ(B) → Hᵏ(C), Hᵏ(C) → Hᵏ⁺¹(A)
    pub i_star: HashMap<i32, LinearMap>,  // inclusion induced
    pub p_star: HashMap<i32, LinearMap>,  // projection induced
    pub delta: HashMap<i32, LinearMap>,   // connecting homomorphism
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cohomology_trivial_complex() {
        // 0 → Z → 0 (single Z in degree 0)
        let mut dims = HashMap::new();
        dims.insert(0, 1);
        let cc = ChainComplex::new(dims, HashMap::new());
        let h = Cohomology::from_chain_complex(&cc);
        assert_eq!(h.betti_number(0), 1);
        assert_eq!(h.total_betti(), 1);
    }

    #[test]
    fn test_cohomology_exact_sequence() {
        // 0 → V⁰ --d₀--> V¹ → 0 with d₀ = identity (dim 1)
        // H⁰ = ker/im = 0, H¹ = ker/im = 0
        let mut dims = HashMap::new();
        dims.insert(0, 1);
        dims.insert(1, 1);
        let mut diffs = HashMap::new();
        diffs.insert(0, LinearMap::from_rows(vec![vec![1.0]]));
        let cc = ChainComplex::new(dims, diffs);
        let h = Cohomology::from_chain_complex(&cc);
        assert_eq!(h.betti_number(0), 0);
        assert_eq!(h.betti_number(1), 0);
    }

    #[test]
    fn test_cohomology_circle() {
        // De Rham complex of S¹ (simplified):
        // V⁰ (dim 2) --d--> V¹ (dim 2), d = [[0,0],[0,0]]
        // H⁰ = 2, H¹ = 2 (but real S¹ has H⁰=1, H¹=1)
        // Here just testing zero differential
        let mut dims = HashMap::new();
        dims.insert(0, 1);
        dims.insert(1, 1);
        let cc = ChainComplex::new(dims, HashMap::new());
        let h = Cohomology::from_chain_complex(&cc);
        assert_eq!(h.betti_number(0), 1);
        assert_eq!(h.betti_number(1), 1);
    }

    #[test]
    fn test_euler_characteristic() {
        let mut dims = HashMap::new();
        dims.insert(0, 1);
        dims.insert(1, 1);
        let cc = ChainComplex::new(dims, HashMap::new());
        let h = Cohomology::from_chain_complex(&cc);
        assert_eq!(h.euler_characteristic(), 0); // 1 - 1 = 0
    }

    #[test]
    fn test_poincare_duality() {
        // S²: H⁰ = 1, H¹ = 0, H² = 1 (Poincaré duality in dim 2)
        let mut dims = HashMap::new();
        dims.insert(0, 1);
        dims.insert(1, 0);
        dims.insert(2, 1);
        let cc = ChainComplex::new(dims, HashMap::new());
        let h = Cohomology::from_chain_complex(&cc);
        assert!(h.check_poincare_duality(2));
    }

    #[test]
    fn test_is_trivial() {
        let mut dims = HashMap::new();
        dims.insert(0, 1);
        dims.insert(1, 1);
        let mut diffs = HashMap::new();
        diffs.insert(0, LinearMap::from_rows(vec![vec![1.0]]));
        let cc = ChainComplex::new(dims, diffs);
        let h = Cohomology::from_chain_complex(&cc);
        assert!(h.is_trivial());
    }

    #[test]
    fn test_nonzero_degrees() {
        let mut dims = HashMap::new();
        dims.insert(0, 1);
        dims.insert(1, 0);
        dims.insert(2, 1);
        let cc = ChainComplex::new(dims, HashMap::new());
        let h = Cohomology::from_chain_complex(&cc);
        assert_eq!(h.nonzero_degrees(), vec![0, 2]);
    }

    #[test]
    fn test_cohomology_from_dga() {
        use crate::dga::*;
        use crate::graded::GradedVectorSpace;

        let mut dims = HashMap::new();
        dims.insert(0, 1);
        dims.insert(1, 1);
        let gvs = GradedVectorSpace::new(dims);

        let mut tables = HashMap::new();
        tables.insert((0, 0), LinearMap::from_rows(vec![vec![1.0]]));
        tables.insert((0, 1), LinearMap::from_rows(vec![vec![1.0]]));
        tables.insert((1, 0), LinearMap::from_rows(vec![vec![1.0]]));
        tables.insert((1, 1), LinearMap::from_rows(vec![vec![0.0]]));

        let dga = DGA::new(gvs, HashMap::new(), MultiplicationTable::new(tables));
        let h = Cohomology::from_dga(&dga);
        assert_eq!(h.betti_number(0), 1);
        assert_eq!(h.betti_number(1), 1);
    }

    #[test]
    fn test_poincare_polynomial() {
        let mut dims = HashMap::new();
        dims.insert(0, 2);
        dims.insert(1, 3);
        dims.insert(2, 1);
        let cc = ChainComplex::new(dims, HashMap::new());
        let h = Cohomology::from_chain_complex(&cc);
        let pp = h.poincare_polynomial();
        assert_eq!(pp[&0], 2);
        assert_eq!(pp[&1], 3);
        assert_eq!(pp[&2], 1);
    }

    #[test]
    fn test_kernel_basis_all_zero() {
        let m = LinearMap::zero(2, 2);
        let basis = find_kernel_basis(&m, 2);
        assert_eq!(basis.len(), 2);
    }

    #[test]
    fn test_kernel_basis_rank_1() {
        let m = LinearMap::from_rows(vec![vec![1.0, 0.0]]);
        let basis = find_kernel_basis(&m, 2);
        assert_eq!(basis.len(), 1);
        assert!((basis[0][0]).abs() < 1e-10);
        assert!((basis[0][1] - 1.0).abs() < 1e-10);
    }
}
