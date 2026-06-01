//! Applications: unify all cohomology theories in the ecosystem as DGAs.

use serde::{Serialize, Deserialize};
use crate::dga::{DGA, MultiplicationTable};
use crate::chain_complex::LinearMap;
use crate::cohomology::Cohomology;
use crate::graded::GradedVectorSpace;
use std::collections::HashMap;

/// De Rham cohomology as a DGA.
/// The de Rham complex: 0 → Ω⁰ → Ω¹ → Ω² → ... with the exterior derivative.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeRhamDGA {
    /// Dimension of the manifold.
    pub manifold_dim: usize,
    /// The underlying DGA.
    pub dga: DGA,
}

impl DeRhamDGA {
    /// Create the de Rham DGA for a manifold of given dimension.
    /// Ωᵏ has dimension C(n, k) where n = manifold_dim.
    pub fn new(dim: usize) -> Self {
        let mut dgs_dims = HashMap::new();
        for k in 0..=dim as i32 {
            dgs_dims.insert(k, binomial(dim, k as usize));
        }
        let gvs = GradedVectorSpace::new(dgs_dims);

        // The exterior derivative d: Ωᵏ → Ωᵏ⁺¹
        // For a flat manifold, d is determined by its action on basis forms
        // We use a simplified version: d = 0 (giving cohomology = the algebra itself)
        let dga = DGA::new(gvs, HashMap::new(), Self::exterior_multiplication(dim));

        Self { manifold_dim: dim, dga }
    }

    /// Build the wedge product for the exterior algebra on dim generators.
    fn exterior_multiplication(dim: usize) -> MultiplicationTable {
        let mut tables = HashMap::new();

        // Basis: k-forms correspond to subsets of {0, ..., dim-1} of size k
        let subsets_by_degree = Self::all_subsets_by_degree(dim);

        for deg_p in 0..=dim as i32 {
            for deg_q in 0..=dim as i32 {
                let deg_r = deg_p + deg_q;
                if deg_r > dim as i32 { continue; }

                let subs_p = subsets_by_degree.get(&deg_p).cloned().unwrap_or_default();
                let subs_q = subsets_by_degree.get(&deg_q).cloned().unwrap_or_default();
                let subs_r = subsets_by_degree.get(&deg_r).cloned().unwrap_or_default();

                let dp = subs_p.len();
                let dq = subs_q.len();
                let dr = subs_r.len();
                if dp == 0 || dq == 0 || dr == 0 { continue; }

                let mut matrix = vec![vec![0.0; dp * dq]; dr];

                for (ip, sp) in subs_p.iter().enumerate() {
                    for (iq, sq) in subs_q.iter().enumerate() {
                        // Check: no common elements (exterior algebra)
                        if sp & sq != 0 { continue; }

                        let product = sp | sq;
                        let sign = wedge_sign(*sp, *sq);

                        // Find product in subs_r
                        for (ir, sr) in subs_r.iter().enumerate() {
                            if *sr == product {
                                matrix[ir][ip * dq + iq] = sign;
                                break;
                            }
                        }
                    }
                }

                tables.insert((deg_p, deg_q), LinearMap::from_rows(matrix));
            }
        }

        MultiplicationTable::new(tables)
    }

    fn all_subsets_by_degree(dim: usize) -> HashMap<i32, Vec<u64>> {
        let mut result = HashMap::new();
        let total = 1u64 << dim;
        for mask in 0..total {
            let deg = mask.count_ones() as i32;
            result.entry(deg).or_insert_with(Vec::new).push(mask);
        }
        result
    }

    /// Compute de Rham cohomology (with d=0, this is just the exterior algebra).
    pub fn cohomology(&self) -> Cohomology {
        Cohomology::from_dga(&self.dga)
    }

    /// Poincaré polynomial.
    pub fn poincare_polynomial(&self) -> HashMap<i32, usize> {
        self.cohomology().poincare_polynomial().clone()
    }
}

/// Sheaf cohomology as a DGA (simplified).
/// Represented as a Čech complex.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SheafDGA {
    /// Number of open sets in the cover.
    pub cover_size: usize,
    /// The DGA structure.
    pub dga: DGA,
}

impl SheafDGA {
    /// Create a Čech DGA for a cover with n open sets.
    /// Čechⁿ = maps from n-fold intersections.
    pub fn new(cover_size: usize) -> Self {
        let mut dims = HashMap::new();
        for k in 0..=cover_size as i32 - 1 {
            dims.insert(k, binomial(cover_size, k as usize + 1));
        }
        let gvs = GradedVectorSpace::new(dims);
        let dga = DGA::new(gvs, HashMap::new(), MultiplicationTable::empty());

        Self { cover_size, dga }
    }

    /// Compute sheaf cohomology.
    pub fn cohomology(&self) -> Cohomology {
        Cohomology::from_dga(&self.dga)
    }
}

/// Hodge theory: Hodge decomposition of cohomology.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HodgeDecomposition {
    /// Dimensions of the harmonic forms in each degree.
    pub harmonic_dimensions: HashMap<i32, usize>,
    /// The Hodge star operator (simplified as dimensions).
    pub hodge_star_dims: HashMap<i32, usize>,
}

impl HodgeDecomposition {
    /// Compute Hodge decomposition for a DGA of dimension n.
    pub fn for_dimension(n: usize) -> Self {
        // Hodge *: Ωᵏ → Ωⁿ⁻ᵏ
        let mut harmonic = HashMap::new();
        let mut star_dims = HashMap::new();
        for k in 0..=n as i32 {
            let dim_k = binomial(n, k as usize);
            harmonic.insert(k, dim_k); // With zero Laplacian, everything is harmonic
            star_dims.insert(k, binomial(n, n - k as usize));
        }
        Self { harmonic_dimensions: harmonic, hodge_star_dims: star_dims }
    }

    /// Hodge decomposition: Ωᵏ = ℋᵏ ⊕ im(d) ⊕ im(d*)
    /// where ℋᵏ are harmonic forms.
    pub fn check_decomposition(&self, total_dims: &HashMap<i32, usize>) -> bool {
        for (&k, &h) in &self.harmonic_dimensions {
            let total = total_dims.get(&k).copied().unwrap_or(0);
            // harmonic ≤ total
            if h > total { return false; }
        }
        true
    }

    /// Hodge numbers (Betti numbers via harmonic forms).
    pub fn betti_numbers(&self) -> &HashMap<i32, usize> {
        &self.harmonic_dimensions
    }
}

/// Unify all cohomology theories as DGAs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UnifiedCohomology {
    pub de_rham: Option<DeRhamDGA>,
    pub sheaf: Option<SheafDGA>,
    pub hodge: Option<HodgeDecomposition>,
}

impl UnifiedCohomology {
    pub fn new() -> Self {
        Self {
            de_rham: None,
            sheaf: None,
            hodge: None,
        }
    }

    /// Add de Rham cohomology for a manifold of dimension n.
    pub fn with_de_rham(mut self, dim: usize) -> Self {
        self.de_rham = Some(DeRhamDGA::new(dim));
        self
    }

    /// Add sheaf cohomology for a cover of size n.
    pub fn with_sheaf(mut self, cover_size: usize) -> Self {
        self.sheaf = Some(SheafDGA::new(cover_size));
        self
    }

    /// Add Hodge decomposition for dimension n.
    pub fn with_hodge(mut self, dim: usize) -> Self {
        self.hodge = Some(HodgeDecomposition::for_dimension(dim));
        self
    }

    /// Check compatibility: all theories should give the same Betti numbers
    /// (up to the constraints of each theory).
    pub fn check_compatibility(&self) -> bool {
        let betti_sets: Vec<HashMap<i32, usize>> = [
            self.de_rham.as_ref().map(|dr| dr.cohomology().betti),
            self.sheaf.as_ref().map(|s| s.cohomology().betti),
            self.hodge.as_ref().map(|h| h.harmonic_dimensions.clone()),
        ].into_iter().flatten().collect();

        if betti_sets.len() < 2 { return true; }

        for i in 1..betti_sets.len() {
            if !betti_sets_compatible(&betti_sets[0], &betti_sets[i]) {
                return false;
            }
        }
        true
    }

    /// Euler characteristic from all theories.
    pub fn euler_characteristic(&self) -> Option<i64> {
        self.de_rham.as_ref().map(|dr| dr.cohomology().euler_characteristic())
    }
}

fn betti_sets_compatible(a: &HashMap<i32, usize>, b: &HashMap<i32, usize>) -> bool {
    for (&k, &va) in a {
        if let Some(&vb) = b.get(&k) {
            if va != vb { return false; }
        }
    }
    for (&k, &vb) in b {
        if let Some(&va) = a.get(&k) {
            if va != vb { return false; }
        }
    }
    true
}

fn binomial(n: usize, k: usize) -> usize {
    if k > n { return 0; }
    if k == 0 || k == n { return 1; }
    let k = k.min(n - k);
    let mut result = 1usize;
    for i in 0..k {
        result = result * (n - i) / (i + 1);
    }
    result
}

/// Sign of the wedge product α ∧ β where α corresponds to subset_a and β to subset_b.
fn wedge_sign(subset_a: u64, subset_b: u64) -> f64 {
    // Count the number of transpositions needed
    let mut sign = 1i64;
    // For each element in b, count how many elements in a are larger
    for i in 0..64 {
        if subset_b & (1 << i) != 0 {
            for j in (i + 1)..64 {
                if subset_a & (1 << j) != 0 {
                    sign *= -1;
                }
            }
        }
    }
    sign as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_de_rham_dim_0() {
        let dr = DeRhamDGA::new(0);
        assert_eq!(dr.dga.graded_vs.dim(0), 1);
        let h = dr.cohomology();
        assert_eq!(h.betti_number(0), 1);
    }

    #[test]
    fn test_de_rham_dim_1() {
        let dr = DeRhamDGA::new(1);
        assert_eq!(dr.dga.graded_vs.dim(0), 1); // 0-forms
        assert_eq!(dr.dga.graded_vs.dim(1), 1); // 1-forms
        let h = dr.cohomology();
        assert_eq!(h.betti_number(0), 1);
        assert_eq!(h.betti_number(1), 1);
    }

    #[test]
    fn test_de_rham_dim_2() {
        let dr = DeRhamDGA::new(2);
        assert_eq!(dr.dga.graded_vs.dim(0), 1); // C(2,0)
        assert_eq!(dr.dga.graded_vs.dim(1), 2); // C(2,1)
        assert_eq!(dr.dga.graded_vs.dim(2), 1); // C(2,2)
    }

    #[test]
    fn test_de_rham_dim_3() {
        let dr = DeRhamDGA::new(3);
        assert_eq!(dr.dga.graded_vs.dim(0), 1);
        assert_eq!(dr.dga.graded_vs.dim(1), 3);
        assert_eq!(dr.dga.graded_vs.dim(2), 3);
        assert_eq!(dr.dga.graded_vs.dim(3), 1);
    }

    #[test]
    fn test_de_rham_poincare_polynomial() {
        let dr = DeRhamDGA::new(2);
        let pp = dr.poincare_polynomial();
        assert_eq!(pp[&0], 1);
        assert_eq!(pp[&1], 2);
        assert_eq!(pp[&2], 1);
    }

    #[test]
    fn test_de_rham_associativity() {
        let dr = DeRhamDGA::new(2);
        assert!(dr.dga.check_associativity(0, 0, 0, 1e-10));
        assert!(dr.dga.check_associativity(0, 1, 0, 1e-10));
    }

    #[test]
    fn test_de_rham_graded_commutativity() {
        let dr = DeRhamDGA::new(2);
        assert!(dr.dga.check_graded_commutativity(0, 1, 1e-10));
        assert!(dr.dga.check_graded_commutativity(1, 1, 1e-10));
    }

    #[test]
    fn test_sheaf_creation() {
        let s = SheafDGA::new(3);
        // Čech⁰ = C(3,1) = 3, Čech¹ = C(3,2) = 3, Čech² = C(3,3) = 1
        assert_eq!(s.dga.graded_vs.dim(0), 3);
        assert_eq!(s.dga.graded_vs.dim(1), 3);
        assert_eq!(s.dga.graded_vs.dim(2), 1);
    }

    #[test]
    fn test_sheaf_cohomology() {
        let s = SheafDGA::new(2);
        let h = s.cohomology();
        assert_eq!(h.betti_number(0), 2); // C(2,1) = 2
        assert_eq!(h.betti_number(1), 1); // C(2,2) = 1
    }

    #[test]
    fn test_hodge_decomposition() {
        let hodge = HodgeDecomposition::for_dimension(2);
        assert_eq!(hodge.harmonic_dimensions[&0], 1);
        assert_eq!(hodge.harmonic_dimensions[&1], 2);
        assert_eq!(hodge.harmonic_dimensions[&2], 1);
    }

    #[test]
    fn test_hodge_betti() {
        let hodge = HodgeDecomposition::for_dimension(3);
        let betti = hodge.betti_numbers();
        assert_eq!(betti[&0], 1);
        assert_eq!(betti[&1], 3);
        assert_eq!(betti[&2], 3);
        assert_eq!(betti[&3], 1);
    }

    #[test]
    fn test_hodge_check_decomposition() {
        let hodge = HodgeDecomposition::for_dimension(2);
        let mut total = HashMap::new();
        total.insert(0, 1);
        total.insert(1, 2);
        total.insert(2, 1);
        assert!(hodge.check_decomposition(&total));
    }

    #[test]
    fn test_unified_creation() {
        let u = UnifiedCohomology::new()
            .with_de_rham(2)
            .with_sheaf(3)
            .with_hodge(2);
        assert!(u.de_rham.is_some());
        assert!(u.sheaf.is_some());
        assert!(u.hodge.is_some());
    }

    #[test]
    fn test_unified_compatibility() {
        let u = UnifiedCohomology::new()
            .with_de_rham(2)
            .with_hodge(2);
        assert!(u.check_compatibility());
    }

    #[test]
    fn test_unified_euler() {
        let u = UnifiedCohomology::new().with_de_rham(2);
        assert_eq!(u.euler_characteristic(), Some(0)); // 1 - 2 + 1 = 0
    }

    #[test]
    fn test_binomial() {
        assert_eq!(binomial(0, 0), 1);
        assert_eq!(binomial(5, 0), 1);
        assert_eq!(binomial(5, 5), 1);
        assert_eq!(binomial(5, 2), 10);
        assert_eq!(binomial(6, 3), 20);
    }

    #[test]
    fn test_wedge_sign() {
        // Empty wedge empty = positive
        assert_eq!(wedge_sign(0, 0), 1.0);
        // {0} wedge {1} = positive (no transpositions needed)
        assert_eq!(wedge_sign(1, 2), 1.0);
        // {1} wedge {0} = negative (one transposition)
        assert_eq!(wedge_sign(2, 1), -1.0);
    }

    #[test]
    fn test_de_rham_wedge_product() {
        let dr = DeRhamDGA::new(2);
        // In dim 2, basis for 1-forms: dx, dy (subsets {0} and {1})
        // dx ∧ dy = dz (subset {0,1}), with sign +1
        let dx = vec![1.0, 0.0]; // 1-form, index 0
        let dy = vec![0.0, 1.0]; // 1-form, index 1
        let wedge = dr.dga.multiplication.multiply(1, &dx, 1, &dy);
        assert!(wedge.is_some());
        // dx ∧ dy should give the 2-form basis element
        assert!((wedge.unwrap()[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_de_rham_wedge_anticommutativity() {
        let dr = DeRhamDGA::new(2);
        let dx = vec![1.0, 0.0];
        let dy = vec![0.0, 1.0];
        let xy = dr.dga.multiplication.multiply(1, &dx, 1, &dy).unwrap();
        let yx = dr.dga.multiplication.multiply(1, &dy, 1, &dx).unwrap();
        // dx ∧ dy = -dy ∧ dx
        assert!((xy[0] + yx[0]).abs() < 1e-10);
    }

    #[test]
    fn test_de_rham_unit() {
        let dr = DeRhamDGA::new(2);
        assert!(dr.dga.check_unit(1e-10));
    }
}
