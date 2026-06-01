//! Differential Graded Algebra (DGA): graded algebra + differential d with d²=0 and Leibniz rule.

use serde::{Serialize, Deserialize};
use crate::graded::{GradedVectorSpace, koszul_sign};
use crate::chain_complex::{ChainComplex, LinearMap};
use std::collections::HashMap;

/// Multiplication table for a graded algebra.
/// For each degree (p, q) → r with p + q = r, stores the structure constants.
/// mult[p][q] = matrix of size (dim Vʳ × (dim Vᵖ * dim Vᵠ))
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiplicationTable {
    /// Map from (degree_p, degree_q) to the multiplication matrix.
    /// Matrix is (dim V^{p+q} × (dim Vᵖ * dim Vᵠ)).
    pub tables: HashMap<(i32, i32), LinearMap>,
}

impl MultiplicationTable {
    pub fn new(tables: HashMap<(i32, i32), LinearMap>) -> Self {
        Self { tables }
    }

    pub fn empty() -> Self {
        Self { tables: HashMap::new() }
    }

    /// Multiply two homogeneous vectors.
    pub fn multiply(&self, a_deg: i32, a: &[f64], b_deg: i32, b: &[f64]) -> Option<Vec<f64>> {
        let key = (a_deg, b_deg);
        let table = self.tables.get(&key)?;
        // Form tensor product a ⊗ b
        let mut tensor = vec![0.0; a.len() * b.len()];
        for i in 0..a.len() {
            for j in 0..b.len() {
                tensor[i * b.len() + j] = a[i] * b[j];
            }
        }
        Some(table.apply(&tensor))
    }
}

/// A differential graded algebra.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DGA {
    /// Underlying graded vector space.
    pub graded_vs: GradedVectorSpace,
    /// The differential d: Vᵏ → Vᵏ⁺¹.
    pub differential: HashMap<i32, LinearMap>,
    /// The multiplication table.
    pub multiplication: MultiplicationTable,
}

impl DGA {
    pub fn new(
        graded_vs: GradedVectorSpace,
        differential: HashMap<i32, LinearMap>,
        multiplication: MultiplicationTable,
    ) -> Self {
        Self { graded_vs, differential, multiplication }
    }

    /// Check d² = 0.
    pub fn check_d_squared_zero(&self, tol: f64) -> bool {
        for (&k, d_k) in &self.differential {
            if let Some(d_k1) = self.differential.get(&(k + 1)) {
                let d2 = d_k1.compose(d_k);
                if !d2.is_zero(tol) {
                    return false;
                }
            }
        }
        true
    }

    /// Apply the differential to a homogeneous vector of degree k.
    pub fn d(&self, k: i32, v: &[f64]) -> Option<Vec<f64>> {
        let d_k = self.differential.get(&k)?;
        Some(d_k.apply(v))
    }

    /// Check the (graded) Leibniz rule: d(a·b) = da·b + (-1)^|a| a·db
    /// for a specific pair of degrees.
    pub fn check_leibniz(&self, deg_a: i32, deg_b: i32, tol: f64) -> bool {
        let dim_a = self.graded_vs.dim(deg_a);
        let dim_b = self.graded_vs.dim(deg_b);
        let deg_ab = deg_a + deg_b;

        // For each basis pair, check Leibniz
        for i in 0..dim_a {
            for j in 0..dim_b {
                let a = {
                    let mut v = vec![0.0; dim_a];
                    v[i] = 1.0;
                    v
                };
                let b = {
                    let mut v = vec![0.0; dim_b];
                    v[j] = 1.0;
                    v
                };

                // a · b
                let ab = match self.multiplication.multiply(deg_a, &a, deg_b, &b) {
                    Some(x) => x,
                    None => continue,
                };

                // d(a · b)
                let d_ab = match self.d(deg_ab, &ab) {
                    Some(x) => x,
                    None => continue,
                };

                // da · b
                let da = match self.d(deg_a, &a) {
                    Some(x) => x,
                    None => continue,
                };
                let da_b = match self.multiplication.multiply(deg_a + 1, &da, deg_b, &b) {
                    Some(x) => x,
                    None => continue,
                };

                // a · db
                let db = match self.d(deg_b, &b) {
                    Some(x) => x,
                    None => continue,
                };
                let a_db = match self.multiplication.multiply(deg_a, &a, deg_b + 1, &db) {
                    Some(x) => x,
                    None => continue,
                };

                // (-1)^|a| a · db
                let sign = koszul_sign(deg_a, 1);
                let rhs_dim = da_b.len().max(a_db.len());
                let mut rhs = vec![0.0; rhs_dim];
                for k in 0..da_b.len() {
                    rhs[k] += da_b[k];
                }
                for k in 0..a_db.len() {
                    rhs[k] += sign * a_db[k];
                }

                if d_ab.len() != rhs.len() { return false; }
                for k in 0..d_ab.len() {
                    if (d_ab[k] - rhs[k]).abs() > tol {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Check associativity: (a·b)·c = a·(b·c) for all degree triples.
    pub fn check_associativity(&self, deg_a: i32, deg_b: i32, deg_c: i32, tol: f64) -> bool {
        let dim_a = self.graded_vs.dim(deg_a);
        let dim_b = self.graded_vs.dim(deg_b);
        let dim_c = self.graded_vs.dim(deg_c);

        for i in 0..dim_a {
            for j in 0..dim_b {
                for k in 0..dim_c {
                    let a = {
                        let mut v = vec![0.0; dim_a]; v[i] = 1.0; v
                    };
                    let b = {
                        let mut v = vec![0.0; dim_b]; v[j] = 1.0; v
                    };
                    let c = {
                        let mut v = vec![0.0; dim_c]; v[k] = 1.0; v
                    };

                    // (a·b)·c
                    let ab = match self.multiplication.multiply(deg_a, &a, deg_b, &b) {
                        Some(x) => x,
                        None => continue,
                    };
                    let ab_c = match self.multiplication.multiply(deg_a + deg_b, &ab, deg_c, &c) {
                        Some(x) => x,
                        None => continue,
                    };

                    // a·(b·c)
                    let bc = match self.multiplication.multiply(deg_b, &b, deg_c, &c) {
                        Some(x) => x,
                        None => continue,
                    };
                    let a_bc = match self.multiplication.multiply(deg_a, &a, deg_b + deg_c, &bc) {
                        Some(x) => x,
                        None => continue,
                    };

                    if ab_c.len() != a_bc.len() { return false; }
                    for m in 0..ab_c.len() {
                        if (ab_c[m] - a_bc[m]).abs() > tol {
                            return false;
                        }
                    }
                }
            }
        }
        true
    }

    /// Extract the underlying chain complex.
    pub fn as_chain_complex(&self) -> ChainComplex {
        ChainComplex::new(self.graded_vs.dimensions.clone(), self.differential.clone())
    }

    /// Check graded commutativity: a·b = (-1)^{|a||b|} b·a
    pub fn check_graded_commutativity(&self, deg_a: i32, deg_b: i32, tol: f64) -> bool {
        let dim_a = self.graded_vs.dim(deg_a);
        let dim_b = self.graded_vs.dim(deg_b);
        let sign = koszul_sign(deg_a, deg_b);

        for i in 0..dim_a {
            for j in 0..dim_b {
                let a = { let mut v = vec![0.0; dim_a]; v[i] = 1.0; v };
                let b = { let mut v = vec![0.0; dim_b]; v[j] = 1.0; v };

                let ab = match self.multiplication.multiply(deg_a, &a, deg_b, &b) {
                    Some(x) => x, None => continue,
                };
                let ba = match self.multiplication.multiply(deg_b, &b, deg_a, &a) {
                    Some(x) => x, None => continue,
                };

                if ab.len() != ba.len() { return false; }
                for k in 0..ab.len() {
                    if (ab[k] - sign * ba[k]).abs() > tol {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// The unit element (if it exists in degree 0).
    pub fn check_unit(&self, tol: f64) -> bool {
        let dim0 = self.graded_vs.dim(0);
        if dim0 == 0 { return false; }

        // Try to find a unit in degree 0
        'outer: for i in 0..dim0 {
            let e = { let mut v = vec![0.0; dim0]; v[i] = 1.0; v };

            for &deg in self.graded_vs.degrees().iter() {
                let dim_k = self.graded_vs.dim(deg);
                for j in 0..dim_k {
                    let x = { let mut v = vec![0.0; dim_k]; v[j] = 1.0; v };

                    // e · x should equal x
                    if let Some(ex) = self.multiplication.multiply(0, &e, deg, &x) {
                        if ex.len() != x.len() { continue 'outer; }
                        for k in 0..x.len() {
                            if (ex[k] - x[k]).abs() > tol { continue 'outer; }
                        }
                    }
                    // x · e should equal x
                    if let Some(xe) = self.multiplication.multiply(deg, &x, 0, &e) {
                        if xe.len() != x.len() { continue 'outer; }
                        for k in 0..x.len() {
                            if (xe[k] - x[k]).abs() > tol { continue 'outer; }
                        }
                    }
                }
            }
            return true;
        }
        false
    }
}

/// A morphism of DGAs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DGAMorphism {
    /// Maps fᵏ: Vᵏ → Wᵏ for each degree.
    pub maps: HashMap<i32, LinearMap>,
}

impl DGAMorphism {
    pub fn new(maps: HashMap<i32, LinearMap>) -> Self {
        Self { maps }
    }

    /// Check that the morphism commutes with differentials.
    pub fn check_chain_map(&self, source: &DGA, target: &DGA, tol: f64) -> bool {
        for (&k, f_k) in &self.maps {
            let d_s = source.differential.get(&k);
            let d_t = target.differential.get(&k);
            let f_k1 = self.maps.get(&(k + 1));

            if let (Some(d_s), Some(f_k1)) = (d_s, f_k1) {
                // f_{k+1} ∘ d_k^s vs d_k^t ∘ f_k
                let lhs = f_k1.compose(d_s);
                if let Some(d_t) = d_t {
                    let rhs = d_t.compose(f_k);
                    let diff = lhs.add(&rhs.scale(-1.0));
                    if let Some(ref d) = diff {
                        if !d.is_zero(tol) { return false; }
                    }
                }
            }
        }
        true
    }

    /// Check that the morphism respects multiplication.
    pub fn check_algebra_map(&self, source: &DGA, target: &DGA, deg_a: i32, deg_b: i32, tol: f64) -> bool {
        let dim_a = source.graded_vs.dim(deg_a);
        let dim_b = source.graded_vs.dim(deg_b);
        let f_a = match self.maps.get(&deg_a) {
            Some(m) => m,
            None => return true, // No map means zero
        };
        let f_b = match self.maps.get(&deg_b) {
            Some(m) => m,
            None => return true,
        };
        let f_ab = match self.maps.get(&(deg_a + deg_b)) {
            Some(m) => m,
            None => return true,
        };

        for i in 0..dim_a {
            for j in 0..dim_b {
                let a = { let mut v = vec![0.0; dim_a]; v[i] = 1.0; v };
                let b = { let mut v = vec![0.0; dim_b]; v[j] = 1.0; v };

                // f(a) · f(b) in target
                let fa = f_a.apply(&a);
                let fb = f_b.apply(&b);
                let fa_fb = match target.multiplication.multiply(deg_a, &fa, deg_b, &fb) {
                    Some(x) => x,
                    None => continue,
                };

                // f(a · b) in target
                let ab = match source.multiplication.multiply(deg_a, &a, deg_b, &b) {
                    Some(x) => x,
                    None => continue,
                };
                let f_ab_val = f_ab.apply(&ab);

                if fa_fb.len() != f_ab_val.len() { return false; }
                for k in 0..fa_fb.len() {
                    if (fa_fb[k] - f_ab_val[k]).abs() > tol {
                        return false;
                    }
                }
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_trivial_dga() -> DGA {
        // Exterior algebra on one generator x of degree 1
        // Basis: {1 (degree 0), x (degree 1)}
        // d = 0
        // 1·1 = 1, 1·x = x, x·1 = x, x·x = 0
        let mut dims = HashMap::new();
        dims.insert(0, 1);
        dims.insert(1, 1);
        let gvs = GradedVectorSpace::new(dims);

        let mult_0_0 = LinearMap::from_rows(vec![vec![1.0]]); // 1·1 = 1
        let mult_0_1 = LinearMap::from_rows(vec![vec![1.0]]); // 1·x = x
        let mult_1_0 = LinearMap::from_rows(vec![vec![1.0]]); // x·1 = x
        let mult_1_1 = LinearMap::from_rows(vec![vec![0.0]]); // x·x = 0

        let mut tables = HashMap::new();
        tables.insert((0, 0), mult_0_0);
        tables.insert((0, 1), mult_0_1);
        tables.insert((1, 0), mult_1_0);
        tables.insert((1, 1), mult_1_1);

        DGA::new(gvs, HashMap::new(), MultiplicationTable::new(tables))
    }

    #[test]
    fn test_trivial_dga_d_squared() {
        let dga = make_trivial_dga();
        assert!(dga.check_d_squared_zero(1e-10));
    }

    #[test]
    fn test_trivial_dga_associativity() {
        let dga = make_trivial_dga();
        assert!(dga.check_associativity(0, 0, 0, 1e-10));
        assert!(dga.check_associativity(0, 0, 1, 1e-10));
        assert!(dga.check_associativity(0, 1, 0, 1e-10));
        assert!(dga.check_associativity(1, 0, 0, 1e-10));
        assert!(dga.check_associativity(1, 1, 0, 1e-10));
        assert!(dga.check_associativity(1, 0, 1, 1e-10));
        assert!(dga.check_associativity(0, 1, 1, 1e-10));
    }

    #[test]
    fn test_trivial_dga_graded_commutativity() {
        let dga = make_trivial_dga();
        assert!(dga.check_graded_commutativity(0, 0, 1e-10));
        assert!(dga.check_graded_commutativity(0, 1, 1e-10));
        assert!(dga.check_graded_commutativity(1, 1, 1e-10));
    }

    #[test]
    fn test_trivial_dga_unit() {
        let dga = make_trivial_dga();
        assert!(dga.check_unit(1e-10));
    }

    #[test]
    fn test_trivial_dga_leibniz() {
        let dga = make_trivial_dga();
        // d = 0, so Leibniz is trivially satisfied
        assert!(dga.check_leibniz(0, 0, 1e-10));
        assert!(dga.check_leibniz(0, 1, 1e-10));
        assert!(dga.check_leibniz(1, 0, 1e-10));
        assert!(dga.check_leibniz(1, 1, 1e-10));
    }

    #[test]
    fn test_trivial_dga_multiply() {
        let dga = make_trivial_dga();
        let one = vec![1.0];
        let x = vec![1.0];

        let ox = dga.multiplication.multiply(0, &one, 1, &x).unwrap();
        assert!((ox[0] - 1.0).abs() < 1e-10);

        let xx = dga.multiplication.multiply(1, &x, 1, &x).unwrap();
        assert!((xx[0]).abs() < 1e-10);
    }

    #[test]
    fn test_dga_as_chain_complex() {
        let dga = make_trivial_dga();
        let cc = dga.as_chain_complex();
        assert_eq!(cc.dim(0), 1);
        assert_eq!(cc.dim(1), 1);
    }

    #[test]
    fn test_dga_morphism_chain_map() {
        let dga = make_trivial_dga();
        let id_map = crate::ChainMap::identity(&dga.as_chain_complex());
        let morphism = DGAMorphism::new(id_map.maps);
        assert!(morphism.check_chain_map(&dga, &dga, 1e-10));
    }

    #[test]
    fn test_dga_morphism_identity_preserves_mult() {
        let dga = make_trivial_dga();
        let id_map = crate::ChainMap::identity(&dga.as_chain_complex());
        let morphism = DGAMorphism::new(id_map.maps);
        assert!(morphism.check_algebra_map(&dga, &dga, 0, 0, 1e-10));
        assert!(morphism.check_algebra_map(&dga, &dga, 0, 1, 1e-10));
        assert!(morphism.check_algebra_map(&dga, &dga, 1, 1, 1e-10));
    }

    #[test]
    fn test_dga_with_nonzero_differential() {
        // V⁰ (dim 1) --d--> V¹ (dim 1), d = [[1]], d² = 0 (no further differential)
        let mut dims = HashMap::new();
        dims.insert(0, 1);
        dims.insert(1, 1);
        let gvs = GradedVectorSpace::new(dims);

        let mut diffs = HashMap::new();
        diffs.insert(0, LinearMap::from_rows(vec![vec![1.0]]));

        // Trivial multiplication for testing d²=0
        let dga = DGA::new(gvs, diffs, MultiplicationTable::empty());
        assert!(dga.check_d_squared_zero(1e-10));
    }
}
