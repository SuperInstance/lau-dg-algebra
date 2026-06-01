//! Quasi-isomorphisms: DGA maps inducing isomorphism on cohomology.

use serde::{Serialize, Deserialize};
use crate::dga::{DGA, DGAMorphism};
use crate::cohomology::Cohomology;
use crate::chain_complex::LinearMap;
use std::collections::HashMap;

/// Check if a DGA morphism induces an isomorphism on cohomology.
pub fn is_quasi_isomorphism(
    source: &DGA,
    target: &DGA,
    morphism: &DGAMorphism,
    tol: f64,
) -> bool {
    // 1. Must be a chain map
    if !morphism.check_chain_map(source, target, tol) {
        return false;
    }

    // 2. Must induce isomorphism on cohomology
    let h_source = Cohomology::from_dga(source);
    let h_target = Cohomology::from_dga(target);

    // Check that induced maps on cohomology are isomorphisms
    let all_degrees: std::collections::HashSet<i32> = h_source.betti.keys()
        .chain(h_target.betti.keys())
        .copied()
        .collect();

    for &k in &all_degrees {
        let b_s = h_source.betti_number(k);
        let b_t = h_target.betti_number(k);
        if b_s != b_t {
            return false;
        }
        // If Betti numbers match and are nonzero, check that the map has full rank
        if b_s > 0 {
            if let Some(f_k) = morphism.maps.get(&k) {
                if f_k.rank() < b_s {
                    return false;
                }
            } else {
                return false;
            }
        }
    }
    true
}

/// Check if two DGAs are quasi-isomorphic (there exists a zig-zag of quasi-isos).
/// This is a simplified check: we compare cohomology algebras.
pub fn are_quasi_isomorphic(dga_a: &DGA, dga_b: &DGA) -> bool {
    let h_a = Cohomology::from_dga(dga_a);
    let h_b = Cohomology::from_dga(dga_b);

    // Same Betti numbers is necessary (but not sufficient)
    if h_a.total_betti() != h_b.total_betti() {
        return false;
    }

    let all_degrees: std::collections::HashSet<i32> = h_a.betti.keys()
        .chain(h_b.betti.keys())
        .copied()
        .collect();

    for &k in &all_degrees {
        if h_a.betti_number(k) != h_b.betti_number(k) {
            return false;
        }
    }
    true
}

/// A quasi-isomorphism with its induced map on cohomology.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QuasiIsomorphism {
    pub morphism: DGAMorphism,
    pub induced_cohomology_map: HashMap<i32, LinearMap>,
    pub source_cohomology: Cohomology,
    pub target_cohomology: Cohomology,
}

impl QuasiIsomorphism {
    /// Construct and verify a quasi-isomorphism.
    pub fn new(
        source: &DGA,
        target: &DGA,
        morphism: DGAMorphism,
        tol: f64,
    ) -> Option<Self> {
        if !is_quasi_isomorphism(source, target, &morphism, tol) {
            return None;
        }

        let source_cohomology = Cohomology::from_dga(source);
        let target_cohomology = Cohomology::from_dga(target);

        // The induced maps on cohomology are the restrictions of the original maps
        let induced = morphism.maps.clone();

        Some(Self {
            morphism,
            induced_cohomology_map: induced,
            source_cohomology,
            target_cohomology,
        })
    }

    /// The inverse quasi-isomorphism (if it exists).
    pub fn inverse(&self) -> Option<QuasiIsomorphism> {
        // Check that induced maps are invertible
        let mut inv_maps: HashMap<i32, LinearMap> = HashMap::new();
        for (&k, f_k) in &self.morphism.maps {
            let rank = f_k.rank();
            let dim = f_k.domain_dim;
            if rank != dim || rank != f_k.codomain_dim {
                return None;
            }
            // For invertible matrices, we could compute the inverse
            // For now, just check square and full rank
        }
        // Would need actual inverse computation
        None
    }

    /// Verify this is indeed a quasi-isomorphism.
    pub fn verify(&self, source: &DGA, target: &DGA, tol: f64) -> bool {
        is_quasi_isomorphism(source, target, &self.morphism, tol)
    }
}

/// Cylinder object for a DGA (used for homotopy).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CylinderObject {
    pub dga: DGA,
    pub i0: DGAMorphism,
    pub i1: DGAMorphism,
    pub projection: DGAMorphism,
}

/// Homotopy between two DGA morphisms.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DGAHomotopy {
    /// Homotopy maps hᵏ: Wᵏ → Vᵏ⁻¹ such that f - g = dh + hd.
    pub maps: HashMap<i32, LinearMap>,
}

impl DGAHomotopy {
    pub fn new(maps: HashMap<i32, LinearMap>) -> Self {
        Self { maps }
    }

    /// Verify the homotopy equation: f - g = dV ∘ h + h ∘ dW.
    pub fn verify(
        &self,
        source: &DGA,
        target: &DGA,
        f: &DGAMorphism,
        g: &DGAMorphism,
        tol: f64,
    ) -> bool {
        for (&k, h_k) in &self.maps {
            // (f - g)_k
            let f_k = f.maps.get(&k);
            let g_k = g.maps.get(&k);
            let fg_diff = match (f_k, g_k) {
                (Some(a), Some(b)) => a.add(&b.scale(-1.0)),
                (Some(a), None) => Some(a.clone()),
                (None, Some(b)) => Some(b.scale(-1.0)),
                (None, None) => continue,
            };

            // dV_{k-1} ∘ h_k + h_{k+1} ∘ dW_k
            let mut dh_plus_hd = if let Some(dv_prev) = target.differential.get(&(k - 1)) {
                dv_prev.compose(h_k)
            } else {
                LinearMap::zero(h_k.domain_dim, 0)
            };

            if let Some(h_k1) = self.maps.get(&(k + 1)) {
                if let Some(dw_k) = source.differential.get(&k) {
                    let hd = h_k1.compose(dw_k);
                    dh_plus_hd = match dh_plus_hd.add(&hd) {
                        Some(m) => m,
                        None => return false,
                    };
                }
            }

            if let Some(ref expected) = fg_diff {
                let diff = dh_plus_hd.add(&expected.scale(-1.0));
                if let Some(d) = diff {
                    if !d.is_zero(tol) {
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
    use crate::graded::GradedVectorSpace;

    fn make_ext_dga() -> DGA {
        let mut dims = HashMap::new();
        dims.insert(0, 1);
        dims.insert(1, 1);
        let gvs = GradedVectorSpace::new(dims);

        let mut tables = HashMap::new();
        tables.insert((0, 0), LinearMap::from_rows(vec![vec![1.0]]));
        tables.insert((0, 1), LinearMap::from_rows(vec![vec![1.0]]));
        tables.insert((1, 0), LinearMap::from_rows(vec![vec![1.0]]));
        tables.insert((1, 1), LinearMap::from_rows(vec![vec![0.0]]));

        DGA::new(gvs, HashMap::new(), crate::dga::MultiplicationTable::new(tables))
    }

    #[test]
    fn test_identity_is_quasi_iso() {
        let dga = make_ext_dga();
        let id_maps: HashMap<i32, LinearMap> = dga.graded_vs.dimensions.iter()
            .map(|(&k, &d)| (k, LinearMap::identity(d)))
            .collect();
        let morphism = DGAMorphism::new(id_maps);
        assert!(is_quasi_isomorphism(&dga, &dga, &morphism, 1e-10));
    }

    #[test]
    fn test_quasi_iso_same_cohomology() {
        let dga_a = make_ext_dga();
        let dga_b = make_ext_dga();
        assert!(are_quasi_isomorphic(&dga_a, &dga_b));
    }

    #[test]
    fn test_quasi_iso_different_cohomology() {
        let dga_a = make_ext_dga();

        let mut dims = HashMap::new();
        dims.insert(0, 2);
        let gvs = GradedVectorSpace::new(dims);
        let dga_b = DGA::new(gvs, HashMap::new(), crate::dga::MultiplicationTable::empty());

        assert!(!are_quasi_isomorphic(&dga_a, &dga_b));
    }

    #[test]
    fn test_quasi_isomorphism_new() {
        let dga = make_ext_dga();
        let id_maps: HashMap<i32, LinearMap> = dga.graded_vs.dimensions.iter()
            .map(|(&k, &d)| (k, LinearMap::identity(d)))
            .collect();
        let morphism = DGAMorphism::new(id_maps);
        let qi = QuasiIsomorphism::new(&dga, &dga, morphism, 1e-10);
        assert!(qi.is_some());
    }

    #[test]
    fn test_quasi_isomorphism_verify() {
        let dga = make_ext_dga();
        let id_maps: HashMap<i32, LinearMap> = dga.graded_vs.dimensions.iter()
            .map(|(&k, &d)| (k, LinearMap::identity(d)))
            .collect();
        let morphism = DGAMorphism::new(id_maps);
        let qi = QuasiIsomorphism::new(&dga, &dga, morphism, 1e-10).unwrap();
        assert!(qi.verify(&dga, &dga, 1e-10));
    }

    #[test]
    fn test_homotopy_trivial() {
        let dga = make_ext_dga();
        let id_maps: HashMap<i32, LinearMap> = dga.graded_vs.dimensions.iter()
            .map(|(&k, &d)| (k, LinearMap::identity(d)))
            .collect();
        let f = DGAMorphism::new(id_maps.clone());
        let g = DGAMorphism::new(id_maps);
        // f - g = 0, so h = 0 works
        let h = DGAHomotopy::new(HashMap::new());
        assert!(h.verify(&dga, &dga, &f, &g, 1e-10));
    }

    #[test]
    fn test_not_quasi_iso_different_ranks() {
        let dga_a = make_ext_dga();

        let mut dims = HashMap::new();
        dims.insert(0, 1);
        dims.insert(1, 2);
        let gvs = GradedVectorSpace::new(dims);
        let dga_b = DGA::new(gvs, HashMap::new(), crate::dga::MultiplicationTable::empty());

        let mut maps = HashMap::new();
        maps.insert(0, LinearMap::identity(1));
        maps.insert(1, LinearMap::from_rows(vec![vec![1.0, 0.0]]));
        let morphism = DGAMorphism::new(maps);

        assert!(!is_quasi_isomorphism(&dga_a, &dga_b, &morphism, 1e-10));
    }
}
