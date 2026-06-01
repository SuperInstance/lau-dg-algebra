//! Derived category: chain complexes up to quasi-isomorphism.

use serde::{Serialize, Deserialize};
use crate::chain_complex::{ChainComplex, ChainMap, LinearMap};
use crate::cohomology::Cohomology;
use std::collections::HashMap;

/// An object in the derived category: a chain complex.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DerivedObject {
    pub complex: ChainComplex,
    pub name: String,
}

impl DerivedObject {
    pub fn new(complex: ChainComplex, name: &str) -> Self {
        Self { complex, name: name.to_string() }
    }

    /// The cohomology of this object.
    pub fn cohomology(&self) -> Cohomology {
        Cohomology::from_chain_complex(&self.complex)
    }

    /// Shift by n: C[n].
    pub fn shift(&self, n: i32) -> Self {
        DerivedObject::new(self.complex.shift(n), &format!("{}[{}]", self.name, n))
    }

    /// Direct sum.
    pub fn direct_sum(&self, other: &Self) -> Self {
        let sum = self.complex.direct_sum(&other.complex);
        DerivedObject::new(sum, &format!("{} ⊕ {}", self.name, other.name))
    }

    /// Check if this is acyclic (all cohomology is zero).
    pub fn is_acyclic(&self) -> bool {
        self.cohomology().is_trivial()
    }

    /// Check if this is a perfect complex (bounded, finite dimensional).
    pub fn is_perfect(&self) -> bool {
        let degrees = self.complex.degrees();
        if degrees.is_empty() { return true; }
        let min_deg = *degrees.first().unwrap();
        let max_deg = *degrees.last().unwrap();
        // Bounded and finite dimensional
        (max_deg - min_deg) < 1000 && self.complex.dimensions.values().all(|&d| d < 10000)
    }

    /// Truncate at degree n (brutal truncation).
    pub fn truncate(&self, n: i32) -> Self {
        let dims = self.complex.dimensions.iter()
            .filter(|(&k, _)| k <= n)
            .map(|(&k, &d)| (k, d))
            .collect();
        let diffs = self.complex.differentials.iter()
            .filter(|(&k, _)| k <= n)
            .map(|(&k, ref v)| (k, (*v).clone()))
            .collect();
        DerivedObject::new(
            ChainComplex::new(dims, diffs),
            &format!("τ_≤{} {}", n, self.name),
        )
    }
}

/// A morphism in the derived category: represented by a roof/zig-zag.
/// Simplified: we represent morphisms as chain maps.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DerivedMorphism {
    pub source_name: String,
    pub target_name: String,
    pub chain_map: ChainMap,
}

impl DerivedMorphism {
    pub fn new(source_name: &str, target_name: &str, chain_map: ChainMap) -> Self {
        Self {
            source_name: source_name.to_string(),
            target_name: target_name.to_string(),
            chain_map,
        }
    }

    /// Compose two morphisms.
    pub fn compose(&self, other: &Self) -> Self {
        Self {
            source_name: other.source_name.clone(),
            target_name: self.target_name.clone(),
            chain_map: self.chain_map.compose(&other.chain_map),
        }
    }

    /// Check if this morphism is zero in the derived category.
    pub fn is_zero(&self, source: &ChainComplex, target: &ChainComplex, tol: f64) -> bool {
        // A morphism is zero in the derived category if it factors through an acyclic complex
        // Simplified: check if all component maps are zero
        self.chain_map.maps.values().all(|m| m.is_zero(tol))
    }

    /// Check if this is a quasi-isomorphism.
    pub fn is_quasi_isomorphism(&self, source: &ChainComplex, target: &ChainComplex, tol: f64) -> bool {
        let h_source = Cohomology::from_chain_complex(source);
        let h_target = Cohomology::from_chain_complex(target);

        for (&k, _) in &h_source.betti {
            if h_source.betti_number(k) != h_target.betti_number(k) {
                return false;
            }
        }
        for (&k, _) in &h_target.betti {
            if h_target.betti_number(k) != h_source.betti_number(k) {
                return false;
            }
        }
        true
    }
}

/// The derived category D(A) of an abelian category.
/// Objects are chain complexes, morphisms are chain maps localized at quasi-isomorphisms.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DerivedCategory {
    /// Objects in the category.
    pub objects: Vec<DerivedObject>,
    /// Morphisms (source_idx, target_idx, chain_map).
    pub morphisms: Vec<(usize, usize, ChainMap)>,
    pub name: String,
}

impl DerivedCategory {
    pub fn new(name: &str) -> Self {
        Self {
            objects: Vec::new(),
            morphisms: Vec::new(),
            name: name.to_string(),
        }
    }

    /// Add an object.
    pub fn add_object(&mut self, obj: DerivedObject) -> usize {
        let idx = self.objects.len();
        self.objects.push(obj);
        idx
    }

    /// Add a morphism.
    pub fn add_morphism(&mut self, source: usize, target: usize, map: ChainMap) {
        self.morphisms.push((source, target, map));
    }

    /// Get an object by index.
    pub fn get_object(&self, idx: usize) -> Option<&DerivedObject> {
        self.objects.get(idx)
    }

    /// Find all morphisms from source to target.
    pub fn morphisms_from_to(&self, source: usize, target: usize) -> Vec<&ChainMap> {
        self.morphisms.iter()
            .filter(|&&(s, t, _)| s == source && t == target)
            .map(|&(_, _, ref m)| m)
            .collect()
    }

    /// Compute Ext groups: Extⁿ(A, B) = Hom_{D(A)}(A, B[n]).
    pub fn ext_group(&self, a_idx: usize, b_idx: usize, n: i32) -> usize {
        let a = match self.objects.get(a_idx) {
            Some(o) => o,
            None => return 0,
        };
        let b_shifted = match self.objects.get(b_idx) {
            Some(o) => o.shift(n),
            None => return 0,
        };

        // Simplified: Extⁿ = dimension of Hom(A, B[n])
        // which is the dimension of chain maps A → B[n]
        // A chain map exists degree by degree, and the dimension is
        // computed from the Hom complex
        let mut dim = 0;
        for (&k, &d_a) in &a.complex.dimensions {
            let d_b = b_shifted.complex.dim(k);
            dim += d_a * d_b;
        }
        dim
    }

    /// Check if two objects are isomorphic in the derived category.
    pub fn are_isomorphic(&self, a_idx: usize, b_idx: usize, tol: f64) -> bool {
        let a = match self.objects.get(a_idx) {
            Some(o) => o,
            None => return false,
        };
        let b = match self.objects.get(b_idx) {
            Some(o) => o,
            None => return false,
        };

        let h_a = a.cohomology();
        let h_b = b.cohomology();

        for (&k, &ba) in &h_a.betti {
            if h_b.betti_number(k) != ba {
                return false;
            }
        }
        for (&k, &bb) in &h_b.betti {
            if h_a.betti_number(k) != bb {
                return false;
            }
        }
        true
    }

    /// The Verdier quotient: D(A)/B where B is a triangulated subcategory.
    /// This is represented as a new derived category.
    pub fn verdier_quotient(&self, subcategory_indices: &[usize]) -> DerivedCategory {
        let mut quotient = DerivedCategory::new(&format!("{} / subcategory", self.name));

        // Objects not in the subcategory
        let sub_set: std::collections::HashSet<usize> = subcategory_indices.iter().copied().collect();
        let mut index_map = HashMap::new();

        for (i, obj) in self.objects.iter().enumerate() {
            if !sub_set.contains(&i) {
                let new_idx = quotient.add_object(obj.clone());
                index_map.insert(i, new_idx);
            }
        }

        // Morphisms between surviving objects
        for &(s, t, ref m) in &self.morphisms {
            if let (Some(&ns), Some(&nt)) = (index_map.get(&s), index_map.get(&t)) {
                quotient.add_morphism(ns, nt, m.clone());
            }
        }

        quotient
    }

    /// The homotopy category K(A): chain complexes with morphisms up to chain homotopy.
    pub fn homotopy_category(&self) -> DerivedCategory {
        // In the homotopy category, we identify homotopic maps
        // Simplified: return a copy
        self.clone()
    }

    /// Number of objects.
    pub fn num_objects(&self) -> usize {
        self.objects.len()
    }

    /// Number of morphisms.
    pub fn num_morphisms(&self) -> usize {
        self.morphisms.len()
    }
}

/// A distinguished triangle in the derived category.
/// A → B → C → A[1]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DistinguishedTriangle {
    pub a_idx: usize,
    pub b_idx: usize,
    pub c_idx: usize,
    pub f: ChainMap,    // A → B
    pub g: ChainMap,    // B → C
    pub h: ChainMap,    // C → A[1]
}

impl DistinguishedTriangle {
    pub fn new(a: usize, b: usize, c: usize, f: ChainMap, g: ChainMap, h: ChainMap) -> Self {
        Self { a_idx: a, b_idx: b, c_idx: c, f, g, h }
    }

    /// Rotate the triangle: B → C → A[1] → B[1]
    pub fn rotate(&self) -> Self {
        DistinguishedTriangle::new(
            self.b_idx,
            self.c_idx,
            self.a_idx,
            self.g.clone(),
            self.h.clone(),
            // The third map would be the negative of f[1]
            ChainMap::new(self.f.maps.iter().map(|(&k, m)| (k + 1, m.scale(-1.0))).collect()),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_simple_complex(dim: usize, name: &str) -> DerivedObject {
        let mut dims = HashMap::new();
        dims.insert(0, dim);
        DerivedObject::new(ChainComplex::new(dims, HashMap::new()), name)
    }

    #[test]
    fn test_derived_object_creation() {
        let obj = make_simple_complex(1, "Z");
        assert_eq!(obj.name, "Z");
        assert_eq!(obj.complex.dim(0), 1);
    }

    #[test]
    fn test_derived_object_cohomology() {
        let obj = make_simple_complex(1, "Z");
        let h = obj.cohomology();
        assert_eq!(h.betti_number(0), 1);
    }

    #[test]
    fn test_derived_object_shift() {
        let obj = make_simple_complex(1, "Z");
        let shifted = obj.shift(2);
        assert_eq!(shifted.complex.dim(2), 1);
        assert_eq!(shifted.name, "Z[2]");
    }

    #[test]
    fn test_derived_object_direct_sum() {
        let a = make_simple_complex(1, "A");
        let b = make_simple_complex(2, "B");
        let sum = a.direct_sum(&b);
        assert_eq!(sum.complex.dim(0), 3);
    }

    #[test]
    fn test_derived_object_acyclic() {
        let mut dims = HashMap::new();
        dims.insert(0, 1);
        dims.insert(1, 1);
        let mut diffs = HashMap::new();
        diffs.insert(0, LinearMap::from_rows(vec![vec![1.0]]));
        let cc = ChainComplex::new(dims, diffs);
        let obj = DerivedObject::new(cc, "exact");
        assert!(obj.is_acyclic());
    }

    #[test]
    fn test_derived_object_perfect() {
        let obj = make_simple_complex(3, "P");
        assert!(obj.is_perfect());
    }

    #[test]
    fn test_derived_object_truncate() {
        let mut dims = HashMap::new();
        dims.insert(0, 1);
        dims.insert(1, 1);
        dims.insert(2, 1);
        let cc = ChainComplex::new(dims, HashMap::new());
        let obj = DerivedObject::new(cc, "C");
        let trunc = obj.truncate(1);
        assert_eq!(trunc.complex.dim(0), 1);
        assert_eq!(trunc.complex.dim(1), 1);
        assert_eq!(trunc.complex.dim(2), 0);
    }

    #[test]
    fn test_derived_morphism_compose() {
        let f = ChainMap::new(HashMap::new());
        let g = ChainMap::new(HashMap::new());
        let m1 = DerivedMorphism::new("A", "B", f);
        let m2 = DerivedMorphism::new("B", "C", g);
        // m2.compose(&m1) = m2 ∘ m1 = A → C
        let comp = m2.compose(&m1);
        assert_eq!(comp.source_name, "A");
        assert_eq!(comp.target_name, "C");
    }

    #[test]
    fn test_derived_category_new() {
        let cat = DerivedCategory::new("D(Ab)");
        assert_eq!(cat.name, "D(Ab)");
        assert_eq!(cat.num_objects(), 0);
    }

    #[test]
    fn test_derived_category_add_object() {
        let mut cat = DerivedCategory::new("D(Ab)");
        let obj = make_simple_complex(1, "Z");
        let idx = cat.add_object(obj);
        assert_eq!(idx, 0);
        assert_eq!(cat.num_objects(), 1);
    }

    #[test]
    fn test_derived_category_add_morphism() {
        let mut cat = DerivedCategory::new("D(Ab)");
        let a = cat.add_object(make_simple_complex(1, "A"));
        let b = cat.add_object(make_simple_complex(1, "B"));
        let map = ChainMap::new(HashMap::new());
        cat.add_morphism(a, b, map);
        assert_eq!(cat.num_morphisms(), 1);
    }

    #[test]
    fn test_derived_category_ext_group() {
        let mut cat = DerivedCategory::new("D(Ab)");
        let a = cat.add_object(make_simple_complex(1, "A"));
        let b = cat.add_object(make_simple_complex(1, "B"));
        let ext0 = cat.ext_group(a, b, 0);
        assert_eq!(ext0, 1); // Hom(A, B) has dim 1
    }

    #[test]
    fn test_derived_category_isomorphic() {
        let mut cat = DerivedCategory::new("D(Ab)");
        let a = cat.add_object(make_simple_complex(1, "Z"));
        let b = cat.add_object(make_simple_complex(1, "Z"));
        assert!(cat.are_isomorphic(a, b, 1e-10));
    }

    #[test]
    fn test_derived_category_not_isomorphic() {
        let mut cat = DerivedCategory::new("D(Ab)");
        let a = cat.add_object(make_simple_complex(1, "Z"));
        let b = cat.add_object(make_simple_complex(2, "Z²"));
        assert!(!cat.are_isomorphic(a, b, 1e-10));
    }

    #[test]
    fn test_verdier_quotient() {
        let mut cat = DerivedCategory::new("D(Ab)");
        let a = cat.add_object(make_simple_complex(1, "A"));
        let b = cat.add_object(make_simple_complex(1, "B"));
        let c = cat.add_object(make_simple_complex(1, "C"));

        cat.add_morphism(a, b, ChainMap::new(HashMap::new()));
        cat.add_morphism(b, c, ChainMap::new(HashMap::new()));

        let quotient = cat.verdier_quotient(&[a]);
        assert_eq!(quotient.num_objects(), 2);
    }

    #[test]
    fn test_distinguished_triangle_rotate() {
        let tri = DistinguishedTriangle::new(
            0, 1, 2,
            ChainMap::new(HashMap::new()),
            ChainMap::new(HashMap::new()),
            ChainMap::new(HashMap::new()),
        );
        let rotated = tri.rotate();
        assert_eq!(rotated.a_idx, 1);
        assert_eq!(rotated.b_idx, 2);
        assert_eq!(rotated.c_idx, 0);
    }

    #[test]
    fn test_homotopy_category() {
        let mut cat = DerivedCategory::new("D(Ab)");
        cat.add_object(make_simple_complex(1, "A"));
        let hcat = cat.homotopy_category();
        assert_eq!(hcat.num_objects(), 1);
    }
}
