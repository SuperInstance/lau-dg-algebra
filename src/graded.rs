//! Graded vector spaces: V = ⊕Vᵏ with degree-k components.

use serde::{Serialize, Deserialize};
use std::collections::HashMap;

/// A homogeneous element of a graded vector space.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GradedElement {
    /// The coefficient (field element), stored as f64 for generality.
    pub coeff: f64,
    /// The degree of this element.
    pub degree: i32,
    /// An index identifying the basis element.
    pub index: usize,
}

impl GradedElement {
    pub fn new(coeff: f64, degree: i32, index: usize) -> Self {
        Self { coeff, degree, index }
    }

    /// The degree |a| of this element.
    pub fn degree(&self) -> i32 {
        self.degree
    }

    /// Scale by a scalar.
    pub fn scale(&self, s: f64) -> Self {
        Self { coeff: self.coeff * s, degree: self.degree, index: self.index }
    }
}

impl std::fmt::Display for GradedElement {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}·e_{}[{}]", self.coeff, self.index, self.degree)
    }
}

/// A graded vector space V = ⊕ Vᵏ, stored as a map from degree to vectors.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GradedVectorSpace {
    /// Map from degree k to the basis dimension of Vᵏ.
    pub dimensions: HashMap<i32, usize>,
}

impl GradedVectorSpace {
    /// Create a new graded vector space with given dimensions per degree.
    pub fn new(dimensions: HashMap<i32, usize>) -> Self {
        Self { dimensions }
    }

    /// The zero graded vector space.
    pub fn zero() -> Self {
        Self { dimensions: HashMap::new() }
    }

    /// Dimension of the degree-k component.
    pub fn dim(&self, k: i32) -> usize {
        *self.dimensions.get(&k).unwrap_or(&0)
    }

    /// Total dimension (sum over all degrees).
    pub fn total_dim(&self) -> usize {
        self.dimensions.values().sum()
    }

    /// Direct sum of two graded vector spaces.
    pub fn direct_sum(&self, other: &Self) -> Self {
        let mut dims = self.dimensions.clone();
        for (&k, &d) in &other.dimensions {
            *dims.entry(k).or_insert(0) += d;
        }
        Self::new(dims)
    }

    /// Tensor product of two graded vector spaces.
    /// (V ⊗ W)ⁿ = ⊕_{p+q=n} Vᵖ ⊗ Wᵠ
    pub fn tensor_product(&self, other: &Self) -> Self {
        let mut dims: HashMap<i32, usize> = HashMap::new();
        for (&p, &dp) in &self.dimensions {
            for (&q, &dq) in &other.dimensions {
                *dims.entry(p + q).or_insert(0) += dp * dq;
            }
        }
        Self::new(dims)
    }

    /// Shift: V[n]ᵏ = Vᵏ⁻ⁿ (suspension/desuspension).
    pub fn shift(&self, n: i32) -> Self {
        let dims = self.dimensions.iter()
            .map(|(&k, &d)| (k + n, d))
            .collect();
        Self::new(dims)
    }

    /// Dual graded vector space.
    pub fn dual(&self) -> Self {
        Self::new(self.dimensions.clone())
    }

    /// Check if the graded vector space is concentrated in a single degree.
    pub fn is_concentrated(&self) -> bool {
        self.dimensions.len() <= 1
    }

    /// All degrees present.
    pub fn degrees(&self) -> Vec<i32> {
        let mut ds: Vec<i32> = self.dimensions.keys().copied().collect();
        ds.sort();
        ds
    }

    /// The Euler characteristic: Σ (-1)ᵏ dim(Vᵏ).
    pub fn euler_characteristic(&self) -> i64 {
        self.dimensions.iter()
            .map(|(&k, &d)| ((-1i64).pow(k as u32)) * d as i64)
            .sum()
    }

    /// Truncate to degrees in [lo, hi].
    pub fn truncate(&self, lo: i32, hi: i32) -> Self {
        let dims = self.dimensions.iter()
            .filter(|(&k, _)| k >= lo && k <= hi)
            .map(|(&k, &d)| (k, d))
            .collect();
        Self::new(dims)
    }
}

/// A homogeneous element of a graded vector space, with value stored
/// as a nalgebra vector for the component in its degree.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GradedVec {
    pub degree: i32,
    pub data: Vec<f64>,
}

impl GradedVec {
    pub fn new(degree: i32, data: Vec<f64>) -> Self {
        Self { degree, data }
    }

    pub fn zero(degree: i32, dim: usize) -> Self {
        Self { degree, data: vec![0.0; dim] }
    }

    pub fn dim(&self) -> usize {
        self.data.len()
    }

    pub fn scale(&self, s: f64) -> Self {
        Self { degree: self.degree, data: self.data.iter().map(|x| x * s).collect() }
    }

    pub fn add(&self, other: &Self) -> Option<Self> {
        if self.degree != other.degree || self.data.len() != other.data.len() {
            return None;
        }
        Some(Self {
            degree: self.degree,
            data: self.data.iter().zip(&other.data).map(|(a, b)| a + b).collect(),
        })
    }

    pub fn is_zero(&self) -> bool {
        self.data.iter().all(|x| x.abs() < 1e-12)
    }

    pub fn norm(&self) -> f64 {
        self.data.iter().map(|x| x * x).sum::<f64>().sqrt()
    }
}

/// A general element of a graded vector space (sum of homogeneous components).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneralGradedVec {
    pub components: HashMap<i32, Vec<f64>>,
}

impl GeneralGradedVec {
    pub fn new(components: HashMap<i32, Vec<f64>>) -> Self {
        Self { components }
    }

    pub fn zero() -> Self {
        Self { components: HashMap::new() }
    }

    pub fn homogeneous(degree: i32, data: Vec<f64>) -> Self {
        let mut c = HashMap::new();
        c.insert(degree, data);
        Self { components: c }
    }

    pub fn add(&self, other: &Self) -> Self {
        let mut result = self.components.clone();
        for (&k, v) in &other.components {
            let entry = result.entry(k).or_insert_with(|| vec![0.0; v.len()]);
            for (i, val) in v.iter().enumerate() {
                if i < entry.len() {
                    entry[i] += val;
                }
            }
        }
        Self { components: result }
    }

    pub fn scale(&self, s: f64) -> Self {
        Self {
            components: self.components.iter()
                .map(|(&k, v)| (k, v.iter().map(|x| x * s).collect()))
                .collect(),
        }
    }

    pub fn degree_part(&self, k: i32) -> Option<&Vec<f64>> {
        self.components.get(&k)
    }

    pub fn is_zero(&self) -> bool {
        self.components.values().all(|v| v.iter().all(|x| x.abs() < 1e-12))
    }

    pub fn degrees(&self) -> Vec<i32> {
        let mut ds: Vec<i32> = self.components.keys().copied().collect();
        ds.sort();
        ds
    }

    /// The top degree with nonzero component.
    pub fn top_degree(&self) -> Option<i32> {
        self.degrees().into_iter()
            .filter(|&k| !self.components[&k].iter().all(|x| x.abs() < 1e-12))
            .last()
    }

    /// The bottom degree with nonzero component.
    pub fn bottom_degree(&self) -> Option<i32> {
        self.degrees().into_iter()
            .filter(|&k| !self.components[&k].iter().all(|x| x.abs() < 1e-12))
            .next()
    }
}

/// Koszul sign: (-1)^{|a||b|}.
pub fn koszul_sign(deg_a: i32, deg_b: i32) -> f64 {
    if (deg_a * deg_b) % 2 == 0 { 1.0 } else { -1.0 }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_graded_element_creation() {
        let e = GradedElement::new(2.0, 3, 0);
        assert_eq!(e.coeff, 2.0);
        assert_eq!(e.degree(), 3);
        assert_eq!(e.index, 0);
    }

    #[test]
    fn test_graded_element_scale() {
        let e = GradedElement::new(3.0, 2, 1);
        let s = e.scale(2.0);
        assert!((s.coeff - 6.0).abs() < 1e-10);
        assert_eq!(s.degree, 2);
    }

    #[test]
    fn test_graded_element_display() {
        let e = GradedElement::new(1.5, 3, 2);
        let s = format!("{}", e);
        assert!(s.contains("1.5"));
        assert!(s.contains("e_2"));
        assert!(s.contains("[3]"));
    }

    #[test]
    fn test_graded_vector_space_new() {
        let mut dims = HashMap::new();
        dims.insert(0, 3);
        dims.insert(1, 2);
        let v = GradedVectorSpace::new(dims);
        assert_eq!(v.dim(0), 3);
        assert_eq!(v.dim(1), 2);
        assert_eq!(v.dim(2), 0);
    }

    #[test]
    fn test_graded_vector_space_zero() {
        let v = GradedVectorSpace::zero();
        assert_eq!(v.total_dim(), 0);
    }

    #[test]
    fn test_total_dim() {
        let mut dims = HashMap::new();
        dims.insert(0, 3);
        dims.insert(1, 2);
        dims.insert(2, 1);
        let v = GradedVectorSpace::new(dims);
        assert_eq!(v.total_dim(), 6);
    }

    #[test]
    fn test_direct_sum() {
        let mut d1 = HashMap::new();
        d1.insert(0, 2);
        d1.insert(1, 1);
        let v1 = GradedVectorSpace::new(d1);

        let mut d2 = HashMap::new();
        d2.insert(1, 3);
        d2.insert(2, 2);
        let v2 = GradedVectorSpace::new(d2);

        let sum = v1.direct_sum(&v2);
        assert_eq!(sum.dim(0), 2);
        assert_eq!(sum.dim(1), 4);
        assert_eq!(sum.dim(2), 2);
    }

    #[test]
    fn test_tensor_product() {
        let mut d1 = HashMap::new();
        d1.insert(0, 2);
        d1.insert(1, 1);
        let v1 = GradedVectorSpace::new(d1);

        let mut d2 = HashMap::new();
        d2.insert(0, 3);
        d2.insert(1, 2);
        let v2 = GradedVectorSpace::new(d2);

        let tensor = v1.tensor_product(&v2);
        assert_eq!(tensor.dim(0), 6);  // 2*3
        assert_eq!(tensor.dim(1), 7);  // 2*2 + 1*3
        assert_eq!(tensor.dim(2), 2);  // 1*2
    }

    #[test]
    fn test_shift() {
        let mut dims = HashMap::new();
        dims.insert(0, 3);
        dims.insert(1, 2);
        let v = GradedVectorSpace::new(dims);

        let shifted = v.shift(2);
        assert_eq!(shifted.dim(2), 3);
        assert_eq!(shifted.dim(3), 2);
        assert_eq!(shifted.dim(0), 0);
    }

    #[test]
    fn test_dual() {
        let mut dims = HashMap::new();
        dims.insert(0, 3);
        let v = GradedVectorSpace::new(dims);
        let dual = v.dual();
        assert_eq!(dual.dim(0), 3);
    }

    #[test]
    fn test_is_concentrated() {
        let mut dims1 = HashMap::new();
        dims1.insert(0, 3);
        let v1 = GradedVectorSpace::new(dims1);
        assert!(v1.is_concentrated());

        let mut dims2 = HashMap::new();
        dims2.insert(0, 3);
        dims2.insert(1, 2);
        let v2 = GradedVectorSpace::new(dims2);
        assert!(!v2.is_concentrated());
    }

    #[test]
    fn test_degrees() {
        let mut dims = HashMap::new();
        dims.insert(2, 1);
        dims.insert(0, 3);
        dims.insert(1, 2);
        let v = GradedVectorSpace::new(dims);
        assert_eq!(v.degrees(), vec![0, 1, 2]);
    }

    #[test]
    fn test_euler_characteristic() {
        let mut dims = HashMap::new();
        dims.insert(0, 3);
        dims.insert(1, 2);
        dims.insert(2, 1);
        let v = GradedVectorSpace::new(dims);
        assert_eq!(v.euler_characteristic(), 3 - 2 + 1); // 2
    }

    #[test]
    fn test_truncate() {
        let mut dims = HashMap::new();
        dims.insert(0, 3);
        dims.insert(1, 2);
        dims.insert(2, 1);
        dims.insert(3, 4);
        let v = GradedVectorSpace::new(dims);
        let trunc = v.truncate(1, 2);
        assert_eq!(trunc.dim(0), 0);
        assert_eq!(trunc.dim(1), 2);
        assert_eq!(trunc.dim(2), 1);
        assert_eq!(trunc.dim(3), 0);
    }

    #[test]
    fn test_graded_vec_add() {
        let v1 = GradedVec::new(0, vec![1.0, 2.0]);
        let v2 = GradedVec::new(0, vec![3.0, 4.0]);
        let sum = v1.add(&v2).unwrap();
        assert!((sum.data[0] - 4.0).abs() < 1e-10);
        assert!((sum.data[1] - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_graded_vec_add_different_degree_fails() {
        let v1 = GradedVec::new(0, vec![1.0]);
        let v2 = GradedVec::new(1, vec![1.0]);
        assert!(v1.add(&v2).is_none());
    }

    #[test]
    fn test_graded_vec_is_zero() {
        let v = GradedVec::new(0, vec![0.0, 0.0]);
        assert!(v.is_zero());
        let v2 = GradedVec::new(0, vec![0.0, 1.0]);
        assert!(!v2.is_zero());
    }

    #[test]
    fn test_graded_vec_norm() {
        let v = GradedVec::new(0, vec![3.0, 4.0]);
        assert!((v.norm() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_general_graded_vec_add() {
        let mut c1 = HashMap::new();
        c1.insert(0, vec![1.0, 2.0]);
        let g1 = GeneralGradedVec::new(c1);

        let mut c2 = HashMap::new();
        c2.insert(0, vec![3.0, 4.0]);
        c2.insert(1, vec![5.0]);
        let g2 = GeneralGradedVec::new(c2);

        let sum = g1.add(&g2);
        assert_eq!(sum.degree_part(0).unwrap()[0], 4.0);
        assert_eq!(sum.degree_part(1).unwrap()[0], 5.0);
    }

    #[test]
    fn test_general_graded_vec_is_zero() {
        let g = GeneralGradedVec::zero();
        assert!(g.is_zero());
        let g2 = GeneralGradedVec::homogeneous(0, vec![0.0]);
        assert!(g2.is_zero());
        let g3 = GeneralGradedVec::homogeneous(0, vec![1.0]);
        assert!(!g3.is_zero());
    }

    #[test]
    fn test_koszul_sign_even_even() {
        assert_eq!(koszul_sign(0, 0), 1.0);
        assert_eq!(koszul_sign(2, 4), 1.0);
    }

    #[test]
    fn test_koszul_sign_odd_odd() {
        assert_eq!(koszul_sign(1, 1), -1.0);
        assert_eq!(koszul_sign(1, 3), -1.0);
    }

    #[test]
    fn test_koszul_sign_even_odd() {
        assert_eq!(koszul_sign(0, 1), 1.0);
        assert_eq!(koszul_sign(2, 1), 1.0);
    }

    #[test]
    fn test_top_bottom_degree() {
        let mut c = HashMap::new();
        c.insert(1, vec![1.0]);
        c.insert(3, vec![2.0]);
        c.insert(5, vec![3.0]);
        let g = GeneralGradedVec::new(c);
        assert_eq!(g.bottom_degree(), Some(1));
        assert_eq!(g.top_degree(), Some(5));
    }
}
