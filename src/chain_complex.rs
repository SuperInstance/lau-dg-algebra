//! Chain complexes: ... → Vⁿ → Vⁿ⁺¹ → ... with d² = 0.

use serde::{Serialize, Deserialize};
use std::collections::HashMap;

/// A linear map from Vᵏ (dim `domain_dim`) to Vᵏ⁺¹ (dim `codomain_dim`),
/// stored as a row-major matrix.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LinearMap {
    pub domain_dim: usize,
    pub codomain_dim: usize,
    pub entries: Vec<Vec<f64>>,
}

impl LinearMap {
    pub fn zero(domain_dim: usize, codomain_dim: usize) -> Self {
        Self {
            domain_dim,
            codomain_dim,
            entries: vec![vec![0.0; domain_dim]; codomain_dim],
        }
    }

    pub fn identity(dim: usize) -> Self {
        let mut m = Self::zero(dim, dim);
        for i in 0..dim {
            m.entries[i][i] = 1.0;
        }
        m
    }

    pub fn from_rows(rows: Vec<Vec<f64>>) -> Self {
        let codomain_dim = rows.len();
        let domain_dim = if codomain_dim > 0 { rows[0].len() } else { 0 };
        Self { domain_dim, codomain_dim, entries: rows }
    }

    /// Apply to a vector.
    pub fn apply(&self, v: &[f64]) -> Vec<f64> {
        let mut result = vec![0.0; self.codomain_dim];
        for i in 0..self.codomain_dim {
            for j in 0..self.domain_dim {
                result[i] += self.entries[i][j] * v[j];
            }
        }
        result
    }

    /// Compose two maps: self ∘ other.
    pub fn compose(&self, other: &Self) -> Self {
        let mut result = Self::zero(other.domain_dim, self.codomain_dim);
        for i in 0..self.codomain_dim {
            for j in 0..other.domain_dim {
                for k in 0..self.domain_dim {
                    result.entries[i][j] += self.entries[i][k] * other.entries[k][j];
                }
            }
        }
        result
    }

    /// Check if this map is (approximately) zero.
    pub fn is_zero(&self, tol: f64) -> bool {
        self.entries.iter().all(|row| row.iter().all(|x| x.abs() < tol))
    }

    /// Scale the map.
    pub fn scale(&self, s: f64) -> Self {
        Self {
            domain_dim: self.domain_dim,
            codomain_dim: self.codomain_dim,
            entries: self.entries.iter()
                .map(|row| row.iter().map(|x| x * s).collect())
                .collect(),
        }
    }

    /// Add two maps.
    pub fn add(&self, other: &Self) -> Option<Self> {
        if self.domain_dim != other.domain_dim || self.codomain_dim != other.codomain_dim {
            return None;
        }
        Some(Self {
            domain_dim: self.domain_dim,
            codomain_dim: self.codomain_dim,
            entries: self.entries.iter().zip(&other.entries)
                .map(|(r1, r2)| r1.iter().zip(r2).map(|(a, b)| a + b).collect())
                .collect(),
        })
    }

    /// Matrix rank (via row reduction, approximate).
    pub fn rank(&self) -> usize {
        let mut mat = self.entries.clone();
        let rows = mat.len();
        if rows == 0 { return 0; }
        let cols = mat[0].len();
        let mut pivot_row = 0;
        for col in 0..cols {
            // Find pivot
            let mut found = None;
            for row in pivot_row..rows {
                if mat[row][col].abs() > 1e-10 {
                    found = Some(row);
                    break;
                }
            }
            if let Some(pr) = found {
                mat.swap(pivot_row, pr);
                let scale = mat[pivot_row][col];
                for j in 0..cols {
                    mat[pivot_row][j] /= scale;
                }
                for row in 0..rows {
                    if row != pivot_row {
                        let factor = mat[row][col];
                        for j in 0..cols {
                            mat[row][j] -= factor * mat[pivot_row][j];
                        }
                    }
                }
                pivot_row += 1;
            }
        }
        pivot_row
    }

    /// Transpose.
    pub fn transpose(&self) -> Self {
        let mut entries = vec![vec![0.0; self.codomain_dim]; self.domain_dim];
        for i in 0..self.codomain_dim {
            for j in 0..self.domain_dim {
                entries[j][i] = self.entries[i][j];
            }
        }
        Self { domain_dim: self.codomain_dim, codomain_dim: self.domain_dim, entries }
    }

    /// Kernel dimension = domain_dim - rank.
    pub fn kernel_dim(&self) -> usize {
        self.domain_dim.saturating_sub(self.rank())
    }

    /// Image dimension = rank.
    pub fn image_dim(&self) -> usize {
        self.rank()
    }
}

/// A chain complex: ... → Vᵏ → Vᵏ⁺¹ → ...
/// Stored as a map from degree k to (dim Vᵏ, differential d: Vᵏ → Vᵏ⁺¹).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChainComplex {
    /// Dimensions: degree → dimension of Vᵏ.
    pub dimensions: HashMap<i32, usize>,
    /// Differentials: degree k → matrix for d: Vᵏ → Vᵏ⁺¹.
    pub differentials: HashMap<i32, LinearMap>,
}

impl ChainComplex {
    pub fn new(
        dimensions: HashMap<i32, usize>,
        differentials: HashMap<i32, LinearMap>,
    ) -> Self {
        Self { dimensions, differentials }
    }

    /// Check that d² = 0 for all degrees.
    pub fn check_d_squared_zero(&self, tol: f64) -> bool {
        for (&k, d_k) in &self.differentials {
            if let Some(d_k1) = self.differentials.get(&(k + 1)) {
                let d2 = d_k1.compose(d_k);
                if !d2.is_zero(tol) {
                    return false;
                }
            }
        }
        true
    }

    /// Dimension of Vᵏ.
    pub fn dim(&self, k: i32) -> usize {
        *self.dimensions.get(&k).unwrap_or(&0)
    }

    /// Compute Betti numbers: βₖ = dim ker(dₖ) - dim im(dₖ₋₁).
    pub fn betti_numbers(&self) -> HashMap<i32, usize> {
        let mut betti = HashMap::new();
        let degrees = self.degrees();
        for &k in &degrees {
            let dk_dim = self.dim(k);
            let ker_d = if let Some(ref d_k) = self.differentials.get(&k) {
                d_k.kernel_dim()
            } else {
                dk_dim
            };
            let im_d_prev = if let Some(ref d_prev) = self.differentials.get(&(k - 1)) {
                d_prev.image_dim()
            } else {
                0
            };
            if ker_d >= im_d_prev {
                betti.insert(k, ker_d - im_d_prev);
            } else {
                betti.insert(k, 0);
            }
        }
        betti
    }

    /// All degrees present.
    pub fn degrees(&self) -> Vec<i32> {
        let mut ds: Vec<i32> = self.dimensions.keys().copied().collect();
        ds.sort();
        ds
    }

    /// Direct sum of two chain complexes.
    pub fn direct_sum(&self, other: &Self) -> Self {
        let mut dims = self.dimensions.clone();
        for (&k, &d) in &other.dimensions {
            *dims.entry(k).or_insert(0) += d;
        }
        let mut diffs = HashMap::new();
        let all_degrees: std::collections::HashSet<i32> = 
            self.differentials.keys().chain(other.differentials.keys()).copied().collect();
        for &k in &all_degrees {
            let d1 = self.differentials.get(&k);
            let d2 = other.differentials.get(&k);
            match (d1, d2) {
                (Some(a), Some(b)) => {
                    // Block diagonal
                    let n = a.domain_dim + b.domain_dim;
                    let m = a.codomain_dim + b.codomain_dim;
                    let mut entries = vec![vec![0.0; n]; m];
                    for i in 0..a.codomain_dim {
                        for j in 0..a.domain_dim {
                            entries[i][j] = a.entries[i][j];
                        }
                    }
                    for i in 0..b.codomain_dim {
                        for j in 0..b.domain_dim {
                            entries[a.codomain_dim + i][a.domain_dim + j] = b.entries[i][j];
                        }
                    }
                    diffs.insert(k, LinearMap { domain_dim: n, codomain_dim: m, entries });
                }
                (Some(a), None) => { diffs.insert(k, a.clone()); }
                (None, Some(b)) => { diffs.insert(k, b.clone()); }
                _ => {}
            }
        }
        Self::new(dims, diffs)
    }

    /// Shift the chain complex by n: (C[n])ᵏ = Cᵏ⁺ⁿ.
    pub fn shift(&self, n: i32) -> Self {
        let dims = self.dimensions.iter().map(|(&k, &d)| (k + n, d)).collect();
        let diffs = self.differentials.iter().map(|(&k, d)| (k + n, d.clone())).collect();
        Self::new(dims, diffs)
    }

    /// The Euler characteristic: Σ (-1)ᵏ dim(Cᵏ).
    pub fn euler_characteristic(&self) -> i64 {
        self.dimensions.iter()
            .map(|(&k, &d)| ((-1i64).pow(k as u32)) * d as i64)
            .sum()
    }

    /// Mapping cone of a chain map f: self → other.
    pub fn mapping_cone(&self, other: &Self, f: &ChainMap) -> Self {
        let mut dims = HashMap::new();
        for (&k, &d) in &self.dimensions {
            *dims.entry(k + 1).or_insert(0) += d;
        }
        for (&k, &d) in &other.dimensions {
            *dims.entry(k).or_insert(0) += d;
        }
        // Simplified: we store direct sum dimensions, differentials would be block matrices
        Self::new(dims, HashMap::new())
    }
}

/// A chain map between two chain complexes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChainMap {
    /// The maps fᵏ: Vᵏ → Wᵏ for each degree.
    pub maps: HashMap<i32, LinearMap>,
}

impl ChainMap {
    pub fn new(maps: HashMap<i32, LinearMap>) -> Self {
        Self { maps }
    }

    /// Check if this is a chain map: f ∘ d_C = d_D ∘ f.
    pub fn check_chain_map(&self, source: &ChainComplex, target: &ChainComplex, tol: f64) -> bool {
        for (&k, f_k) in &self.maps {
            let d_source = source.differentials.get(&k);
            let d_target = target.differentials.get(&k);
            // f_{k+1} ∘ d_k^source vs d_k^target ∘ f_k
            if let (Some(d_s), Some(f_k1)) = (d_source, self.maps.get(&(k + 1))) {
                let lhs = f_k1.compose(d_s);
                if let Some(d_t) = d_target {
                    let rhs = d_t.compose(f_k);
                    let rhs_neg = rhs.scale(-1.0);
                    if let Some(diff) = lhs.add(&rhs_neg) {
                        if !diff.is_zero(tol) {
                            return false;
                        }
                    }
                }
            }
        }
        true
    }

    /// Zero chain map.
    pub fn zero() -> Self {
        Self { maps: HashMap::new() }
    }

    /// Identity chain map.
    pub fn identity(complex: &ChainComplex) -> Self {
        let maps = complex.dimensions.iter()
            .map(|(&k, &d)| (k, LinearMap::identity(d)))
            .collect();
        Self { maps }
    }

    /// Compose two chain maps.
    pub fn compose(&self, other: &Self) -> Self {
        let mut maps = HashMap::new();
        for (&k, f) in &self.maps {
            if let Some(g) = other.maps.get(&k) {
                maps.insert(k, f.compose(g));
            }
        }
        Self { maps }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linear_map_zero() {
        let m = LinearMap::zero(3, 2);
        assert!(m.is_zero(1e-10));
        assert_eq!(m.domain_dim, 3);
        assert_eq!(m.codomain_dim, 2);
    }

    #[test]
    fn test_linear_map_identity() {
        let m = LinearMap::identity(3);
        let v = vec![1.0, 2.0, 3.0];
        assert_eq!(m.apply(&v), v);
    }

    #[test]
    fn test_linear_map_apply() {
        let m = LinearMap::from_rows(vec![
            vec![1.0, 2.0],
            vec![3.0, 4.0],
        ]);
        let result = m.apply(&[1.0, 0.0]);
        assert!((result[0] - 1.0).abs() < 1e-10);
        assert!((result[1] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_linear_map_compose() {
        let a = LinearMap::from_rows(vec![vec![1.0, 0.0], vec![0.0, 1.0]]);
        let b = LinearMap::from_rows(vec![vec![2.0, 0.0], vec![0.0, 3.0]]);
        let c = a.compose(&b);
        assert!((c.entries[0][0] - 2.0).abs() < 1e-10);
        assert!((c.entries[1][1] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_linear_map_scale() {
        let m = LinearMap::from_rows(vec![vec![1.0, 2.0]]);
        let s = m.scale(3.0);
        assert!((s.entries[0][0] - 3.0).abs() < 1e-10);
        assert!((s.entries[0][1] - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_linear_map_rank() {
        let m = LinearMap::from_rows(vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
        ]);
        assert_eq!(m.rank(), 2);

        let m2 = LinearMap::from_rows(vec![
            vec![1.0, 2.0],
            vec![2.0, 4.0],
        ]);
        assert_eq!(m2.rank(), 1);
    }

    #[test]
    fn test_linear_map_transpose() {
        let m = LinearMap::from_rows(vec![vec![1.0, 2.0, 3.0]]);
        let t = m.transpose();
        assert_eq!(t.codomain_dim, 3);
        assert_eq!(t.domain_dim, 1);
        assert!((t.entries[0][0] - 1.0).abs() < 1e-10);
        assert!((t.entries[2][0] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_chain_complex_d_squared() {
        // Simple: V⁰ --d₀--> V¹ --d₁--> V²
        // d₀ = [[1], [-1]], d₁ = [[1, 1]]
        // d₁ ∘ d₀ = [[1*1 + 1*(-1)]] = [[0]] ✓
        let mut dims = HashMap::new();
        dims.insert(0, 1);
        dims.insert(1, 2);
        dims.insert(2, 1);

        let mut diffs = HashMap::new();
        diffs.insert(0, LinearMap::from_rows(vec![vec![1.0], vec![-1.0]]));
        diffs.insert(1, LinearMap::from_rows(vec![vec![1.0, 1.0]]));

        let cc = ChainComplex::new(dims, diffs);
        assert!(cc.check_d_squared_zero(1e-10));
    }

    #[test]
    fn test_chain_complex_betti() {
        // 0 → Z → 0 (single Z in degree 0, no differentials)
        let mut dims = HashMap::new();
        dims.insert(0, 1);
        let cc = ChainComplex::new(dims, HashMap::new());
        let betti = cc.betti_numbers();
        assert_eq!(betti.get(&0), Some(&1));
    }

    #[test]
    fn test_chain_complex_euler() {
        let mut dims = HashMap::new();
        dims.insert(0, 3);
        dims.insert(1, 2);
        dims.insert(2, 1);
        let cc = ChainComplex::new(dims, HashMap::new());
        assert_eq!(cc.euler_characteristic(), 3 - 2 + 1);
    }

    #[test]
    fn test_chain_complex_shift() {
        let mut dims = HashMap::new();
        dims.insert(0, 1);
        let cc = ChainComplex::new(dims, HashMap::new());
        let shifted = cc.shift(2);
        assert_eq!(shifted.dim(2), 1);
    }

    #[test]
    fn test_chain_map_identity() {
        let mut dims = HashMap::new();
        dims.insert(0, 2);
        let cc = ChainComplex::new(dims, HashMap::new());
        let id = ChainMap::identity(&cc);
        let v = vec![1.0, 2.0];
        assert_eq!(id.maps[&0].apply(&v), v);
    }

    #[test]
    fn test_chain_complex_direct_sum() {
        let mut dims1 = HashMap::new();
        dims1.insert(0, 1);
        let cc1 = ChainComplex::new(dims1, HashMap::new());

        let mut dims2 = HashMap::new();
        dims2.insert(0, 2);
        let cc2 = ChainComplex::new(dims2, HashMap::new());

        let sum = cc1.direct_sum(&cc2);
        assert_eq!(sum.dim(0), 3);
    }
}
