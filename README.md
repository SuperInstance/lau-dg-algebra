# lau-dg-algebra

> **Differential graded algebras (DGAs) for agents** — the algebraic structure underlying every cohomology theory in the PLATO ecosystem.

## What This Does

This crate implements the full machinery of **differential graded algebras** from scratch: graded vector spaces, chain complexes, differentials satisfying d² = 0, the graded Leibniz rule, cohomology computation, quasi-isomorphisms, minimal models, Massey products, A∞-algebras, and derived categories. It also includes concrete applications — de Rham cohomology as a DGA, the exterior algebra, and Poincaré duality verification.

Every cohomology theory in the PLATO/LAU ecosystem (de Rham, sheaf, Hodge) is a DGA. This crate provides the shared foundation.

Part of the **PLATO/LAU ecosystem** — a mathematically rigorous framework for building educational agents that learn, teach, and evolve.

## Key Idea

A **differential graded algebra** is a graded vector space V = ⊕Vᵏ equipped with:
1. A **differential** d: Vᵏ → Vᵏ⁺¹ with d² = 0
2. A **graded multiplication** Vᵖ ⊗ Vᵠ → Vᵖ⁺ᵠ
3. The **graded Leibniz rule**: d(a·b) = da·b + (−1)^|a| a·db

From d² = 0 alone, you get **cohomology**: H(d) = ker(d)/im(d), a graded algebra. The Betti numbers βₖ = dim Hᵏ are topological invariants. The product structure on cohomology captures information that Betti numbers alone miss — and Massey products capture even more.

This crate doesn't just *use* these structures — it **builds** them from linear algebra primitives, checking every axiom along the way.

## Install

```bash
cargo add lau-dg-algebra
```

Dependencies: `nalgebra` (linear algebra), `serde` (serialization).

## Quick Start

```rust
use lau_dg_algebra::*;
use std::collections::HashMap;

// Build a chain complex: V⁰ → V¹ → V² with d² = 0
let mut dims = HashMap::new();
dims.insert(0, 2);
dims.insert(1, 3);
dims.insert(2, 1);

let mut diffs = HashMap::new();
diffs.insert(0, LinearMap::from_rows(vec![
    vec![1.0, 0.0],   // d₀: V⁰ → V¹
    vec![-1.0, 1.0],
    vec![0.0, -1.0],
]));
diffs.insert(1, LinearMap::from_rows(vec![
    vec![0.0, 1.0, 1.0],  // d₁: V¹ → V²
]));

let cc = ChainComplex::new(dims, diffs);
assert!(cc.check_d_squared_zero(1e-10));

// Compute cohomology and Betti numbers
let cohom = Cohomology::from_chain_complex(&cc);
println!("Betti numbers: {:?}", cohom.betti);
println!("Total Betti: {}", cohom.total_betti());

// Graded vector space operations
let mut d1 = HashMap::new();
d1.insert(0, 2);
d1.insert(1, 3);
let v1 = GradedVectorSpace::new(d1);

let tensor = v1.tensor_product(&v1);
println!("Tensor product dims: {:?}", tensor.dimensions);

// Koszul sign rule
let sign = koszul_sign(1, 1); // (-1)^{1·1} = -1
assert_eq!(sign, -1.0);

// De Rham DGA for a 3-manifold
let de_rham = DeRhamDGA::new(3);
println!("Manifold dim: {}", de_rham.manifold_dim);
// Ω⁰ has dim 1, Ω¹ has dim 3, Ω² has dim 3, Ω³ has dim 1
```

## API Reference

### Graded Vector Spaces (`graded`)

```rust
// A homogeneous element
pub struct GradedElement {
    pub coeff: f64,
    pub degree: i32,
    pub index: usize,
}

// A graded vector space V = ⊕Vᵏ
pub struct GradedVectorSpace {
    pub dimensions: HashMap<i32, usize>,  // degree → dimension
}

impl GradedVectorSpace {
    pub fn new(dimensions: HashMap<i32, usize>) -> Self;
    pub fn zero() -> Self;
    pub fn dim(&self, k: i32) -> usize;
    pub fn total_dim(&self) -> usize;
    pub fn direct_sum(&self, other: &Self) -> Self;
    pub fn tensor_product(&self, other: &Self) -> Self;  // (V⊗W)ⁿ = ⊕_{p+q=n} Vᵖ⊗Wᵠ
    pub fn shift(&self, n: i32) -> Self;                  // V[n]ᵏ = Vᵏ⁻ⁿ
    pub fn dual(&self) -> Self;
    pub fn euler_characteristic(&self) -> i64;
    pub fn truncate(&self, lo: i32, hi: i32) -> Self;
    pub fn degrees(&self) -> Vec<i32>;
}

// Concrete vectors
pub struct GradedVec { pub degree: i32, pub data: Vec<f64> }
pub struct GeneralGradedVec { pub components: HashMap<i32, Vec<f64>> }

// Koszul sign: (-1)^{|a|·|b|}
pub fn koszul_sign(deg_a: i32, deg_b: i32) -> f64;
```

### Chain Complexes (`chain_complex`)

```rust
// A linear map stored as a row-major matrix
pub struct LinearMap {
    pub domain_dim: usize,
    pub codomain_dim: usize,
    pub entries: Vec<Vec<f64>>,
}

impl LinearMap {
    pub fn zero(domain_dim: usize, codomain_dim: usize) -> Self;
    pub fn identity(dim: usize) -> Self;
    pub fn from_rows(rows: Vec<Vec<f64>>) -> Self;
    pub fn apply(&self, v: &[f64]) -> Vec<f64>;
    pub fn compose(&self, other: &Self) -> Self;
    pub fn rank(&self) -> usize;              // via row reduction
    pub fn kernel_dim(&self) -> usize;        // domain_dim - rank
    pub fn image_dim(&self) -> usize;         // rank
    pub fn transpose(&self) -> Self;
    pub fn is_zero(&self, tol: f64) -> bool;
}

// A chain complex: ... → Vᵏ → Vᵏ⁺¹ → ...
pub struct ChainComplex {
    pub dimensions: HashMap<i32, usize>,
    pub differentials: HashMap<i32, LinearMap>,
}

impl ChainComplex {
    pub fn new(dimensions, differentials) -> Self;
    pub fn check_d_squared_zero(&self, tol: f64) -> bool;
    pub fn betti_numbers(&self) -> HashMap<i32, usize>;  // βₖ = dim ker(dₖ) - dim im(dₖ₋₁)
    pub fn direct_sum(&self, other: &Self) -> Self;
    pub fn shift(&self, n: i32) -> Self;
    pub fn euler_characteristic(&self) -> i64;
    pub fn mapping_cone(&self, other: &Self, f: &ChainMap) -> Self;
}

// A chain map f: C → D (commutes with differentials)
pub struct ChainMap {
    pub maps: HashMap<i32, LinearMap>,
}
```

### Differential Graded Algebras (`dga`)

```rust
pub struct MultiplicationTable {
    pub tables: HashMap<(i32, i32), LinearMap>,  // (deg_p, deg_q) → multiplication matrix
}

impl MultiplicationTable {
    pub fn multiply(&self, a_deg: i32, a: &[f64], b_deg: i32, b: &[f64]) -> Option<Vec<f64>>;
}

pub struct DGA {
    pub graded_vs: GradedVectorSpace,
    pub differential: HashMap<i32, LinearMap>,
    pub multiplication: MultiplicationTable,
}

impl DGA {
    pub fn new(graded_vs, differential, multiplication) -> Self;
    pub fn check_d_squared_zero(&self, tol: f64) -> bool;
    pub fn check_leibniz(&self, deg_a: i32, deg_b: i32, tol: f64) -> bool;
    pub fn check_associativity(&self, deg_a: i32, deg_b: i32, deg_c: i32, tol: f64) -> bool;
    pub fn d(&self, k: i32, v: &[f64]) -> Option<Vec<f64>>;
    pub fn as_chain_complex(&self) -> ChainComplex;
}

// A DGA morphism (preserves degree, multiplication, and commutes with d)
pub struct DGAMorphism {
    pub maps: HashMap<i32, LinearMap>,
}
```

### Cohomology (`cohomology`)

```rust
pub struct Cohomology {
    pub betti: HashMap<i32, usize>,
    pub representatives: HashMap<i32, Vec<Vec<f64>>>,
    pub graded_dims: HashMap<i32, usize>,
}

impl Cohomology {
    pub fn from_chain_complex(cc: &ChainComplex) -> Self;
    pub fn from_dga(dga: &DGA) -> Self;
    pub fn betti_number(&self, k: i32) -> usize;
    pub fn total_betti(&self) -> usize;
    pub fn euler_characteristic(&self) -> i64;
    pub fn is_trivial(&self) -> bool;  // all Betti numbers zero
    pub fn poincare_polynomial(&self) -> String;  // Σ βₖ · tᵏ
}
```

### Quasi-Isomorphisms (`quasi_iso`)

```rust
// Check if a morphism induces isomorphism on cohomology
pub fn is_quasi_isomorphism(source: &DGA, target: &DGA, morphism: &DGAMorphism, tol: f64) -> bool;

// Check if two DGAs have isomorphic cohomology (necessary condition)
pub fn are_quasi_isomorphic(dga_a: &DGA, dga_b: &DGA) -> bool;

pub struct QuasiIsomorphism {
    pub morphism: DGAMorphism,
    pub induced_cohomology_map: HashMap<i32, LinearMap>,
    pub source_cohomology: Cohomology,
    pub target_cohomology: Cohomology,
}
```

### Minimal Models (`minimal_model`)

```rust
// Free graded commutative algebra on generators
pub struct FreeGradedAlgebra {
    pub generator_degrees: Vec<i32>,
    pub generator_parity: Vec<bool>,  // true = odd (exterior), false = even (polynomial)
}

impl FreeGradedAlgebra {
    pub fn exterior(degrees: Vec<i32>) -> Self;
    pub fn polynomial(degrees: Vec<i32>) -> Self;
    pub fn dimensions(&self, max_degree: i32) -> HashMap<i32, usize>;
}

// A minimal model: free DGA quasi-isomorphic to the original
pub struct MinimalModel {
    pub free_algebra: FreeGradedAlgebra,
    pub differential: HashMap<i32, LinearMap>,
    pub quasi_iso: DGAMorphism,
}

pub fn compute_minimal_model(dga: &DGA, max_degree: i32) -> Option<MinimalModel>;
```

### Massey Products (`massey`)

```rust
// Triple Massey product ⟨α, β, γ⟩
pub struct TripleMasseyProduct {
    pub deg_a: i32, pub deg_b: i32, pub deg_c: i32,
    pub alpha: Vec<f64>, pub beta: Vec<f64>, pub gamma: Vec<f64>,
    pub result: Vec<f64>,
    pub result_degree: i32,
    pub indeterminacy: Vec<Vec<f64>>,
}

impl TripleMasseyProduct {
    pub fn compute(dga: &DGA, deg_a: i32, alpha: &[f64], deg_b: i32, beta: &[f64],
                   deg_c: i32, gamma: &[f64], tol: f64) -> Option<Self>;
}

// Higher-order Massey products ⟨α₁, ..., αₙ⟩
pub struct HigherMasseyProduct { /* ... */ }
```

### A∞-Algebras (`a_infinity`)

```rust
pub struct AInfinityAlgebra {
    pub dimensions: HashMap<i32, usize>,
    pub operations: Vec<AInfinityOperation>,  // m₁, m₂, m₃, ...
}

impl AInfinityAlgebra {
    pub fn m(&self, n: usize) -> Option<&AInfinityOperation>;
    pub fn check_relation(&self, n: usize, tol: f64) -> bool;  // Stasheff identities
    pub fn from_dga(dga: &DGA) -> Self;  // any DGA is an A∞-algebra with mₙ=0 for n≥3
}
```

### Derived Category (`derived`)

```rust
pub struct DerivedObject {
    pub complex: ChainComplex,
    pub name: String,
}

impl DerivedObject {
    pub fn cohomology(&self) -> Cohomology;
    pub fn shift(&self, n: i32) -> Self;       // C[n]
    pub fn direct_sum(&self, other: &Self) -> Self;
    pub fn is_acyclic(&self) -> bool;
    pub fn is_perfect(&self) -> bool;
    pub fn truncate(&self, n: i32) -> Self;
}

pub struct DerivedMorphism { pub chain_map: ChainMap, /* ... */ }
```

### Applications (`applications`)

```rust
// De Rham cohomology as a DGA
pub struct DeRhamDGA {
    pub manifold_dim: usize,
    pub dga: DGA,
}

impl DeRhamDGA {
    pub fn new(dim: usize) -> Self;  // Ωᵏ has dimension C(n,k)
    pub fn cohomology(&self) -> Cohomology;
}

// Poincaré duality check
pub fn check_poincare_duality(cohom: &Cohomology, manifold_dim: usize) -> bool;
// Exterior algebra construction
pub fn exterior_algebra_dga(n_generators: usize, generator_degrees: &[i32]) -> DGA;
```

## How It Works

### Architecture

The crate is layered: each module builds on the previous ones.

```
graded.rs           →  Graded vector spaces, Koszul signs
     ↓
chain_complex.rs    →  Chain complexes, differentials, Betti numbers
     ↓
dga.rs              →  Differential graded algebras (graded algebra + d + Leibniz)
     ↓
cohomology.rs       →  H(d) = ker(d)/im(d), representatives, Poincaré polynomial
     ↓
quasi_iso.rs        →  Quasi-isomorphisms (maps inducing cohomology isomorphisms)
     ↓
minimal_model.rs    →  Minimal models (free DGAs with trivial generator differential)
massey.rs           →  Massey products (higher-order cohomology operations)
a_infinity.rs       →  A∞-algebras (associativity up to homotopy)
derived.rs          →  Derived category (chain complexes up to quasi-isomorphism)
applications.rs     →  Concrete: de Rham DGA, Poincaré duality, exterior algebra
```

### Linear Algebra Primitives

All linear maps are stored as row-major `Vec<Vec<f64>>` matrices (the `LinearMap` type). Key operations:
- **Composition**: standard matrix multiplication
- **Rank**: Gaussian elimination with partial pivoting
- **Kernel/Image dimension**: rank-nullity theorem
- **d² = 0 check**: compose consecutive differentials and verify the result is zero

### Cohomology Computation

Betti numbers are computed via the rank-nullity theorem applied to each differential:

```
βₖ = dim ker(dₖ) - dim im(dₖ₋₁)
   = (dim Vᵏ - rank(dₖ)) - rank(dₖ₋₁)
```

Representative cocycles are found by computing a kernel basis for each differential using row reduction.

### Leibniz Rule Verification

For each pair of basis elements (a, b), the crate checks:

```
d(a·b) = da·b + (-1)^|a| · a·db
```

This is verified numerically with a configurable tolerance, ensuring the DGA structure is correct.

## The Math

### Graded Vector Spaces

A **graded vector space** is a direct sum V = ⊕ᵏ Vᵏ where each Vᵏ is a finite-dimensional vector space. Elements of Vᵏ have **degree** k. The **Koszul sign convention** introduces signs when swapping graded elements: a ⊗ b = (−1)^|a||b| b ⊗ a.

### Chain Complexes and d² = 0

A **chain complex** is a sequence of vector spaces and linear maps:

```
... → Vᵏ⁻¹ →ᵈ Vᵏ →ᵈ Vᵏ⁺¹ → ...
```

with d² = 0. This means im(dₖ₋₁) ⊆ ker(dₖ), so we can define **cohomology**:

```
Hᵏ = ker(dₖ) / im(dₖ₋₁)
```

The **Betti number** βₖ = dim Hᵏ counts the "holes" of dimension k. The **Euler characteristic** is χ = Σ (−1)ᵏ βₖ = Σ (−1)ᵏ dim Vᵏ (these are equal by the Hopf trace formula).

### Differential Graded Algebras

A **DGA** is a graded vector space with:
1. A differential d: Vᵏ → Vᵏ⁺¹ with d² = 0
2. A graded multiplication ·: Vᵖ × Vᵠ → Vᵖ⁺ᵠ
3. The **graded Leibniz rule**: d(a·b) = da·b + (−1)^|a| a·db

The Leibniz rule makes cohomology H(d) into a **graded algebra** — the cup product in topology, the wedge product in de Rham theory.

### Quasi-Isomorphisms

A **quasi-isomorphism** is a DGA morphism f: A → B that induces an isomorphism on cohomology. Quasi-isomorphic DGAs have the same "homotopy type" even if they look different algebraically. Two DGAs are quasi-isomorphic if and only if they have isomorphic minimal models.

### Minimal Models

A **minimal model** of a DGA A is a free graded commutative algebra M with a differential that is decomposable (d(M) ⊆ M⁺ · M⁺, where M⁺ is the augmentation ideal), together with a quasi-isomorphism M → A. Minimal models are unique up to isomorphism and classify DGAs up to quasi-isomorphism.

### Massey Products

When the cup product of two cohomology classes vanishes (α · β = 0 in cohomology), the **triple Massey product** ⟨α, β, γ⟩ can be nonzero. It detects higher-order linking:

```
Given: d(a₁₂) = α·β, d(a₂₃) = β·γ
Then: ⟨α, β, γ⟩ = [a₁₂·γ + (-1)^|α|·α·a₂₃] ∈ H^{|α|+|β|+|γ|-1}
```

Massey products are the simplest higher-order cohomology operations. They can distinguish spaces that have isomorphic cohomology rings but different homotopy types (e.g., the Borromean rings).

### A∞-Algebras

An **A∞-algebra** generalizes a DGA by replacing strict associativity with a family of operations mₙ: V^⊗ⁿ → V satisfying the **Stasheff identities**:

```
m₁ = d (differential)
m₂ = multiplication (associative up to m₃)
m₃ = homotopy for associativity
m₄, m₅, ... = higher coherence data
```

Every DGA is an A∞-algebra with mₙ = 0 for n ≥ 3. But A∞-algebras are strictly more general — they capture "associativity up to homotopy."

### Derived Categories

The **derived category** D(A) is the category of chain complexes over A, localized at quasi-isomorphisms. Objects are chain complexes; morphisms are "roofs" X ← Z → Y where the backward arrow is a quasi-isomorphism. Key constructions:
- **Shift**: C[n]ᵏ = Cᵏ⁺ⁿ
- **Mapping cone**: cone(f) for a chain map f
- **Truncation**: τ_≤ₙ C (kill everything above degree n)

## Testing

**129 tests** across 10 modules:

| Module | Tests | Coverage |
|---|---|---|
| `graded` | 24 | Elements, vector spaces, tensor products, shifts, Koszul signs |
| `chain_complex` | 13 | Linear maps, composition, rank, d²=0, Betti numbers |
| `dga` | 10 | Leibniz rule, associativity, d²=0, morphisms |
| `cohomology` | 11 | Betti numbers, representatives, Euler characteristic, Poincaré polynomial |
| `quasi_iso` | 7 | Quasi-isomorphism detection, cohomology comparison |
| `minimal_model` | 8 | Free algebras, exterior/polynomial, minimal model computation |
| `massey` | 7 | Triple products, indeterminacy, cocycle verification |
| `a_infinity` | 12 | Stasheff identities, DGA conversion, operations |
| `derived` | 17 | Objects, shifts, direct sums, truncation, acyclicity |
| `applications` | 20 | De Rham DGA, Poincaré duality, exterior algebra |

Run with:

```bash
cargo test
```

## License

MIT
