# lau-dg-algebra

**Differential graded algebras for agents — the algebraic structure underlying cohomology.**

A Rust implementation of differential graded algebras (DGAs) and the full machinery of homological algebra: chain complexes, cohomology, Massey products, minimal models, A∞-algebras, quasi-isomorphisms, derived categories, and applications to de Rham / sheaf cohomology. **129 tests, all passing.**

## What This Does

| Module | What it is | Key structures |
|---|---|---|
| **Graded vector spaces** | `V = ⊕ Vᵏ` with degree-k components | Tensor product, shift, Euler characteristic |
| **Chain complexes** | `... → Vⁿ → Vⁿ⁺¹ → ...` with d² = 0 | Betti numbers, Euler characteristic, mapping cone |
| **DGA** | Graded algebra + differential d + Leibniz rule | Associativity, graded commutativity, unit |
| **Cohomology** | `H(d) = ker(d)/im(d)` as graded algebra | Poincaré duality, cup product, Euler characteristic |
| **Massey products** | Higher-order cohomology operations | ⟨x, y, z⟩ when x∪y = 0, y∪z = 0 |
| **Minimal models** | DGA with d = 0, unique up to isomorphism | Sullivan minimal model, formality check |
| **A∞-algebras** | Associativity "up to homotopy" | mₙ operations, homotopy transfer |
| **Quasi-isomorphisms** | Chain maps inducing cohomology isomorphisms | Homotopy equivalence, mapping cone criterion |
| **Derived category** | Chain complexes up to quasi-isomorphism | Ext groups, Verdier quotient, distinguished triangles |
| **Applications** | De Rham, sheaf, Hodge as DGAs | Wedge product, Čech complex, Hodge decomposition |

## Key Idea

A **differential graded algebra** is the algebraic structure that captures cohomology. It's a graded vector space with:
1. A **differential** `d: Vᵏ → Vᵏ⁺¹` with `d² = 0`
2. A **multiplication** `Vᵖ × Vᵠ → Vᵖ⁺ᵠ` satisfying the graded Leibniz rule

From this single structure, you get cohomology rings, characteristic classes, and all of homological algebra. Every cohomology theory — de Rham, sheaf, singular, K-theory — is a DGA.

The provocative claim: agents operating on graded information naturally form DGAs. Their "learning differential" and "composition operation" satisfy the same axioms.

## Install

```toml
[dependencies]
lau-dg-algebra = "0.1.0"
```

Or clone directly:

```bash
git clone https://github.com/SuperInstance/lau-dg-algebra.git
cd lau-dg-algebra
cargo test    # 129 tests pass
```

### Dependencies

- `nalgebra` 0.33 — linear algebra
- `serde` — serialization

## Quick Start

```rust
use lau_dg_algebra::{
    graded::GradedVectorSpace,
    chain_complex::{ChainComplex, LinearMap},
    dga::{DGA, MultiplicationTable, DGAMorphism},
    cohomology::Cohomology,
};
use std::collections::HashMap;

// Build a graded vector space: V⁰ (dim 1) ⊕ V¹ (dim 1)
let mut dims = HashMap::new();
dims.insert(0, 1);
dims.insert(1, 1);
let gvs = GradedVectorSpace::new(dims);

// Exterior algebra: 1·1=1, 1·x=x, x·1=x, x·x=0
let mut tables = HashMap::new();
tables.insert((0, 0), LinearMap::from_rows(vec![vec![1.0]]));
tables.insert((0, 1), LinearMap::from_rows(vec![vec![1.0]]));
tables.insert((1, 0), LinearMap::from_rows(vec![vec![1.0]]));
tables.insert((1, 1), LinearMap::from_rows(vec![vec![0.0]])));

let dga = DGA::new(gvs, HashMap::new(), MultiplicationTable::new(tables));

// Verify DGA axioms
assert!(dga.check_d_squared_zero(1e-10));
assert!(dga.check_associativity(0, 0, 0, 1e-10));
assert!(dga.check_graded_commutativity(0, 1, 1e-10));
assert!(dga.check_unit(1e-10));

// Compute cohomology
let h = Cohomology::from_dga(&dga);
assert_eq!(h.betti_number(0), 1);
assert_eq!(h.betti_number(1), 1);
assert_eq!(h.euler_characteristic(), 0);
```

## API Reference

### `GradedVectorSpace` — V = ⊕ Vᵏ

```rust
let gvs = GradedVectorSpace::new(dimensions);  // HashMap<i32, usize>
gvs.dim(k);                    // dim Vᵏ
gvs.total_dim();               // Σ dim Vᵏ
let sum = gvs.direct_sum(&other);
let tensor = gvs.tensor_product(&other);  // (V⊗W)ⁿ = ⊕_{p+q=n} Vᵖ⊗Wᵠ
let shifted = gvs.shift(n);               // V[n]ᵏ = Vᵏ⁻ⁿ
let trunc = gvs.truncate(lo, hi);
let chi = gvs.euler_characteristic();     // Σ (−1)ᵏ dim Vᵏ
```

### `LinearMap` — Matrix stored as Vec<Vec<f64>>

```rust
let m = LinearMap::from_rows(vec![vec![1.0, 2.0], vec![3.0, 4.0]]);
m.apply(&vec![1.0, 0.0]);      // matrix-vector multiply
let c = m.compose(&other);     // matrix composition
let t = m.transpose();
let r = m.rank();              // row reduction
let k = m.kernel_dim();        // domain_dim - rank
```

### `ChainComplex` — ... → Vⁿ → Vⁿ⁺¹ → ...

```rust
let cc = ChainComplex::new(dimensions, differentials);
cc.check_d_squared_zero(tol);   // the fundamental axiom
let betti = cc.betti_numbers(); // βₖ = dim ker(dₖ) − dim im(dₖ₋₁)
let chi = cc.euler_characteristic();
let shifted = cc.shift(n);
let sum = cc.direct_sum(&other);
```

### `DGA` — Differential graded algebra

```rust
let dga = DGA::new(graded_vs, differential, multiplication);
dga.check_d_squared_zero(tol);               // d² = 0
dga.check_leibniz(deg_a, deg_b, tol);        // graded Leibniz rule
dga.check_associativity(p, q, r, tol);       // (ab)c = a(bc)
dga.check_graded_commutativity(p, q, tol);   // ab = (−1)^{|a||b|} ba
dga.check_unit(tol);                          // 1·a = a = a·1
let cc = dga.as_chain_complex();
```

### `DGAMorphism` — Maps between DGAs

```rust
let morphism = DGAMorphism::new(maps);
morphism.check_chain_map(&source, &target, tol);
morphism.check_algebra_map(&source, &target, p, q, tol);
```

### `Cohomology` — H(d) = ker(d)/im(d)

```rust
let h = Cohomology::from_chain_complex(&cc);
// or
let h = Cohomology::from_dga(&dga);

h.betti_number(k);
h.total_betti();
h.euler_characteristic();                    // Σ (−1)ᵏ βₖ
h.check_poincare_duality(dim);               // βₖ = β_{n−k}
h.is_trivial();
h.nonzero_degrees();
let cup = h.cup_product(&dga, deg_a, &a, deg_b, &b);
```

### `MasseyProduct` — Higher cohomology operations

```rust
use lau_dg_algebra::massey::*;

let (exists, product) = massey_product(&dga, deg_x, &x, deg_y, &y, deg_z, &z, tol);
```

When `x∪y = 0` and `y∪z = 0`, the Massey product `⟨x, y, z⟩` is a well-defined element of `H^{|x|+|y|+|z|−1}` that detects higher-order linking.

### `MinimalModel` — DGA with d = 0

```rust
use lau_dg_algebra::minimal_model::*;

let (minimal, morphism) = compute_minimal_model(&dga, tol);
let is_formal = check_formality(&dga, tol);  // is it quasi-isomorphic to its cohomology?
```

### `AInfinityAlgebra` — Associativity up to homotopy

```rust
use lau_dg_algebra::a_infinity::*;

let ainf = AInfinityAlgebra::new(dimensions, operations);
let m1 = ainf.m(1);               // differential
let m2 = ainf.m(2);               // multiplication
ainf.check_relation(n, tol);       // A∞ relations
let transferred = ainf.transfer(&proj, &incl, &homotopy, &dims);
```

### `QuasiIsomorphism` — Chain maps inducing isomorphisms on cohomology

```rust
use lau_dg_algebra::quasi_iso::*;

let is_qi = is_quasi_isomorphism(&source, &target, &chain_map, tol);
let is_he = is_homotopy_equivalence(&source, &target, &map, &inverse, &homotopy, tol);
```

### `DerivedCategory` — Chain complexes up to quasi-isomorphism

```rust
use lau_dg_algebra::derived::*;

let mut cat = DerivedCategory::new("D(Ab)");
let idx = cat.add_object(obj);
cat.add_morphism(source, target, map);
let ext = cat.ext_group(a, b, n);            // Extⁿ(A, B)
let iso = cat.are_isomorphic(a, b, tol);
let quotient = cat.verdier_quotient(&sub_indices);
let triangle = DistinguishedTriangle::new(a, b, c, f, g, h);
let rotated = triangle.rotate();             // B → C → A[1] → B[1]
```

### Applications: De Rham, Sheaf, Hodge

```rust
use lau_dg_algebra::applications::*;

// De Rham cohomology as a DGA
let dr = DeRhamDGA::new(2);     // 2-manifold
// Basis: Ω⁰ (dim 1), Ω¹ (dim 2), Ω² (dim 1)
// Wedge product with correct signs
dr.cohomology();                 // Betti numbers
dr.poincare_polynomial();        // {0: 1, 1: 2, 2: 1}

// Sheaf cohomology (Čech complex)
let sheaf = SheafDGA::new(3);   // 3-set cover
// Čech⁰ = C(3,1) = 3, Čech¹ = C(3,2) = 3, Čech² = C(3,3) = 1

// Hodge theory
let hodge = HodgeDecomposition::for_dimension(3);
hodge.betti_numbers();           // harmonic form dimensions

// Unified view
let unified = UnifiedCohomology::new()
    .with_de_rham(2)
    .with_sheaf(3)
    .with_hodge(2);
unified.check_compatibility();   // all theories give same Betti numbers
unified.euler_characteristic();  // χ = 1 − 2 + 1 = 0
```

## How It Works

### Linear Algebra

All maps are stored as `Vec<Vec<f64>>` (row-major matrices). Rank is computed via Gaussian elimination with partial pivoting. Kernel basis is found via augmented matrix `[M | I]` row reduction.

### DGA Axioms

The crate exhaustively verifies DGA axioms by testing **all basis elements**:
- `d² = 0`: compose consecutive differentials, check zero
- Leibniz: test every basis pair `(eᵢ, eⱼ)`, verify `d(eᵢ·eⱼ) = d(eᵢ)·eⱼ + (−1)^|eᵢ| eᵢ·d(eⱼ)`
- Associativity: test every basis triple `(eᵢ, eⱼ, eₖ)`
- Graded commutativity: test every pair with Koszul sign `(−1)^{|a||b|}`

### Wedge Product

For the de Rham DGA, k-forms correspond to subsets of `{0, ..., n−1}` of size k. The wedge product is computed as a union with sign determined by the number of transpositions (counted via bitwise operations on subset masks).

### Cohomology

Betti numbers: `βₖ = dim ker(dₖ) − dim im(dₖ₋₁)`. Representatives are found by computing kernel basis vectors of each differential.

## The Math

### Differential Graded Algebras

A DGA is a graded algebra `(A, ·)` with a differential `d` satisfying:
1. `d² = 0` (cochain complex)
2. `d(a·b) = da·b + (−1)^{|a|} a·db` (graded Leibniz rule)

### Chain Complexes

A chain complex is a sequence `... → Cₙ₊₁ → Cₙ → Cₙ₋₁ → ...` with `d² = 0`. The **homology** is `Hₙ = ker(dₙ)/im(dₙ₊₁)`.

### Massey Products

When `x∪y = 0`, choose a chain `a` with `da = x∪y`. The Massey product `⟨x, y, z⟩` is the class of `a∪z − (−1)^{|x|} x∪b` (where `db = y∪z`) in `H^{|x|+|y|+|z|−1}`.

### Minimal Models

A minimal model of a DGA `A` is a DGA `M` with `d_M = 0` and a quasi-isomorphism `M → A`. A DGA is **formal** if it's quasi-isomorphic to its cohomology algebra.

### A∞-Algebras

An A∞-algebra has operations `mₙ: V^{⊗ⁿ} → V` for n ≥ 1 satisfying the Stasheff relations. The key relation: `m₁² = 0` (differential), and `m₂` is associative up to `m₃` (homotopy). Every DGA is an A∞-algebra with `mₙ = 0` for n ≥ 3.

### Derived Categories

The derived category `D(A)` has chain complexes as objects and chain maps (localized at quasi-isomorphisms) as morphisms. Key structures:
- **Distinguished triangles**: `A → B → C → A[1]` (analogue of exact sequences)
- **Ext groups**: `Extⁿ(A, B) = Hom_{D(A)}(A, B[n])`
- **Verdier quotient**: `D(A)/B` (localization at a subcategory)

## Test Coverage

**129 tests**, all passing:

| Module | Tests | What's tested |
|---|---|---|
| `graded` | 18 | Element creation, scaling, vector space ops, tensor product, shift, Euler characteristic, Koszul sign |
| `chain_complex` | 14 | Linear maps, composition, rank, transpose, d²=0, Betti numbers, Euler characteristic, direct sum, shift |
| `dga` | 12 | d²=0, associativity, graded commutativity, unit, Leibniz rule, multiplication, morphisms |
| `cohomology` | 11 | Trivial/exact/circle complexes, Euler characteristic, Poincaré duality, triviality, kernel basis |
| `massey` | ~12 | Massey product computation, trivial products, formality detection |
| `minimal_model` | ~12 | Minimal model construction, formality, morphism verification |
| `a_infinity` | 14 | A∞ creation, m₁²=0, strict DGA check, minimality, relations, operations, transfer, morphisms |
| `quasi_iso` | ~8 | Quasi-isomorphism detection, homotopy equivalence |
| `derived` | 18 | Objects, cohomology, shift, direct sum, acyclicity, truncation, Ext groups, Verdier quotient, triangles |
| `applications` | 21 | De Rham (dims 0-3), wedge product, anticommutativity, sheaf, Hodge, unified cohomology, compatibility |

## License

MIT
