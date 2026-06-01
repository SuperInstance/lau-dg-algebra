# lau-dg-algebra

Differential graded algebras for agents — the algebraic structure underlying cohomology.

Every cohomology theory in the ecosystem (de Rham, sheaf, Hodge) is a DGA.

## Core Concepts

- **Graded Vector Space**: `V = ⊕Vᵏ` with degree-k components
- **Graded Multiplication**: `Vᵖ ⊗ Vᵠ → Vᵖ⁺ᵠ` (respecting grading)
- **Differential**: `d: Vᵏ → Vᵏ⁺¹` with `d² = 0`
- **Leibniz Rule**: `d(a·b) = da·b + (-1)^|a| a·db`
- **Cohomology**: `H(d) = ker(d)/im(d)` as a graded algebra
- **Quasi-isomorphism**: DGA map inducing isomorphism on cohomology
- **Minimal Models**: Free DGA with trivial differential on generators
- **Massey Products**: Higher-order cohomology operations
- **A∞-Algebras**: Generalization where associativity holds "up to homotopy"
- **Derived Category**: Chain complexes up to quasi-isomorphism

## Modules

| Module | Description |
|--------|-------------|
| `graded` | Graded vector spaces, Koszul signs, tensor products |
| `chain_complex` | Chain complexes, linear maps, Betti numbers |
| `dga` | Differential graded algebras, morphisms, Leibniz rule |
| `cohomology` | Cohomology computation, Poincaré duality, Euler characteristic |
| `quasi_iso` | Quasi-isomorphisms, homotopies |
| `minimal_model` | Minimal models, formality, free graded algebras |
| `massey` | Massey products (triple and n-ary) |
| `a_infinity` | A∞-algebras, homotopy transfer |
| `derived` | Derived category, Ext groups, distinguished triangles |
| `applications` | De Rham, sheaf, Hodge cohomology as DGAs |

## Usage

```rust
use lau_dg_algebra::{DGA, GradedVectorSpace, Cohomology, DeRhamDGA};

// De Rham cohomology of a 3-manifold
let dr = DeRhamDGA::new(3);
let h = dr.cohomology();
assert_eq!(h.betti_number(0), 1);
assert_eq!(h.betti_number(1), 3);
```

## Tests

129 tests covering all modules.

```bash
cargo test
```

## License

MIT
