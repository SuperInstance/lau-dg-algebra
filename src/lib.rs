//! # lau-dg-algebra
//!
//! Differential graded algebras (DGAs) for agents — the algebraic structure
//! underlying cohomology. Every cohomology theory in the ecosystem (de Rham,
//! sheaf, Hodge) is a DGA.
//!
//! ## Core Concepts
//!
//! - **Graded Vector Space**: `V = ⊕Vᵏ` with degree-k components
//! - **Graded Multiplication**: `Vᵖ ⊗ Vᵠ → Vᵖ⁺ᵠ`
//! - **Differential**: `d: Vᵏ → Vᵏ⁺¹` with `d² = 0`
//! - **Leibniz Rule**: `d(a·b) = da·b + (-1)^|a| a·db`
//! - **Cohomology**: `H(d) = ker(d)/im(d)` as a graded algebra
//! - **Quasi-isomorphism**: DGA map inducing isomorphism on cohomology
//! - **Minimal Models**: Free DGA with trivial differential on generators
//! - **Massey Products**: Higher-order cohomology operations
//! - **A∞-Algebras**: Associativity holds "up to homotopy"
//! - **Derived Category**: Chain complexes up to quasi-isomorphism

pub mod graded;
pub mod chain_complex;
pub mod dga;
pub mod cohomology;
pub mod quasi_iso;
pub mod minimal_model;
pub mod massey;
pub mod a_infinity;
pub mod derived;
pub mod applications;

pub use graded::*;
pub use chain_complex::*;
pub use dga::*;
pub use cohomology::*;
pub use quasi_iso::*;
pub use minimal_model::*;
pub use massey::*;
pub use a_infinity::*;
pub use derived::*;
pub use applications::*;
