// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Matrix-free sum-factorized high-order operators (Kronbichler-Kormann 2012).
//!
//! Provides element-local tensor-product evaluation of high-order Poisson
//! stiffness operators on regular hexahedral meshes, replacing the global
//! assembled sparse mat-vec with sum-factorized contractions.

pub mod operator;
pub mod sum_factorization;

pub use operator::{MatrixFreeError, MatrixFreeOperator, MatrixFreePcg};
pub use sum_factorization::SumFactPoisson;
