// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Surface reconstruction from oriented point clouds.
//!
//! This module implements screened-Poisson surface reconstruction following
//! Kazhdan & Hoppe, "Screened Poisson Surface Reconstruction" (ACM TOG 2013).
//!
//! The implementation uses the regular-grid variant: oriented normals are
//! splatted onto a uniform voxel grid via trilinear weights to form a vector
//! field `V`, the screened-Poisson system `(-Δ + α·S) χ = ∇·V` is assembled in
//! compressed-sparse-row form, and solved with a Jacobi-preconditioned
//! conjugate-gradient iteration. The resulting implicit indicator field `χ` is
//! polygonised with marching cubes at an iso-value taken from the average of
//! `χ` evaluated at the input samples.

pub mod poisson;

pub use poisson::{OrientedPoint, PoissonError, PoissonSurface, screened_poisson_reconstruct};
