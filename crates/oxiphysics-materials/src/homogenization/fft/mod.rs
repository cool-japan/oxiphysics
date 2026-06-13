// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! FFT spectral homogenization kernels (Moulinec–Suquet basic scheme and the
//! Eyre–Milton accelerated scheme), built on OxiFFT split-complex transforms.

pub mod eyre_milton;
pub mod green_operator;
pub mod lippmann_schwinger;

pub use eyre_milton::eyre_milton;
pub use green_operator::apply_gamma0;
pub use lippmann_schwinger::lippmann_schwinger;

pub(crate) use eyre_milton::eyre_milton_inner;
pub(crate) use lippmann_schwinger::lippmann_schwinger_inner;
