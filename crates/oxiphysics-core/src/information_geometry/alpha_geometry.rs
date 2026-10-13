// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Amari's α-geometry on a [`StatisticalManifold`].

use super::statistical_manifold::StatisticalManifold;
use crate::error::PhysicsResult;

/// Relative finite-difference step used to differentiate the connection
/// coefficients when assembling the curvature tensor.
const CURVATURE_STEP_SCALE: f64 = 1.0e-3;

/// α-geometry: Amari's one-parameter family of affine connections
/// `Γ^(α) = Γ^(0) − (α/2) T` on a statistical manifold, where `Γ^(0)` is the
/// Levi-Civita connection of the Fisher metric and `T` the Amari–Chentsov
/// tensor.
///
/// `α = 0` is the Levi-Civita connection, `α = 1` the exponential connection
/// and `α = −1` the mixture connection; `α` and `−α` are dual with respect to
/// the Fisher metric.  Exponential families are flat for `α = ±1`.
#[derive(Debug, Clone)]
pub struct AlphaGeometry {
    /// The α parameter controlling the connection.
    pub alpha: f64,
    manifold: StatisticalManifold,
}

impl AlphaGeometry {
    /// Create the α-geometry with parameter `alpha` on `manifold`.
    pub fn new(alpha: f64, manifold: StatisticalManifold) -> Self {
        Self { alpha, manifold }
    }

    /// The underlying statistical manifold.
    pub fn manifold(&self) -> &StatisticalManifold {
        &self.manifold
    }

    /// α-connection coefficients `Γ^(α)k_{ij}` at `params`, indexed
    /// `[k][i][j]`.
    pub fn alpha_connection(&self, params: &[f64]) -> PhysicsResult<Vec<Vec<Vec<f64>>>> {
        self.manifold.alpha_christoffel_symbols(params, self.alpha)
    }

    /// Coefficients of the dual `(−α)`-connection at `params`, indexed
    /// `[k][i][j]`.
    pub fn dual_connection(&self, params: &[f64]) -> PhysicsResult<Vec<Vec<Vec<f64>>>> {
        self.manifold.alpha_christoffel_symbols(params, -self.alpha)
    }

    /// Riemann curvature tensor `R^l_{kij}` of the α-connection, indexed
    /// `[l][k][i][j]`:
    ///
    /// `R^l_{kij} = ∂_i Γ^l_{jk} − ∂_j Γ^l_{ik} + Γ^l_{im} Γ^m_{jk} − Γ^l_{jm} Γ^m_{ik}`,
    ///
    /// i.e. the `l`-component of `R(∂_i, ∂_j) ∂_k`.  Derivatives of the
    /// connection coefficients are taken by central differences.
    pub fn curvature_tensor(&self, params: &[f64]) -> PhysicsResult<Vec<Vec<Vec<Vec<f64>>>>> {
        let n = params.len();
        let gamma = self.alpha_connection(params)?;
        // dgamma[m][l][i][j] = ∂_m Γ^l_{ij}
        let mut dgamma = vec![vec![vec![vec![0.0f64; n]; n]; n]; n];
        for (m, dgamma_m) in dgamma.iter_mut().enumerate() {
            let h = CURVATURE_STEP_SCALE * params[m].abs().max(1.0);
            let mut pm = params.to_vec();
            let mut mm = params.to_vec();
            pm[m] += h;
            mm[m] -= h;
            let gp = self.alpha_connection(&pm)?;
            let gm = self.alpha_connection(&mm)?;
            for l in 0..n {
                for i in 0..n {
                    for j in 0..n {
                        dgamma_m[l][i][j] = (gp[l][i][j] - gm[l][i][j]) / (2.0 * h);
                    }
                }
            }
        }
        let mut r = vec![vec![vec![vec![0.0f64; n]; n]; n]; n];
        for l in 0..n {
            for k in 0..n {
                for i in 0..n {
                    for j in 0..n {
                        let term1 = dgamma[i][l][j][k];
                        let term2 = dgamma[j][l][i][k];
                        let mut term3 = 0.0f64;
                        let mut term4 = 0.0f64;
                        for (mm, &g_lim) in gamma[l][i].iter().enumerate() {
                            term3 += g_lim * gamma[mm][j][k];
                            term4 += gamma[l][j][mm] * gamma[mm][i][k];
                        }
                        r[l][k][i][j] = term1 - term2 + term3 - term4;
                    }
                }
            }
        }
        Ok(r)
    }

    /// Sectional curvature `⟨R(∂_i, ∂_j)∂_j, ∂_i⟩ / (g_ii g_jj − g_ij²)` of the
    /// α-connection in the coordinate plane `(i, j)` at `params`, with the
    /// first index lowered by the Fisher metric.
    pub fn sectional_curvature(&self, params: &[f64], i: usize, j: usize) -> PhysicsResult<f64> {
        let n = params.len();
        if i >= n || j >= n || i == j {
            return Err(crate::error::PhysicsError::InvalidInput {
                field: "i, j".to_string(),
                reason: "need two distinct coordinate indices below dim".to_string(),
            });
        }
        let r = self.curvature_tensor(params)?;
        let g = self.manifold.fisher_metric(params)?;
        let num: f64 = (0..n).map(|m| g[i][m] * r[m][j][i][j]).sum();
        let den = g[i][i] * g[j][j] - g[i][j] * g[i][j];
        if den <= 0.0 {
            return Err(crate::error::PhysicsError::InvalidInput {
                field: "params".to_string(),
                reason: "Fisher metric is degenerate on the coordinate plane".to_string(),
            });
        }
        Ok(num / den)
    }
}

#[cfg(test)]
#[path = "alpha_geometry_tests.rs"]
mod tests;
