// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Buoyancy-driven convection: Boussinesq mantle convection and double-diffusive
//! thermohaline circulation.
//!
//! In both models the temperature (and salinity) are passive-scalar
//! advection–diffusion populations that are advected by the *fluid* velocity
//! obtained from the momentum populations, and the fluid feels the scalars
//! through a Boussinesq body force applied with the Guo scheme.

use super::{CS2, NQ, OPP, W, feq, geq_linear, guo_force, macroscopic};
use super::{stream_periodic, stream_walled_y};

/// Temperature of the hot bottom wall (the layer temperature difference is 1).
const T_BOTTOM: f64 = 1.0;
/// Temperature of the cold top wall.
const T_TOP: f64 = 0.0;
/// Reference (mean) temperature of the Boussinesq buoyancy.
const T_REF: f64 = 0.5;

/// Zou–He no-slip condition on a bottom wall node: reconstruct the unknown
/// populations 2, 5, 6 (pointing into the fluid) for zero wall velocity.
fn zou_he_noslip_bottom(f: &mut [f64; NQ]) {
    let shear = 0.5 * (f[1] - f[3]);
    f[2] = f[4];
    f[5] = f[7] - shear;
    f[6] = f[8] + shear;
}

/// Zou–He no-slip condition on a top wall node: reconstruct the unknown
/// populations 4, 7, 8 for zero wall velocity.
fn zou_he_noslip_top(f: &mut [f64; NQ]) {
    let shear = 0.5 * (f[1] - f[3]);
    f[4] = f[2];
    f[8] = f[6] - shear;
    f[7] = f[5] + shear;
}

/// Dirichlet temperature `t_wall` on a wall node whose populations `unknown`
/// point into the fluid: anti-bounce-back of the populations that streamed in,
/// `g_q = 2 w_q T_w − g_opp(q)`, followed by a weight-proportional correction of
/// the unknown ones so that the zeroth moment equals `t_wall` exactly.
fn dirichlet_wall(g: &mut [f64; NQ], unknown: [usize; 3], t_wall: f64) {
    for &q in &unknown {
        g[q] = 2.0 * W[q] * t_wall - g[OPP[q]];
    }
    let missing = t_wall - g.iter().sum::<f64>();
    let weight_sum: f64 = unknown.iter().map(|&q| W[q]).sum();
    for &q in &unknown {
        g[q] += missing * W[q] / weight_sum;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MantelConvection
// ─────────────────────────────────────────────────────────────────────────────

/// Rayleigh–Bénard (Boussinesq) convection model for Earth's mantle.
///
/// Uses a double-distribution LBM: D2Q9 populations `f` for momentum and
/// D2Q9 populations `g` (advection–diffusion equilibrium, linear in the
/// velocity) for the temperature `T ∈ [0, 1]`.  The scalar is advected by the
/// fluid velocity computed from `f`; the fluid feels the temperature through
/// the buoyancy force density `F_y = ρ g β (T − 1/2)` applied with the Guo
/// scheme.
///
/// The layer is bounded by the wall nodes `iy = 0` (hot, `T = 1`) and
/// `iy = ny − 1` (cold, `T = 0`), both no-slip (Zou–He) and isothermal, and is
/// periodic in x.  The layer thickness is therefore `H = ny − 1` lattice
/// units.  The Rayleigh number `Ra = g β ΔT H³ / (ν κ)` with `ΔT = 1`, the
/// kinematic viscosity `ν = (1/ω_f − 1/2)/3` and the thermal diffusivity
/// `κ = (1/ω_t − 1/2)/3` fixes the lattice buoyancy
/// `g β = Ra ν κ / H³` ([`buoyancy_coefficient`](Self::buoyancy_coefficient)).
pub struct MantelConvection {
    /// Grid width.
    pub nx: usize,
    /// Grid height (`H = ny − 1` is the layer thickness).
    pub ny: usize,
    /// Momentum relaxation frequency.
    pub omega_f: f64,
    /// Temperature relaxation frequency.
    pub omega_t: f64,
    /// Rayleigh number.
    pub ra: f64,
    /// Nominal Prandtl number supplied at construction (informational).  The
    /// Prandtl number of the simulation is `ν/κ`, fixed by `omega_f` and
    /// `omega_t`; see [`prandtl_number`](Self::prandtl_number).
    pub pr: f64,
    /// Momentum distribution functions.
    pub f: Vec<Vec<[f64; NQ]>>,
    /// Temperature distribution functions (advection–diffusion D2Q9).
    pub g: Vec<Vec<[f64; NQ]>>,
    /// Temperature field (zeroth moment of `g`).
    pub temp: Vec<Vec<f64>>,
}

impl MantelConvection {
    /// Create a Boussinesq mantle convection simulation at rest with
    /// `T = 1/2` everywhere.
    ///
    /// The Rayleigh number `Ra` controls convective vigour; for Earth's mantle
    /// `Ra ~ 10^6 – 10^8`.  The effective Prandtl number of the mantle is
    /// astronomically large, but an LBM run uses `Pr ~ 1 – 100`, which is set
    /// through the ratio of the two relaxation times (`pr` merely records the
    /// intended value).
    pub fn new(nx: usize, ny: usize, omega_f: f64, omega_t: f64, ra: f64, pr: f64) -> Self {
        let f = vec![vec![feq(1.0, 0.0, 0.0); ny]; nx];
        let g = vec![vec![geq_linear(T_REF, 0.0, 0.0); ny]; nx];
        let temp = vec![vec![T_REF; ny]; nx];
        Self {
            nx,
            ny,
            omega_f,
            omega_t,
            ra,
            pr,
            f,
            g,
            temp,
        }
    }

    /// Set a conductive profile (linear temperature 1 at the bottom wall to 0
    /// at the top wall).
    pub fn conductive_ic(&mut self) {
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let t = T_BOTTOM - (T_BOTTOM - T_TOP) * iy as f64 / self.layer_height();
                self.temp[ix][iy] = t;
                self.g[ix][iy] = geq_linear(t, 0.0, 0.0);
            }
        }
    }

    /// Kinematic viscosity `ν = (1/ω_f − 1/2)/3` (lattice units).
    pub fn viscosity(&self) -> f64 {
        (1.0 / self.omega_f - 0.5) * CS2
    }

    /// Thermal diffusivity `κ = (1/ω_t − 1/2)/3` (lattice units).
    pub fn diffusivity(&self) -> f64 {
        (1.0 / self.omega_t - 0.5) * CS2
    }

    /// Prandtl number `ν/κ` implied by the two relaxation frequencies.
    pub fn prandtl_number(&self) -> f64 {
        self.viscosity() / self.diffusivity()
    }

    /// Layer thickness `H = ny − 1` between the two wall nodes.
    pub fn layer_height(&self) -> f64 {
        self.ny.saturating_sub(1).max(1) as f64
    }

    /// Lattice buoyancy `g β = Ra ν κ / H³` per unit temperature difference.
    pub fn buoyancy_coefficient(&self) -> f64 {
        let h = self.layer_height();
        self.ra * self.viscosity() * self.diffusivity() / (h * h * h)
    }

    /// Recompute the cached temperature field from the zeroth moment of `g`.
    fn refresh_temperature(&mut self) {
        for (temp_row, g_row) in self.temp.iter_mut().zip(self.g.iter()) {
            for (t, g_cell) in temp_row.iter_mut().zip(g_row.iter()) {
                *t = g_cell.iter().sum();
            }
        }
    }

    /// Advance one coupled timestep: collision (momentum with buoyancy,
    /// temperature advected by the fluid), periodic streaming, then the wall
    /// boundary conditions.
    pub fn step(&mut self) {
        let (nx, ny) = (self.nx, self.ny);
        let (omega_f, omega_t) = (self.omega_f, self.omega_t);
        let gbeta = self.buoyancy_coefficient();
        self.refresh_temperature();
        for ix in 0..nx {
            for iy in 0..ny {
                let wall = iy == 0 || iy + 1 == ny;
                let (rho, ux0, uy0) = macroscopic(&self.f[ix][iy]);
                let t_local = self.temp[ix][iy];
                // Boussinesq buoyancy; wall nodes are no-slip and carry no force.
                let fy = if wall {
                    0.0
                } else {
                    rho * gbeta * (t_local - T_REF)
                };
                // Force-corrected fluid velocity u = (Σ c f + F/2) / ρ.
                let (ux, uy) = if wall || rho <= 1e-12 {
                    (0.0, 0.0)
                } else {
                    (ux0, uy0 + 0.5 * fy / rho)
                };
                let feq_arr = feq(rho, ux, uy);
                let cell = &mut self.f[ix][iy];
                for (f_q, feq_q) in cell.iter_mut().zip(feq_arr.iter()) {
                    *f_q += omega_f * (*feq_q - *f_q);
                }
                guo_force(cell, 0.0, fy, ux, uy, omega_f);
                // The temperature is advected by the fluid velocity.
                let geq_arr = geq_linear(t_local, ux, uy);
                let g_cell = &mut self.g[ix][iy];
                for (g_q, geq_q) in g_cell.iter_mut().zip(geq_arr.iter()) {
                    *g_q += omega_t * (*geq_q - *g_q);
                }
            }
        }
        self.f = stream_periodic(&self.f, nx, ny);
        self.g = stream_periodic(&self.g, nx, ny);
        // Isothermal no-slip walls.
        if ny >= 3 {
            for ix in 0..nx {
                zou_he_noslip_bottom(&mut self.f[ix][0]);
                zou_he_noslip_top(&mut self.f[ix][ny - 1]);
                dirichlet_wall(&mut self.g[ix][0], [2, 5, 6], T_BOTTOM);
                dirichlet_wall(&mut self.g[ix][ny - 1], [4, 7, 8], T_TOP);
            }
        }
        self.refresh_temperature();
    }

    /// Mean kinetic energy density `½ ρ |u|²` over the fluid nodes (the
    /// interior rows between the two wall nodes).
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        let mut count = 0usize;
        for ix in 0..self.nx {
            for iy in 1..self.ny.saturating_sub(1) {
                let (rho, ux, uy) = macroscopic(&self.f[ix][iy]);
                ke += 0.5 * rho * (ux * ux + uy * uy);
                count += 1;
            }
        }
        if count == 0 { 0.0 } else { ke / count as f64 }
    }

    /// Nusselt number from the conductive heat flux through the bottom wall,
    /// `Nu = −(∂T/∂y)_wall H / ΔT` with a one-sided second-order gradient and
    /// `ΔT = 1`; pure conduction gives `Nu = 1`.
    pub fn nusselt_number(&self) -> f64 {
        let ny = self.ny;
        if ny < 2 {
            return 0.0;
        }
        let flux: f64 = (0..self.nx)
            .map(|ix| {
                let t0 = self.temp[ix][0];
                let t1 = self.temp[ix][1];
                if ny >= 3 {
                    let t2 = self.temp[ix][2];
                    0.5 * (3.0 * t0 - 4.0 * t1 + t2)
                } else {
                    t0 - t1
                }
            })
            .sum::<f64>()
            / self.nx as f64;
        flux * self.layer_height() / (T_BOTTOM - T_TOP)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ThermohalineCirculation
// ─────────────────────────────────────────────────────────────────────────────

/// Double-diffusive thermohaline circulation model.
///
/// Two scalar fields (temperature T, salinity S) drive density via the
/// relative anomaly `ρ' = −α (T − 1/2) + β (S − 1/2)`.  Uses three coupled
/// D2Q9 distributions: momentum `f`, temperature `g_t` and salinity `g_s`.
/// Both scalars are advected by the fluid velocity of `f`, and the fluid
/// feels them through the Boussinesq body force `F_y = −g ρ ρ'` (`g` is the
/// field [`gravity`](Self::gravity)).
///
/// The basin is periodic in x and closed at the top and bottom by half-way
/// bounce-back walls: no-slip for the flow and no-flux for temperature and
/// salinity.
pub struct ThermohalineCirculation {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Momentum relaxation rate.
    pub omega_f: f64,
    /// Thermal diffusion relaxation rate.
    pub omega_t: f64,
    /// Saline diffusion relaxation rate.
    pub omega_s: f64,
    /// Thermal expansion coefficient.
    pub alpha_t: f64,
    /// Haline contraction coefficient.
    pub beta_s: f64,
    /// Gravitational acceleration (lattice units) scaling the density-anomaly
    /// buoyancy; `0.01` unless changed (see [`with_gravity`](Self::with_gravity)).
    pub gravity: f64,
    /// Momentum distribution.
    pub f: Vec<Vec<[f64; NQ]>>,
    /// Temperature distribution.
    pub g_t: Vec<Vec<[f64; NQ]>>,
    /// Salinity distribution.
    pub g_s: Vec<Vec<[f64; NQ]>>,
}

impl ThermohalineCirculation {
    /// Construct a thermohaline circulation model at rest with
    /// `T = S = 1/2` and the default gravity `0.01`.
    pub fn new(
        nx: usize,
        ny: usize,
        omega_f: f64,
        omega_t: f64,
        omega_s: f64,
        alpha_t: f64,
        beta_s: f64,
    ) -> Self {
        let f = vec![vec![feq(1.0, 0.0, 0.0); ny]; nx];
        let g_t = vec![vec![geq_linear(0.5, 0.0, 0.0); ny]; nx];
        let g_s = vec![vec![geq_linear(0.5, 0.0, 0.0); ny]; nx];
        Self {
            nx,
            ny,
            omega_f,
            omega_t,
            omega_s,
            alpha_t,
            beta_s,
            gravity: 0.01,
            f,
            g_t,
            g_s,
        }
    }

    /// Builder: set the gravitational acceleration of the buoyancy force.
    pub fn with_gravity(mut self, gravity: f64) -> Self {
        self.gravity = gravity;
        self
    }

    /// Advance one step: collision (momentum with density buoyancy, both
    /// scalars advected by the fluid), then streaming with closed walls at the
    /// top and bottom.
    pub fn step(&mut self) {
        let (nx, ny) = (self.nx, self.ny);
        let (omega_f, omega_t, omega_s) = (self.omega_f, self.omega_t, self.omega_s);
        let (alpha_t, beta_s, gravity) = (self.alpha_t, self.beta_s, self.gravity);
        for ix in 0..nx {
            for iy in 0..ny {
                let (rho, ux0, uy0) = macroscopic(&self.f[ix][iy]);
                let (t, _, _) = macroscopic(&self.g_t[ix][iy]);
                let (s, _, _) = macroscopic(&self.g_s[ix][iy]);
                // Boussinesq buoyancy from the relative density anomaly.
                let rho_prime = -alpha_t * (t - T_REF) + beta_s * (s - T_REF);
                let fy = -gravity * rho * rho_prime;
                // Force-corrected fluid velocity u = (Σ c f + F/2) / ρ.
                let (ux, uy) = if rho > 1e-12 {
                    (ux0, uy0 + 0.5 * fy / rho)
                } else {
                    (0.0, 0.0)
                };
                let feq_arr = feq(rho, ux, uy);
                let cell = &mut self.f[ix][iy];
                for (f_q, feq_q) in cell.iter_mut().zip(feq_arr.iter()) {
                    *f_q += omega_f * (*feq_q - *f_q);
                }
                guo_force(cell, 0.0, fy, ux, uy, omega_f);
                // Temperature and salinity are advected by the fluid velocity.
                let geq_t = geq_linear(t, ux, uy);
                for (g_q, geq_q) in self.g_t[ix][iy].iter_mut().zip(geq_t.iter()) {
                    *g_q += omega_t * (*geq_q - *g_q);
                }
                let geq_s = geq_linear(s, ux, uy);
                for (g_q, geq_q) in self.g_s[ix][iy].iter_mut().zip(geq_s.iter()) {
                    *g_q += omega_s * (*geq_q - *g_q);
                }
            }
        }
        self.f = stream_walled_y(&self.f, nx, ny, &OPP);
        self.g_t = stream_walled_y(&self.g_t, nx, ny, &OPP);
        self.g_s = stream_walled_y(&self.g_s, nx, ny, &OPP);
    }

    /// Return density anomaly field `ρ' = -αT' + βS'`.
    pub fn density_anomaly(&self) -> Vec<Vec<f64>> {
        let mut out = vec![vec![0.0f64; self.ny]; self.nx];
        for (out_row, (g_t_row, g_s_row)) in
            out.iter_mut().zip(self.g_t.iter().zip(self.g_s.iter()))
        {
            for (out_val, (g_t_cell, g_s_cell)) in
                out_row.iter_mut().zip(g_t_row.iter().zip(g_s_row.iter()))
            {
                let (t, _a, _b) = macroscopic(g_t_cell);
                let (s, _c, _d) = macroscopic(g_s_cell);
                *out_val = -self.alpha_t * (t - T_REF) + self.beta_s * (s - T_REF);
            }
        }
        out
    }

    /// Diagnose thermohaline overturning circulation strength.
    pub fn overturning_strength(&self) -> f64 {
        let ny = self.ny;
        let ny2 = ny / 2;
        let psi: f64 = (0..self.nx)
            .map(|ix| {
                let (_rho, ux, _uy) = macroscopic(&self.f[ix][ny2]);
                ux
            })
            .sum();
        psi.abs()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Install a uniform flow `(ux, uy)` in the momentum populations.
    fn set_uniform_flow(f: &mut [Vec<[f64; NQ]>], ux: f64, uy: f64) {
        for row in f.iter_mut() {
            for cell in row.iter_mut() {
                *cell = feq(1.0, ux, uy);
            }
        }
    }

    /// x-centroid of `T - 1/2` over the rows `rows` of a temperature field.
    fn blob_centroid(temp: &[Vec<f64>], rows: std::ops::Range<usize>) -> f64 {
        let (mut mass, mut moment) = (0.0, 0.0);
        for (ix, column) in temp.iter().enumerate() {
            for iy in rows.clone() {
                let excess = column[iy] - 0.5;
                mass += excess;
                moment += excess * ix as f64;
            }
        }
        moment / mass
    }

    fn gaussian_blob(nx: usize, ny: usize, amplitude: f64) -> Vec<Vec<f64>> {
        let (x0, y0, sigma) = (nx as f64 * 0.5, ny as f64 * 0.5, 4.0);
        (0..nx)
            .map(|ix| {
                (0..ny)
                    .map(|iy| {
                        let r2 = (ix as f64 - x0).powi(2) + (iy as f64 - y0).powi(2);
                        0.5 + amplitude * (-r2 / (2.0 * sigma * sigma)).exp()
                    })
                    .collect()
            })
            .collect()
    }

    // ── MantelConvection ──────────────────────────────────────────────────────

    #[test]
    fn test_mantel_construction() {
        let m = MantelConvection::new(16, 16, 1.0, 1.0, 1e4, 7.0);
        assert_eq!(m.nx, 16);
    }

    #[test]
    fn test_mantel_conductive_ic() {
        let mut m = MantelConvection::new(16, 16, 1.0, 1.0, 1e4, 7.0);
        m.conductive_ic();
        assert!((m.temp[0][0] - 1.0).abs() < 1e-12, "bottom should be hot");
        assert!(
            (m.temp[0][m.ny - 1] - 0.0).abs() < 1e-12,
            "top should be cold"
        );
    }

    /// Pure conduction has Nu = 1 exactly.
    #[test]
    fn test_mantel_nusselt_conductive() {
        let mut m = MantelConvection::new(16, 16, 1.0, 1.0, 1e3, 7.0);
        m.conductive_ic();
        let nu = m.nusselt_number();
        assert!((nu - 1.0).abs() < 1e-12, "conductive Nu = {nu}");
    }

    #[test]
    fn test_mantel_step_runs() {
        let mut m = MantelConvection::new(8, 8, 1.0, 1.0, 1e3, 7.0);
        m.conductive_ic();
        for _ in 0..3 {
            m.step();
        }
        assert!(m.temp.iter().flatten().all(|&t| t.is_finite()));
    }

    #[test]
    fn mantel_transport_coefficients_follow_the_relaxation_frequencies() {
        // tau_f = 1 -> nu = 1/6; tau_t = 0.8 -> kappa = 0.1
        let m = MantelConvection::new(20, 18, 1.0, 1.0 / 0.8, 5000.0, 7.0);
        assert!((m.viscosity() - 1.0 / 6.0).abs() < 1e-15);
        assert!((m.diffusivity() - 0.1).abs() < 1e-15);
        assert!((m.prandtl_number() - 5.0 / 3.0).abs() < 1e-14);
        assert_eq!(m.layer_height(), 17.0);
        let expected = 5000.0 * (1.0 / 6.0) * 0.1 / (17.0_f64).powi(3);
        assert!((m.buoyancy_coefficient() - expected).abs() < 1e-18);
        // The Rayleigh number is recovered from the lattice buoyancy.
        let ra = m.buoyancy_coefficient() * 17.0_f64.powi(3) / (m.viscosity() * m.diffusivity());
        assert!((ra - 5000.0).abs() < 1e-8, "Ra = {ra}");
    }

    /// The temperature is advected by the *fluid* velocity: a uniform flow
    /// u = 0.05 carries a temperature blob 0.05 cells per step.
    #[test]
    fn mantel_temperature_is_advected_by_the_fluid() {
        let (nx, ny) = (64, 40);
        let mut m = MantelConvection::new(nx, ny, 1.0, 1.0, 1e-6, 1.0);
        set_uniform_flow(&mut m.f, 0.05, 0.0);
        let blob = gaussian_blob(nx, ny, 0.2);
        for ix in 0..nx {
            for iy in 0..ny {
                m.g[ix][iy] = geq_linear(blob[ix][iy], 0.0, 0.0);
            }
        }
        m.refresh_temperature();
        let start = blob_centroid(&m.temp, 12..28);
        let steps = 100;
        for _ in 0..steps {
            m.step();
        }
        let shift = blob_centroid(&m.temp, 12..28) - start;
        let expected = 0.05 * steps as f64;
        assert!(
            (shift - expected).abs() < 0.05 * expected,
            "blob moved {shift} cells, expected {expected}"
        );
    }

    /// A zero-velocity fluid does not advect the scalar, whatever the scalar's
    /// own first moment (the old code advected T with its own moment).
    #[test]
    fn mantel_scalar_first_moment_does_not_drive_the_flow() {
        let (nx, ny) = (32, 24);
        let mut m = MantelConvection::new(nx, ny, 1.0, 1.0, 1e-6, 1.0);
        // Give the scalar populations a net x-momentum while the fluid is at rest.
        for ix in 0..nx {
            for iy in 0..ny {
                m.g[ix][iy] = geq_linear(0.5, 0.05, 0.0);
            }
        }
        m.refresh_temperature();
        for _ in 0..50 {
            m.step();
        }
        // The scalar relaxes to the (zero) fluid velocity instead of keeping 0.05.
        let (_, ux, _) = macroscopic(&m.g[nx / 2][ny / 2]);
        assert!(ux.abs() < 1e-3, "scalar first moment persisted: {ux}");
    }

    /// Hot fluid (T > 1/2) is pushed up, cold fluid down, with the force
    /// F_y = rho g beta (T - 1/2): after n steps the total momentum is the
    /// accumulated impulse.
    #[test]
    fn mantel_buoyancy_impulse_matches_the_boussinesq_force() {
        let (nx, ny) = (24, 41);
        let mut m = MantelConvection::new(nx, ny, 1.0, 1.0, 5.0e4, 1.0);
        let blob = gaussian_blob(nx, ny, 0.2);
        for ix in 0..nx {
            for iy in 0..ny {
                m.g[ix][iy] = geq_linear(blob[ix][iy], 0.0, 0.0);
            }
        }
        m.refresh_temperature();
        let gbeta = m.buoyancy_coefficient();
        let excess: f64 = m.temp.iter().flatten().map(|t| t - 0.5).sum();
        let steps = 5;
        for _ in 0..steps {
            m.step();
        }
        let mut momentum_y = 0.0;
        for column in &m.f {
            for cell in column {
                for q in 0..NQ {
                    momentum_y += super::super::CY[q] * cell[q];
                }
            }
        }
        let expected = steps as f64 * gbeta * excess;
        assert!(expected > 0.0, "test setup");
        assert!(
            (momentum_y - expected).abs() < 1e-3 * expected,
            "momentum {momentum_y} vs impulse {expected}"
        );
    }

    #[test]
    fn mantel_cold_anomaly_sinks() {
        let (nx, ny) = (24, 41);
        let mut m = MantelConvection::new(nx, ny, 1.0, 1.0, 5.0e4, 1.0);
        let blob = gaussian_blob(nx, ny, -0.2);
        for ix in 0..nx {
            for iy in 0..ny {
                m.g[ix][iy] = geq_linear(blob[ix][iy], 0.0, 0.0);
            }
        }
        m.refresh_temperature();
        for _ in 0..10 {
            m.step();
        }
        let (_, _, uy) = macroscopic(&m.f[nx / 2][ny / 2]);
        assert!(uy < 0.0, "cold fluid must sink, u_y = {uy}");
    }

    /// Both walls are no-slip and hold their temperatures.
    #[test]
    fn mantel_walls_are_isothermal_and_no_slip() {
        let mut m = MantelConvection::new(16, 12, 1.2, 1.2, 2000.0, 1.0);
        m.conductive_ic();
        for _ in 0..40 {
            m.step();
        }
        for ix in 0..m.nx {
            assert!(
                (m.temp[ix][0] - 1.0).abs() < 1e-12,
                "bottom T {}",
                m.temp[ix][0]
            );
            assert!(
                m.temp[ix][m.ny - 1].abs() < 1e-12,
                "top T {}",
                m.temp[ix][m.ny - 1]
            );
            for iy in [0, m.ny - 1] {
                let (_, ux, uy) = macroscopic(&m.f[ix][iy]);
                assert!(
                    ux.abs() < 1e-12 && uy.abs() < 1e-12,
                    "wall velocity ({ux}, {uy})"
                );
            }
        }
    }

    /// Kinetic energy after `steps` for a conductive state with a critical-mode
    /// temperature perturbation, at Prandtl number 1 in a box of aspect ratio 2.
    fn convection_kinetic_energy(ra: f64, steps: &[usize]) -> Vec<f64> {
        let (nx, ny) = (34, 18);
        let mut m = MantelConvection::new(nx, ny, 1.0, 1.0, ra, 1.0);
        m.conductive_ic();
        let h = m.layer_height();
        for ix in 0..nx {
            for iy in 1..ny - 1 {
                let perturbation = 1e-3
                    * (2.0 * std::f64::consts::PI * ix as f64 / nx as f64).cos()
                    * (std::f64::consts::PI * iy as f64 / h).sin();
                let t = m.temp[ix][iy] + perturbation;
                m.g[ix][iy] = geq_linear(t, 0.0, 0.0);
            }
        }
        m.refresh_temperature();
        let mut done = 0;
        let mut out = Vec::new();
        for &target in steps {
            while done < target {
                m.step();
                done += 1;
            }
            out.push(m.kinetic_energy());
        }
        out
    }

    /// Rayleigh–Bénard onset between no-slip isothermal plates is at
    /// Ra_c = 1708: perturbations decay below it and grow above it.  This
    /// checks the Ra nu kappa / H^3 normalisation of the lattice buoyancy.
    #[test]
    fn mantel_convection_sets_in_above_the_critical_rayleigh_number() {
        let sub = convection_kinetic_energy(1000.0, &[400, 800]);
        assert!(
            sub[1] < 0.5 * sub[0],
            "Ra = 1000 < Ra_c: kinetic energy must decay, {sub:?}"
        );
        let sup = convection_kinetic_energy(4000.0, &[400, 800]);
        assert!(
            sup[1] > 4.0 * sup[0],
            "Ra = 4000 > Ra_c: kinetic energy must grow, {sup:?}"
        );
    }

    // ── ThermohalineCirculation ───────────────────────────────────────────────

    #[test]
    fn test_thermohaline_construction() {
        let thc = ThermohalineCirculation::new(16, 16, 1.0, 1.0, 0.8, 0.1, 0.08);
        assert_eq!(thc.nx, 16);
        assert_eq!(thc.gravity, 0.01);
    }

    #[test]
    fn test_thermohaline_density_anomaly() {
        let thc = ThermohalineCirculation::new(8, 8, 1.0, 1.0, 0.8, 0.1, 0.08);
        let da = thc.density_anomaly();
        assert_eq!(da.len(), 8);
        assert_eq!(da[0].len(), 8);
    }

    #[test]
    fn test_thermohaline_step_finite() {
        let mut thc = ThermohalineCirculation::new(8, 8, 1.0, 1.0, 0.8, 0.1, 0.08);
        for _ in 0..5 {
            thc.step();
        }
        let da = thc.density_anomaly();
        assert!(da.iter().flatten().all(|&v| v.is_finite()));
    }

    #[test]
    fn thermohaline_with_gravity_sets_the_field() {
        let thc = ThermohalineCirculation::new(8, 8, 1.0, 1.0, 1.0, 0.1, 0.08).with_gravity(0.03);
        assert_eq!(thc.gravity, 0.03);
    }

    /// Temperature and salinity are both advected by the fluid velocity.
    #[test]
    fn thermohaline_scalars_are_advected_by_the_fluid() {
        let (nx, ny) = (64, 40);
        let mut thc = ThermohalineCirculation::new(nx, ny, 1.0, 1.0, 1.0, 1e-9, 1e-9);
        set_uniform_flow(&mut thc.f, 0.05, 0.0);
        let temp_blob = gaussian_blob(nx, ny, 0.2);
        let sal_blob = gaussian_blob(nx, ny, 0.15);
        for ix in 0..nx {
            for iy in 0..ny {
                thc.g_t[ix][iy] = geq_linear(temp_blob[ix][iy], 0.0, 0.0);
                thc.g_s[ix][iy] = geq_linear(sal_blob[ix][iy], 0.0, 0.0);
            }
        }
        let scalar = |g: &[Vec<[f64; NQ]>]| -> Vec<Vec<f64>> {
            g.iter()
                .map(|col| col.iter().map(|cell| cell.iter().sum()).collect())
                .collect()
        };
        let (t0, s0) = (
            blob_centroid(&scalar(&thc.g_t), 12..28),
            blob_centroid(&scalar(&thc.g_s), 12..28),
        );
        let steps = 100;
        for _ in 0..steps {
            thc.step();
        }
        let expected = 0.05 * steps as f64;
        let dt = blob_centroid(&scalar(&thc.g_t), 12..28) - t0;
        let ds = blob_centroid(&scalar(&thc.g_s), 12..28) - s0;
        assert!(
            (dt - expected).abs() < 0.05 * expected,
            "temperature blob moved {dt}, expected {expected}"
        );
        assert!(
            (ds - expected).abs() < 0.05 * expected,
            "salinity blob moved {ds}, expected {expected}"
        );
    }

    /// Total y-momentum after n steps equals the accumulated buoyancy impulse
    /// n g alpha sum(T - 1/2); the field `gravity` (not a hard-coded 0.01)
    /// sets its magnitude.
    #[test]
    fn thermohaline_buoyancy_impulse_scales_with_gravity() {
        let (nx, ny) = (24, 41);
        let blob = gaussian_blob(nx, ny, 0.1);
        let mut impulses = Vec::new();
        for &gravity in &[0.01, 0.03] {
            let mut thc =
                ThermohalineCirculation::new(nx, ny, 1.0, 1.0, 1.0, 0.5, 0.2).with_gravity(gravity);
            let mut excess = 0.0;
            for ix in 0..nx {
                for iy in 0..ny {
                    thc.g_t[ix][iy] = geq_linear(blob[ix][iy], 0.0, 0.0);
                    excess += blob[ix][iy] - 0.5;
                }
            }
            let steps = 5;
            for _ in 0..steps {
                thc.step();
            }
            let mut momentum_y = 0.0;
            for column in &thc.f {
                for cell in column {
                    for q in 0..NQ {
                        momentum_y += super::super::CY[q] * cell[q];
                    }
                }
            }
            let expected = steps as f64 * gravity * 0.5 * excess;
            assert!(
                (momentum_y - expected).abs() < 1e-3 * expected,
                "g = {gravity}: momentum {momentum_y} vs impulse {expected}"
            );
            impulses.push(momentum_y);
        }
        assert!(
            (impulses[1] / impulses[0] - 3.0).abs() < 0.01,
            "impulse ratio {} should be 3",
            impulses[1] / impulses[0]
        );
    }

    /// Warm water rises, salty water sinks.
    #[test]
    fn thermohaline_warm_rises_and_salty_sinks() {
        let (nx, ny) = (24, 41);
        let blob = gaussian_blob(nx, ny, 0.2);
        let mut warm = ThermohalineCirculation::new(nx, ny, 1.0, 1.0, 1.0, 0.5, 0.5);
        let mut salty = ThermohalineCirculation::new(nx, ny, 1.0, 1.0, 1.0, 0.5, 0.5);
        for ix in 0..nx {
            for iy in 0..ny {
                warm.g_t[ix][iy] = geq_linear(blob[ix][iy], 0.0, 0.0);
                salty.g_s[ix][iy] = geq_linear(blob[ix][iy], 0.0, 0.0);
            }
        }
        for _ in 0..10 {
            warm.step();
            salty.step();
        }
        let (_, _, uy_warm) = macroscopic(&warm.f[nx / 2][ny / 2]);
        let (_, _, uy_salty) = macroscopic(&salty.f[nx / 2][ny / 2]);
        assert!(uy_warm > 0.0, "warm anomaly must rise, u_y = {uy_warm}");
        assert!(uy_salty < 0.0, "salty anomaly must sink, u_y = {uy_salty}");
    }

    /// The basin is closed: total mass, heat and salt are conserved and a
    /// vertical through-flow is stopped at the walls (periodic wrapping would
    /// keep it forever).
    #[test]
    fn thermohaline_basin_is_closed() {
        let (nx, ny) = (16, 40);
        let mut thc = ThermohalineCirculation::new(nx, ny, 1.0, 1.0, 0.8, 0.1, 0.08);
        set_uniform_flow(&mut thc.f, 0.0, 0.05);
        let blob = gaussian_blob(nx, ny, 0.1);
        for ix in 0..nx {
            for iy in 0..ny {
                thc.g_t[ix][iy] = geq_linear(blob[ix][iy], 0.0, 0.0);
                thc.g_s[ix][iy] = geq_linear(1.0 - blob[ix][iy], 0.0, 0.0);
            }
        }
        let total = |g: &[Vec<[f64; NQ]>]| -> f64 {
            g.iter()
                .flatten()
                .map(|cell| cell.iter().sum::<f64>())
                .sum()
        };
        let (m0, t0, s0) = (total(&thc.f), total(&thc.g_t), total(&thc.g_s));
        for _ in 0..20 {
            thc.step();
        }
        let (m1, t1, s1) = (total(&thc.f), total(&thc.g_t), total(&thc.g_s));
        assert!((m1 - m0).abs() < 1e-10 * m0, "mass {m0} -> {m1}");
        assert!((t1 - t0).abs() < 1e-10 * t0, "heat {t0} -> {t1}");
        assert!((s1 - s0).abs() < 1e-10 * s0, "salt {s0} -> {s1}");
        let (_, _, uy_wall) = macroscopic(&thc.f[nx / 2][0]);
        let (_, _, uy_mid) = macroscopic(&thc.f[nx / 2][ny / 2]);
        assert!(
            uy_wall.abs() < 0.5 * 0.05,
            "the wall row must have lost its through-flow, u_y = {uy_wall}"
        );
        assert!(
            (uy_mid - 0.05).abs() < 0.1 * 0.05,
            "the mid-depth flow is not yet affected, u_y = {uy_mid}"
        );
    }
}
