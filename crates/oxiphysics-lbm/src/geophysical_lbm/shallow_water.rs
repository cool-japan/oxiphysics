// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Lattice Boltzmann shallow-water equations (Zhou 2004).

use super::{CX, CY, NQ, guo_force, macroscopic, stream_periodic};

/// Shallow-water LBM following Zhou (2004), *Lattice Boltzmann Methods for
/// Shallow Water Flows*.
///
/// The populations carry the water depth `h = Σ f_q` and the discharge
/// `h u = Σ c_q f_q` on the D2Q9 lattice with unit lattice speed
/// (`Δx = Δt = 1`).  The equilibrium of [`ShallowWaterLbm`] has the second
/// moment `½ g h² δ_ab + h u_a u_b`, i.e. the shallow-water momentum flux, so
/// the long-wave limit recovers
///
/// ```text
/// ∂t h + ∇·(h u) = 0
/// ∂t (h u) + ∇·(h u ⊗ u) = −g h ∇(h + b) + ν ∇²(h u)
/// ```
///
/// with the bed elevation `b` (field `bathymetry`), the lattice viscosity
/// `ν = (2/ω − 1)/6` and long-wave gravity speed `c = √(g h)`.  The bed
/// slope force `−g h ∇b` (central differences, periodic) is applied with the
/// Guo scheme using the force-corrected velocity, so total mass is conserved
/// to rounding for any bathymetry.  For stability `g h` has to stay well below
/// `1/3` (wave speed far below the lattice sound speed).  The domain is
/// periodic in both directions.
pub struct ShallowWaterLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// BGK relaxation frequency.
    pub omega: f64,
    /// Gravitational acceleration in lattice units.
    pub gravity: f64,
    /// Mean water depth.
    pub h_mean: f64,
    /// Distribution functions representing the height field.
    pub f: Vec<Vec<[f64; NQ]>>,
    /// Bed elevation `b` above the reference level; the free surface is
    /// `h + b`.
    pub bathymetry: Vec<Vec<f64>>,
}

impl ShallowWaterLbm {
    /// Create a new shallow-water LBM domain at rest with uniform depth.
    ///
    /// # Arguments
    /// * `nx`, `ny` — domain size.
    /// * `omega` — relaxation frequency.
    /// * `gravity` — gravitational acceleration (lattice units).
    /// * `h_mean` — mean water depth.
    pub fn new(nx: usize, ny: usize, omega: f64, gravity: f64, h_mean: f64) -> Self {
        let f = vec![vec![Self::feq_sw(h_mean, 0.0, 0.0, gravity); ny]; nx];
        let bathymetry = vec![vec![0.0f64; ny]; nx];
        Self {
            nx,
            ny,
            omega,
            gravity,
            h_mean,
            f,
            bathymetry,
        }
    }

    /// Set a Gaussian dam-break initial condition.
    ///
    /// Elevates the water depth by `amplitude` at the grid centre, at rest.
    pub fn dam_break_ic(&mut self, amplitude: f64) {
        let cx = self.nx as f64 * 0.5;
        let cy = self.ny as f64 * 0.5;
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let r2 = (ix as f64 - cx).powi(2) + (iy as f64 - cy).powi(2);
                let h = self.h_mean + amplitude * (-r2 / 50.0).exp();
                self.f[ix][iy] = Self::feq_sw(h, 0.0, 0.0, self.gravity);
            }
        }
    }

    /// Shallow-water equilibrium of Zhou (2004) for lattice speed 1.
    ///
    /// ```text
    /// f_0     = h − 5 g h²/6 − 2 h u²/3
    /// f_{1-4} = g h²/6  + h (c·u)/3  + h (c·u)²/2 − h u²/6
    /// f_{5-8} = g h²/24 + h (c·u)/12 + h (c·u)²/8 − h u²/24
    /// ```
    ///
    /// Its moments are `Σ f = h`, `Σ c f = h u` and
    /// `Σ c_a c_b f = ½ g h² δ_ab + h u_a u_b`.
    fn feq_sw(h: f64, ux: f64, uy: f64, gravity: f64) -> [f64; NQ] {
        let u2 = ux * ux + uy * uy;
        let gh2 = gravity * h * h;
        let mut f = [0.0f64; NQ];
        f[0] = h - 5.0 * gh2 / 6.0 - 2.0 * h * u2 / 3.0;
        for q in 1..NQ {
            let cu = CX[q] * ux + CY[q] * uy;
            f[q] = if q <= 4 {
                gh2 / 6.0 + h * cu / 3.0 + h * cu * cu / 2.0 - h * u2 / 6.0
            } else {
                gh2 / 24.0 + h * cu / 12.0 + h * cu * cu / 8.0 - h * u2 / 24.0
            };
        }
        f
    }

    /// Long-wave gravity wave speed `√(g h_mean)` (lattice units).
    pub fn wave_speed(&self) -> f64 {
        (self.gravity * self.h_mean).sqrt()
    }

    /// Advance one timestep: BGK collision with the bed-slope force
    /// `F = −g h ∇b`, then periodic streaming.
    pub fn step(&mut self) {
        let omega = self.omega;
        let g = self.gravity;
        let (nx, ny) = (self.nx, self.ny);
        for ix in 0..nx {
            let ixp = (ix + 1) % nx;
            let ixm = (ix + nx - 1) % nx;
            for iy in 0..ny {
                let iyp = (iy + 1) % ny;
                let iym = (iy + ny - 1) % ny;
                let (h, ux0, uy0) = macroscopic(&self.f[ix][iy]);
                // Bed-slope force density F = −g h ∇b (central differences).
                let bx = 0.5 * (self.bathymetry[ixp][iy] - self.bathymetry[ixm][iy]);
                let by = 0.5 * (self.bathymetry[ix][iyp] - self.bathymetry[ix][iym]);
                let fx = -g * h * bx;
                let fy = -g * h * by;
                // Force-corrected velocity u = (Σ c f + F/2) / h.
                let (ux, uy) = if h > 1e-12 {
                    (ux0 + 0.5 * fx / h, uy0 + 0.5 * fy / h)
                } else {
                    (0.0, 0.0)
                };
                let feq_arr = Self::feq_sw(h, ux, uy, g);
                let cell = &mut self.f[ix][iy];
                for (f_q, feq_q) in cell.iter_mut().zip(feq_arr.iter()) {
                    *f_q += omega * (*feq_q - *f_q);
                }
                guo_force(cell, fx, fy, ux, uy, omega);
            }
        }
        self.f = stream_periodic(&self.f, nx, ny);
    }

    /// Return the water depth `h` as a flat vector.
    pub fn height_field(&self) -> Vec<f64> {
        let mut h = Vec::with_capacity(self.nx * self.ny);
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (rho, _ux, _uy) = macroscopic(&self.f[ix][iy]);
                h.push(rho);
            }
        }
        h
    }

    /// Return the free-surface elevation `h + b` as a flat vector.
    pub fn surface_elevation(&self) -> Vec<f64> {
        let mut eta = Vec::with_capacity(self.nx * self.ny);
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (h, _ux, _uy) = macroscopic(&self.f[ix][iy]);
                eta.push(h + self.bathymetry[ix][iy]);
            }
        }
        eta
    }

    /// Total mass `Σ h` over the grid (conserved to rounding).
    pub fn total_mass(&self) -> f64 {
        self.f
            .iter()
            .flatten()
            .map(|cell| cell.iter().sum::<f64>())
            .sum()
    }

    /// Compute total energy: kinetic `½ h |u|²` plus potential
    /// `½ g h² + g h b` (free-surface and bed contributions).
    pub fn total_energy(&self) -> f64 {
        let g = self.gravity;
        let mut e = 0.0;
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (h, ux, uy) = macroscopic(&self.f[ix][iy]);
                e += 0.5 * h * (ux * ux + uy * uy)
                    + 0.5 * g * h * h
                    + g * h * self.bathymetry[ix][iy];
            }
        }
        e
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Moments of a population vector: (h, j_x, j_y, Π_xx, Π_xy, Π_yy).
    fn moments(f: &[f64; NQ]) -> (f64, f64, f64, f64, f64, f64) {
        let mut m = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        for q in 0..NQ {
            m.0 += f[q];
            m.1 += CX[q] * f[q];
            m.2 += CY[q] * f[q];
            m.3 += CX[q] * CX[q] * f[q];
            m.4 += CX[q] * CY[q] * f[q];
            m.5 += CY[q] * CY[q] * f[q];
        }
        m
    }

    #[test]
    fn equilibrium_has_shallow_water_moments() {
        let (h, ux, uy, g) = (1.3, 0.04, -0.03, 0.05);
        let f = ShallowWaterLbm::feq_sw(h, ux, uy, g);
        let (m0, jx, jy, pxx, pxy, pyy) = moments(&f);
        assert!((m0 - h).abs() < 1e-14, "depth {m0}");
        assert!((jx - h * ux).abs() < 1e-14, "j_x {jx}");
        assert!((jy - h * uy).abs() < 1e-14, "j_y {jy}");
        assert!(
            (pxx - (0.5 * g * h * h + h * ux * ux)).abs() < 1e-14,
            "Pxx {pxx}"
        );
        assert!((pxy - h * ux * uy).abs() < 1e-14, "Pxy {pxy}");
        assert!(
            (pyy - (0.5 * g * h * h + h * uy * uy)).abs() < 1e-14,
            "Pyy {pyy}"
        );
    }

    #[test]
    fn equilibrium_depends_on_gravity() {
        let rest_low = ShallowWaterLbm::feq_sw(1.0, 0.0, 0.0, 0.01);
        let rest_high = ShallowWaterLbm::feq_sw(1.0, 0.0, 0.0, 0.04);
        assert!(rest_low[0] > rest_high[0], "rest population shrinks with g");
        assert!(rest_high[1] > rest_low[1], "axis populations grow with g");
        assert!(
            rest_high[5] > rest_low[5],
            "diagonal populations grow with g"
        );
    }

    /// Place a Gaussian hump `h = h_mean + A exp(-(x - x0)^2 / (2 s^2))`,
    /// uniform along y, at rest.
    fn gaussian_hump(sw: &mut ShallowWaterLbm, x0: f64, sigma: f64, amplitude: f64) {
        for ix in 0..sw.nx {
            let x = ix as f64 - x0;
            let h = sw.h_mean + amplitude * (-x * x / (2.0 * sigma * sigma)).exp();
            for iy in 0..sw.ny {
                sw.f[ix][iy] = ShallowWaterLbm::feq_sw(h, 0.0, 0.0, sw.gravity);
            }
        }
    }

    /// Distance travelled by the right-going half of a split hump.
    fn right_pulse_travel(gravity: f64, steps: usize) -> f64 {
        let (nx, ny) = (200, 4);
        let x0 = 100.0;
        let mut sw = ShallowWaterLbm::new(nx, ny, 1.0, gravity, 1.0);
        gaussian_hump(&mut sw, x0, 6.0, 0.005);
        for _ in 0..steps {
            sw.step();
        }
        let h = sw.height_field();
        let (mut mass, mut moment) = (0.0, 0.0);
        for ix in (x0 as usize)..nx {
            let dh = h[ix * ny + 1] - sw.h_mean;
            mass += dh;
            moment += dh * ix as f64;
        }
        moment / mass - x0
    }

    #[test]
    fn long_waves_travel_at_sqrt_g_h() {
        let steps = 150;
        for &g in &[0.02, 0.04] {
            let expected = (g * 1.0_f64).sqrt() * steps as f64;
            let travel = right_pulse_travel(g, steps);
            assert!(
                (travel - expected).abs() < 0.03 * expected,
                "g = {g}: pulse travelled {travel}, expected sqrt(g h) t = {expected}"
            );
        }
    }

    #[test]
    fn wave_speed_scales_with_the_square_root_of_gravity() {
        let slow = right_pulse_travel(0.01, 150);
        let fast = right_pulse_travel(0.04, 150);
        assert!(
            (fast / slow - 2.0).abs() < 0.06,
            "quadrupling g must double the wave speed, ratio {}",
            fast / slow
        );
    }

    #[test]
    fn wave_speed_accessor_matches_definition() {
        let sw = ShallowWaterLbm::new(8, 8, 1.0, 0.04, 2.0);
        assert!((sw.wave_speed() - (0.08_f64).sqrt()).abs() < 1e-15);
    }

    fn bump_and_hump(sw: &mut ShallowWaterLbm) {
        let (cx, cy) = (sw.nx as f64 * 0.5, sw.ny as f64 * 0.5);
        for ix in 0..sw.nx {
            for iy in 0..sw.ny {
                let r2 = (ix as f64 - cx).powi(2) + (iy as f64 - cy).powi(2);
                sw.bathymetry[ix][iy] = 0.3 * (-r2 / 32.0).exp();
                let h = sw.h_mean + 0.1 * (-r2 / 50.0).exp();
                sw.f[ix][iy] = ShallowWaterLbm::feq_sw(h, 0.0, 0.0, sw.gravity);
            }
        }
    }

    #[test]
    fn mass_is_conserved_with_nonzero_bathymetry() {
        let mut sw = ShallowWaterLbm::new(32, 32, 1.2, 0.02, 1.0);
        bump_and_hump(&mut sw);
        let mass0 = sw.total_mass();
        for _ in 0..300 {
            sw.step();
        }
        let mass1 = sw.total_mass();
        assert!(
            (mass1 - mass0).abs() < 1e-10 * mass0,
            "mass changed from {mass0} to {mass1}"
        );
        let h = sw.height_field();
        assert!(h.iter().all(|&v| v.is_finite() && v > 0.0));
        // The bathymetry must actually act: compare with a flat-bed run.
        let mut flat = ShallowWaterLbm::new(32, 32, 1.2, 0.02, 1.0);
        bump_and_hump(&mut flat);
        for row in flat.bathymetry.iter_mut() {
            row.fill(0.0);
        }
        for _ in 0..300 {
            flat.step();
        }
        let hf = flat.height_field();
        let max_diff = h
            .iter()
            .zip(hf.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_diff > 1e-3,
            "bathymetry has no effect (diff {max_diff})"
        );
    }

    #[test]
    fn bed_slope_force_pushes_water_downhill() {
        let (nx, ny) = (48, 8);
        let mut sw = ShallowWaterLbm::new(nx, ny, 1.0, 0.02, 1.0);
        let x0 = 24.0;
        let mut expected_right = 0.0;
        for ix in 0..nx {
            let x = ix as f64 - x0;
            for iy in 0..ny {
                sw.bathymetry[ix][iy] = 0.2 * (-x * x / (2.0 * 16.0)).exp();
            }
        }
        for ix in 0..nx {
            let ixp = (ix + 1) % nx;
            let ixm = (ix + nx - 1) % nx;
            let slope = 0.5 * (sw.bathymetry[ixp][0] - sw.bathymetry[ixm][0]);
            if ix as f64 > x0 {
                expected_right += -sw.gravity * sw.h_mean * slope * ny as f64;
            }
        }
        assert!(expected_right > 0.0, "test setup");
        sw.step();
        let mut momentum_right = 0.0;
        let mut momentum_left = 0.0;
        for ix in 0..nx {
            for iy in 0..ny {
                let (h, ux, _) = macroscopic(&sw.f[ix][iy]);
                if ix as f64 > x0 {
                    momentum_right += h * ux;
                } else {
                    momentum_left += h * ux;
                }
            }
        }
        assert!(
            momentum_right > 0.0 && momentum_left < 0.0,
            "water must flow away from the bump: right {momentum_right}, left {momentum_left}"
        );
        assert!(
            (momentum_right - expected_right).abs() < 0.05 * expected_right,
            "momentum {momentum_right} vs force impulse {expected_right}"
        );
    }

    #[test]
    fn energy_includes_the_bed_potential() {
        let mut sw = ShallowWaterLbm::new(8, 8, 1.0, 0.02, 1.0);
        let flat = sw.total_energy();
        for row in sw.bathymetry.iter_mut() {
            row.fill(0.5);
        }
        let raised = sw.total_energy();
        assert!((raised - flat - 0.02 * 1.0 * 0.5 * 64.0).abs() < 1e-12);
    }

    #[test]
    fn surface_elevation_adds_the_bed() {
        let mut sw = ShallowWaterLbm::new(4, 4, 1.0, 0.02, 1.0);
        sw.bathymetry[1][2] = 0.25;
        let eta = sw.surface_elevation();
        assert!((eta[1 * 4 + 2] - 1.25).abs() < 1e-14);
        assert!((eta[0] - 1.0).abs() < 1e-14);
    }
}
