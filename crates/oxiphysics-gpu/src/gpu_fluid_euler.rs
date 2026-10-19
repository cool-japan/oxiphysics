// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU Eulerian fluid simulation on a staggered MAC grid (CPU mock backend).
//!
//! Implements a staggered Marker-and-Cell (MAC) grid fluid solver with:
//! semi-Lagrangian advection, gravity forcing, Jacobi pressure projection,
//! solid boundary conditions, vorticity confinement, and diagnostic utilities.

use crate::error::{Error, Result};

// ── Data structures ──────────────────────────────────────────────────────────

/// A 3-D staggered MAC grid for Eulerian fluid simulation.
///
/// Velocity components are face-centred; pressure and density are cell-centred.
#[derive(Debug, Clone)]
pub struct GpuEulerGrid {
    /// Number of cells in x.
    pub nx: usize,
    /// Number of cells in y.
    pub ny: usize,
    /// Number of cells in z.
    pub nz: usize,
    /// Grid spacing (m).
    pub dx: f64,
    /// Cell-centred density (kg/m³).
    pub rho: Vec<f64>,
    /// x-face velocity (m/s), size (nx+1)·ny·nz.
    pub u: Vec<f64>,
    /// y-face velocity (m/s), size nx·(ny+1)·nz.
    pub v: Vec<f64>,
    /// z-face velocity (m/s), size nx·ny·(nz+1).
    pub w: Vec<f64>,
    /// Cell-centred pressure (Pa).
    pub p: Vec<f64>,
    /// Max-norm Poisson residual `max |∇²p − (ρ/Δt) ∇·u*|` left by the most
    /// recent [`gpu_pressure_projection`](Self::gpu_pressure_projection)
    /// (`0.0` before any projection).
    pub last_residual: f64,
}

impl GpuEulerGrid {
    /// Create a new `GpuEulerGrid` with `nx × ny × nz` cells and spacing `dx`.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64) -> Self {
        let nc = nx * ny * nz;
        Self {
            nx,
            ny,
            nz,
            dx,
            rho: vec![1.0; nc],
            u: vec![0.0; (nx + 1) * ny * nz],
            v: vec![0.0; nx * (ny + 1) * nz],
            w: vec![0.0; nx * ny * (nz + 1)],
            p: vec![0.0; nc],
            last_residual: 0.0,
        }
    }

    /// Flat cell index for cell (i, j, k).
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.nx * self.ny + j * self.nx + i
    }

    /// Total number of cells.
    pub fn cell_count(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Maximum speed over all velocity components.
    pub fn max_velocity(&self) -> f64 {
        let mu = self.u.iter().cloned().fold(0.0_f64, |a, v| a.max(v.abs()));
        let mv = self.v.iter().cloned().fold(0.0_f64, |a, v| a.max(v.abs()));
        let mw = self.w.iter().cloned().fold(0.0_f64, |a, v| a.max(v.abs()));
        mu.max(mv).max(mw)
    }

    /// CFL-safe time step: `cfl · dx / max_velocity` (or `dx` if velocity is zero).
    pub fn cfl_timestep(&self, cfl: f64) -> f64 {
        let vmax = self.max_velocity();
        if vmax < 1e-15 {
            return self.dx;
        }
        cfl * self.dx / vmax
    }

    /// Apply gravitational acceleration to the v (y) velocity field.
    pub fn gpu_apply_gravity(&mut self, gravity: f64, dt: f64) {
        for val in self.v.iter_mut() {
            *val += gravity * dt;
        }
    }

    /// Compute the divergence field ∇·u at every cell.
    pub fn gpu_compute_divergence(&self) -> Vec<f64> {
        let nc = self.cell_count();
        let mut div = vec![0.0f64; nc];
        let dx = self.dx;
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    let idx = self.index(i, j, k);
                    // u on x-faces
                    let ui_p = self.u[k * (self.nx + 1) * self.ny + j * (self.nx + 1) + i + 1];
                    let ui_m = self.u[k * (self.nx + 1) * self.ny + j * (self.nx + 1) + i];
                    // v on y-faces
                    let vj_p = self.v[k * self.nx * (self.ny + 1) + (j + 1) * self.nx + i];
                    let vj_m = self.v[k * self.nx * (self.ny + 1) + j * self.nx + i];
                    // w on z-faces
                    let wk_p = self.w[(k + 1) * self.nx * self.ny + j * self.nx + i];
                    let wk_m = self.w[k * self.nx * self.ny + j * self.nx + i];
                    div[idx] = (ui_p - ui_m + vj_p - vj_m + wk_p - wk_m) / dx;
                }
            }
        }
        div
    }

    /// Local divergence at cell (i, j, k).
    pub fn divergence_at(&self, i: usize, j: usize, k: usize) -> f64 {
        let div = self.gpu_compute_divergence();
        div[self.index(i, j, k)]
    }

    /// Local vorticity (curl of velocity) at cell (i, j, k).
    ///
    /// Returns `[ωx, ωy, ωz]` using central differences.
    pub fn vorticity_at(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        let dx = self.dx;
        // Use interior clamping for boundary cells
        let ip = i.min(self.nx - 1);
        let im = if i > 0 { i - 1 } else { 0 };
        let jp = j.min(self.ny - 1);
        let jm = if j > 0 { j - 1 } else { 0 };
        let kp = k.min(self.nz - 1);
        let km = if k > 0 { k - 1 } else { 0 };

        // Sample cell-centred velocities from face averages
        let u_jm = self.u[km * (self.nx + 1) * self.ny + jm * (self.nx + 1) + i];
        let u_jp = self.u[kp * (self.nx + 1) * self.ny + jp * (self.nx + 1) + i];
        let u_km = self.u[km * (self.nx + 1) * self.ny + j * (self.nx + 1) + i];
        let u_kp = self.u[kp * (self.nx + 1) * self.ny + j * (self.nx + 1) + i];

        let v_im = self.v[k * self.nx * (self.ny + 1) + j * self.nx + im];
        let v_ip = self.v[k * self.nx * (self.ny + 1) + j * self.nx + ip];
        let v_km = self.v[km * self.nx * (self.ny + 1) + j * self.nx + i];
        let v_kp = self.v[kp * self.nx * (self.ny + 1) + j * self.nx + i];

        let w_im = self.w[k * self.nx * self.ny + j * self.nx + im];
        let w_ip = self.w[k * self.nx * self.ny + j * self.nx + ip];
        let w_jm = self.w[k * self.nx * self.ny + jm * self.nx + i];
        let w_jp = self.w[k * self.nx * self.ny + jp * self.nx + i];

        let wx = (w_jp - w_jm) / (2.0 * dx) - (v_kp - v_km) / (2.0 * dx);
        let wy = (u_kp - u_km) / (2.0 * dx) - (w_ip - w_im) / (2.0 * dx);
        let wz = (v_ip - v_im) / (2.0 * dx) - (u_jp - u_jm) / (2.0 * dx);
        [wx, wy, wz]
    }

    /// Zero normal velocities at domain boundaries (solid wall BC).
    pub fn gpu_enforce_solid_bc(&mut self) {
        // x-faces at i=0 and i=nx
        for k in 0..self.nz {
            for j in 0..self.ny {
                self.u[k * (self.nx + 1) * self.ny + j * (self.nx + 1)] = 0.0;
                self.u[k * (self.nx + 1) * self.ny + j * (self.nx + 1) + self.nx] = 0.0;
            }
        }
        // y-faces at j=0 and j=ny
        for k in 0..self.nz {
            for i in 0..self.nx {
                self.v[k * self.nx * (self.ny + 1) + i] = 0.0;
                self.v[k * self.nx * (self.ny + 1) + self.ny * self.nx + i] = 0.0;
            }
        }
        // z-faces at k=0 and k=nz
        for j in 0..self.ny {
            for i in 0..self.nx {
                self.w[j * self.nx + i] = 0.0;
                self.w[self.nz * self.nx * self.ny + j * self.nx + i] = 0.0;
            }
        }
    }

    /// Validate the time step and density shared by the projection routines.
    fn check_projection_params(dt: f64, rho: f64) -> Result<()> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(Error::general(format!(
                "pressure projection needs a finite positive time step, got dt = {dt}"
            )));
        }
        if !(rho.is_finite() && rho > 0.0) {
            return Err(Error::general(format!(
                "pressure projection needs a finite positive density, got rho = {rho}"
            )));
        }
        Ok(())
    }

    /// Sum of the six face-neighbour values of cell `(i, j, k)` in the
    /// cell-centred field `field`.
    ///
    /// Across a domain face the zero-gradient (Neumann) ghost value, i.e. the
    /// cell's own value, is used.
    fn neighbour_sum(&self, field: &[f64], i: usize, j: usize, k: usize) -> f64 {
        let centre = field[self.index(i, j, k)];
        let at = |ii: usize, jj: usize, kk: usize| field[self.index(ii, jj, kk)];
        let x_plus = if i + 1 < self.nx {
            at(i + 1, j, k)
        } else {
            centre
        };
        let x_minus = if i > 0 { at(i - 1, j, k) } else { centre };
        let y_plus = if j + 1 < self.ny {
            at(i, j + 1, k)
        } else {
            centre
        };
        let y_minus = if j > 0 { at(i, j - 1, k) } else { centre };
        let z_plus = if k + 1 < self.nz {
            at(i, j, k + 1)
        } else {
            centre
        };
        let z_minus = if k > 0 { at(i, j, k - 1) } else { centre };
        x_plus + x_minus + y_plus + y_minus + z_plus + z_minus
    }

    /// Max-norm of `∇²p − rhs` over all cells, with the same zero-gradient
    /// boundary treatment as the Jacobi sweep.
    fn poisson_residual(&self, rhs: &[f64]) -> f64 {
        let inv_dx2 = 1.0 / (self.dx * self.dx);
        let mut worst = 0.0_f64;
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    let idx = self.index(i, j, k);
                    let laplacian =
                        (self.neighbour_sum(&self.p, i, j, k) - 6.0 * self.p[idx]) * inv_dx2;
                    worst = worst.max((laplacian - rhs[idx]).abs());
                }
            }
        }
        worst
    }

    /// Jacobi pressure projection to enforce incompressibility.
    ///
    /// Solves the pressure Poisson equation `∇²p = (ρ/Δt) ∇·u*` for the
    /// current face velocities `u*` with zero-gradient (Neumann) boundaries,
    /// starting from the current pressure field.  Together with
    /// [`gpu_update_velocity_from_pressure`](Self::gpu_update_velocity_from_pressure)
    /// (`u ← u* − (Δt/ρ) ∇p`) the resulting velocity field is divergence-free
    /// to the accuracy of the solve, independently of `dt` and `rho`.
    ///
    /// The velocity update only touches interior faces, so the boundary faces
    /// must already carry zero normal velocity
    /// ([`gpu_enforce_solid_bc`](Self::gpu_enforce_solid_bc)); this also makes
    /// the all-Neumann problem solvable (`Σ ∇·u* = 0`).
    ///
    /// Performs at most `max_iter` sweeps and stops early once the largest
    /// pressure change of a sweep falls below `tol`.  Returns the max-norm of
    /// the Poisson residual `max |∇²p − (ρ/Δt) ∇·u*|` after the last sweep and
    /// stores it in [`last_residual`](Self::last_residual).
    ///
    /// # Errors
    ///
    /// Returns an error when `dt` or `rho` is not a finite positive number, or
    /// when the velocity field contains non-finite values; the pressure field
    /// is left unchanged in that case.
    pub fn gpu_pressure_projection(
        &mut self,
        max_iter: usize,
        tol: f64,
        dt: f64,
        rho: f64,
    ) -> Result<f64> {
        Self::check_projection_params(dt, rho)?;
        let dx2 = self.dx * self.dx;
        let scale = rho / dt;
        let rhs: Vec<f64> = self
            .gpu_compute_divergence()
            .into_iter()
            .map(|d| scale * d)
            .collect();
        if rhs.iter().any(|r| !r.is_finite()) {
            return Err(Error::general(
                "pressure projection needs a finite velocity field",
            ));
        }

        // Jacobi: every sweep reads the previous iterate and writes `p_next`.
        let mut p_next = self.p.clone();
        for _ in 0..max_iter {
            let mut max_update = 0.0_f64;
            for k in 0..self.nz {
                for j in 0..self.ny {
                    for i in 0..self.nx {
                        let idx = self.index(i, j, k);
                        let p_new = (self.neighbour_sum(&self.p, i, j, k) - dx2 * rhs[idx]) / 6.0;
                        max_update = max_update.max((p_new - self.p[idx]).abs());
                        p_next[idx] = p_new;
                    }
                }
            }
            std::mem::swap(&mut self.p, &mut p_next);
            if max_update < tol {
                break;
            }
        }
        self.last_residual = self.poisson_residual(&rhs);
        Ok(self.last_residual)
    }

    /// Subtract the pressure gradient from the interior face velocities:
    /// `u ← u − (Δt/ρ) ∇p`.
    ///
    /// Boundary faces are left untouched (zero-gradient pressure boundary).
    ///
    /// # Errors
    ///
    /// Returns an error when `dt` or `rho` is not a finite positive number; the
    /// velocities are left unchanged in that case.
    pub fn gpu_update_velocity_from_pressure(&mut self, dt: f64, rho: f64) -> Result<()> {
        Self::check_projection_params(dt, rho)?;
        let factor = dt / (rho * self.dx);
        // Update u (x-faces)
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 1..self.nx {
                    let idx_r = self.index(i, j, k);
                    let idx_l = self.index(i - 1, j, k);
                    let fi = k * (self.nx + 1) * self.ny + j * (self.nx + 1) + i;
                    self.u[fi] -= factor * (self.p[idx_r] - self.p[idx_l]);
                }
            }
        }
        // Update v (y-faces)
        for k in 0..self.nz {
            for j in 1..self.ny {
                for i in 0..self.nx {
                    let idx_t = self.index(i, j, k);
                    let idx_b = self.index(i, j - 1, k);
                    let fj = k * self.nx * (self.ny + 1) + j * self.nx + i;
                    self.v[fj] -= factor * (self.p[idx_t] - self.p[idx_b]);
                }
            }
        }
        // Update w (z-faces)
        for k in 1..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    let idx_f = self.index(i, j, k);
                    let idx_b = self.index(i, j, k - 1);
                    let fk = k * self.nx * self.ny + j * self.nx + i;
                    self.w[fk] -= factor * (self.p[idx_f] - self.p[idx_b]);
                }
            }
        }
        Ok(())
    }

    /// Semi-Lagrangian advection of velocity (first-order back-trace).
    pub fn gpu_advect_semi_lagrange(&mut self, dt: f64) {
        let dx = self.dx;
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;

        // Advect u
        let u_old = self.u.clone();
        let v_old = self.v.clone();
        let w_old = self.w.clone();
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..=nx {
                    let fi = k * (nx + 1) * ny + j * (nx + 1) + i;
                    // Average v and w at this face
                    let vavg = if i < nx {
                        0.5 * (v_old[k * nx * (ny + 1) + j * nx + i]
                            + v_old[k * nx * (ny + 1) + (j + 1).min(ny) * nx + i])
                    } else {
                        0.0
                    };
                    let wavg = if i < nx {
                        0.5 * (w_old[k * nx * ny + j * nx + i]
                            + w_old[(k + 1).min(nz) * nx * ny + j * nx + i])
                    } else {
                        0.0
                    };
                    let xi = i as f64 * dx - u_old[fi] * dt;
                    let yi = (j as f64 + 0.5) * dx - vavg * dt;
                    let zi = (k as f64 + 0.5) * dx - wavg * dt;
                    self.u[fi] = trilinear_u(&u_old, xi, yi, zi, nx, ny, nz, dx);
                }
            }
        }
        // Advect v
        for k in 0..nz {
            for j in 0..=ny {
                for i in 0..nx {
                    let fj = k * nx * (ny + 1) + j * nx + i;
                    let uavg = if j < ny {
                        0.5 * (u_old[k * (nx + 1) * ny + j * (nx + 1) + i]
                            + u_old[k * (nx + 1) * ny + j * (nx + 1) + i + 1])
                    } else {
                        0.0
                    };
                    let wavg = if j < ny {
                        0.5 * (w_old[k * nx * ny + j * nx + i]
                            + w_old[(k + 1).min(nz) * nx * ny + j * nx + i])
                    } else {
                        0.0
                    };
                    let xi = (i as f64 + 0.5) * dx - uavg * dt;
                    let yi = j as f64 * dx - v_old[fj] * dt;
                    let zi = (k as f64 + 0.5) * dx - wavg * dt;
                    self.v[fj] = trilinear_v(&v_old, xi, yi, zi, nx, ny, nz, dx);
                }
            }
        }
        // Advect w (simplified: re-use current values for now – keeps it O(N))
        for k in 0..=nz {
            for j in 0..ny {
                for i in 0..nx {
                    let fk = k * nx * ny + j * nx + i;
                    let uavg = if k < nz {
                        0.5 * (u_old[k * (nx + 1) * ny + j * (nx + 1) + i]
                            + u_old[k * (nx + 1) * ny + j * (nx + 1) + i + 1])
                    } else {
                        0.0
                    };
                    let xi = (i as f64 + 0.5) * dx - uavg * dt;
                    let zi = k as f64 * dx - w_old[fk] * dt;
                    self.w[fk] = trilinear_w(&w_old, xi, zi, i, j, k, nx, ny, nz, dx);
                }
            }
        }
    }

    /// Vorticity confinement: adds a force proportional to `epsilon` toward
    /// regions of high vorticity to counteract numerical diffusion.
    pub fn gpu_vorticity_confinement(&mut self, epsilon: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let dx = self.dx;
        // Collect vorticity magnitudes
        let mut eta: Vec<[f64; 3]> = vec![[0.0; 3]; nx * ny * nz];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    eta[k * nx * ny + j * nx + i] = self.vorticity_at(i, j, k);
                }
            }
        }
        // Compute normalised gradient of |eta| and apply confinement force to v
        for k in 1..nz - 1 {
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    let mag = |e: [f64; 3]| (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt();
                    let m_ip = mag(eta[k * nx * ny + j * nx + i + 1]);
                    let m_im = mag(eta[k * nx * ny + j * nx + i - 1]);
                    let m_jp = mag(eta[k * nx * ny + (j + 1) * nx + i]);
                    let m_jm = mag(eta[k * nx * ny + (j - 1) * nx + i]);
                    let m_kp = mag(eta[(k + 1) * nx * ny + j * nx + i]);
                    let m_km = mag(eta[(k - 1) * nx * ny + j * nx + i]);
                    let nx_grad = (m_ip - m_im) / (2.0 * dx);
                    let ny_grad = (m_jp - m_jm) / (2.0 * dx);
                    let nz_grad = (m_kp - m_km) / (2.0 * dx);
                    let norm = (nx_grad * nx_grad + ny_grad * ny_grad + nz_grad * nz_grad).sqrt();
                    if norm > 1e-15 {
                        let omega = eta[k * nx * ny + j * nx + i];
                        let nx_n = nx_grad / norm;
                        let ny_n = ny_grad / norm;
                        let nz_n = nz_grad / norm;
                        // f = epsilon * (N × omega)
                        let fx = epsilon * (ny_n * omega[2] - nz_n * omega[1]);
                        let fy = epsilon * (nz_n * omega[0] - nx_n * omega[2]);
                        let _fz = epsilon * (nx_n * omega[1] - ny_n * omega[0]);
                        // Apply to y-face velocity (simplified: affect cell-centre via faces)
                        let fj = k * nx * (ny + 1) + j * nx + i;
                        self.v[fj] += fy;
                        let fi = k * (nx + 1) * ny + j * (nx + 1) + i;
                        self.u[fi] += fx;
                    }
                }
            }
        }
    }

    /// Compute simulation diagnostics.
    pub fn compute_stats(&self) -> FluidSimStats {
        let max_velocity = self.max_velocity();
        let total_kinetic_energy = self
            .u
            .iter()
            .chain(self.v.iter())
            .chain(self.w.iter())
            .map(|&vi| 0.5 * vi * vi)
            .sum();
        let div = self.gpu_compute_divergence();
        let max_divergence = div.iter().cloned().fold(0.0_f64, |a, v| a.max(v.abs()));
        FluidSimStats {
            max_velocity,
            total_kinetic_energy,
            max_divergence,
            pressure_residual: self.last_residual,
        }
    }
}

// ── Trilinear interpolation helpers (private) ─────────────────────────────────

fn clamp_idx(v: f64, n: usize) -> (usize, usize, f64) {
    let scaled = v.max(0.0).min((n as f64) - 1.0 - 1e-9);
    let i0 = scaled as usize;
    let i1 = (i0 + 1).min(n - 1);
    (i0, i1, scaled - i0 as f64)
}

fn trilinear_u(u: &[f64], x: f64, y: f64, z: f64, nx: usize, ny: usize, nz: usize, dx: f64) -> f64 {
    let (i0, i1, tx) = clamp_idx(x / dx, nx + 1);
    let (j0, j1, ty) = clamp_idx(y / dx - 0.5, ny);
    let (k0, k1, tz) = clamp_idx(z / dx - 0.5, nz);
    let idx = |i: usize, j: usize, k: usize| k * (nx + 1) * ny + j * (nx + 1) + i;
    let u000 = u[idx(i0, j0, k0)];
    let u100 = u[idx(i1, j0, k0)];
    let u010 = u[idx(i0, j1, k0)];
    let u110 = u[idx(i1, j1, k0)];
    let u001 = u[idx(i0, j0, k1)];
    let u101 = u[idx(i1, j0, k1)];
    let u011 = u[idx(i0, j1, k1)];
    let u111 = u[idx(i1, j1, k1)];
    let lerp = |a: f64, b: f64, t: f64| a + t * (b - a);
    lerp(
        lerp(lerp(u000, u100, tx), lerp(u010, u110, tx), ty),
        lerp(lerp(u001, u101, tx), lerp(u011, u111, tx), ty),
        tz,
    )
}

fn trilinear_v(v: &[f64], x: f64, y: f64, z: f64, nx: usize, ny: usize, nz: usize, dx: f64) -> f64 {
    let (i0, i1, tx) = clamp_idx(x / dx - 0.5, nx);
    let (j0, j1, ty) = clamp_idx(y / dx, ny + 1);
    let (k0, k1, tz) = clamp_idx(z / dx - 0.5, nz);
    let idx = |i: usize, j: usize, k: usize| k * nx * (ny + 1) + j * nx + i;
    let v000 = v[idx(i0, j0, k0)];
    let v100 = v[idx(i1, j0, k0)];
    let v010 = v[idx(i0, j1, k0)];
    let v110 = v[idx(i1, j1, k0)];
    let v001 = v[idx(i0, j0, k1)];
    let v101 = v[idx(i1, j0, k1)];
    let v011 = v[idx(i0, j1, k1)];
    let v111 = v[idx(i1, j1, k1)];
    let lerp = |a: f64, b: f64, t: f64| a + t * (b - a);
    lerp(
        lerp(lerp(v000, v100, tx), lerp(v010, v110, tx), ty),
        lerp(lerp(v001, v101, tx), lerp(v011, v111, tx), ty),
        tz,
    )
}

fn trilinear_w(
    w: &[f64],
    x: f64,
    z: f64,
    _i: usize,
    j: usize,
    k: usize,
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
) -> f64 {
    // Simplified: clamp back-trace to same j row
    let (i0, i1, tx) = clamp_idx(x / dx - 0.5, nx);
    let (k0, k1, tz) = clamp_idx(z / dx, nz + 1);
    let idx = |ii: usize, jj: usize, kk: usize| kk * nx * ny + jj * nx + ii;
    let j_c = j.min(ny - 1);
    let w000 = w[idx(i0, j_c, k0)];
    let w100 = w[idx(i1, j_c, k0)];
    let w001 = w[idx(i0, j_c, k1)];
    let w101 = w[idx(i1, j_c, k1)];
    let lerp = |a: f64, b: f64, t: f64| a + t * (b - a);
    let _ = (k, dx);
    lerp(lerp(w000, w100, tx), lerp(w001, w101, tx), tz)
}

// ── Free functions ────────────────────────────────────────────────────────────

/// Statistics for a fluid simulation step.
#[derive(Debug, Clone)]
pub struct FluidSimStats {
    /// Maximum speed across all velocity components.
    pub max_velocity: f64,
    /// Total kinetic energy (½ Σ v²).
    pub total_kinetic_energy: f64,
    /// Maximum absolute divergence across all cells.
    pub max_divergence: f64,
    /// Max-norm Poisson residual left by the last pressure projection
    /// (`0.0` before any projection), see
    /// [`GpuEulerGrid::gpu_pressure_projection`].
    pub pressure_residual: f64,
}

/// Initialise fluid velocity in a rectangular region of the grid.
///
/// Sets u, v, w to `(u0, v0, w0)` for all cells with indices in `region = [i0, i1, j0, j1, k0, k1]`.
pub fn marker_and_cell_init(
    grid: &mut GpuEulerGrid,
    region: [usize; 6],
    u0: f64,
    v0: f64,
    w0: f64,
) {
    let [ri0, ri1, rj0, rj1, rk0, rk1] = region;
    for k in rk0..rk1.min(grid.nz) {
        for j in rj0..rj1.min(grid.ny) {
            for i in ri0..ri1.min(grid.nx) {
                // Set x-face velocities
                let fi = k * (grid.nx + 1) * grid.ny + j * (grid.nx + 1) + i;
                grid.u[fi] = u0;
                // Set y-face velocities
                let fj = k * grid.nx * (grid.ny + 1) + j * grid.nx + i;
                grid.v[fj] = v0;
                // Set z-face velocities
                let fk = k * grid.nx * grid.ny + j * grid.nx + i;
                grid.w[fk] = w0;
            }
        }
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_grid() -> GpuEulerGrid {
        GpuEulerGrid::new(4, 4, 4, 0.1)
    }

    #[test]
    fn test_cell_count() {
        let g = make_grid();
        assert_eq!(g.cell_count(), 64);
    }

    #[test]
    fn test_cell_count_custom() {
        let g = GpuEulerGrid::new(3, 5, 7, 0.1);
        assert_eq!(g.cell_count(), 3 * 5 * 7);
    }

    #[test]
    fn test_index_origin() {
        let g = make_grid();
        assert_eq!(g.index(0, 0, 0), 0);
    }

    #[test]
    fn test_index_last_cell() {
        let g = make_grid();
        assert_eq!(g.index(3, 3, 3), 63);
    }

    #[test]
    fn test_index_roundtrip() {
        let g = make_grid();
        for k in 0..g.nz {
            for j in 0..g.ny {
                for i in 0..g.nx {
                    let flat = g.index(i, j, k);
                    let kk = flat / (g.nx * g.ny);
                    let jj = (flat % (g.nx * g.ny)) / g.nx;
                    let ii = flat % g.nx;
                    assert_eq!((ii, jj, kk), (i, j, k));
                }
            }
        }
    }

    #[test]
    fn test_initial_max_velocity_zero() {
        let g = make_grid();
        assert!(g.max_velocity() < 1e-15);
    }

    #[test]
    fn test_cfl_timestep_positive() {
        let mut g = make_grid();
        g.u[10] = 1.0;
        let dt = g.cfl_timestep(0.5);
        assert!(dt > 0.0);
    }

    #[test]
    fn test_cfl_timestep_no_velocity() {
        let g = make_grid();
        let dt = g.cfl_timestep(0.5);
        assert_eq!(dt, g.dx);
    }

    #[test]
    fn test_gravity_increases_v() {
        let mut g = make_grid();
        let v_before: Vec<f64> = g.v.clone();
        g.gpu_apply_gravity(-9.81, 0.01);
        let v_after = &g.v;
        // At least one v should have decreased (gravity = -9.81)
        let changed = v_before
            .iter()
            .zip(v_after.iter())
            .any(|(&a, &b)| (b - a).abs() > 1e-12);
        assert!(changed);
    }

    #[test]
    fn test_gravity_value() {
        let mut g = GpuEulerGrid::new(2, 2, 2, 0.1);
        g.gpu_apply_gravity(-10.0, 0.1);
        // All v values should be -1.0
        assert!(g.v.iter().all(|&vi| (vi + 1.0).abs() < 1e-12));
    }

    #[test]
    fn test_solid_bc_zeroes_boundary_u() {
        let mut g = make_grid();
        for v in g.u.iter_mut() {
            *v = 5.0;
        }
        g.gpu_enforce_solid_bc();
        // Check i=0 face
        for k in 0..g.nz {
            for j in 0..g.ny {
                assert_eq!(g.u[k * (g.nx + 1) * g.ny + j * (g.nx + 1)], 0.0);
            }
        }
    }

    #[test]
    fn test_solid_bc_zeroes_boundary_v() {
        let mut g = make_grid();
        for v in g.v.iter_mut() {
            *v = 3.0;
        }
        g.gpu_enforce_solid_bc();
        // Check j=0 face
        for k in 0..g.nz {
            for i in 0..g.nx {
                assert_eq!(g.v[k * g.nx * (g.ny + 1) + i], 0.0);
            }
        }
    }

    #[test]
    fn test_solid_bc_zeroes_boundary_w() {
        let mut g = make_grid();
        for w in g.w.iter_mut() {
            *w = 2.0;
        }
        g.gpu_enforce_solid_bc();
        // Check k=0 face
        for j in 0..g.ny {
            for i in 0..g.nx {
                assert_eq!(g.w[j * g.nx + i], 0.0);
            }
        }
    }

    #[test]
    fn test_divergence_zero_initially() {
        let g = make_grid();
        let div = g.gpu_compute_divergence();
        assert!(div.iter().all(|&d| d.abs() < 1e-12));
    }

    #[test]
    fn test_pressure_projection_reduces_divergence() {
        let mut g = GpuEulerGrid::new(4, 4, 4, 0.1);
        // Introduce divergence
        marker_and_cell_init(&mut g, [1, 3, 1, 3, 1, 3], 1.0, 0.5, 0.25);
        let div_before: f64 = g.gpu_compute_divergence().iter().map(|v| v.abs()).sum();
        assert!(div_before > 1.0, "the seeded field must be divergent");
        g.gpu_pressure_projection(100, 1e-6, 0.01, 1.0)
            .expect("valid projection parameters");
        g.gpu_update_velocity_from_pressure(0.01, 1.0)
            .expect("valid projection parameters");
        let div_after: f64 = g.gpu_compute_divergence().iter().map(|v| v.abs()).sum();
        // The old scaling removed only dt (1 %) of the divergence.
        assert!(
            div_after < 1e-2 * div_before,
            "divergence {div_before} -> {div_after}"
        );
    }

    #[test]
    fn test_pressure_projection_returns_residual() {
        let mut g = make_grid();
        let res = g
            .gpu_pressure_projection(10, 1e-6, 0.01, 1.0)
            .expect("valid projection parameters");
        // A field at rest has zero divergence: the zero pressure is exact.
        assert_eq!(res, 0.0);
        assert!(g.p.iter().all(|&p| p == 0.0));
    }

    #[test]
    fn test_update_velocity_from_pressure_finite() {
        let mut g = make_grid();
        g.p[0] = 10.0;
        g.gpu_update_velocity_from_pressure(0.01, 1.0)
            .expect("valid projection parameters");
        assert!(g.u.iter().all(|v| v.is_finite()));
        assert!(g.v.iter().all(|v| v.is_finite()));
    }

    // ── Pressure projection (regression for the dt-scaled divergence) ────────

    /// Deterministic pseudo-random value in `[-1, 1)`.
    fn next_unit(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((*state >> 11) as f64) / ((1_u64 << 53) as f64) * 2.0 - 1.0
    }

    /// An `n³` grid (spacing 0.1) with random face velocities and zero normal
    /// velocity on the solid boundary.
    fn random_grid(n: usize, seed: u64) -> GpuEulerGrid {
        let mut g = GpuEulerGrid::new(n, n, n, 0.1);
        let mut state = seed;
        for component in [&mut g.u, &mut g.v, &mut g.w] {
            for value in component.iter_mut() {
                *value = next_unit(&mut state);
            }
        }
        g.gpu_enforce_solid_bc();
        g
    }

    fn max_abs(values: &[f64]) -> f64 {
        values.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
    }

    /// Project `g` and return `(max|div| before, max|div| after)`.
    fn project(g: &mut GpuEulerGrid, sweeps: usize, dt: f64, rho: f64) -> (f64, f64) {
        let before = max_abs(&g.gpu_compute_divergence());
        g.gpu_pressure_projection(sweeps, 0.0, dt, rho)
            .expect("valid projection parameters");
        g.gpu_update_velocity_from_pressure(dt, rho)
            .expect("valid projection parameters");
        (before, max_abs(&g.gpu_compute_divergence()))
    }

    /// With dt = 0.01 the old projection left `(1 - dt)` of the divergence;
    /// the corrected one removes it to the accuracy of the Jacobi solve: 500
    /// sweeps shrink the slowest Jacobi mode of an 8^3 Neumann grid by
    /// `((2 + cos(pi/8)) / 3)^500`, about 3e-6, so the divergence of this
    /// random field falls by about seven orders of magnitude.
    #[test]
    fn test_projection_removes_divergence_of_random_field() {
        let mut g = random_grid(8, 0x00C0_FFEE);
        let (before, after) = project(&mut g, 500, 0.01, 1.0);
        assert!(before > 1.0, "the random field must be strongly divergent");
        assert!(
            after < 1e-3 * before,
            "max|div| only fell from {before} to {after}"
        );
    }

    /// `dt` and `rho` cancel between the Poisson right-hand side and the
    /// gradient update, so the projected velocity does not depend on them.
    #[test]
    fn test_projection_is_independent_of_dt_and_rho() {
        let reference = {
            let mut g = random_grid(6, 7);
            project(&mut g, 400, 0.01, 1.0);
            g
        };
        for &(dt, rho) in &[(0.5, 1.0), (0.001, 997.0), (2.0, 0.25)] {
            let mut g = random_grid(6, 7);
            let (before, after) = project(&mut g, 400, dt, rho);
            assert!(before / after > 100.0, "dt={dt} rho={rho}");
            for (name, got, want) in [
                ("u", &g.u, &reference.u),
                ("v", &g.v, &reference.v),
                ("w", &g.w, &reference.w),
            ] {
                for (a, b) in got.iter().zip(want.iter()) {
                    assert!(
                        (a - b).abs() < 1e-8,
                        "{name} differs for dt={dt} rho={rho}: {a} vs {b}"
                    );
                }
            }
        }
    }

    /// The velocity update is `u -= (dt / rho) * dp / dx` on interior faces
    /// and leaves boundary faces alone.
    #[test]
    fn test_velocity_update_scales_with_dt_over_rho() {
        let (dt, rho) = (0.2, 4.0);
        let mut g = GpuEulerGrid::new(3, 3, 3, 0.5);
        // p = 2 i + 3 j + 5 k  => constant gradients (2, 3, 5) / dx.
        for k in 0..3 {
            for j in 0..3 {
                for i in 0..3 {
                    let idx = g.index(i, j, k);
                    g.p[idx] = 2.0 * i as f64 + 3.0 * j as f64 + 5.0 * k as f64;
                }
            }
        }
        g.gpu_update_velocity_from_pressure(dt, rho)
            .expect("valid projection parameters");
        let factor = dt / rho / 0.5;
        // x-face (i = 1, j = 1, k = 1): index k (nx+1) ny + j (nx+1) + i.
        assert!((g.u[4 * 3 + 4 + 1] + factor * 2.0).abs() < 1e-14);
        assert!((g.v[3 * 4 + 3 + 1] + factor * 3.0).abs() < 1e-14);
        assert!((g.w[3 * 3 + 3 + 1] + factor * 5.0).abs() < 1e-14);
        // Boundary faces stay untouched.
        assert_eq!(g.u[0], 0.0);
        assert_eq!(g.v[0], 0.0);
        assert_eq!(g.w[0], 0.0);
        assert_eq!(g.u[3], 0.0);
    }

    /// The reported residual is the max-norm of `∇²p - (ρ/Δt) ∇·u`,
    /// recomputed here with an independent stencil.
    #[test]
    fn test_projection_returns_the_poisson_residual() {
        let (dt, rho) = (0.02, 3.0);
        let mut g = random_grid(5, 99);
        let rhs: Vec<f64> = g
            .gpu_compute_divergence()
            .iter()
            .map(|d| rho / dt * d)
            .collect();
        let returned = g
            .gpu_pressure_projection(7, 0.0, dt, rho)
            .expect("valid projection parameters");
        let mut worst = 0.0_f64;
        for k in 0..5_usize {
            for j in 0..5_usize {
                for i in 0..5_usize {
                    let idx = g.index(i, j, k);
                    let mut laplacian = 0.0;
                    for (di, dj, dk) in [
                        (1_i32, 0_i32, 0_i32),
                        (-1, 0, 0),
                        (0, 1, 0),
                        (0, -1, 0),
                        (0, 0, 1),
                        (0, 0, -1),
                    ] {
                        let (ni, nj, nk) = (i as i32 + di, j as i32 + dj, k as i32 + dk);
                        if (0..5).contains(&ni) && (0..5).contains(&nj) && (0..5).contains(&nk) {
                            let nb = g.index(ni as usize, nj as usize, nk as usize);
                            laplacian += g.p[nb] - g.p[idx];
                        }
                    }
                    laplacian /= g.dx * g.dx;
                    worst = worst.max((laplacian - rhs[idx]).abs());
                }
            }
        }
        assert!(returned > 0.0, "seven sweeps cannot solve the system");
        assert!(
            (returned - worst).abs() <= 1e-9 * worst,
            "returned {returned} vs independent {worst}"
        );
        assert_eq!(g.last_residual, returned);
    }

    /// More sweeps drive the Poisson residual down and `compute_stats`
    /// reports it instead of a hard-coded zero.
    #[test]
    fn test_compute_stats_reports_pressure_residual() {
        let (dt, rho) = (0.01, 1.0);
        let mut g = random_grid(6, 5);
        assert_eq!(g.compute_stats().pressure_residual, 0.0);
        let rhs_scale = rho / dt * max_abs(&g.gpu_compute_divergence());
        let coarse = g
            .gpu_pressure_projection(3, 0.0, dt, rho)
            .expect("valid projection parameters");
        assert_eq!(g.compute_stats().pressure_residual, coarse);
        let fine = g
            .gpu_pressure_projection(600, 0.0, dt, rho)
            .expect("valid projection parameters");
        let stats = g.compute_stats();
        assert_eq!(stats.pressure_residual, fine);
        assert!(fine < coarse, "residual must fall: {coarse} -> {fine}");
        assert!(
            fine < 1e-3 * rhs_scale,
            "converged residual {fine} vs rhs scale {rhs_scale}"
        );
    }

    /// `tol` stops the iteration after the first sweep when the update is
    /// already below it.
    #[test]
    fn test_projection_tolerance_stops_early() {
        let mut one = random_grid(4, 11);
        let mut many = one.clone();
        one.gpu_pressure_projection(1, 0.0, 0.01, 1.0)
            .expect("valid projection parameters");
        many.gpu_pressure_projection(50, 1e30, 0.01, 1.0)
            .expect("valid projection parameters");
        assert_eq!(one.p, many.p);
        assert_eq!(one.last_residual, many.last_residual);
    }

    /// Invalid time steps / densities and non-finite velocities are rejected
    /// without touching the grid.
    #[test]
    fn test_projection_rejects_invalid_parameters() {
        let mut g = random_grid(4, 3);
        let u_before = g.u.clone();
        for &(dt, rho) in &[
            (0.0, 1.0),
            (-0.01, 1.0),
            (f64::NAN, 1.0),
            (f64::INFINITY, 1.0),
            (0.01, 0.0),
            (0.01, -1.0),
            (0.01, f64::NAN),
        ] {
            assert!(g.gpu_pressure_projection(10, 0.0, dt, rho).is_err());
            assert!(g.gpu_update_velocity_from_pressure(dt, rho).is_err());
        }
        assert_eq!(g.u, u_before);
        assert!(g.p.iter().all(|&p| p == 0.0));
        g.u[5] = f64::NAN;
        assert!(g.gpu_pressure_projection(10, 0.0, 0.01, 1.0).is_err());
    }

    #[test]
    fn test_marker_and_cell_init_sets_u() {
        let mut g = make_grid();
        marker_and_cell_init(&mut g, [0, 2, 0, 2, 0, 2], 3.0, 0.0, 0.0);
        // At least some u values should be 3.0
        assert!(g.u.iter().any(|&v| (v - 3.0).abs() < 1e-12));
    }

    #[test]
    fn test_vorticity_at_zero_field() {
        let g = make_grid();
        let vort = g.vorticity_at(2, 2, 2);
        assert!(vort.iter().all(|v| v.abs() < 1e-12));
    }

    #[test]
    fn test_divergence_at_cell() {
        let g = make_grid();
        let d = g.divergence_at(0, 0, 0);
        assert!(d.abs() < 1e-12);
    }

    #[test]
    fn test_stats_max_velocity_nonneg() {
        let g = make_grid();
        let stats = g.compute_stats();
        assert!(stats.max_velocity >= 0.0);
    }

    #[test]
    fn test_stats_kinetic_energy_nonneg() {
        let g = make_grid();
        let stats = g.compute_stats();
        assert!(stats.total_kinetic_energy >= 0.0);
    }

    #[test]
    fn test_stats_max_divergence_nonneg() {
        let g = make_grid();
        let stats = g.compute_stats();
        assert!(stats.max_divergence >= 0.0);
    }

    #[test]
    fn test_advect_does_not_panic() {
        let mut g = GpuEulerGrid::new(4, 4, 4, 0.1);
        g.u[5] = 0.5;
        g.gpu_advect_semi_lagrange(0.01);
        // Should complete without panic
    }

    #[test]
    fn test_vorticity_confinement_does_not_panic() {
        let mut g = GpuEulerGrid::new(4, 4, 4, 0.1);
        g.u[10] = 1.0;
        g.gpu_vorticity_confinement(0.01);
    }

    #[test]
    fn test_grid_clone() {
        let g = make_grid();
        let g2 = g.clone();
        assert_eq!(g2.cell_count(), g.cell_count());
    }

    #[test]
    fn test_fluid_sim_stats_debug() {
        let stats = FluidSimStats {
            max_velocity: 1.0,
            total_kinetic_energy: 2.0,
            max_divergence: 0.1,
            pressure_residual: 0.01,
        };
        let s = format!("{stats:?}");
        assert!(s.contains("FluidSimStats"));
    }

    #[test]
    fn test_cfl_formula() {
        let mut g = GpuEulerGrid::new(4, 4, 4, 0.2);
        g.u[0] = 2.0;
        let dt = g.cfl_timestep(0.5);
        assert!((dt - 0.5 * 0.2 / 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_rho_initialised_to_one() {
        let g = make_grid();
        assert!(g.rho.iter().all(|&r| (r - 1.0).abs() < 1e-12));
    }
}
