// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the D3Q19 MRT matrix and collision, the D2Q9 TRT magic
//! parameter and the cascaded central-moment collision.

use super::*;
use crate::collision::compute_equilibrium_d3q19;
use crate::mrt3d;

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

/// D2Q9 velocity x-components (ordering of `mrt.rs`).
const CX9: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
/// D2Q9 velocity y-components.
const CY9: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

/// Indices of the conserved moments (ρ, jx, jy, jz) in the D3Q19 moment basis.
const CONSERVED: [usize; 4] = [0, 3, 5, 7];

fn matmul19(a: &[[f64; 19]; 19], b: &[[f64; 19]; 19]) -> [[f64; 19]; 19] {
    let mut out = [[0.0_f64; 19]; 19];
    for (out_row, a_row) in out.iter_mut().zip(a.iter()) {
        for (c, out_rc) in out_row.iter_mut().enumerate() {
            *out_rc = a_row
                .iter()
                .zip(b.iter())
                .map(|(&a_rk, b_row)| a_rk * b_row[c])
                .sum();
        }
    }
    out
}

fn transpose19(a: &[[f64; 19]; 19]) -> [[f64; 19]; 19] {
    let mut out = [[0.0_f64; 19]; 19];
    for (r, a_row) in a.iter().enumerate() {
        for (c, &a_rc) in a_row.iter().enumerate() {
            out[c][r] = a_rc;
        }
    }
    out
}

/// Density and velocity of a D3Q19 population.
fn macros19(f: &[f64; 19]) -> (f64, [f64; 3]) {
    let rho: f64 = f.iter().sum();
    let j = momentum19(f);
    (rho, [j[0] / rho, j[1] / rho, j[2] / rho])
}

/// Momentum `Σ f_i c_i` of a D3Q19 population.
fn momentum19(f: &[f64; 19]) -> [f64; 3] {
    let mut j = [0.0_f64; 3];
    for (i, &fi) in f.iter().enumerate() {
        j[0] += fi * f64::from(MRT_EX[i]);
        j[1] += fi * f64::from(MRT_EY[i]);
        j[2] += fi * f64::from(MRT_EZ[i]);
    }
    j
}

/// Equilibrium plus the populations that carry the moment perturbation `dm`
/// (given in the row order of `mrt3d::M`; the conserved rows must be zero).
fn state_with_moment_perturbation(rho: f64, u: [f64; 3], dm: &[f64; 19]) -> [f64; 19] {
    let mut f = compute_equilibrium_d3q19(rho, u);
    let delta_f = mrt3d::inverse_transform(dm);
    for (fi, d) in f.iter_mut().zip(delta_f.iter()) {
        *fi += d;
    }
    f
}

/// A moment perturbation with a distinct non-zero value in every
/// non-conserved row of `mrt3d::M`.
fn full_moment_perturbation() -> [f64; 19] {
    let mut dm = [0.0_f64; 19];
    for (k, d) in dm.iter_mut().enumerate() {
        if !CONSERVED.contains(&k) {
            *d = 1.0e-3 * (1.0 + 0.1 * k as f64);
        }
    }
    dm
}

/// Rates (in `MrtRelaxation` order) that differ for every non-conserved mode.
fn distinct_rates() -> MrtRelaxation {
    let mut s = [0.0_f64; 19];
    for (k, s_k) in s.iter_mut().enumerate() {
        if !CONSERVED.contains(&k) {
            *s_k = 0.30 + 0.08 * k as f64;
        }
    }
    MrtRelaxation { s }
}

/// Moments `M f` in the row order of `mrt3d::M`.
fn moments(f: &[f64; 19]) -> [f64; 19] {
    mrt3d::transform_to_moment_space(f)
}

// ---------------------------------------------------------------------------
// F2608 — D3Q19 transformation matrix
// ---------------------------------------------------------------------------

/// `M · M⁻¹ = M⁻¹ · M = I` to 1e-12.  The previous matrix had rank 18, so its
/// "inverse" was garbage.
#[test]
fn transform_matrix_times_inverse_is_identity() {
    let m = mrt_transform_matrix();
    let inv = mrt_inverse_matrix();
    let left = matmul19(&m, &inv);
    let right = matmul19(&inv, &m);
    for r in 0..19 {
        for c in 0..19 {
            let expected = if r == c { 1.0 } else { 0.0 };
            assert!(
                (left[r][c] - expected).abs() < 1e-12,
                "(M·M⁻¹)[{r}][{c}] = {}",
                left[r][c]
            );
            assert!(
                (right[r][c] - expected).abs() < 1e-12,
                "(M⁻¹·M)[{r}][{c}] = {}",
                right[r][c]
            );
        }
    }
}

/// The inverse is a bounded matrix: the largest entry of the exact inverse is
/// 1/4, while the inverse of the singular matrix blew up.
#[test]
fn inverse_matrix_is_finite_and_bounded() {
    let inv = mrt_inverse_matrix();
    let largest = inv
        .iter()
        .flatten()
        .fold(0.0_f64, |acc, &x| acc.max(x.abs()));
    assert!(inv.iter().flatten().all(|x| x.is_finite()));
    assert!(
        (largest - 0.25).abs() < 1e-12,
        "largest |M⁻¹| entry = {largest}"
    );
}

/// The rows of `M` are mutually orthogonal with the d'Humières norms.
#[test]
fn transform_matrix_rows_are_orthogonal_with_known_norms() {
    let m = mrt_transform_matrix();
    let gram = matmul19(&m, &transpose19(&m));
    let norms: [f64; 19] = [
        19.0, 2394.0, 252.0, 10.0, 40.0, 10.0, 40.0, 10.0, 40.0, 36.0, 72.0, 12.0, 24.0, 4.0, 4.0,
        4.0, 8.0, 8.0, 8.0,
    ];
    for (r, (gram_row, &norm)) in gram.iter().zip(norms.iter()).enumerate() {
        for (c, &entry) in gram_row.iter().enumerate() {
            let expected = if r == c { norm } else { 0.0 };
            assert!(
                (entry - expected).abs() < 1e-12,
                "Gram[{r}][{c}] = {entry}, expected {expected}"
            );
        }
    }
}

/// Every row equals the polynomial of the documentation, evaluated at the
/// lattice velocities.
#[test]
fn transform_matrix_matches_documented_polynomials() {
    let m = mrt_transform_matrix();
    for i in 0..19 {
        let x = f64::from(MRT_EX[i]);
        let y = f64::from(MRT_EY[i]);
        let z = f64::from(MRT_EZ[i]);
        let r2 = x * x + y * y + z * z;
        let expected: [f64; 19] = [
            1.0,
            19.0 * r2 - 30.0,
            (21.0 * r2 * r2 - 53.0 * r2 + 24.0) / 2.0,
            x,
            (5.0 * r2 - 9.0) * x,
            y,
            (5.0 * r2 - 9.0) * y,
            z,
            (5.0 * r2 - 9.0) * z,
            3.0 * x * x - r2,
            (3.0 * r2 - 5.0) * (3.0 * x * x - r2),
            y * y - z * z,
            (3.0 * r2 - 5.0) * (y * y - z * z),
            x * y,
            x * z,
            y * z,
            (y * y - z * z) * x,
            (z * z - x * x) * y,
            (x * x - y * y) * z,
        ];
        for (row, &e) in expected.iter().enumerate() {
            assert!(
                (m[row][i] - e).abs() < 1e-12,
                "M[{row}][{i}] = {}, expected {e}",
                m[row][i]
            );
        }
    }
}

/// Literature values: ε = 12, −4, 1 and the heat flux ±4, ±1 on the rest,
/// face and edge velocities (the previous rows gave 12, 1.5, −9 and −1, 2).
#[test]
fn transform_matrix_has_literature_energy_and_heat_flux_rows() {
    let m = mrt_transform_matrix();
    assert_eq!(m[1][0], -30.0);
    assert_eq!(m[1][1], -11.0);
    assert_eq!(m[1][7], 8.0);
    assert_eq!(m[2][0], 12.0);
    assert_eq!(m[2][1], -4.0);
    assert_eq!(m[2][7], 1.0);
    // q_x on (1,0,0), (-1,0,0) and (1,1,0).
    assert_eq!(m[4][1], -4.0);
    assert_eq!(m[4][2], 4.0);
    assert_eq!(m[4][7], 1.0);
}

/// `mrt_transform_matrix` is the `mrt3d` matrix with the stress rows 14/15
/// exchanged, and the inverse is the `mrt3d` inverse with the columns exchanged.
#[test]
fn transform_matrix_agrees_with_mrt3d_up_to_the_stress_swap() {
    let m = mrt_transform_matrix();
    let inv = mrt_inverse_matrix();
    let perm = |k: usize| match k {
        14 => 15,
        15 => 14,
        other => other,
    };
    for r in 0..19 {
        for c in 0..19 {
            assert_eq!(m[r][c], mrt3d::M[perm(r)][c], "M[{r}][{c}]");
            assert!(
                (inv[r][c] - mrt3d::M_INV[r][perm(c)]).abs() < 1e-12,
                "M⁻¹[{r}][{c}] = {}, mrt3d: {}",
                inv[r][c],
                mrt3d::M_INV[r][perm(c)]
            );
        }
    }
}

/// `mrt3d::M · mrt3d::M_INV = I` to near machine precision (the collision
/// below relies on it).
#[test]
fn mrt3d_matrix_pair_is_consistent() {
    let left = matmul19(&mrt3d::M, &mrt3d::M_INV);
    for (r, row) in left.iter().enumerate() {
        for (c, &value) in row.iter().enumerate() {
            let expected = if r == c { 1.0 } else { 0.0 };
            assert!(
                (value - expected).abs() < 1e-14,
                "(M·M_INV)[{r}][{c}] = {value}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// F2609 — MrtD3Q19::collide_cell is a real MRT operator
// ---------------------------------------------------------------------------

#[test]
fn relaxation_rates_convert_to_mrt3d_order() {
    let relax = distinct_rates();
    let converted = relax.to_mrt3d_order();
    for (k, &rate) in converted.iter().enumerate() {
        let src = match k {
            14 => 15,
            15 => 14,
            other => other,
        };
        assert_eq!(rate, relax.s[src], "rate {k}");
    }
}

#[test]
fn bgk_equivalent_has_zero_conserved_and_unit_ghost_rates() {
    let nu = 0.05;
    let relax = MrtRelaxation::bgk_equivalent(nu);
    let s_nu = 1.0 / (3.0 * nu + 0.5);
    for (k, &s_k) in relax.s.iter().enumerate() {
        let expected = if CONSERVED.contains(&k) {
            0.0
        } else if [9, 11, 13, 14, 15].contains(&k) {
            s_nu
        } else {
            1.0
        };
        assert!((s_k - expected).abs() < 1e-14, "s[{k}] = {s_k}");
    }
}

/// Every non-conserved moment relaxes toward its equilibrium at its own rate:
/// `m*_k = m_eq,k + (1 − s_k) (m_k − m_eq,k)`.  A BGK operator driven by the
/// single rate `s[9]` cannot satisfy this for distinct rates.
#[test]
fn collide_cell_relaxes_each_moment_at_its_own_rate() {
    let relax = distinct_rates();
    let sim = MrtD3Q19::new(2, 2, 2, MrtRelaxation { s: relax.s });
    let rho = 1.05;
    let u = [0.04, -0.02, 0.03];
    let dm = full_moment_perturbation();
    let f = state_with_moment_perturbation(rho, u, &dm);
    let (rho_f, u_f) = macros19(&f);
    assert!((rho_f - rho).abs() < 1e-13);

    let f_star = sim.collide_cell(f, rho_f, u_f);
    let m_eq = moments(&compute_equilibrium_d3q19(rho_f, u_f));
    let m_star = moments(&f_star);
    let s_mrt3d = relax.to_mrt3d_order();
    for k in 0..19 {
        let expected = m_eq[k] + (1.0 - s_mrt3d[k]) * dm[k];
        assert!(
            (m_star[k] - expected).abs() < 1e-12,
            "moment {k}: got {}, expected {expected}",
            m_star[k]
        );
    }
}

/// Changing only the ghost rates leaves the conserved and the stress moments
/// untouched but does change the ghost moments (so the operator is not BGK).
#[test]
fn ghost_rates_change_only_the_ghost_moments() {
    let dm = full_moment_perturbation();
    let f = state_with_moment_perturbation(1.0, [0.03, 0.02, -0.01], &dm);
    let (rho, u) = macros19(&f);

    let base = MrtRelaxation::bgk_equivalent(0.1);
    let mut changed = MrtRelaxation { s: base.s };
    for k in [10, 12, 16, 17, 18] {
        changed.s[k] = 0.4;
    }
    let sim_a = MrtD3Q19::new(2, 2, 2, base);
    let sim_b = MrtD3Q19::new(2, 2, 2, changed);
    let m_a = moments(&sim_a.collide_cell(f, rho, u));
    let m_b = moments(&sim_b.collide_cell(f, rho, u));

    for k in [0, 3, 5, 7, 9, 11, 13, 14, 15] {
        assert!(
            (m_a[k] - m_b[k]).abs() < 1e-13,
            "moment {k} changed with the ghost rates: {} vs {}",
            m_a[k],
            m_b[k]
        );
    }
    for k in [10, 12, 16, 17, 18] {
        assert!(
            (m_a[k] - m_b[k]).abs() > 1e-4,
            "ghost moment {k} did not respond to its rate: {} vs {}",
            m_a[k],
            m_b[k]
        );
    }
}

/// The shear rates act on their own rows only, and `s[14]`/`s[15]` hit the
/// `p_xz`/`p_yz` rows even though `mrt3d` stores them in the opposite order.
#[test]
fn each_shear_rate_changes_only_its_own_moment() {
    let dm = full_moment_perturbation();
    let f = state_with_moment_perturbation(1.0, [0.02, -0.03, 0.01], &dm);
    let (rho, u) = macros19(&f);
    let base = MrtRelaxation::bgk_equivalent(0.1);
    let m_base =
        moments(&MrtD3Q19::new(2, 2, 2, MrtRelaxation { s: base.s }).collide_cell(f, rho, u));

    // (index in `MrtRelaxation`, row of `mrt3d::M` that must react)
    for (index, row) in [(9, 9), (11, 11), (13, 13), (14, 15), (15, 14)] {
        let mut s = base.s;
        s[index] = 0.7;
        let m_changed =
            moments(&MrtD3Q19::new(2, 2, 2, MrtRelaxation { s }).collide_cell(f, rho, u));
        for k in 0..19 {
            let diff = (m_changed[k] - m_base[k]).abs();
            if k == row {
                assert!(diff > 1e-5, "s[{index}] left moment {row} unchanged");
            } else {
                assert!(diff < 1e-13, "s[{index}] changed moment {k} by {diff}");
            }
        }
    }
}

/// Mass and momentum are conserved to 1e-13 with the full, distinct rate set —
/// also when `(rho, u)` is not the pair of moments of `f`.
#[test]
fn collide_cell_conserves_mass_and_momentum() {
    let sim = MrtD3Q19::new(2, 2, 2, distinct_rates());
    let dm = full_moment_perturbation();
    let f = state_with_moment_perturbation(1.1, [0.05, 0.03, -0.04], &dm);
    let (rho, u) = macros19(&f);
    let j_in = momentum19(&f);
    for (rho_arg, u_arg) in [(rho, u), (0.7, [0.2, -0.1, 0.05])] {
        let f_star = sim.collide_cell(f, rho_arg, u_arg);
        let mass_out: f64 = f_star.iter().sum();
        let j_out = momentum19(&f_star);
        assert!((mass_out - rho).abs() < 1e-13, "mass {rho} -> {mass_out}");
        for d in 0..3 {
            assert!(
                (j_out[d] - j_in[d]).abs() < 1e-13,
                "momentum[{d}] {} -> {}",
                j_in[d],
                j_out[d]
            );
        }
    }
}

/// With one rate for every mode the MRT operator is the BGK operator.
#[test]
fn uniform_rates_reduce_to_bgk() {
    let omega = 0.8;
    let sim = MrtD3Q19::new(2, 2, 2, MrtRelaxation { s: [omega; 19] });
    let dm = full_moment_perturbation();
    let f = state_with_moment_perturbation(0.95, [0.02, 0.04, -0.03], &dm);
    let (rho, u) = macros19(&f);
    let feq = compute_equilibrium_d3q19(rho, u);
    let f_star = sim.collide_cell(f, rho, u);
    for i in 0..19 {
        let bgk = f[i] - omega * (f[i] - feq[i]);
        assert!(
            (f_star[i] - bgk).abs() < 1e-13,
            "population {i}: MRT {}, BGK {bgk}",
            f_star[i]
        );
    }
}

/// A different rate for the energy mode makes MRT differ from BGK at the
/// viscosity-matched rate, which the previous `collide_cell` could not do.
#[test]
fn bgk_equivalent_with_unit_ghost_rates_is_not_bgk() {
    let nu = 0.05;
    let relax = MrtRelaxation::bgk_equivalent(nu);
    let omega = relax.s[9];
    let sim = MrtD3Q19::new(2, 2, 2, relax);
    let dm = full_moment_perturbation();
    let f = state_with_moment_perturbation(1.0, [0.03, 0.0, 0.01], &dm);
    let (rho, u) = macros19(&f);
    let feq = compute_equilibrium_d3q19(rho, u);
    let f_star = sim.collide_cell(f, rho, u);
    let deviation = (0..19)
        .map(|i| (f_star[i] - (f[i] - omega * (f[i] - feq[i]))).abs())
        .fold(0.0_f64, f64::max);
    assert!(deviation > 1e-5, "MRT is indistinguishable from BGK");
}

/// Shear waves of the two out-of-plane stress modes decay with the viscosity
/// of the rate that belongs to them: `u_x(z)` with `s[14]` (p_xz) and `u_y(z)`
/// with `s[15]` (p_yz), `ν = (1/s − 1/2)/3`.
#[test]
fn shear_wave_decay_follows_the_matching_stress_rate() {
    let nz = 32;
    let steps = 100;
    let amplitude = 0.01;
    let s_xz = 1.2;
    let s_yz = 0.8;
    let k2 = (2.0 * std::f64::consts::PI / nz as f64).powi(2);

    for (component, s_own, s_other) in [(0usize, s_xz, s_yz), (1, s_yz, s_xz)] {
        let mut relax = MrtRelaxation::bgk_equivalent(0.05);
        relax.s[14] = s_xz;
        relax.s[15] = s_yz;
        let mut sim = MrtD3Q19::new(2, 2, nz, relax);
        for z in 0..nz {
            let phase = 2.0 * std::f64::consts::PI * z as f64 / nz as f64;
            let mut u = [0.0_f64; 3];
            u[component] = amplitude * phase.sin();
            for y in 0..2 {
                for x in 0..2 {
                    let idx = sim.idx(x, y, z);
                    sim.f[idx] = MrtD3Q19::equilibrium(1.0, u);
                }
            }
        }
        for _ in 0..steps {
            sim.step();
        }
        let mut projection = 0.0;
        let mut cells = 0.0;
        for z in 0..nz {
            let phase = 2.0 * std::f64::consts::PI * z as f64 / nz as f64;
            for y in 0..2 {
                for x in 0..2 {
                    let (_, u) = sim.macros(sim.idx(x, y, z));
                    projection += u[component] * phase.sin();
                    cells += 1.0;
                }
            }
        }
        let measured = 2.0 * projection / cells / amplitude;
        let nu_own = (1.0 / s_own - 0.5) / 3.0;
        let nu_other = (1.0 / s_other - 0.5) / 3.0;
        let expected = (-nu_own * k2 * steps as f64).exp();
        let wrong = (-nu_other * k2 * steps as f64).exp();
        assert!(
            (measured / expected - 1.0).abs() < 0.01,
            "component {component}: measured {measured}, expected {expected}"
        );
        assert!(
            (measured / wrong - 1.0).abs() > 0.2,
            "component {component}: decay is indistinguishable from the other stress rate"
        );
    }
}

/// A full MRT step conserves the total mass and momentum of the periodic box.
#[test]
fn mrt_step_conserves_global_mass_and_momentum() {
    let mut sim = MrtD3Q19::new(4, 4, 4, distinct_rates());
    for k in 0..sim.f.len() {
        let phase = k as f64 * 0.37;
        let u = [
            0.02 * phase.sin(),
            0.015 * phase.cos(),
            0.01 * (2.0 * phase).sin(),
        ];
        sim.f[k] = MrtD3Q19::equilibrium(1.0 + 0.05 * phase.cos(), u);
    }
    let totals = |sim: &MrtD3Q19| {
        let mass: f64 = sim.f.iter().map(|c| c.iter().sum::<f64>()).sum();
        let mut j = [0.0_f64; 3];
        for cell in &sim.f {
            let jc = momentum19(cell);
            for d in 0..3 {
                j[d] += jc[d];
            }
        }
        (mass, j)
    };
    let (mass0, j0) = totals(&sim);
    for _ in 0..20 {
        sim.step();
    }
    let (mass1, j1) = totals(&sim);
    assert!((mass1 - mass0).abs() < 1e-10, "mass {mass0} -> {mass1}");
    for d in 0..3 {
        assert!(
            (j1[d] - j0[d]).abs() < 1e-10,
            "momentum[{d}] {} -> {}",
            j0[d],
            j1[d]
        );
    }
}

// ---------------------------------------------------------------------------
// D2Q9 helpers shared by the TRT and cascaded tests
// ---------------------------------------------------------------------------

/// Density and velocity of a D2Q9 population.
fn macros9(f: &[f64; 9]) -> (f64, f64, f64) {
    let rho: f64 = f.iter().sum();
    let jx: f64 = f.iter().zip(CX9.iter()).map(|(a, c)| a * c).sum();
    let jy: f64 = f.iter().zip(CY9.iter()).map(|(a, c)| a * c).sum();
    (rho, jx / rho, jy / rho)
}

/// Momentum `(Σ f cx, Σ f cy)` of a D2Q9 population.
fn momentum9(f: &[f64; 9]) -> (f64, f64) {
    let jx: f64 = f.iter().zip(CX9.iter()).map(|(a, c)| a * c).sum();
    let jy: f64 = f.iter().zip(CY9.iter()).map(|(a, c)| a * c).sum();
    (jx, jy)
}

/// Periodic streaming of a node-major D2Q9 field.
fn stream_periodic(f: &[[f64; 9]], nx: usize, ny: usize) -> Vec<[f64; 9]> {
    let mut out = vec![[0.0_f64; 9]; nx * ny];
    for y in 0..ny {
        for x in 0..nx {
            for q in 0..9 {
                let xd = (x as i64 + CX9[q] as i64).rem_euclid(nx as i64) as usize;
                let yd = (y as i64 + CY9[q] as i64).rem_euclid(ny as i64) as usize;
                out[yd * nx + xd][q] = f[y * nx + x][q];
            }
        }
    }
    out
}

/// Orientation of the periodic shear wave used to probe the viscosity.
#[derive(Clone, Copy)]
enum Wave {
    /// `u_x = U sin(2π y / N)`: excites the shear stress `Π_xy` only.
    Vertical,
    /// `u = U sin(2π (x + y) / N) (1, −1)`: excites the normal-stress
    /// difference `Π_xx − Π_yy` only.
    Diagonal,
}

impl Wave {
    fn velocity(self, x: usize, y: usize, n: usize, amplitude: f64) -> (f64, f64) {
        let tau = 2.0 * std::f64::consts::PI / n as f64;
        match self {
            Wave::Vertical => (amplitude * (tau * y as f64).sin(), 0.0),
            Wave::Diagonal => {
                let u = amplitude * (tau * (x + y) as f64).sin();
                (u, -u)
            }
        }
    }

    fn phase(self, x: usize, y: usize, n: usize) -> f64 {
        let tau = 2.0 * std::f64::consts::PI / n as f64;
        match self {
            Wave::Vertical => (tau * y as f64).sin(),
            Wave::Diagonal => (tau * (x + y) as f64).sin(),
        }
    }

    /// Squared wave number `|k|²` on an `n × n` periodic box.
    fn k_squared(self, n: usize) -> f64 {
        let k = 2.0 * std::f64::consts::PI / n as f64;
        match self {
            Wave::Vertical => k * k,
            Wave::Diagonal => 2.0 * k * k,
        }
    }
}

/// Run a decaying periodic shear wave and return `(measured, exact_theory)`
/// amplitude ratios after `steps` steps for a fluid of viscosity `nu`.
fn shear_wave_ratio(
    wave: Wave,
    nu: f64,
    steps: usize,
    equilibrium: impl Fn(f64, f64, f64) -> [f64; 9],
    collide: impl Fn(&[f64; 9]) -> [f64; 9],
) -> (f64, f64) {
    let n = 32;
    let amplitude = 0.01;
    let mut field = vec![[0.0_f64; 9]; n * n];
    for y in 0..n {
        for x in 0..n {
            let (ux, uy) = wave.velocity(x, y, n, amplitude);
            field[y * n + x] = equilibrium(1.0, ux, uy);
        }
    }
    for _ in 0..steps {
        for cell in field.iter_mut() {
            *cell = collide(cell);
        }
        field = stream_periodic(&field, n, n);
    }
    let mut projection = 0.0;
    for y in 0..n {
        for x in 0..n {
            let (_, ux, _) = macros9(&field[y * n + x]);
            projection += ux * wave.phase(x, y, n);
        }
    }
    let measured = 2.0 * projection / (n * n) as f64 / amplitude;
    let theory = (-nu * wave.k_squared(n) * steps as f64).exp();
    (measured, theory)
}

fn deterministic_perturbation(base: &[f64; 9], scale: f64) -> [f64; 9] {
    let mut f = *base;
    for (i, fi) in f.iter_mut().enumerate() {
        *fi += scale * ((i as f64 + 1.0) * 1.7).sin();
    }
    f
}

// ---------------------------------------------------------------------------
// F2611 — TRT magic parameter
// ---------------------------------------------------------------------------

#[test]
fn trt_magic_parameter_is_stable_for_all_viscosities() {
    for &nu in &[0.001_f64, 0.01, 0.1, 1.0 / 6.0, 0.5, 1.0, 5.0] {
        let trt = TrtCollision2D::new(nu);
        let lambda = (trt.tau_plus - 0.5) * (trt.tau_minus - 0.5);
        assert!(
            (lambda - 3.0 / 16.0).abs() < 1e-12 * lambda.max(1.0),
            "ν={nu}"
        );
        assert!((trt.lambda() - lambda).abs() < 1e-14);
        let s_minus = 1.0 / trt.tau_minus;
        assert!(s_minus > 0.0 && s_minus < 2.0, "ν={nu}: s⁻ = {s_minus}");
        assert!((trt.viscosity() - nu).abs() < 1e-14 * nu.max(1.0));
    }
}

/// The 2D operator uses the same antisymmetric rate as `mrt3d::TrtCollision3D`.
#[test]
fn trt2d_antisymmetric_rate_matches_trt3d() {
    for &nu in &[0.01_f64, 0.1, 1.0 / 6.0, 0.4] {
        let trt2 = TrtCollision2D::new(nu);
        let trt3 = mrt3d::TrtCollision3D::new(nu);
        assert!(
            (1.0 / trt2.tau_minus - trt3.s_minus).abs() < 1e-12,
            "ν={nu}: 2D s⁻ = {}, 3D s⁻ = {}",
            1.0 / trt2.tau_minus,
            trt3.s_minus
        );
        assert!((1.0 / trt2.tau_plus - trt3.s_plus).abs() < 1e-14);
    }
}

#[test]
fn trt_with_magic_parameter_realises_the_requested_lambda() {
    for &(nu, lambda) in &[(0.05_f64, 0.25_f64), (0.2, 1.0 / 12.0), (0.5, 3.0 / 16.0)] {
        let trt = TrtCollision2D::with_magic_parameter(nu, lambda);
        assert!((trt.lambda() - lambda).abs() < 1e-14, "ν={nu}, Λ={lambda}");
    }
}

#[test]
fn trt_rejects_unusable_viscosity_or_magic_parameter() {
    for &nu in &[0.0_f64, -0.1, f64::NAN, f64::INFINITY] {
        let trt = TrtCollision2D::new(nu);
        assert!(trt.tau_minus.is_nan(), "ν={nu}: τ⁻ = {}", trt.tau_minus);
        let feq = d2q9_equilibrium_trt(1.0, 0.02, 0.01);
        let out = trt.collide(&feq, 1.0, 0.02, 0.01);
        assert!(out.iter().all(|x| x.is_nan()), "ν={nu}: {out:?}");
        assert!(TrtCollision2D::try_new(nu).is_err(), "ν={nu}");
    }
    assert!(
        TrtCollision2D::with_magic_parameter(0.1, 0.0)
            .tau_minus
            .is_nan()
    );
    assert!(
        TrtCollision2D::with_magic_parameter(0.1, -1.0)
            .tau_minus
            .is_nan()
    );
    assert!(TrtCollision2D::try_with_magic_parameter(0.1, 0.0).is_err());
    assert!(TrtCollision2D::try_with_magic_parameter(0.1, f64::NAN).is_err());
    let ok = TrtCollision2D::try_new(0.1).expect("ν = 0.1 is valid");
    assert!((ok.lambda() - 3.0 / 16.0).abs() < 1e-14);
}

/// With `τ⁻ = τ⁺` the TRT operator is BGK.
#[test]
fn trt_with_equal_times_is_bgk() {
    let tau = 0.8;
    let trt = TrtCollision2D::with_rates(tau, tau);
    let f = deterministic_perturbation(&d2q9_equilibrium_trt(1.1, 0.03, -0.02), 0.01);
    let (rho, ux, uy) = macros9(&f);
    let feq = d2q9_equilibrium_trt(rho, ux, uy);
    let out = trt.collide(&f, rho, ux, uy);
    for i in 0..9 {
        let bgk = f[i] - (f[i] - feq[i]) / tau;
        assert!((out[i] - bgk).abs() < 1e-13, "population {i}");
    }
}

#[test]
fn trt_conserves_mass_and_momentum() {
    let trt = TrtCollision2D::new(0.05);
    let f = deterministic_perturbation(&d2q9_equilibrium_trt(1.2, 0.04, 0.01), 0.01);
    let (rho, ux, uy) = macros9(&f);
    let out = trt.collide(&f, rho, ux, uy);
    let (jx_in, jy_in) = momentum9(&f);
    let (jx_out, jy_out) = momentum9(&out);
    assert!((out.iter().sum::<f64>() - rho).abs() < 1e-13);
    assert!((jx_out - jx_in).abs() < 1e-13);
    assert!((jy_out - jy_in).abs() < 1e-13);
}

/// A shear wave decays at the viscosity `ν = (τ⁺ − 1/2)/3` for low and high
/// viscosity alike — the antisymmetric rate must stay stable.
#[test]
fn trt_shear_wave_decays_with_the_set_viscosity() {
    for &(nu, steps) in &[(1.0_f64 / 6.0, 200usize), (0.01, 200)] {
        let trt = TrtCollision2D::new(nu);
        let (measured, theory) =
            shear_wave_ratio(Wave::Vertical, nu, steps, d2q9_equilibrium_trt, |f| {
                let (rho, ux, uy) = macros9(f);
                trt.collide(f, rho, ux, uy)
            });
        assert!(
            (measured / theory - 1.0).abs() < 0.01,
            "ν={nu}: measured {measured}, theory {theory}"
        );
    }
}

/// The previous parameterisation `τ⁻ = 3 / (16 τ⁺)` gives `s⁻ = 16/3 > 2` for
/// ν = 1/6; the shear wave then blows up instead of decaying.
#[test]
fn trt_with_the_old_antisymmetric_time_is_unstable() {
    let nu = 1.0 / 6.0;
    let tau_plus = 3.0 * nu + 0.5;
    let old = TrtCollision2D::with_rates(tau_plus, 3.0 / (16.0 * tau_plus));
    assert!(1.0 / old.tau_minus > 2.0);
    let (measured, _) = shear_wave_ratio(Wave::Vertical, nu, 200, d2q9_equilibrium_trt, |f| {
        let (rho, ux, uy) = macros9(f);
        old.collide(f, rho, ux, uy)
    });
    assert!(
        !(measured.is_finite() && measured.abs() <= 1.0),
        "the unstable rate unexpectedly produced a decaying wave: {measured}"
    );
}

// ---------------------------------------------------------------------------
// F2612 — cascaded (central-moment) collision
// ---------------------------------------------------------------------------

const CENTRAL_STATES: [(f64, f64, f64); 4] = [
    (1.0, 0.0, 0.0),
    (1.0, 0.05, 0.02),
    (1.3, -0.08, 0.06),
    (0.9, 0.1, -0.1),
];

#[test]
fn cascaded_rates_layout_and_accessors() {
    let nu = 0.05;
    let s_nu = 1.0 / (3.0 * nu + 0.5);
    let clbm = CascadedCollision2D::new(nu);
    assert_eq!(clbm.viscosity(), nu);
    let expected = [0.0, 0.0, 0.0, 1.0, s_nu, s_nu, 1.0, 1.0, 1.0];
    for (k, (&rate, &want)) in clbm.rates.iter().zip(expected.iter()).enumerate() {
        assert!((rate - want).abs() < 1e-14, "rates[{k}] = {rate}");
    }
    let custom = CascadedCollision2D::with_bulk_rate(nu, 1.4);
    assert_eq!(custom.rates[3], 1.4);
    assert!((custom.rates[5] - s_nu).abs() < 1e-14);
}

/// The Maxwellian-central-moment equilibrium has exactly the central moments
/// `(ρ, 0, 0, ρ/3, ρ/3, 0, 0, 0, ρ/9)` and the product form `ρ ψ(a) ψ(b)`.
#[test]
fn cascaded_equilibrium_has_maxwellian_central_moments() {
    for &(rho, ux, uy) in &CENTRAL_STATES {
        let feq = CascadedCollision2D::equilibrium(rho, ux, uy);
        let kappa = CascadedCollision2D::central_moments(&feq, ux, uy);
        let expected = [
            rho,
            0.0,
            0.0,
            rho / 3.0,
            rho / 3.0,
            0.0,
            0.0,
            0.0,
            rho / 9.0,
        ];
        for (k, (&got, &want)) in kappa.iter().zip(expected.iter()).enumerate() {
            assert!(
                (got - want).abs() < 1e-13,
                "ρ={rho}, u=({ux},{uy}): κ[{k}] = {got}, expected {want}"
            );
        }
        let psi = |c: f64, u: f64| {
            if c == 0.0 {
                2.0 / 3.0 - u * u
            } else {
                0.5 * (1.0 / 3.0 + u * u + c * u)
            }
        };
        for i in 0..9 {
            let product = rho * psi(CX9[i], ux) * psi(CY9[i], uy);
            assert!((feq[i] - product).abs() < 1e-13, "population {i}");
        }
    }
}

/// The cascaded equilibrium is a fixed point of the collision for any rates.
#[test]
fn cascaded_collision_leaves_its_equilibrium_unchanged() {
    let operators = [
        CascadedCollision2D::new(1.0 / 6.0),
        CascadedCollision2D::new(0.02),
        CascadedCollision2D {
            nu: 0.1,
            rates: [0.0, 0.0, 0.0, 1.3, 0.7, 1.9, 0.4, 0.9, 1.6],
        },
    ];
    for op in &operators {
        for &(rho, ux, uy) in &CENTRAL_STATES {
            let feq = CascadedCollision2D::equilibrium(rho, ux, uy);
            let out = op.collide(&feq, rho, ux, uy);
            for i in 0..9 {
                assert!(
                    (out[i] - feq[i]).abs() < 1e-13,
                    "ρ={rho}, u=({ux},{uy}), population {i}: {} vs {}",
                    out[i],
                    feq[i]
                );
            }
        }
    }
}

/// At rest the cascaded and the polynomial equilibrium coincide, so the
/// polynomial equilibrium is an exact fixed point there.
#[test]
fn cascaded_collision_leaves_the_rest_state_unchanged() {
    let clbm = CascadedCollision2D::new(0.05);
    for &rho in &[1.0_f64, 0.8, 1.4] {
        let feq = d2q9_equilibrium_trt(rho, 0.0, 0.0);
        let out = clbm.collide(&feq, rho, 0.0, 0.0);
        for i in 0..9 {
            assert!((out[i] - feq[i]).abs() < 1e-13, "ρ={rho}, population {i}");
        }
    }
}

/// The input of the audit finding: the previous operator reversed the
/// momentum of `f_eq(1, 0.05, 0.02)` and moved the populations by 0.44.
/// Now the momentum and every moment up to second order are preserved, and the
/// populations change only by the `O(u³)` difference of the two equilibria.
#[test]
fn cascaded_collision_of_polynomial_equilibrium_keeps_hydrodynamic_moments() {
    let (rho, ux, uy) = (1.0, 0.05, 0.02);
    let clbm = CascadedCollision2D::new(1.0 / 6.0);
    let feq = d2q9_equilibrium_trt(rho, ux, uy);
    let out = clbm.collide(&feq, rho, ux, uy);

    let (jx, jy) = momentum9(&out);
    assert!((jx - rho * ux).abs() < 1e-13, "jx = {jx}");
    assert!((jy - rho * uy).abs() < 1e-13, "jy = {jy}");

    let raw_in = CascadedCollision2D::central_moments(&feq, 0.0, 0.0);
    let raw_out = CascadedCollision2D::central_moments(&out, 0.0, 0.0);
    for k in 0..6 {
        assert!(
            (raw_in[k] - raw_out[k]).abs() < 1e-13,
            "raw moment {k}: {} -> {}",
            raw_in[k],
            raw_out[k]
        );
    }
    let deviation = (0..9)
        .map(|i| (out[i] - feq[i]).abs())
        .fold(0.0_f64, f64::max);
    assert!(deviation < 1e-4, "max |f* − f_eq| = {deviation}");
}

/// The cascaded and the polynomial equilibrium differ by `O(|u|³)`: halving
/// the velocity divides the difference by 8.
#[test]
fn cascaded_equilibrium_differs_from_polynomial_one_by_third_order_terms() {
    let gap = |scale: f64| {
        let (ux, uy) = (0.08 * scale, -0.05 * scale);
        let cascaded = CascadedCollision2D::equilibrium(1.0, ux, uy);
        let polynomial = d2q9_equilibrium_trt(1.0, ux, uy);
        (0..9)
            .map(|i| (cascaded[i] - polynomial[i]).abs())
            .fold(0.0_f64, f64::max)
    };
    let full = gap(1.0);
    let half = gap(0.5);
    assert!(full > 1e-6, "the equilibria coincide: {full}");
    let ratio = full / half;
    assert!((ratio - 8.0).abs() < 0.5, "gap ratio {ratio}, expected 8");
}

/// Mass and momentum are conserved for arbitrary input — including states
/// with negative populations (no clamping) and arguments that are not the
/// moments of `f` (no renormalisation).
#[test]
fn cascaded_collision_conserves_mass_and_momentum_without_clipping() {
    let clbm = CascadedCollision2D::new(0.02);
    let mut negative = d2q9_equilibrium_trt(1.0, 0.1, 0.0);
    negative[5] -= 0.1;
    assert!(
        negative.iter().any(|&x| x < 0.0),
        "test input must be negative somewhere"
    );
    let inputs = [
        negative,
        deterministic_perturbation(&d2q9_equilibrium_trt(1.1, 0.04, -0.03), 0.02),
    ];
    for f in &inputs {
        let (rho, ux, uy) = macros9(f);
        let (jx_in, jy_in) = momentum9(f);
        for (rho_arg, ux_arg, uy_arg) in [(rho, ux, uy), (0.7, 0.15, -0.1)] {
            let out = clbm.collide(f, rho_arg, ux_arg, uy_arg);
            let (jx_out, jy_out) = momentum9(&out);
            assert!((out.iter().sum::<f64>() - rho).abs() < 1e-13);
            assert!((jx_out - jx_in).abs() < 1e-13);
            assert!((jy_out - jy_in).abs() < 1e-13);
        }
    }
}

/// The defining property of a cascaded operator: in the frame of the fluid
/// each central moment relaxes at its own rate toward the Maxwellian value.
#[test]
fn cascaded_collision_relaxes_central_moments_at_their_rates() {
    let mut clbm = CascadedCollision2D::new(0.04);
    clbm.rates[6] = 0.5;
    clbm.rates[8] = 0.8;
    let s_bulk = clbm.rates[3];
    let s_nu = clbm.rates[5];
    let f = deterministic_perturbation(&d2q9_equilibrium_trt(1.1, 0.06, -0.04), 0.015);
    let (rho, ux, uy) = macros9(&f);
    let out = clbm.collide(&f, rho, ux, uy);
    let k_in = CascadedCollision2D::central_moments(&f, ux, uy);
    let k_out = CascadedCollision2D::central_moments(&out, ux, uy);

    // Conserved moments are untouched (κ10 and κ01 vanish in the fluid frame).
    for k in 0..3 {
        assert!((k_out[k] - k_in[k]).abs() < 1e-13, "κ[{k}]");
    }
    // Bulk mode: trace of the second-order block.
    let trace_in = k_in[3] + k_in[4];
    let trace_out = k_out[3] + k_out[4];
    let trace_expected = trace_in - s_bulk * (trace_in - 2.0 * rho / 3.0);
    assert!((trace_out - trace_expected).abs() < 1e-13);
    // Shear modes: the deviator and κ11 relax toward zero.
    let dev_in = k_in[3] - k_in[4];
    let dev_out = k_out[3] - k_out[4];
    assert!((dev_out - (1.0 - s_nu) * dev_in).abs() < 1e-13);
    assert!((k_out[5] - (1.0 - s_nu) * k_in[5]).abs() < 1e-13);
    // Higher-order moments relax toward (0, 0, ρ/9) at their own rates.
    assert!((k_out[6] - 0.5 * k_in[6]).abs() < 1e-13);
    assert!(k_out[7].abs() < 1e-13, "κ12 is relaxed at rate 1");
    assert!((k_out[8] - (0.2 * k_in[8] + 0.8 * rho / 9.0)).abs() < 1e-13);
}

/// With one rate for all non-conserved modes the cascaded operator is the BGK
/// operator toward the cascaded equilibrium.
#[test]
fn cascaded_collision_with_uniform_rates_is_bgk() {
    let omega = 0.9;
    let clbm = CascadedCollision2D {
        nu: 0.0,
        rates: [0.0, 0.0, 0.0, omega, omega, omega, omega, omega, omega],
    };
    let f = deterministic_perturbation(&d2q9_equilibrium_trt(0.95, -0.03, 0.05), 0.01);
    let (rho, ux, uy) = macros9(&f);
    let feq = CascadedCollision2D::equilibrium(rho, ux, uy);
    let out = clbm.collide(&f, rho, ux, uy);
    for i in 0..9 {
        let bgk = f[i] - omega * (f[i] - feq[i]);
        assert!((out[i] - bgk).abs() < 1e-13, "population {i}");
    }
}

/// A vertical shear wave probes `κ11` (rate `rates[5]`), a diagonal wave
/// `u ∥ (1, −1)`, `k ∥ (1, 1)` probes `κ20 − κ02` (rate `rates[4]`): both
/// decay with the set viscosity `ν`.
#[test]
fn cascaded_shear_waves_decay_with_the_set_viscosity() {
    for &(wave, nu, steps) in &[
        (Wave::Vertical, 0.05_f64, 100usize),
        (Wave::Vertical, 1.0 / 6.0, 100),
        (Wave::Diagonal, 0.05, 100),
        (Wave::Diagonal, 1.0 / 6.0, 100),
    ] {
        let clbm = CascadedCollision2D::new(nu);
        let (measured, theory) =
            shear_wave_ratio(wave, nu, steps, CascadedCollision2D::equilibrium, |f| {
                let (rho, ux, uy) = macros9(f);
                clbm.collide(f, rho, ux, uy)
            });
        assert!(
            (measured / theory - 1.0).abs() < 0.015,
            "ν={nu}: measured {measured}, theory {theory}"
        );
    }
}

/// The previous rate table relaxed `κ11` at 1 whatever `ν` was; with the
/// right `rates[5]` and a wrong `rates[4]` the diagonal wave would not decay
/// at `ν` — the test above is therefore sensitive to both shear rates.
#[test]
fn cascaded_shear_wave_detects_a_wrong_deviator_rate() {
    let nu = 0.05;
    let mut clbm = CascadedCollision2D::new(nu);
    clbm.rates[4] = 1.0;
    let (measured, theory) = shear_wave_ratio(
        Wave::Diagonal,
        nu,
        100,
        CascadedCollision2D::equilibrium,
        |f| {
            let (rho, ux, uy) = macros9(f);
            clbm.collide(f, rho, ux, uy)
        },
    );
    assert!(
        (measured / theory - 1.0).abs() > 0.2,
        "{measured} vs {theory}"
    );
}
