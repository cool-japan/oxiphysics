//! Regression tests for the upper-convected conformation-tensor terms and for
//! the polymer LBM channel grid: stretching by the velocity gradient, FENE-P
//! rest state, halfway bounce-back walls and upwind transport of the
//! conformation tensor.

use super::*;

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

/// Reference `L·A + A·Lᵀ` by explicit 2x2 matrix products.
fn reference_stretch_2d(
    a: &ConformationTensor2D,
    dudx: f64,
    dudy: f64,
    dvdx: f64,
    dvdy: f64,
) -> [[f64; 2]; 2] {
    let l = [[dudx, dudy], [dvdx, dvdy]];
    let m = [[a.xx, a.xy], [a.xy, a.yy]];
    let mut s = [[0.0; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                // (L A)_ij + (A L^T)_ij
                s[i][j] += l[i][k] * m[k][j] + m[i][k] * l[j][k];
            }
        }
    }
    s
}

// ---------------------------------------------------------------------------
// Upper-convected stretching (2D / 3D)
// ---------------------------------------------------------------------------

#[test]
fn stretch_2d_matches_matrix_products() {
    let a = ConformationTensor2D {
        xx: 2.0,
        yy: 3.0,
        xy: 0.5,
    };
    let (dudx, dudy, dvdx, dvdy) = (0.1, 0.2, 0.3, 0.4);
    let s = upper_convected_stretch_2d(&a, dudx, dudy, dvdx, dvdy);
    let r = reference_stretch_2d(&a, dudx, dudy, dvdx, dvdy);
    assert!(close(s.xx, r[0][0], 1e-14), "xx: {} vs {}", s.xx, r[0][0]);
    assert!(close(s.yy, r[1][1], 1e-14), "yy: {} vs {}", s.yy, r[1][1]);
    assert!(close(s.xy, r[0][1], 1e-14), "xy: {} vs {}", s.xy, r[0][1]);
    assert!(close(r[0][1], r[1][0], 1e-14));
}

#[test]
fn stretch_3d_matches_matrix_products() {
    let a = ConformationTensor3D {
        xx: 2.0,
        yy: 1.5,
        zz: 1.2,
        xy: 0.3,
        xz: -0.2,
        yz: 0.4,
    };
    let l: [f64; 9] = [0.10, 0.25, -0.15, 0.30, -0.05, 0.20, 0.12, -0.22, 0.07];
    let s = upper_convected_stretch_3d(&a, &l);
    let m = [[a.xx, a.xy, a.xz], [a.xy, a.yy, a.yz], [a.xz, a.yz, a.zz]];
    let mut r = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i][j] += l[3 * i + k] * m[k][j] + m[i][k] * l[3 * j + k];
            }
        }
    }
    assert!(close(s.xx, r[0][0], 1e-14));
    assert!(close(s.yy, r[1][1], 1e-14));
    assert!(close(s.zz, r[2][2], 1e-14));
    assert!(close(s.xy, r[0][1], 1e-14));
    assert!(close(s.xz, r[0][2], 1e-14));
    assert!(close(s.yz, r[1][2], 1e-14));
}

#[test]
fn oldroyd_b_simple_shear_has_the_analytic_steady_state() {
    // du/dy = gamma, lambda = 1: Axx = 1 + 2 W^2, Axy = W, Ayy = 1 with
    // W = lambda gamma.  The explicit-Euler fixed point is exact, whatever dt.
    for &w in &[0.5_f64, 1.0, 2.0] {
        let lambda = 1.0;
        let gamma = w / lambda;
        let mut a = ConformationTensor2D::identity();
        for _ in 0..20_000 {
            a = oldroyd_b_evolve_2d(&a, 0.0, gamma, 0.0, 0.0, lambda, 0.005);
        }
        assert!(
            close(a.xx, 1.0 + 2.0 * w * w, 1e-6),
            "W = {w}: Axx = {}",
            a.xx
        );
        assert!(close(a.xy, w, 1e-6), "W = {w}: Axy = {}", a.xy);
        assert!(close(a.yy, 1.0, 1e-6), "W = {w}: Ayy = {}", a.yy);

        // First normal stress difference is positive in shear flow.
        let tau = PolymerStress2D::from_conformation(&a, 0.1, lambda);
        assert!(tau.n1() > 0.0, "N1 = {}", tau.n1());
        assert!(close(tau.n1(), 0.1 / lambda * 2.0 * w * w, 1e-6));
    }
}

#[test]
fn oldroyd_b_shear_along_the_other_axis_stretches_ayy() {
    // dv/dx = gamma: the roles of x and y are swapped.
    let (lambda, gamma) = (1.0, 1.0);
    let mut a = ConformationTensor2D::identity();
    for _ in 0..20_000 {
        a = oldroyd_b_evolve_2d(&a, 0.0, 0.0, gamma, 0.0, lambda, 0.005);
    }
    assert!(close(a.yy, 3.0, 1e-6), "Ayy = {}", a.yy);
    assert!(close(a.xy, 1.0, 1e-6), "Axy = {}", a.xy);
    assert!(close(a.xx, 1.0, 1e-6), "Axx = {}", a.xx);
}

#[test]
fn oldroyd_b_planar_extension_has_the_analytic_steady_state() {
    // u = e x, v = -e y with W = lambda e = 0.2 < 0.5:
    // Axx = 1/(1 - 2W), Ayy = 1/(1 + 2W).
    let (lambda, e) = (1.0, 0.2);
    let mut a = ConformationTensor2D::identity();
    for _ in 0..20_000 {
        a = oldroyd_b_evolve_2d(&a, e, 0.0, 0.0, -e, lambda, 0.005);
    }
    assert!(close(a.xx, 1.0 / (1.0 - 0.4), 1e-6), "Axx = {}", a.xx);
    assert!(close(a.yy, 1.0 / (1.0 + 0.4), 1e-6), "Ayy = {}", a.yy);
    assert!(close(a.xy, 0.0, 1e-9));
}

/// Steady simple shear of the FENE-P model with Peterlin factor
/// `f = 1/(1 - tr(A)/b)`:
/// `Ayy = 1/f`, `Axy = W/f^2`, `Axx = (1 + 2 W^2/f^2)/f`, with `f` the root of
/// `1 - (Axx + Ayy)/b - 1/f = 0`.
fn fene_p_steady_shear(w: f64, b: f64) -> ConformationTensor2D {
    let residual = |f: f64| {
        let ayy = 1.0 / f;
        let axx = (1.0 + 2.0 * w * w / (f * f)) / f;
        1.0 - (axx + ayy) / b - 1.0 / f
    };
    let (mut lo, mut hi) = (1.0_f64, 1.0e6_f64);
    for _ in 0..300 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let f = 0.5 * (lo + hi);
    ConformationTensor2D {
        xx: (1.0 + 2.0 * w * w / (f * f)) / f,
        yy: 1.0 / f,
        xy: w / (f * f),
    }
}

#[test]
fn fene_p_simple_shear_has_the_analytic_steady_state() {
    let (lambda, b) = (1.0, 100.0);
    for &w in &[0.5_f64, 1.0, 3.0] {
        let reference = fene_p_steady_shear(w, b);
        let mut a = fene_p_equilibrium_conformation(b);
        for _ in 0..40_000 {
            a = fene_p_evolve_2d(&a, 0.0, w / lambda, 0.0, 0.0, lambda, b, 0.002);
        }
        assert!(
            close(a.xx, reference.xx, 1e-6),
            "W = {w}: Axx = {} vs {}",
            a.xx,
            reference.xx
        );
        assert!(
            close(a.xy, reference.xy, 1e-6),
            "W = {w}: Axy = {} vs {}",
            a.xy,
            reference.xy
        );
        assert!(
            close(a.yy, reference.yy, 1e-6),
            "W = {w}: Ayy = {} vs {}",
            a.yy,
            reference.yy
        );
        // The chain stays below full extension.
        assert!(a.trace() < b);
        // Positive N1 for the FENE-P stress, too.
        let tau = fene_p_stress_2d(&a, 0.1, lambda, b);
        assert!(tau.n1() > 0.0, "N1 = {}", tau.n1());
    }
}

#[test]
fn fene_p_rest_state_is_stationary_and_stress_free() {
    let b = 50.0;
    let a = fene_p_equilibrium_conformation(b);
    assert!(close(a.xx, b / (b + 2.0), 1e-15));
    assert!(close(a.yy, b / (b + 2.0), 1e-15));
    let tau = fene_p_stress_2d(&a, 0.3, 2.0, b);
    assert!(tau.txx.abs() < 1e-14 && tau.tyy.abs() < 1e-14 && tau.txy.abs() < 1e-14);
    let next = fene_p_evolve_2d(&a, 0.0, 0.0, 0.0, 0.0, 2.0, b, 1.0);
    assert!(close(next.xx, a.xx, 1e-14));
    assert!(close(next.yy, a.yy, 1e-14));
    assert!(close(next.xy, 0.0, 1e-14));
    // b -> infinity is the Oldroyd-B identity.
    let id = fene_p_equilibrium_conformation(f64::INFINITY);
    assert_eq!((id.xx, id.yy, id.xy), (1.0, 1.0, 0.0));
}

#[test]
fn oldroyd_b_3d_simple_shear_has_the_analytic_steady_state_in_every_plane() {
    let (lambda, w) = (1.0, 1.0);
    let gamma = w / lambda;
    for i in 0..3 {
        for j in 0..3 {
            if i == j {
                continue;
            }
            // u_i = gamma x_j  ->  L_ij = gamma
            let mut l = [0.0_f64; 9];
            l[3 * i + j] = gamma;
            let mut a = ConformationTensor3D::identity();
            for _ in 0..4_000 {
                a = oldroyd_b_evolve_3d(&a, &l, lambda, 0.01);
            }
            let full = [[a.xx, a.xy, a.xz], [a.xy, a.yy, a.yz], [a.xz, a.yz, a.zz]];
            for r in 0..3 {
                for c in 0..3 {
                    let expected = if r == c {
                        if r == i { 1.0 + 2.0 * w * w } else { 1.0 }
                    } else if (r == i && c == j) || (r == j && c == i) {
                        w
                    } else {
                        0.0
                    };
                    assert!(
                        close(full[r][c], expected, 1e-6),
                        "shear u_{i}(x_{j}): A[{r}][{c}] = {} expected {expected}",
                        full[r][c]
                    );
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// PolymerLbmGrid2D: walls, gradients, transport
// ---------------------------------------------------------------------------

/// BGK frequency for which the halfway bounce-back wall location is exact for
/// Poiseuille flow, `(tau - 1/2)^2 = 3/16`.
fn magic_omega() -> f64 {
    1.0 / (0.5 + 3.0_f64.sqrt() / 4.0)
}

fn fluid_mass(grid: &PolymerLbmGrid2D) -> f64 {
    grid.fluid_rows()
        .flat_map(|iy| (0..grid.nx).map(move |ix| (ix, iy)))
        .map(|(ix, iy)| grid.cells[grid.idx(ix, iy)].f.iter().sum::<f64>())
        .sum()
}

#[test]
fn grid_streaming_reflects_populations_at_the_walls() {
    let mut grid = PolymerLbmGrid2D::new(4, 6, 1.0, 0.0, 10.0, f64::INFINITY);
    // Give every population of every cell a unique value.
    for (c, cell) in grid.cells.iter_mut().enumerate() {
        for (q, f) in cell.f.iter_mut().enumerate() {
            *f = 1.0 + c as f64 * 0.01 + q as f64 * 1e-4;
        }
    }
    let before: Vec<[f64; 9]> = grid.cells.iter().map(|c| c.f).collect();
    grid.stream();
    let nx = grid.nx;
    let top = grid.ny - 2;
    for ix in 0..nx {
        // Upper fluid row: populations moving up (2, 5, 6) hit the wall and come
        // back as 4, 7, 8 from the same node.
        let dst = grid.idx(ix, top);
        assert_eq!(grid.cells[dst].f[4], before[dst][2]);
        assert_eq!(grid.cells[dst].f[7], before[dst][5]);
        assert_eq!(grid.cells[dst].f[8], before[dst][6]);
        // Lower fluid row: 4, 7, 8 come back as 2, 5, 6.
        let dst = grid.idx(ix, 1);
        assert_eq!(grid.cells[dst].f[2], before[dst][4]);
        assert_eq!(grid.cells[dst].f[5], before[dst][7]);
        assert_eq!(grid.cells[dst].f[6], before[dst][8]);
    }
    // The wall rows are left untouched: nothing wraps around in y.
    for ix in 0..nx {
        let lo = grid.idx(ix, 0);
        let hi = grid.idx(ix, grid.ny - 1);
        assert_eq!(grid.cells[lo].f, before[lo]);
        assert_eq!(grid.cells[hi].f, before[hi]);
    }
}

#[test]
fn grid_conserves_mass_with_active_polymer_stress() {
    let mut grid = PolymerLbmGrid2D::new(6, 12, 1.2, 0.05, 20.0, f64::INFINITY);
    grid.set_poiseuille_profile(0.02);
    let m0 = fluid_mass(&grid);
    for _ in 0..200 {
        grid.step();
    }
    let m1 = fluid_mass(&grid);
    assert!(
        (m1 - m0).abs() / m0 < 1e-12,
        "mass drifted from {m0} to {m1}"
    );
}

#[test]
fn grid_at_rest_stays_at_rest_for_both_models() {
    for &b in &[f64::INFINITY, 40.0] {
        let mut grid = PolymerLbmGrid2D::new(5, 8, 1.0, 0.1, 10.0, b);
        let a0 = fene_p_equilibrium_conformation(b);
        for _ in 0..50 {
            grid.step();
        }
        for cell in &grid.cells {
            assert!(close(cell.conformation.xx, a0.xx, 1e-12), "b = {b}");
            assert!(close(cell.conformation.yy, a0.yy, 1e-12), "b = {b}");
            assert!(close(cell.conformation.xy, 0.0, 1e-12), "b = {b}");
            assert!(cell.stress.txx.abs() < 1e-12 && cell.stress.txy.abs() < 1e-12);
            assert!(cell.ux.abs() < 1e-12 && cell.uy.abs() < 1e-12);
        }
        assert!(close(grid.max_n1(), 0.0, 1e-12));
    }
}

#[test]
fn grid_newtonian_channel_reproduces_the_poiseuille_profile_between_halfway_walls() {
    let (nx, ny) = (4, 16);
    let omega = magic_omega();
    let nu = (1.0 / omega - 0.5) / 3.0;
    let force = 1.0e-5;
    // No polymer: pure BGK channel flow.
    let mut grid = PolymerLbmGrid2D::new(nx, ny, omega, 0.0, 10.0, f64::INFINITY);
    let ex = D2Q9_EX;
    let w = D2Q9_WEIGHTS;
    for _ in 0..6_000 {
        for iy in grid.fluid_rows() {
            for ix in 0..nx {
                let idx = grid.idx(ix, iy);
                for i in 0..9 {
                    grid.cells[idx].f[i] += 3.0 * w[i] * force * ex[i];
                }
            }
        }
        grid.step();
    }
    let h = (ny - 2) as f64;
    let u_max = force * h * h / (8.0 * nu);
    for iy in grid.fluid_rows() {
        let y = iy as f64 - 0.5;
        let analytic = force / (2.0 * nu) * y * (h - y);
        let idx = grid.idx(1, iy);
        let measured = grid.cells[idx].ux;
        assert!(
            (measured - analytic).abs() < 0.03 * u_max,
            "row {iy}: u = {measured}, analytic {analytic}"
        );
        // The profile is uniform along the periodic direction and has no
        // cross-channel velocity.
        let other = grid.cells[grid.idx(3, iy)].ux;
        assert!((measured - other).abs() < 1e-12);
        assert!(grid.cells[idx].uy.abs() < 1e-12);
    }
    let centre = grid.cells[grid.idx(1, ny / 2)].ux;
    assert!(
        (centre - u_max).abs() < 0.03 * u_max,
        "centre velocity {centre} vs {u_max}"
    );
    // Walls are at rest.
    assert_eq!(grid.cells[grid.idx(1, 0)].ux, 0.0);
    assert_eq!(grid.cells[grid.idx(1, ny - 1)].ux, 0.0);
}

#[test]
fn grid_velocity_gradient_is_exact_for_the_parabolic_profile() {
    let mut grid = PolymerLbmGrid2D::new(4, 14, 1.0, 0.0, 10.0, f64::INFINITY);
    let u_max = 0.02;
    grid.set_poiseuille_profile(u_max);
    let h = (grid.ny - 2) as f64;
    for iy in grid.fluid_rows() {
        let y = iy as f64 - 0.5;
        let expected = 4.0 * u_max * (1.0 - 2.0 * y / h) / h;
        let [dudx, dudy, dvdx, dvdy] = grid.velocity_gradient(2, iy);
        assert!(
            close(dudy, expected, 1e-12),
            "row {iy}: du/dy = {dudy}, expected {expected}"
        );
        assert!(dudx.abs() < 1e-15 && dvdx.abs() < 1e-15 && dvdy.abs() < 1e-15);
    }
}

#[test]
fn grid_upwind_advection_translates_the_conformation_at_the_flow_speed() {
    let (nx, ny) = (16, 7);
    // lambda = infinity: no relaxation, pure transport.
    let mut grid = PolymerLbmGrid2D::new(nx, ny, 1.0, 0.0, f64::INFINITY, f64::INFINITY);
    let u = 0.25;
    for cell in grid.cells.iter_mut() {
        cell.ux = u;
    }
    // A conformation bump on the middle row (whose neighbouring rows move at the
    // same speed, so the velocity gradient there vanishes).
    let iy = 3;
    for ix in 5..8 {
        let idx = grid.idx(ix, iy);
        grid.cells[idx].conformation.xx = 2.0;
    }
    let moment = |g: &PolymerLbmGrid2D| -> (f64, f64) {
        let mut mass = 0.0;
        let mut first = 0.0;
        for ix in 0..nx {
            let e = g.cells[g.idx(ix, iy)].conformation.xx - 1.0;
            mass += e;
            first += ix as f64 * e;
        }
        (mass, first / mass)
    };
    let (m0, x0) = moment(&grid);
    let steps = 8;
    for _ in 0..steps {
        grid.advance_conformation(1.0);
    }
    let (m1, x1) = moment(&grid);
    assert!(
        close(m1, m0, 1e-12),
        "transported amount changed: {m0} -> {m1}"
    );
    assert!(
        close(x1 - x0, u * steps as f64, 1e-9),
        "centre of mass moved by {} instead of {}",
        x1 - x0,
        u * steps as f64
    );
    // Upwind transport is monotone: the bump never overshoots.
    for ix in 0..nx {
        let a = grid.cells[grid.idx(ix, iy)].conformation.xx;
        assert!((1.0 - 1e-12..=2.0 + 1e-12).contains(&a), "A_xx = {a}");
    }
}

#[test]
fn polymer_channel_with_weissenberg_one_stretches_the_polymer_near_the_walls() {
    // omega = 1 -> nu_s = 1/6; the body force sustains the initial profile with
    // u_max = 0.01 (u_max = F H^2 / (8 nu), H = 14).  The wall shear rate is
    // about 2.65e-3, so lambda = 350 gives Wi ~ 0.93 at the first fluid row.
    let cfg = PolymerChannelConfig {
        nx: 4,
        ny: 16,
        nu_s: 1.0 / 6.0,
        eta_p: 0.02,
        lambda: 350.0,
        fene_b: f64::INFINITY,
        body_force: 6.8e-5,
        n_steps: 3_000,
    };
    let grid = run_polymer_channel_grid(&cfg);
    let mean_trace = run_polymer_channel(&cfg);
    assert!(
        mean_trace > 2.0 + 1e-3,
        "the polymer is not stretched: mean tr(A) = {mean_trace}"
    );
    assert!(
        close(mean_trace, grid.mean_polymer_trace(), 1e-9),
        "run_polymer_channel and run_polymer_channel_grid disagree"
    );

    let wall_row = 1;
    let mid_row = cfg.ny / 2;
    let a_wall = grid.cells[grid.idx(0, wall_row)].conformation;
    let a_mid = grid.cells[grid.idx(0, mid_row)].conformation;
    assert!(
        a_wall.xx > a_mid.xx + 0.5,
        "A_xx near the wall ({}) must exceed the centreline value ({})",
        a_wall.xx,
        a_mid.xx
    );
    assert!(a_mid.xx < 1.1, "centreline A_xx = {}", a_mid.xx);
    // Positive first normal stress difference in the sheared region.
    assert!(grid.cells[grid.idx(0, wall_row)].stress.n1() > 0.0);
    // The shear component has opposite signs at the two walls and the
    // normal component is mirror symmetric.
    let upper_row = cfg.ny - 2;
    let a_upper = grid.cells[grid.idx(0, upper_row)].conformation;
    assert!(a_wall.xy > 0.0 && a_upper.xy < 0.0);
    assert!(close(a_wall.xx, a_upper.xx, 1e-6 * a_wall.xx));
    assert!(close(a_wall.xy, -a_upper.xy, 1e-6 * a_wall.xy.abs()));
}

#[test]
fn fene_p_channel_limits_the_chain_extension() {
    let b = 30.0;
    let cfg = PolymerChannelConfig {
        nx: 4,
        ny: 16,
        nu_s: 1.0 / 6.0,
        eta_p: 0.02,
        lambda: 700.0,
        fene_b: b,
        body_force: 6.8e-5,
        n_steps: 3_000,
    };
    let grid = run_polymer_channel_grid(&cfg);
    for iy in grid.fluid_rows() {
        let a = grid.cells[grid.idx(0, iy)].conformation;
        assert!(
            a.trace().is_finite() && a.trace() < b,
            "tr(A) = {}",
            a.trace()
        );
    }
    let a_wall = grid.cells[grid.idx(0, 1)].conformation;
    let a_mid = grid.cells[grid.idx(0, cfg.ny / 2)].conformation;
    assert!(
        a_wall.xx > a_mid.xx + 0.5,
        "FENE-P stretching near the wall: {} vs centre {}",
        a_wall.xx,
        a_mid.xx
    );
    // The centreline is close to the FENE-P rest state b/(b+2) (its small
    // residual shear rate stretches it only slightly).
    assert!(close(a_mid.xx, b / (b + 2.0), 0.15));
}

#[test]
fn polymer_channel_conformation_matches_the_steady_shear_solution_row_by_row() {
    // Long enough (many relaxation times) for the x-uniform channel to reach a
    // quasi-steady state: every fluid row then carries the steady simple-shear
    // conformation of its own local shear rate.
    let (lambda, b) = (350.0, 40.0);
    for &fene_b in &[f64::INFINITY, b] {
        let cfg = PolymerChannelConfig {
            nx: 4,
            ny: 16,
            nu_s: 1.0 / 6.0,
            eta_p: 0.02,
            lambda,
            fene_b,
            body_force: 6.8e-5,
            n_steps: 8_000,
        };
        let grid = run_polymer_channel_grid(&cfg);
        for iy in grid.fluid_rows() {
            let [_, dudy, _, _] = grid.velocity_gradient(0, iy);
            let w = lambda * dudy.abs();
            let steady = if fene_b.is_finite() {
                fene_p_steady_shear(w, fene_b)
            } else {
                ConformationTensor2D {
                    xx: 1.0 + 2.0 * w * w,
                    yy: 1.0,
                    xy: w,
                }
            };
            let a = grid.cells[grid.idx(0, iy)].conformation;
            let sign = dudy.signum();
            assert!(
                close(a.xx, steady.xx, 0.02 * steady.xx),
                "b = {fene_b}, row {iy}: Axx = {} vs steady {}",
                a.xx,
                steady.xx
            );
            assert!(
                close(a.xy, sign * steady.xy, 0.02 * steady.xy.abs() + 1e-4),
                "b = {fene_b}, row {iy}: Axy = {} vs steady {}",
                a.xy,
                sign * steady.xy
            );
            assert!(
                close(a.yy, steady.yy, 0.02 * steady.yy),
                "b = {fene_b}, row {iy}: Ayy = {} vs steady {}",
                a.yy,
                steady.yy
            );
            if fene_b.is_finite() {
                assert!(a.trace() < fene_b);
            }
        }
    }
}

#[test]
fn grids_without_fluid_rows_do_not_panic() {
    for ny in 0..3 {
        let mut grid = PolymerLbmGrid2D::new(4, ny, 1.0, 0.1, 10.0, f64::INFINITY);
        grid.set_poiseuille_profile(0.01);
        grid.step();
        assert!(grid.fluid_rows().is_empty());
        assert_eq!(grid.step, 1);
        assert!(grid.mean_polymer_trace().is_finite());
    }
}
