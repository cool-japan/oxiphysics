//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
use super::types::OcOptimizer;
use super::types::{BesoParams, SimpParams, TopologyGrid};

/// Fraction of 2×2 patches exhibiting a checkerboard pattern.
///
/// A patch at (col, row) is checkerboard when the four corners alternate
/// between ρ > 0.5 and ρ ≤ 0.5 in a checker fashion.
pub fn checkerboard_measure(grid: &TopologyGrid) -> f64 {
    if grid.nx < 2 || grid.ny < 2 {
        return 0.0;
    }
    let mut checker = 0u32;
    let mut total = 0u32;
    for row in 0..(grid.ny - 1) {
        for col in 0..(grid.nx - 1) {
            let idx = |r: usize, c: usize| r * grid.nx + c;
            let a = grid.elements[idx(row, col)].rho > 0.5;
            let b = grid.elements[idx(row, col + 1)].rho > 0.5;
            let c = grid.elements[idx(row + 1, col)].rho > 0.5;
            let d = grid.elements[idx(row + 1, col + 1)].rho > 0.5;
            if (a && !b && !c && d) || (!a && b && c && !d) {
                checker += 1;
            }
            total += 1;
        }
    }
    if total == 0 {
        0.0
    } else {
        checker as f64 / total as f64
    }
}
/// Fraction of elements in the gray zone 0.1 < ρ < 0.9.
pub fn gray_scale_measure(grid: &TopologyGrid) -> f64 {
    let gray = grid
        .elements
        .iter()
        .filter(|e| e.rho > 0.1 && e.rho < 0.9)
        .count();
    gray as f64 / grid.elements.len() as f64
}
/// Heaviside projection: smooth step from 0 → 1 around `eta` with sharpness `beta`.
///
/// Formula: (tanh(β·η) + tanh(β·(ρ − η))) / (tanh(β·η) + tanh(β·(1 − η)))
pub fn heaviside_projection(rho: f64, beta: f64, eta: f64) -> f64 {
    let num = beta.mul_add(eta, 0.0).tanh() + (beta * (rho - eta)).tanh();
    let den = (beta * eta).tanh() + (beta * (1.0 - eta)).tanh();
    if den.abs() < 1e-30 { rho } else { num / den }
}
/// Derivative of Heaviside projection with respect to rho.
///
/// ∂H/∂ρ = β · (1 - tanh²(β·(ρ - η))) / (tanh(β·η) + tanh(β·(1 − η)))
pub fn heaviside_derivative(rho: f64, beta: f64, eta: f64) -> f64 {
    let t = (beta * (rho - eta)).tanh();
    let den = (beta * eta).tanh() + (beta * (1.0 - eta)).tanh();
    if den.abs() < 1e-30 {
        1.0
    } else {
        beta * (1.0 - t * t) / den
    }
}
/// Sum of |ρ_i − ρ_j| for all horizontally or vertically adjacent element pairs.
pub fn perimeter_measure(grid: &TopologyGrid) -> f64 {
    let mut perim = 0.0_f64;
    for row in 0..grid.ny {
        for col in 0..grid.nx {
            let i = row * grid.nx + col;
            if col + 1 < grid.nx {
                let j = row * grid.nx + (col + 1);
                perim += (grid.elements[i].rho - grid.elements[j].rho).abs();
            }
            if row + 1 < grid.ny {
                let j = (row + 1) * grid.nx + col;
                perim += (grid.elements[i].rho - grid.elements[j].rho).abs();
            }
        }
    }
    perim
}
/// Discreteness measure: 4 * ρ * (1 - ρ) averaged over all elements.
///
/// Equals 0 for a fully discrete (0/1) design and 1 for ρ = 0.5 everywhere.
pub fn discreteness_measure(grid: &TopologyGrid) -> f64 {
    let n = grid.elements.len() as f64;
    if n < 1.0 {
        return 0.0;
    }
    let sum: f64 = grid
        .elements
        .iter()
        .map(|e| 4.0 * e.rho * (1.0 - e.rho))
        .sum();
    sum / n
}
/// Compute the compliance objective function value.
///
/// c = Σ E(ρ_i) · sensitivity_i (proxy when FEA is not available).
pub fn compliance_objective(grid: &TopologyGrid, params: &SimpParams) -> f64 {
    grid.elements
        .iter()
        .map(|e| {
            let ee = params.emin + e.rho.powf(params.penalty) * (params.e0 - params.emin);
            ee * e.sensitivity.abs()
        })
        .sum()
}
/// Volume constraint value: V(ρ) / V_target - 1.
///
/// Negative means under budget, positive means over budget.
pub fn volume_constraint(grid: &TopologyGrid, target_vf: f64) -> f64 {
    grid.volume_fraction() / target_vf - 1.0
}
/// Continuation strategy for penalty parameter.
///
/// Returns updated penalty: p_new = min(p_current + increment, p_max).
pub fn penalty_continuation(current_penalty: f64, increment: f64, max_penalty: f64) -> f64 {
    (current_penalty + increment).min(max_penalty)
}
/// Continuation strategy for Heaviside beta parameter.
///
/// Returns updated beta: β_new = min(β_current * factor, β_max).
pub fn beta_continuation(current_beta: f64, factor: f64, max_beta: f64) -> f64 {
    (current_beta * factor).min(max_beta)
}
/// Compute a distance-weighted sensitivity filter with normalization.
///
/// This is an improved version of the hat-function filter that uses a
/// Gaussian kernel instead of a tent function for smoother sensitivity fields.
pub fn filter_sensitivities_gaussian(
    grid: &TopologyGrid,
    sensitivities: &[f64],
    radius: f64,
) -> Vec<f64> {
    let n = grid.elements.len();
    let sigma = radius / 3.0;
    let mut filtered = vec![0.0f64; n];
    for (i, filt_i) in filtered.iter_mut().enumerate() {
        let ei = &grid.elements[i];
        let mut num = 0.0f64;
        let mut den = 0.0f64;
        for (j, ej) in grid.elements.iter().enumerate() {
            let dist_sq = (ei.x0 - ej.x0).powi(2) + (ei.y0 - ej.y0).powi(2);
            if dist_sq < (radius * radius) {
                let w = (-dist_sq / (2.0 * sigma * sigma)).exp();
                num += w * ej.rho * sensitivities[j];
                den += w * ej.rho;
            }
        }
        *filt_i = if den.abs() > 1e-30 {
            num / den
        } else {
            sensitivities[i]
        };
    }
    filtered
}
/// Compute a PDE-based (Helmholtz) smoothed sensitivity field.
///
/// Solves the PDE: `(1 - r²∇²) φ = ψ` approximately via a single
/// explicit Euler step on a uniform grid.  This is the PDE filter
/// of Lazarov & Sigmund (2016).
pub fn filter_sensitivities_pde(
    grid: &TopologyGrid,
    sensitivities: &[f64],
    r_pde: f64,
) -> Vec<f64> {
    let nx = grid.nx;
    let ny = grid.ny;
    let _n = nx * ny;
    let dx = if nx > 0 { grid.elements[0].dx } else { 1.0 };
    let dy = if ny > 0 { grid.elements[0].dy } else { 1.0 };
    let r2 = r_pde * r_pde;
    let mut phi = sensitivities.to_vec();
    let mut phi_new = phi.clone();
    for row in 0..ny {
        for col in 0..nx {
            let idx = row * nx + col;
            let lapl = {
                let east = if col + 1 < nx {
                    phi[row * nx + col + 1]
                } else {
                    phi[idx]
                };
                let west = if col > 0 {
                    phi[row * nx + col - 1]
                } else {
                    phi[idx]
                };
                let north = if row + 1 < ny {
                    phi[(row + 1) * nx + col]
                } else {
                    phi[idx]
                };
                let south = if row > 0 {
                    phi[(row - 1) * nx + col]
                } else {
                    phi[idx]
                };
                (east - 2.0 * phi[idx] + west) / (dx * dx)
                    + (north - 2.0 * phi[idx] + south) / (dy * dy)
            };
            phi_new[idx] = (sensitivities[idx] + r2 * lapl)
                / (1.0 + 2.0 * r2 / (dx * dx) + 2.0 * r2 / (dy * dy));
        }
    }
    phi = phi_new;
    for v in phi.iter_mut() {
        if !v.is_finite() {
            *v = 0.0;
        }
    }
    phi
}
/// BESO binary density update.
///
/// Sorts elements by their sensitivity number and switches the bottom `er`
/// fraction from solid to void (or vice-versa) to move toward the volume
/// target.
pub fn beso_update(densities: &[f64], sensitivities: &[f64], params: &BesoParams) -> Vec<f64> {
    let n = densities.len();
    let target_solid = (params.target_vf * n as f64).round() as usize;
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        sensitivities[a]
            .partial_cmp(&sensitivities[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut new_dens = densities.to_vec();
    let current_solid = densities.iter().filter(|&&d| d > 0.5).count();
    if current_solid > target_solid {
        let remove = (current_solid.saturating_sub(target_solid))
            .min((params.er * n as f64).ceil() as usize);
        let mut removed = 0;
        for &idx in &order {
            if removed >= remove {
                break;
            }
            if densities[idx] > 0.5 {
                new_dens[idx] = params.x_min;
                removed += 1;
            }
        }
    } else if current_solid < target_solid {
        let add = (target_solid.saturating_sub(current_solid))
            .min((params.er * n as f64).ceil() as usize);
        let mut added = 0;
        for &idx in order.iter().rev() {
            if added >= add {
                break;
            }
            if densities[idx] <= 0.5 {
                new_dens[idx] = params.x_max;
                added += 1;
            }
        }
    }
    new_dens
}
/// Run multiple BESO iterations.
///
/// Returns the final density vector and the volume fraction history.
pub fn beso_optimize(
    initial_densities: &[f64],
    sensitivities_fn: &dyn Fn(&[f64]) -> Vec<f64>,
    params: &BesoParams,
    max_iter: usize,
) -> (Vec<f64>, Vec<f64>) {
    let mut densities = initial_densities.to_vec();
    let mut vf_history = Vec::with_capacity(max_iter);
    for _ in 0..max_iter {
        let sens = sensitivities_fn(&densities);
        densities = beso_update(&densities, &sens, params);
        let vf = densities.iter().filter(|&&d| d > 0.5).count() as f64 / densities.len() as f64;
        vf_history.push(vf);
    }
    (densities, vf_history)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology_opt::*;
    #[test]
    fn test_default_steel_e0() {
        let p = SimpParams::default_steel();
        assert!((p.e0 - 210e9).abs() < 1.0);
    }
    #[test]
    fn test_default_steel_emin() {
        let p = SimpParams::default_steel();
        assert!((p.emin - 1e-9 * p.e0).abs() < 1e-3);
    }
    #[test]
    fn test_default_steel_penalty() {
        let p = SimpParams::default_steel();
        assert_eq!(p.penalty, 3.0);
    }
    #[test]
    fn test_default_steel_vf() {
        let p = SimpParams::default_steel();
        assert_eq!(p.volume_fraction, 0.5);
    }
    #[test]
    fn test_default_steel_move_limit() {
        let p = SimpParams::default_steel();
        assert_eq!(p.move_limit, 0.2);
    }
    #[test]
    fn test_default_aluminum() {
        let p = SimpParams::default_aluminum();
        assert!((p.e0 - 70e9).abs() < 1.0);
        assert_eq!(p.volume_fraction, 0.4);
    }
    fn make_elem(rho: f64) -> SimpElement {
        SimpElement {
            rho,
            stiffness: 0.0,
            sensitivity: 0.0,
            x0: 0.5,
            y0: 0.5,
            dx: 1.0,
            dy: 1.0,
        }
    }
    #[test]
    fn test_effective_modulus_solid() {
        let p = SimpParams::default_steel();
        let e = make_elem(1.0);
        assert!((e.effective_modulus(&p) - p.e0).abs() < 1.0);
    }
    #[test]
    fn test_effective_modulus_void() {
        let p = SimpParams::default_steel();
        let e = make_elem(0.0);
        assert!((e.effective_modulus(&p) - p.emin).abs() < 1e-6 * p.e0);
    }
    #[test]
    fn test_effective_modulus_half() {
        let p = SimpParams::default_steel();
        let e = make_elem(0.5);
        let expected = p.emin + 0.5_f64.powf(p.penalty) * (p.e0 - p.emin);
        assert!((e.effective_modulus(&p) - expected).abs() < 1.0);
    }
    #[test]
    fn test_compliance_sensitivity_solid() {
        let p = SimpParams::default_steel();
        let e = make_elem(1.0);
        let s = e.compliance_sensitivity(1.0, &p);
        let expected = -p.penalty * (p.e0 - p.emin);
        assert!((s - expected).abs() < 1.0);
    }
    #[test]
    fn test_compliance_sensitivity_negative() {
        let p = SimpParams::default_steel();
        let e = make_elem(0.5);
        assert!(e.compliance_sensitivity(1.0, &p) <= 0.0);
    }
    #[test]
    fn test_ramp_modulus_solid() {
        let p = SimpParams::default_steel();
        let e = make_elem(1.0);
        assert!((e.ramp_modulus(&p) - p.e0).abs() < 1.0);
    }
    #[test]
    fn test_ramp_modulus_void() {
        let p = SimpParams::default_steel();
        let e = make_elem(0.0);
        assert!((e.ramp_modulus(&p) - p.emin).abs() < 1e-6 * p.e0);
    }
    #[test]
    fn test_ramp_sensitivity_negative() {
        let p = SimpParams::default_steel();
        let e = make_elem(0.5);
        assert!(e.ramp_sensitivity(1.0, &p) <= 0.0);
    }
    #[test]
    fn test_element_area() {
        let e = make_elem(0.5);
        assert!((e.area() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_grid_uniform_size() {
        let g = TopologyGrid::new_uniform(4, 3, 0.5, 210e9, 1.0, 1.0);
        assert_eq!(g.elements.len(), 12);
        assert_eq!(g.nx, 4);
        assert_eq!(g.ny, 3);
    }
    #[test]
    fn test_grid_volume_fraction() {
        let g = TopologyGrid::new_uniform(5, 5, 0.4, 210e9, 1.0, 1.0);
        assert!((g.volume_fraction() - 0.4).abs() < 1e-12);
    }
    #[test]
    fn test_grid_volume_fraction_weighted() {
        let g = TopologyGrid::new_uniform(5, 5, 0.4, 210e9, 1.0, 1.0);
        assert!((g.volume_fraction_weighted() - 0.4).abs() < 1e-12);
    }
    #[test]
    fn test_grid_centroids() {
        let g = TopologyGrid::new_uniform(2, 2, 0.5, 210e9, 1.0, 1.0);
        assert!((g.elements[0].x0 - 0.5).abs() < 1e-12);
        assert!((g.elements[0].y0 - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_filter_sensitivities_length() {
        let mut g = TopologyGrid::new_uniform(4, 4, 0.5, 210e9, 1.0, 1.0);
        for e in g.elements.iter_mut() {
            e.sensitivity = -1.0;
        }
        let f = g.filter_sensitivities(1.5);
        assert_eq!(f.len(), 16);
    }
    #[test]
    fn test_filter_densities_uniform() {
        let g = TopologyGrid::new_uniform(4, 4, 0.5, 210e9, 1.0, 1.0);
        let filtered = g.filter_densities(1.5);
        for &d in &filtered {
            assert!(
                (d - 0.5).abs() < 1e-10,
                "Uniform density should be unchanged by filter"
            );
        }
    }
    #[test]
    fn test_oc_update_volume_constraint() {
        let params = SimpParams {
            e0: 1.0,
            emin: 1e-9,
            penalty: 3.0,
            volume_fraction: 0.4,
            filter_radius: 1.5,
            move_limit: 0.2,
        };
        let mut g = TopologyGrid::new_uniform(6, 6, 0.5, 1.0, 1.0, 1.0);
        let sens = vec![-1.0; 36];
        g.update_densities_oc(&sens, &params);
        let vf = g.volume_fraction();
        assert!((vf - 0.4).abs() < 0.05, "vf={vf}");
    }
    #[test]
    fn test_mma_update_volume_constraint() {
        let params = SimpParams {
            e0: 1.0,
            emin: 1e-9,
            penalty: 3.0,
            volume_fraction: 0.4,
            filter_radius: 1.5,
            move_limit: 0.2,
        };
        let mut g = TopologyGrid::new_uniform(6, 6, 0.5, 1.0, 1.0, 1.0);
        let sens = vec![-1.0; 36];
        let lo = vec![0.0; 36];
        let hi = vec![1.0; 36];
        g.update_densities_mma(&sens, &params, &lo, &hi);
        let vf = g.volume_fraction();
        assert!((vf - 0.4).abs() < 0.05, "vf={vf}");
    }
    #[test]
    fn test_element_index_roundtrip() {
        let g = TopologyGrid::new_uniform(5, 4, 0.5, 1.0, 1.0, 1.0);
        for row in 0..4 {
            for col in 0..5 {
                let idx = g.element_index(row, col);
                let (r, c) = g.row_col(idx);
                assert_eq!(r, row);
                assert_eq!(c, col);
            }
        }
    }
    #[test]
    fn test_density_gradient_uniform() {
        let g = TopologyGrid::new_uniform(4, 4, 0.5, 1.0, 1.0, 1.0);
        let grad = g.density_gradient_magnitude();
        for &v in &grad {
            assert!(v.abs() < 1e-12, "Uniform field should have zero gradient");
        }
    }
    #[test]
    fn test_oc_optimizer_new() {
        let p = SimpParams::default_steel();
        let opt = OcOptimizer::new(4, 4, p);
        assert_eq!(opt.iteration, 0);
        assert_eq!(opt.grid.elements.len(), 16);
    }
    #[test]
    fn test_oc_optimizer_step_increments_iter() {
        let p = SimpParams::default_steel();
        let mut opt = OcOptimizer::new(4, 4, p);
        opt.step();
        assert_eq!(opt.iteration, 1);
    }
    #[test]
    fn test_oc_optimizer_records_history() {
        let p = SimpParams::default_steel();
        let mut opt = OcOptimizer::new(4, 4, p);
        opt.step();
        opt.step();
        assert_eq!(opt.compliance_history.len(), 2);
        assert_eq!(opt.volume_history.len(), 2);
    }
    #[test]
    fn test_compliance_change_none_initially() {
        let p = SimpParams::default_steel();
        let opt = OcOptimizer::new(4, 4, p);
        assert!(opt.compliance_change().is_none());
    }
    #[test]
    fn test_density_field_length() {
        let p = SimpParams::default_steel();
        let opt = OcOptimizer::new(3, 5, p);
        assert_eq!(opt.density_field().len(), 15);
    }
    #[test]
    fn test_checkerboard_full() {
        let mut g = TopologyGrid::new_uniform(2, 2, 0.5, 1.0, 1.0, 1.0);
        g.elements[0].rho = 1.0;
        g.elements[1].rho = 0.0;
        g.elements[2].rho = 0.0;
        g.elements[3].rho = 1.0;
        assert!((checkerboard_measure(&g) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_checkerboard_none() {
        let g = TopologyGrid::new_uniform(4, 4, 1.0, 1.0, 1.0, 1.0);
        assert_eq!(checkerboard_measure(&g), 0.0);
    }
    #[test]
    fn test_gray_scale_all_gray() {
        let g = TopologyGrid::new_uniform(4, 4, 0.5, 1.0, 1.0, 1.0);
        assert!((gray_scale_measure(&g) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_gray_scale_none() {
        let g = TopologyGrid::new_uniform(4, 4, 1.0, 1.0, 1.0, 1.0);
        assert_eq!(gray_scale_measure(&g), 0.0);
    }
    #[test]
    fn test_heaviside_mid() {
        let h = heaviside_projection(0.5, 10.0, 0.5);
        assert!((h - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_heaviside_zero() {
        let h = heaviside_projection(0.0, 50.0, 0.5);
        assert!(h < 0.05, "h={h}");
    }
    #[test]
    fn test_heaviside_one() {
        let h = heaviside_projection(1.0, 50.0, 0.5);
        assert!(h > 0.95, "h={h}");
    }
    #[test]
    fn test_heaviside_derivative_positive() {
        let dh = heaviside_derivative(0.5, 10.0, 0.5);
        assert!(dh > 0.0, "Derivative should be positive at rho=eta");
    }
    #[test]
    fn test_heaviside_derivative_sharp() {
        let dh_small = heaviside_derivative(0.5, 1.0, 0.5);
        let dh_large = heaviside_derivative(0.5, 50.0, 0.5);
        assert!(
            dh_large > dh_small,
            "Larger beta should give sharper derivative"
        );
    }
    #[test]
    fn test_perimeter_uniform() {
        let g = TopologyGrid::new_uniform(4, 4, 0.7, 1.0, 1.0, 1.0);
        assert!(perimeter_measure(&g).abs() < 1e-12);
    }
    #[test]
    fn test_perimeter_nonzero() {
        let mut g = TopologyGrid::new_uniform(2, 1, 0.5, 1.0, 1.0, 1.0);
        g.elements[0].rho = 0.0;
        g.elements[1].rho = 1.0;
        assert!((perimeter_measure(&g) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_discreteness_measure_half() {
        let g = TopologyGrid::new_uniform(4, 4, 0.5, 1.0, 1.0, 1.0);
        let dm = discreteness_measure(&g);
        assert!((dm - 1.0).abs() < 1e-12, "All ρ=0.5 → discreteness=1");
    }
    #[test]
    fn test_discreteness_measure_binary() {
        let mut g = TopologyGrid::new_uniform(4, 4, 1.0, 1.0, 1.0, 1.0);
        for i in 0..8 {
            g.elements[i].rho = 0.0;
        }
        let dm = discreteness_measure(&g);
        assert!(dm.abs() < 1e-12, "Binary design → discreteness=0");
    }
    #[test]
    fn test_volume_constraint_at_target() {
        let g = TopologyGrid::new_uniform(4, 4, 0.5, 1.0, 1.0, 1.0);
        let vc = volume_constraint(&g, 0.5);
        assert!(vc.abs() < 1e-12, "At target → constraint = 0");
    }
    #[test]
    fn test_volume_constraint_over_budget() {
        let g = TopologyGrid::new_uniform(4, 4, 0.8, 1.0, 1.0, 1.0);
        let vc = volume_constraint(&g, 0.5);
        assert!(vc > 0.0, "Over budget → positive constraint");
    }
    #[test]
    fn test_penalty_continuation() {
        let p = penalty_continuation(3.0, 0.5, 5.0);
        assert!((p - 3.5).abs() < 1e-12);
        let p2 = penalty_continuation(4.8, 0.5, 5.0);
        assert!((p2 - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_beta_continuation() {
        let b = beta_continuation(1.0, 2.0, 64.0);
        assert!((b - 2.0).abs() < 1e-12);
        let b2 = beta_continuation(50.0, 2.0, 64.0);
        assert!((b2 - 64.0).abs() < 1e-12);
    }
    #[test]
    fn test_compliance_objective_positive() {
        let p = SimpParams::default_steel();
        let mut g = TopologyGrid::new_uniform(4, 4, 0.5, p.e0, 1.0, 1.0);
        for e in g.elements.iter_mut() {
            e.sensitivity = -1.0;
        }
        let c = compliance_objective(&g, &p);
        assert!(c > 0.0, "Compliance should be positive");
    }
    #[test]
    fn test_projection_filter_apply_dimensions() {
        let g = TopologyGrid::new_uniform(4, 4, 0.5, 1.0, 1.0, 1.0);
        let pf = ProjectionFilter::new(1.5, 10.0, 0.5);
        let (rho_bar, rho_phys) = pf.apply(&g);
        assert_eq!(rho_bar.len(), 16);
        assert_eq!(rho_phys.len(), 16);
    }
    #[test]
    fn test_projection_filter_uniform_density() {
        let g = TopologyGrid::new_uniform(4, 4, 0.5, 1.0, 1.0, 1.0);
        let pf = ProjectionFilter::new(1.5, 10.0, 0.5);
        let (_, rho_phys) = pf.apply(&g);
        for &rp in &rho_phys {
            assert!((rp - 0.5).abs() < 0.05, "rho_phys = {rp}");
        }
    }
    #[test]
    fn test_projection_filter_chain_sensitivity_length() {
        let g = TopologyGrid::new_uniform(4, 4, 0.5, 1.0, 1.0, 1.0);
        let pf = ProjectionFilter::new(1.5, 10.0, 0.5);
        let (rho_bar, _) = pf.apply(&g);
        let dc_phys = vec![-1.0; 16];
        let dc = pf.chain_sensitivity(&g, &dc_phys, &rho_bar);
        assert_eq!(dc.len(), 16);
    }
    #[test]
    fn test_projection_filter_sharp_projection() {
        let g = TopologyGrid::new_uniform(4, 4, 0.8, 1.0, 1.0, 1.0);
        let pf = ProjectionFilter::new(1.5, 100.0, 0.5);
        let (_, rho_phys) = pf.apply(&g);
        for &rp in &rho_phys {
            assert!(rp > 0.9, "Sharp projection: rho_phys={rp}");
        }
    }
    #[test]
    fn test_gaussian_sensitivity_filter_length() {
        let g = TopologyGrid::new_uniform(4, 4, 0.5, 1.0, 1.0, 1.0);
        let sens = vec![-1.0; 16];
        let filtered = filter_sensitivities_gaussian(&g, &sens, 1.5);
        assert_eq!(filtered.len(), 16);
    }
    #[test]
    fn test_gaussian_sensitivity_filter_uniform() {
        let g = TopologyGrid::new_uniform(4, 4, 0.5, 1.0, 1.0, 1.0);
        let sens = vec![-2.0; 16];
        let filtered = filter_sensitivities_gaussian(&g, &sens, 1.5);
        for &f in &filtered {
            assert!(
                (f - (-2.0)).abs() < 1e-6,
                "Gaussian filter changed uniform: {f}"
            );
        }
    }
    #[test]
    fn test_pde_sensitivity_filter_length() {
        let g = TopologyGrid::new_uniform(4, 4, 0.5, 1.0, 1.0, 1.0);
        let sens = vec![-1.0; 16];
        let filtered = filter_sensitivities_pde(&g, &sens, 0.5);
        assert_eq!(filtered.len(), 16);
    }
    #[test]
    fn test_pde_sensitivity_filter_uniform() {
        let g = TopologyGrid::new_uniform(4, 4, 0.5, 1.0, 1.0, 1.0);
        let sens = vec![-1.5; 16];
        let filtered = filter_sensitivities_pde(&g, &sens, 0.5);
        for &f in &filtered {
            assert!(f.is_finite(), "PDE filter produced non-finite: {f}");
        }
    }
    #[test]
    fn test_beso_params_default() {
        let p = BesoParams::default_structural();
        assert!((p.target_vf - 0.5).abs() < 1e-12);
        assert!((p.er - 0.02).abs() < 1e-12);
    }
    #[test]
    fn test_beso_update_length_preserved() {
        let dens = vec![1.0; 10];
        let sens = vec![-1.0; 10];
        let params = BesoParams::default_structural();
        let new_dens = beso_update(&dens, &sens, &params);
        assert_eq!(new_dens.len(), 10);
    }
    #[test]
    fn test_beso_update_removes_material() {
        let dens = vec![1.0; 20];
        let sens: Vec<f64> = (0..20).map(|i| -(i as f64)).collect();
        let params = BesoParams::default_structural();
        let new_dens = beso_update(&dens, &sens, &params);
        let solid = new_dens.iter().filter(|&&d| d > 0.5).count();
        assert!(
            solid < 20,
            "Should have removed some elements; solid={solid}"
        );
    }
    #[test]
    fn test_beso_update_adds_material() {
        let dens = vec![1e-3; 20];
        let sens: Vec<f64> = (0..20).map(|i| -(i as f64)).collect();
        let mut params = BesoParams::default_structural();
        params.target_vf = 0.5;
        let new_dens = beso_update(&dens, &sens, &params);
        let solid = new_dens.iter().filter(|&&d| d > 0.5).count();
        assert!(solid > 0, "Should have added some elements; solid={solid}");
    }
    #[test]
    fn test_beso_optimize_converges_toward_vf() {
        let n = 20;
        let dens_init = vec![1.0; n];
        let params = BesoParams {
            target_vf: 0.5,
            er: 0.05,
            penalty: 3.0,
            x_max: 1.0,
            x_min: 1e-3,
        };
        let (final_dens, vf_history) = beso_optimize(
            &dens_init,
            &|d: &[f64]| d.iter().map(|&x| -x).collect(),
            &params,
            15,
        );
        assert_eq!(final_dens.len(), n);
        assert!(!vf_history.is_empty());
        let final_vf = vf_history.last().unwrap();
        assert!((0.0..=1.0).contains(final_vf));
    }
    #[test]
    fn test_level_set_new_uniform_volume_fraction() {
        let ls = LevelSetField::new_uniform(4, 4, 1.0, 1.0, 0.6);
        let vf = ls.volume_fraction();
        assert!((vf - 1.0).abs() < 1e-12, "vf = {vf}");
    }
    #[test]
    fn test_level_set_half_domain() {
        let mut phi = vec![-1.0; 16];
        for p in phi.iter_mut().take(8) {
            *p = 1.0;
        }
        let ls = LevelSetField::from_sdf(4, 4, 1.0, 1.0, phi);
        let vf = ls.volume_fraction();
        assert!((vf - 0.5).abs() < 1e-12, "vf = {vf}");
    }
    #[test]
    fn test_level_set_heaviside_smooth() {
        let h = LevelSetField::heaviside_smooth(1.0, 0.1);
        assert!((h - 1.0).abs() < 1e-10);
        let h2 = LevelSetField::heaviside_smooth(-1.0, 0.1);
        assert!(h2.abs() < 1e-10);
        let h3 = LevelSetField::heaviside_smooth(0.0, 0.1);
        assert!((h3 - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_level_set_dirac_smooth() {
        let d = LevelSetField::dirac_smooth(1.0, 0.1);
        assert!(d.abs() < 1e-12, "dirac outside bandwidth = {d}");
        let d0 = LevelSetField::dirac_smooth(0.0, 0.1);
        assert!(d0 > 0.0, "dirac at zero should be positive: {d0}");
    }
    #[test]
    fn test_level_set_material_indicator() {
        let phi = vec![1.0, -1.0, 0.0, 0.5];
        let ls = LevelSetField::from_sdf(4, 1, 1.0, 1.0, phi);
        let ind = ls.material_indicator(0.1);
        assert!((ind[0] - 1.0).abs() < 1e-10);
        assert!(ind[1].abs() < 1e-10);
        assert!((ind[2] - 0.5).abs() < 0.1);
    }
    #[test]
    fn test_level_set_advect_changes_phi() {
        let mut ls = LevelSetField::from_sdf(4, 1, 1.0, 1.0, vec![1.0, 0.5, -0.5, -1.0]);
        let vel = vec![1.0, 1.0, -1.0, -1.0];
        let phi_before = ls.phi.clone();
        ls.advect(&vel, 0.01);
        let changed = ls
            .phi
            .iter()
            .zip(phi_before.iter())
            .any(|(a, b)| (a - b).abs() > 1e-15);
        assert!(changed, "Advection should change phi");
    }
    #[test]
    fn test_level_set_optimizer_step_increments_iter() {
        let mut opt = LevelSetOptimizer::new(4, 4, 1.0, 1.0, 0.5, 0.01, 0.1);
        let sens = vec![-1.0; 16];
        opt.step(&sens);
        assert_eq!(opt.iteration, 1);
    }
    #[test]
    fn test_level_set_optimizer_material_indicator_length() {
        let opt = LevelSetOptimizer::new(4, 4, 1.0, 1.0, 0.5, 0.01, 0.1);
        let ind = opt.material_indicator();
        assert_eq!(ind.len(), 16);
    }
    #[test]
    fn test_multi_load_new_valid() {
        let p = SimpParams::default_steel();
        let opt = MultiLoadOptimizer::new(3, vec![0.5, 0.3, 0.2], p, 16);
        assert_eq!(opt.weights.len(), 3);
        assert_eq!(opt.accumulated_sensitivity.len(), 16);
    }
    #[test]
    fn test_multi_load_accumulate_sensitivity() {
        let p = SimpParams::default_steel();
        let mut opt = MultiLoadOptimizer::new(2, vec![0.5, 0.5], p, 4);
        opt.accumulate_sensitivity(0, &[-1.0, -2.0, -3.0, -4.0]);
        opt.accumulate_sensitivity(1, &[-2.0, -4.0, -6.0, -8.0]);
        assert!((opt.accumulated_sensitivity[0] - (-1.5)).abs() < 1e-12);
        assert!((opt.accumulated_sensitivity[1] - (-3.0)).abs() < 1e-12);
    }
    #[test]
    fn test_multi_load_reset_sensitivity() {
        let p = SimpParams::default_steel();
        let mut opt = MultiLoadOptimizer::new(2, vec![0.5, 0.5], p, 4);
        opt.accumulate_sensitivity(0, &[-1.0, -2.0, -3.0, -4.0]);
        opt.reset_sensitivity();
        for &s in &opt.accumulated_sensitivity {
            assert_eq!(s, 0.0);
        }
    }
    #[test]
    fn test_multi_load_total_compliance() {
        let p = SimpParams::default_steel();
        let opt = MultiLoadOptimizer::new(2, vec![0.6, 0.4], p, 4);
        let c = opt.total_compliance(&[10.0, 20.0]);
        assert!((c - 14.0).abs() < 1e-12, "Total compliance = {c}");
    }
    #[test]
    fn test_multi_load_update_densities() {
        let p = SimpParams {
            e0: 1.0,
            emin: 1e-9,
            penalty: 3.0,
            volume_fraction: 0.4,
            filter_radius: 1.5,
            move_limit: 0.2,
        };
        let mut opt = MultiLoadOptimizer::new(2, vec![0.5, 0.5], p, 36);
        opt.accumulate_sensitivity(0, &vec![-1.0; 36]);
        opt.accumulate_sensitivity(1, &vec![-1.0; 36]);
        let mut g = TopologyGrid::new_uniform(6, 6, 0.5, 1.0, 1.0, 1.0);
        opt.update_densities_oc(&mut g);
        let vf = g.volume_fraction();
        assert!((vf - 0.4).abs() < 0.1, "Multi-load OC update: vf={vf}");
    }
}
/// Filter sensitivities using a weighted average over elements within radius r_min.
///
/// `mesh` is a list of element centroid (x, y) coordinates.
pub fn filter_sensitivities_mesh(
    sensitivities: &[f64],
    mesh: &[(f64, f64)],
    r_min: f64,
) -> Vec<f64> {
    assert_eq!(sensitivities.len(), mesh.len());
    let n = sensitivities.len();
    let mut filtered = vec![0.0_f64; n];
    for i in 0..n {
        let (xi, yi) = mesh[i];
        let mut weight_sum = 0.0_f64;
        let mut sens_sum = 0.0_f64;
        for j in 0..n {
            let (xj, yj) = mesh[j];
            let dx = xi - xj;
            let dy = yi - yj;
            let dist = (dx * dx + dy * dy).sqrt();
            let w = (r_min - dist).max(0.0);
            weight_sum += w;
            sens_sum += w * sensitivities[j];
        }
        filtered[i] = if weight_sum > 0.0 {
            sens_sum / weight_sum
        } else {
            sensitivities[i]
        };
    }
    filtered
}
/// Compliance objective C = u^T * K * u computed from sparse triplets.
///
/// `K_global` is a list of `(row, col, value)` sparse matrix triplets.
/// `u` is the displacement vector; `_f` is accepted for API consistency.
pub fn compliance_objective_sparse(k_global: &[(usize, usize, f64)], u: &[f64], _f: &[f64]) -> f64 {
    let n = u.len();
    let mut c = 0.0_f64;
    for &(row, col, val) in k_global {
        if row < n && col < n {
            c += u[row] * val * u[col];
        }
    }
    c
}
/// Total (weighted) volume of the design domain.
///
/// V = sum_i (densities\[i\] * element_volumes\[i\])
pub fn total_volume(densities: &[f64], element_volumes: &[f64]) -> f64 {
    assert_eq!(densities.len(), element_volumes.len());
    densities
        .iter()
        .zip(element_volumes.iter())
        .map(|(&r, &v)| r * v)
        .sum()
}
/// Compute the gradient of the volume constraint V(ρ) = Σ vₑ·ρₑ / V_total
/// with respect to density.
///
/// Returns `∂V/∂ρₑ = vₑ / V_total` for each element.
pub fn volume_constraint_gradient(element_volumes: &[f64], total_volume: f64) -> Vec<f64> {
    assert!(total_volume > 0.0);
    element_volumes.iter().map(|&v| v / total_volume).collect()
}
/// Compute the augmented Lagrangian penalty contribution for a volume constraint.
///
/// `g(ρ) = V(ρ) - V_target` is the constraint violation.
/// Returns `λ * g + (ρ_aug / 2) * g²`.
pub fn augmented_lagrangian_volume_penalty(
    current_vf: f64,
    target_vf: f64,
    lambda: f64,
    rho_aug: f64,
) -> f64 {
    let g = current_vf - target_vf;
    lambda * g + 0.5 * rho_aug * g * g
}
/// Normalize sensitivities by their maximum absolute value.
///
/// Prevents numerical issues when sensitivities span many orders of magnitude.
/// Returns a zero vector if all sensitivities are zero.
pub fn normalize_sensitivities(sensitivities: &[f64]) -> Vec<f64> {
    let max_abs = sensitivities
        .iter()
        .map(|&s| s.abs())
        .fold(0.0_f64, f64::max);
    if max_abs < 1e-30 {
        return vec![0.0; sensitivities.len()];
    }
    sensitivities.iter().map(|&s| s / max_abs).collect()
}
/// Scale sensitivities so their sum equals `n_elements`.
///
/// Useful when the OC update requires sensitivities summing to the
/// number of design variables.
pub fn scale_sensitivities_to_sum(sensitivities: &[f64]) -> Vec<f64> {
    let n = sensitivities.len();
    if n == 0 {
        return vec![];
    }
    let sum: f64 = sensitivities.iter().sum::<f64>().abs();
    if sum < 1e-30 {
        return vec![1.0; n];
    }
    sensitivities.iter().map(|&s| s * n as f64 / sum).collect()
}
#[cfg(test)]
mod tests_topopt_new {
    use super::*;
    use crate::topology_opt::*;
    #[test]
    fn test_simp_model_at_rho1_gives_e0() {
        let m = SimpModel {
            penalty: 3.0,
            e_min: 1e-9,
            e_0: 200e9,
        };
        let e = m.young_modulus(1.0);
        assert!(
            (e - 200e9).abs() < 1.0,
            "E(rho=1) should equal E_0, got {e:.3e}"
        );
    }
    #[test]
    fn test_simp_model_at_rho0_gives_emin() {
        let e_min = 1e-3_f64;
        let m = SimpModel {
            penalty: 3.0,
            e_min,
            e_0: 200e9,
        };
        let e = m.young_modulus(0.0);
        assert!(
            (e - e_min).abs() < 1e-10,
            "E(rho=0) should equal E_min, got {e:.3e}"
        );
    }
    #[test]
    fn test_oc_solver_moves_toward_target() {
        let n = 100;
        let initial = 0.8_f64;
        let target = 0.5_f64;
        let mut solver = OcSolver::new(n, initial, 0.2);
        let sens = vec![1.0_f64; n];
        solver.update(&sens, target);
        let vf = solver.densities.iter().sum::<f64>() / n as f64;
        assert!(
            vf < initial + 1e-6,
            "OC update should reduce density from {initial} toward {target}, got {vf}"
        );
        assert!(vf > 0.0, "density must remain positive, got {vf}");
    }
    #[test]
    fn test_filter_sensitivities_mesh_smooth() {
        let n = 9usize;
        let mesh: Vec<(f64, f64)> = (0..n).map(|i| (i as f64, 0.0)).collect();
        let sens = vec![2.0_f64; n];
        let filtered = filter_sensitivities_mesh(&sens, &mesh, 2.5);
        for (i, &f) in filtered.iter().enumerate() {
            assert!((f - 2.0).abs() < 1e-10, "uniform input: filtered[{i}]={f}");
        }
    }
    #[test]
    fn test_total_volume_basic() {
        let densities = vec![0.5, 0.5, 1.0];
        let volumes = vec![1.0, 2.0, 1.0];
        let v = total_volume(&densities, &volumes);
        assert!((v - 2.5).abs() < 1e-12, "total volume = {v}");
    }
    #[test]
    fn test_compliance_objective_sparse_basic() {
        let k: Vec<(usize, usize, f64)> = vec![(0, 0, 2.0), (1, 1, 2.0)];
        let u = vec![1.0, 1.0];
        let f = vec![2.0, 2.0];
        let c = compliance_objective_sparse(&k, &u, &f);
        assert!((c - 4.0).abs() < 1e-12, "compliance = {c}");
    }
    #[test]
    fn test_ramp_at_rho1_gives_e0() {
        let m = RampModel {
            q: 3.0,
            e_min: 1e-6,
            e_0: 200e9,
        };
        assert!((m.young_modulus(1.0) - 200e9).abs() < 1.0);
    }
    #[test]
    fn test_ramp_at_rho0_gives_emin() {
        let m = RampModel {
            q: 3.0,
            e_min: 1e-6,
            e_0: 200e9,
        };
        assert!((m.young_modulus(0.0) - 1e-6).abs() < 1e-15);
    }
    #[test]
    fn test_ramp_derivative_positive() {
        let m = RampModel {
            q: 3.0,
            e_min: 1e-6,
            e_0: 200e9,
        };
        let d = m.derivative(0.5);
        assert!(d > 0.0, "RAMP derivative should be positive: {d}");
    }
    #[test]
    fn test_sensitivity_history_average() {
        let mut hist = SensitivityHistory::new(4, 2);
        hist.add_iteration(&[1.0, 2.0, 3.0, 4.0]);
        hist.add_iteration(&[3.0, 4.0, 5.0, 6.0]);
        let avg = hist.average_sensitivity();
        assert!((avg[0] - 2.0).abs() < 1e-12, "avg[0]={}", avg[0]);
        assert!((avg[1] - 3.0).abs() < 1e-12, "avg[1]={}", avg[1]);
    }
    #[test]
    fn test_sensitivity_history_single_iter() {
        let mut hist = SensitivityHistory::new(3, 5);
        let sens = [0.5, 1.0, 1.5];
        hist.add_iteration(&sens);
        let avg = hist.average_sensitivity();
        for (i, &v) in avg.iter().enumerate() {
            assert!((v - sens[i]).abs() < 1e-12, "avg[{i}]={v}");
        }
    }
    #[test]
    fn test_volume_gradient_uniform_elements() {
        let vols = vec![2.0, 2.0, 2.0];
        let total = vols.iter().sum::<f64>();
        let grad = volume_constraint_gradient(&vols, total);
        for &g in &grad {
            assert!((g - 1.0 / 3.0).abs() < 1e-12, "grad={g}");
        }
    }
    #[test]
    fn test_convergence_monitor_detects_convergence() {
        let mut mon = TopOptConvergenceMonitor::new(5, 1e-3);
        for _ in 0..5 {
            mon.record_change(1e-5);
        }
        assert!(mon.is_converged(), "should detect convergence");
    }
    #[test]
    fn test_convergence_monitor_no_convergence() {
        let mut mon = TopOptConvergenceMonitor::new(5, 1e-3);
        mon.record_change(0.5);
        assert!(
            !mon.is_converged(),
            "should not converge with one large change"
        );
    }
    #[test]
    fn test_augmented_lagrangian_update() {
        let mut al = AugmentedLagrangianVolumeConstraint::new(0.5, 1.0);
        let vf = 0.6;
        al.update_multiplier(vf);
        assert!(
            al.lambda != 0.0,
            "multiplier should be nonzero after update"
        );
    }
    #[test]
    fn test_density_projection_threshold() {
        let proj = DensityProjection {
            beta: 8.0,
            eta: 0.5,
        };
        let v = proj.project(0.5);
        assert!(
            v > 0.0 && v < 1.0,
            "projected density should be in (0,1): {v}"
        );
    }
    #[test]
    fn test_density_projection_monotone() {
        let proj = DensityProjection {
            beta: 4.0,
            eta: 0.5,
        };
        let v1 = proj.project(0.3);
        let v2 = proj.project(0.7);
        assert!(v1 < v2, "projection should be monotone: {v1} < {v2}");
    }
}
