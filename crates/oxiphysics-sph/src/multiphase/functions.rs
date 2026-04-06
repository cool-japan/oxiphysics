//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{DropletInfo, PhaseMorphologyAnalyzer};

/// Cubic-spline kernel value W(r, h).
///
/// Uses the 3-D normalisation factor `3/(2π h³)`.
pub fn cubic_spline_w(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 3.0 / (2.0 * std::f64::consts::PI * h * h * h);
    if q < 1.0 {
        sigma * (2.0 / 3.0 - q * q + 0.5 * q * q * q)
    } else if q < 2.0 {
        sigma * (1.0 / 6.0) * (2.0 - q).powi(3)
    } else {
        0.0
    }
}
/// Scalar derivative dW/dr of the cubic-spline kernel.
pub fn cubic_spline_dw(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 3.0 / (2.0 * std::f64::consts::PI * h * h * h);
    if q < 1.0 {
        sigma / h * (-2.0 * q + 1.5 * q * q)
    } else if q < 2.0 {
        sigma / h * (-0.5 * (2.0 - q).powi(2))
    } else {
        0.0
    }
}
/// Kernel gradient vector ∇W(x_i − x_j, h).
pub fn kernel_grad(dx: [f64; 3], r: f64, h: f64) -> [f64; 3] {
    if r < 1e-14 {
        return [0.0; 3];
    }
    let dw_dr = cubic_spline_dw(r, h);
    let scale = dw_dr / r;
    [scale * dx[0], scale * dx[1], scale * dx[2]]
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::multiphase::types::*;
    #[test]
    fn test_multiphase_pressure_correction() {
        let rho_i = 1000.0_f64;
        let rho_j = 1.2_f64;
        let p_i = 1000.0_f64;
        let p_j = 100.0_f64;
        let result = DensityRatioMultiphase::pressure_correction(rho_i, rho_j, p_i, p_j);
        let expected = p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j);
        assert!(
            (result - expected).abs() < 1e-14,
            "pressure_correction mismatch: {result} vs {expected}"
        );
    }
    #[test]
    fn test_density_ratio_same_phase() {
        let rho = 1000.0_f64;
        let p = 5000.0_f64;
        let result = DensityRatioMultiphase::pressure_correction(rho, rho, p, p);
        let expected = 2.0 * p / (rho * rho);
        assert!(
            (result - expected).abs() < 1e-14,
            "Same-phase correction should be 2p/ρ²: {result} vs {expected}"
        );
    }
    #[test]
    fn test_viscous_force_direction() {
        let mut mp = DensityRatioMultiphase::new();
        let i = mp.add_particle(Phase::Liquid, 1000.0, 1e-3);
        let j = mp.add_particle(Phase::Gas, 1.2, 1.8e-5);
        let r_ij = [0.05_f64, 0.0, 0.0];
        let v_ij = [1.0_f64, 0.0, 0.0];
        let rho_i = 1000.0_f64;
        let rho_j = 1.2_f64;
        let h = 0.1_f64;
        let f = mp.viscous_force(i, j, r_ij, v_ij, rho_i, rho_j, h);
        assert!(
            f[0] < 0.0,
            "Viscous force should oppose relative motion: f[0] = {}",
            f[0]
        );
        assert!(f[1].abs() < 1e-30, "f[1] should be zero: {}", f[1]);
        assert!(f[2].abs() < 1e-30, "f[2] should be zero: {}", f[2]);
    }
    #[test]
    fn test_harmonic_mean_viscosity() {
        let mu_i = 1e-3_f64;
        let mu_j = 1e-5_f64;
        let expected_mu = 2.0 * mu_i * mu_j / (mu_i + mu_j);
        let mut mp = DensityRatioMultiphase::new();
        let i = mp.add_particle(Phase::Liquid, 1000.0, mu_i);
        let j = mp.add_particle(Phase::Gas, 1.2, mu_j);
        let mu_i2 = 2.0 * mu_i;
        let mu_j2 = 2.0 * mu_j;
        let expected_mu2 = 2.0 * mu_i2 * mu_j2 / (mu_i2 + mu_j2);
        let mut mp2 = DensityRatioMultiphase::new();
        let i2 = mp2.add_particle(Phase::Liquid, 1000.0, mu_i2);
        let j2 = mp2.add_particle(Phase::Gas, 1.2, mu_j2);
        let r_ij = [0.05_f64, 0.0, 0.0];
        let v_ij = [1.0_f64, 0.0, 0.0];
        let rho_i = 1000.0_f64;
        let rho_j = 1.2_f64;
        let h = 0.1_f64;
        let f1 = mp.viscous_force(i, j, r_ij, v_ij, rho_i, rho_j, h);
        let f2 = mp2.viscous_force(i2, j2, r_ij, v_ij, rho_i, rho_j, h);
        let ratio_force = f2[0] / f1[0];
        let ratio_mu = expected_mu2 / expected_mu;
        assert!(
            (ratio_force - ratio_mu).abs() < 1e-10,
            "Force ratio {ratio_force} should equal μ ratio {ratio_mu}"
        );
    }
    #[test]
    fn test_phase_change_temperature() {
        let pc = PhaseChange::new(373.15, 273.15, 2_260_000.0, 334_000.0);
        assert_eq!(pc.phase_from_temperature(250.0), Phase::Solid);
        assert_eq!(pc.phase_from_temperature(273.15), Phase::Liquid);
        assert_eq!(pc.phase_from_temperature(300.0), Phase::Liquid);
        assert_eq!(pc.phase_from_temperature(373.15), Phase::Gas);
        assert_eq!(pc.phase_from_temperature(400.0), Phase::Gas);
    }
    #[test]
    fn test_latent_heat_values() {
        let pc = PhaseChange::new(373.15, 273.15, 2_260_000.0, 334_000.0);
        assert!(
            pc.latent_heat_fusion < pc.latent_heat_vaporization,
            "Fusion latent heat ({}) should be less than vaporisation ({})",
            pc.latent_heat_fusion,
            pc.latent_heat_vaporization,
        );
        assert!(pc.latent_heat(Phase::Liquid, Phase::Gas) > 0.0);
        assert!(pc.latent_heat(Phase::Solid, Phase::Liquid) > 0.0);
        assert!(pc.latent_heat(Phase::Gas, Phase::Liquid) < 0.0);
        assert!(pc.latent_heat(Phase::Liquid, Phase::Solid) < 0.0);
    }
    #[test]
    fn test_mixture_density() {
        let mut mix = MixtureSph::new(1);
        mix.set_fraction(0, 0.5, 0.5, 0.0);
        let rho_l = 1000.0_f64;
        let rho_g = 1.2_f64;
        let rho_s = 2700.0_f64;
        let rho_eff = mix.effective_density(0, rho_l, rho_g, rho_s);
        let expected = (rho_l + rho_g) / 2.0;
        assert!(
            (rho_eff - expected).abs() < 1e-10,
            "effective density = {rho_eff}, expected {expected}"
        );
    }
    #[test]
    fn test_mixture_normalize() {
        let mut mix = MixtureSph::new(1);
        mix.set_fraction(0, 3.0, 1.0, 2.0);
        let [fl, fg, fs] = mix.phase_fractions[0];
        let sum = fl + fg + fs;
        assert!(
            (sum - 1.0).abs() < 1e-14,
            "Phase fractions should sum to 1.0, got {sum}"
        );
        assert!((fl - 0.5).abs() < 1e-14, "f_liquid should be 0.5, got {fl}");
        assert!(
            (fg - 1.0 / 6.0).abs() < 1e-14,
            "f_gas should be 1/6, got {fg}"
        );
        assert!(
            (fs - 1.0 / 3.0).abs() < 1e-14,
            "f_solid should be 1/3, got {fs}"
        );
    }
    #[test]
    fn test_interphase_tension_set_get() {
        let mut t = InterphaseTension::new(3);
        t.set(0, 1, 0.07);
        t.set(0, 2, 0.05);
        assert!(
            (t.get(0, 1) - 0.07).abs() < 1e-14,
            "sigma(0,1) should be 0.07"
        );
        assert!(
            (t.get(1, 0) - 0.07).abs() < 1e-14,
            "sigma(1,0) should be 0.07 (symmetric)"
        );
        assert!(
            (t.get(0, 2) - 0.05).abs() < 1e-14,
            "sigma(0,2) should be 0.05"
        );
        assert!(
            (t.get(1, 2) - 0.0).abs() < 1e-14,
            "sigma(1,2) should be 0 (unset)"
        );
    }
    #[test]
    fn test_cubic_spline_w_at_origin() {
        let w = cubic_spline_w(0.0, 1.0);
        assert!(w > 0.0, "W(0, h) must be positive, got {w}");
    }
    #[test]
    fn test_cubic_spline_w_outside_support() {
        let w = cubic_spline_w(2.0, 1.0);
        assert_eq!(w, 0.0, "W(2h, h) should be 0 (outside support), got {w}");
    }
    #[test]
    fn test_multiphase_system_add_phase_particle() {
        let mut sys = MultiphaseSystem::new(0.1, [0.0, -9.81, 0.0]);
        let water = FluidPhase {
            id: 0,
            density_ref: 1000.0,
            dynamic_viscosity: 1e-3,
            surface_tension: 0.07,
            color: [0.0, 0.0, 1.0],
        };
        let oil = FluidPhase {
            id: 1,
            density_ref: 800.0,
            dynamic_viscosity: 1e-2,
            surface_tension: 0.03,
            color: [1.0, 0.8, 0.0],
        };
        let id_w = sys.add_phase(water);
        let id_o = sys.add_phase(oil);
        assert_eq!(id_w, 0);
        assert_eq!(id_o, 1);
        assert_eq!(sys.phases.len(), 2);
        let pi = sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], id_w);
        let pj = sys.add_particle([0.05, 0.0, 0.0], [0.0; 3], id_o);
        assert_eq!(pi, 0);
        assert_eq!(pj, 1);
        assert_eq!(sys.particles.len(), 2);
        assert_eq!(sys.particles[0].phase_id, id_w);
        assert_eq!(sys.particles[1].phase_id, id_o);
    }
    #[test]
    fn test_compute_density_single_particle() {
        let mut sys = MultiphaseSystem::new(0.1, [0.0, -9.81, 0.0]);
        let phase = FluidPhase {
            id: 0,
            density_ref: 1000.0,
            dynamic_viscosity: 1e-3,
            surface_tension: 0.07,
            color: [0.0, 0.0, 1.0],
        };
        sys.add_phase(phase);
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 0);
        sys.particles[0].mass = 1.0;
        sys.compute_density();
        let h = 0.1;
        let expected = cubic_spline_w(0.0, h);
        let got = sys.particles[0].density;
        assert!(
            (got - expected).abs() < 1e-12,
            "single-particle density: expected {expected}, got {got}"
        );
    }
    #[test]
    fn test_phase_separation_index_same_phase() {
        let mut sys = MultiphaseSystem::new(1.0, [0.0, -9.81, 0.0]);
        let phase = FluidPhase {
            id: 0,
            density_ref: 1000.0,
            dynamic_viscosity: 1e-3,
            surface_tension: 0.07,
            color: [0.0, 0.0, 1.0],
        };
        sys.add_phase(phase);
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 0);
        sys.add_particle([0.5, 0.0, 0.0], [0.0; 3], 0);
        sys.particles[0].density = 1000.0;
        sys.particles[1].density = 1000.0;
        let idx = sys.phase_separation_index();
        assert!(
            (idx - 1.0).abs() < 1e-10,
            "Same-phase particles should give separation index 1.0, got {idx}"
        );
    }
    #[test]
    fn test_phase_separation_index_different_phases() {
        let mut sys = MultiphaseSystem::new(1.0, [0.0, -9.81, 0.0]);
        let water = FluidPhase {
            id: 0,
            density_ref: 1000.0,
            dynamic_viscosity: 1e-3,
            surface_tension: 0.07,
            color: [0.0, 0.0, 1.0],
        };
        let oil = FluidPhase {
            id: 1,
            density_ref: 800.0,
            dynamic_viscosity: 1e-2,
            surface_tension: 0.03,
            color: [1.0, 0.8, 0.0],
        };
        sys.add_phase(water);
        sys.add_phase(oil);
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 0);
        sys.add_particle([0.5, 0.0, 0.0], [0.0; 3], 1);
        sys.particles[0].density = 1000.0;
        sys.particles[1].density = 800.0;
        let idx = sys.phase_separation_index();
        assert!(
            idx < 0.1,
            "Different-phase neighbour pair should give separation index near 0, got {idx}"
        );
    }
}
#[allow(dead_code)]
/// Extract droplet descriptors from the current particle state.
///
/// Returns one [`DropletInfo`] per connected component of `phase_id` particles.
pub fn extract_droplets(
    positions: &[[f64; 3]],
    phase_ids: &[usize],
    h: f64,
    phase_id: usize,
) -> Vec<DropletInfo> {
    let analyzer = PhaseMorphologyAnalyzer::new(h);
    let n = positions.len();
    let (filtered_pos, orig_idx): (Vec<[f64; 3]>, Vec<usize>) = (0..n)
        .filter(|&i| phase_ids[i] == phase_id)
        .map(|i| (positions[i], i))
        .unzip();
    if filtered_pos.is_empty() {
        return Vec::new();
    }
    let filtered_phase: Vec<usize> = vec![phase_id; filtered_pos.len()];
    let labels = analyzer.label_components(&filtered_pos, &filtered_phase);
    let max_label = labels.iter().cloned().max().unwrap_or(0);
    let n_comps = max_label + 1;
    let mut counts = vec![0usize; n_comps];
    let mut sum_pos = vec![[0.0_f64; 3]; n_comps];
    for (k, &label) in labels.iter().enumerate() {
        counts[label] += 1;
        let p = filtered_pos[k];
        sum_pos[label][0] += p[0];
        sum_pos[label][1] += p[1];
        sum_pos[label][2] += p[2];
    }
    let mut droplets = Vec::new();
    for comp in 0..n_comps {
        let nc = counts[comp];
        if nc == 0 {
            continue;
        }
        let center = [
            sum_pos[comp][0] / nc as f64,
            sum_pos[comp][1] / nc as f64,
            sum_pos[comp][2] / nc as f64,
        ];
        let rg2 = labels
            .iter()
            .zip(filtered_pos.iter())
            .filter(|(l, _)| **l == comp)
            .map(|(_, p)| {
                let d = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
                d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
            })
            .sum::<f64>()
            / nc as f64;
        let _ = orig_idx[0];
        droplets.push(DropletInfo {
            label: comp,
            particle_count: nc,
            center,
            radius_of_gyration: rg2.sqrt(),
        });
    }
    droplets
}
#[cfg(test)]
mod tests_ext {
    use super::*;
    use crate::multiphase::types::*;
    #[test]
    fn contact_angle_water_on_glass() {
        let m = ContactAngleModel::new(0.073, 0.049, 0.072);
        let theta = m.contact_angle_deg().expect("should have valid angle");
        assert!(
            (theta - 70.5_f64).abs() < 1.0,
            "water on glass θ ≈ 70.5°, got {theta:.2}°"
        );
    }
    #[test]
    fn contact_angle_complete_wetting() {
        let m = ContactAngleModel::new(1.0, 0.001, 0.072);
        assert!(
            m.contact_angle_rad().is_none(),
            "cos θ > 1 should return None"
        );
    }
    #[test]
    fn spreading_coefficient_positive_for_wetting() {
        let m = ContactAngleModel::new(0.1, 0.01, 0.07);
        assert!(
            m.spreading_coefficient() > 0.0,
            "S should be positive for wetting"
        );
    }
    #[test]
    fn wetting_regime_classification() {
        let m_hydrophobic = ContactAngleModel::new(0.03, 0.08, 0.07);
        assert_eq!(m_hydrophobic.wetting_regime(), WettingRegime::Hydrophobic);
        let m_hydrophilic = ContactAngleModel::new(0.08, 0.03, 0.07);
        assert_eq!(m_hydrophilic.wetting_regime(), WettingRegime::Hydrophilic);
    }
    #[test]
    fn adhesion_work_positive_for_partial_wetting() {
        let m = ContactAngleModel::new(0.073, 0.049, 0.072);
        let w = m.adhesion_work();
        assert!(w > 0.0, "adhesion work should be positive, got {w}");
    }
    #[test]
    fn capillary_pressure_sign_for_wetting_fluid() {
        let m = ContactAngleModel::new(0.08, 0.03, 0.07);
        let r = 1e-3;
        let pc = m.capillary_pressure(r).expect("valid angle");
        assert!(
            pc > 0.0,
            "capillary pressure should be positive for wetting, got {pc}"
        );
    }
    #[test]
    fn capillary_pressure_zero_radius_returns_none() {
        let m = ContactAngleModel::new(0.073, 0.049, 0.072);
        assert!(m.capillary_pressure(0.0).is_none());
    }
    #[test]
    fn atwood_number_equal_densities_is_zero() {
        let rt = RayleighTaylorAnalyzer::new(1000.0, 1000.0, 9.81);
        assert!((rt.atwood_number()).abs() < 1e-14);
    }
    #[test]
    fn atwood_number_range() {
        let rt = RayleighTaylorAnalyzer::new(2000.0, 1000.0, 9.81);
        let a = rt.atwood_number();
        assert!(
            a > 0.0 && a < 1.0,
            "Atwood number should be in (0,1), got {a}"
        );
    }
    #[test]
    fn linear_growth_rate_positive() {
        let rt = RayleighTaylorAnalyzer::new(2000.0, 1000.0, 9.81);
        let sigma = rt.linear_growth_rate(10.0);
        assert!(sigma > 0.0, "growth rate should be positive, got {sigma}");
    }
    #[test]
    fn mixing_layer_width_increases_with_time() {
        let rt = RayleighTaylorAnalyzer::new(2000.0, 1000.0, 9.81);
        let h1 = rt.mixing_layer_width(1.0, 0.06);
        let h2 = rt.mixing_layer_width(2.0, 0.06);
        assert!(h2 > h1, "mixing layer should grow with time");
    }
    #[test]
    fn phase_centers_of_mass_simple() {
        let positions = vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let phase_ids = vec![0, 1];
        let masses = vec![1.0, 1.0];
        let (y_h, y_l) =
            RayleighTaylorAnalyzer::phase_centers_of_mass(&positions, &phase_ids, &masses);
        assert!((y_h - 0.0).abs() < 1e-14);
        assert!((y_l - 1.0).abs() < 1e-14);
    }
    #[test]
    fn morphology_single_cluster_same_phase() {
        let positions = vec![
            [0.0_f64, 0.0, 0.0],
            [0.05, 0.0, 0.0],
            [0.10, 0.0, 0.0],
            [0.15, 0.0, 0.0],
        ];
        let phase_ids = vec![0usize; 4];
        let analyzer = PhaseMorphologyAnalyzer::new(0.12);
        let n_comps = analyzer.component_count(&positions, &phase_ids, 0);
        assert_eq!(n_comps, 1, "all within h → single component");
    }
    #[test]
    fn morphology_two_separate_clusters() {
        let positions = vec![
            [0.0_f64, 0.0, 0.0],
            [0.05, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [10.05, 0.0, 0.0],
        ];
        let phase_ids = vec![0usize; 4];
        let analyzer = PhaseMorphologyAnalyzer::new(0.1);
        let n_comps = analyzer.component_count(&positions, &phase_ids, 0);
        assert_eq!(n_comps, 2, "two separated clusters → two components");
    }
    #[test]
    fn morphology_labels_length_equals_particle_count() {
        let positions: Vec<[f64; 3]> = (0..6).map(|i| [i as f64, 0.0, 0.0]).collect();
        let phase_ids = vec![0, 1, 0, 1, 0, 1];
        let analyzer = PhaseMorphologyAnalyzer::new(0.5);
        let labels = analyzer.label_components(&positions, &phase_ids);
        assert_eq!(labels.len(), 6);
    }
    #[test]
    fn morphology_component_sizes_sum_to_n() {
        let positions: Vec<[f64; 3]> = (0..8).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let phase_ids = vec![0usize; 8];
        let analyzer = PhaseMorphologyAnalyzer::new(0.15);
        let sizes = analyzer.component_sizes(&positions, &phase_ids);
        assert_eq!(sizes.iter().sum::<usize>(), 8);
    }
    #[test]
    fn radius_of_gyration_symmetric_dumbbell() {
        let d = 0.5_f64;
        let positions = vec![[-d, 0.0, 0.0], [d, 0.0, 0.0]];
        let phase_ids = vec![0usize; 2];
        let analyzer = PhaseMorphologyAnalyzer::new(2.0 * d + 0.1);
        let rg = analyzer
            .radius_of_gyration(&positions, &phase_ids, 0)
            .expect("should compute Rg");
        assert!(
            (rg - d).abs() < 1e-12,
            "Rg of symmetric dumbbell should be {d}, got {rg}"
        );
    }
    #[test]
    fn extract_droplets_single_droplet() {
        let positions = vec![[0.0_f64, 0.0, 0.0], [0.05, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let phase_ids = vec![0, 0, 1];
        let droplets = extract_droplets(&positions, &phase_ids, 0.1, 0);
        assert_eq!(droplets.len(), 1, "should find 1 liquid droplet");
        assert_eq!(droplets[0].particle_count, 2);
    }
    #[test]
    fn extract_droplets_empty_phase() {
        let positions = vec![[0.0_f64; 3], [0.1, 0.0, 0.0]];
        let phase_ids = vec![1usize, 1];
        let droplets = extract_droplets(&positions, &phase_ids, 0.2, 0);
        assert!(droplets.is_empty(), "no phase-0 particles → no droplets");
    }
}
/// Laplacian of the cubic spline kernel W(r, h) in 3D.
#[allow(dead_code)]
pub(super) fn cubic_spline_laplacian(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 3.0 / (2.0 * std::f64::consts::PI * h * h * h);
    if q >= 2.0 {
        0.0
    } else if q >= 1.0 {
        let t = 2.0 - q;
        sigma / (h * h) * t
    } else {
        sigma / (h * h) * (-2.0 + 3.0 * q)
    }
}
/// Compute the volume fraction of each phase within a uniform 3D grid of cells.
///
/// Divides the domain `[x_min, x_max] × [y_min, y_max] × [z_min, z_max]` into
/// `nx × ny × nz` cells and counts how many particles of each phase fall into
/// each cell.  The returned array has shape `[n_cells][n_phases]` where
/// `n_cells = nx * ny * nz`.
#[allow(dead_code)]
pub fn volume_fraction_in_cells(
    positions: &[[f64; 3]],
    phase_ids: &[usize],
    n_phases: usize,
    bounds: ([f64; 3], [f64; 3]),
    grid: (usize, usize, usize),
) -> Vec<Vec<f64>> {
    let (lo, hi) = bounds;
    let (nx, ny, nz) = grid;
    let n_cells = nx * ny * nz;
    let mut counts = vec![vec![0usize; n_phases]; n_cells];
    for (pos, &pid) in positions.iter().zip(phase_ids) {
        if pid >= n_phases {
            continue;
        }
        let ix = ((pos[0] - lo[0]) / (hi[0] - lo[0]) * nx as f64).floor() as isize;
        let iy = ((pos[1] - lo[1]) / (hi[1] - lo[1]) * ny as f64).floor() as isize;
        let iz = ((pos[2] - lo[2]) / (hi[2] - lo[2]) * nz as f64).floor() as isize;
        if ix < 0 || iy < 0 || iz < 0 || ix >= nx as isize || iy >= ny as isize || iz >= nz as isize
        {
            continue;
        }
        let cell = ix as usize * ny * nz + iy as usize * nz + iz as usize;
        counts[cell][pid] += 1;
    }
    counts
        .iter()
        .map(|cell_counts| {
            let total: usize = cell_counts.iter().sum();
            if total == 0 {
                vec![0.0; n_phases]
            } else {
                cell_counts
                    .iter()
                    .map(|&c| c as f64 / total as f64)
                    .collect()
            }
        })
        .collect()
}
#[cfg(test)]
mod tests_new_multiphase {
    use super::*;
    use crate::multiphase::types::*;
    #[test]
    fn sharpening_step_clamped_to_unit_interval() {
        let positions = vec![[0.0_f64, 0.0, 0.0], [0.05, 0.0, 0.0], [0.1, 0.0, 0.0]];
        let volumes = vec![1e-4_f64; 3];
        let color = vec![0.1, 0.5, 0.9];
        let sharp = InterfaceSharpening::new(0.01, 0.1);
        let new_c = sharp.step(&positions, &volumes, &color, 0.001);
        for &c in &new_c {
            assert!(
                (0.0..=1.0).contains(&c),
                "color must stay in [0,1], got {c}"
            );
        }
    }
    #[test]
    fn sharpening_step_preserves_pure_phases() {
        let positions = vec![[0.0_f64, 0.0, 0.0], [0.05, 0.0, 0.0], [0.1, 0.0, 0.0]];
        let volumes = vec![1e-4_f64; 3];
        let color = vec![0.0, 0.0, 0.0];
        let sharp = InterfaceSharpening::new(0.01, 0.1);
        let new_c = sharp.step(&positions, &volumes, &color, 0.01);
        for &c in &new_c {
            assert!(c.abs() < 1e-12, "pure-phase c=0 must not change: {c}");
        }
    }
    #[test]
    fn sharpening_step_output_length_matches_input() {
        let n = 10;
        let positions: Vec<[f64; 3]> = (0..n).map(|i| [i as f64 * 0.05, 0.0, 0.0]).collect();
        let volumes = vec![5e-5_f64; n];
        let color: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        let sharp = InterfaceSharpening::new(0.01, 0.08);
        let new_c = sharp.step(&positions, &volumes, &color, 0.001);
        assert_eq!(new_c.len(), n);
    }
    #[test]
    fn vdw_eos_pressure_positive_for_gas() {
        let eos = VanDerWaalsEos::from_critical(647.0, 22.1e6, 461.5);
        let p = eos.pressure(0.6, 700.0);
        assert!(p.is_some(), "should return a pressure");
        assert!(p.unwrap() > 0.0, "pressure should be positive");
    }
    #[test]
    fn vdw_eos_speed_of_sound_positive_in_stable_region() {
        let eos = VanDerWaalsEos::from_critical(647.0, 22.1e6, 461.5);
        let c = eos.speed_of_sound(1000.0, 300.0);
        if let Some(c_val) = c {
            assert!(c_val > 0.0, "speed of sound must be positive: {c_val}");
        }
    }
    #[test]
    fn vdw_eos_spinodal_detection() {
        let eos = VanDerWaalsEos::from_critical(647.0, 22.1e6, 461.5);
        let in_spinodal = eos.in_spinodal_region(0.1, 800.0);
        assert!(
            !in_spinodal,
            "superheated gas should not be in spinodal region"
        );
    }
    #[test]
    fn vdw_eos_reduced_temperature_at_critical_point() {
        let eos = VanDerWaalsEos::from_critical(647.0, 22.1e6, 461.5);
        let tr = eos.reduced_temperature(647.0);
        assert!(
            (tr - 1.0).abs() < 1e-10,
            "T_r at critical point should be 1.0, got {tr}"
        );
    }
    #[test]
    fn vdw_eos_unphysical_state_returns_none() {
        let eos = VanDerWaalsEos::from_critical(647.0, 22.1e6, 461.5);
        let p = eos.pressure(1.0 / eos.b + 100.0, 300.0);
        assert!(p.is_none(), "beyond close-packing should return None");
    }
    #[test]
    fn volume_fraction_single_phase_in_single_cell() {
        let positions = vec![[0.5_f64, 0.5, 0.5]];
        let phase_ids = vec![0_usize];
        let vf =
            volume_fraction_in_cells(&positions, &phase_ids, 2, ([0.0; 3], [1.0; 3]), (1, 1, 1));
        assert_eq!(vf.len(), 1);
        assert!((vf[0][0] - 1.0).abs() < 1e-14, "single phase → fraction 1");
        assert!((vf[0][1] - 0.0).abs() < 1e-14, "other phase → fraction 0");
    }
    #[test]
    fn volume_fraction_two_phases_equal_split() {
        let positions = vec![[0.5_f64, 0.5, 0.5], [0.5, 0.5, 0.5]];
        let phase_ids = vec![0_usize, 1];
        let vf =
            volume_fraction_in_cells(&positions, &phase_ids, 2, ([0.0; 3], [1.0; 3]), (1, 1, 1));
        assert!((vf[0][0] - 0.5).abs() < 1e-14, "50 % phase 0");
        assert!((vf[0][1] - 0.5).abs() < 1e-14, "50 % phase 1");
    }
    #[test]
    fn volume_fraction_out_of_bounds_particle_ignored() {
        let positions = vec![[2.0_f64, 0.5, 0.5]];
        let phase_ids = vec![0_usize];
        let vf =
            volume_fraction_in_cells(&positions, &phase_ids, 1, ([0.0; 3], [1.0; 3]), (1, 1, 1));
        let total: f64 = vf[0].iter().sum();
        assert!(
            total.abs() < 1e-14,
            "out-of-bounds particle should be ignored"
        );
    }
    #[test]
    fn volume_fraction_fractions_sum_to_one() {
        let positions: Vec<[f64; 3]> = (0..8)
            .map(|i| [(i % 2) as f64 * 0.4 + 0.1, 0.5, 0.5])
            .collect();
        let phase_ids: Vec<usize> = (0..8).map(|i| i % 2).collect();
        let vf =
            volume_fraction_in_cells(&positions, &phase_ids, 2, ([0.0; 3], [1.0; 3]), (2, 1, 1));
        for cell in &vf {
            let s: f64 = cell.iter().sum();
            if s > 1e-14 {
                assert!((s - 1.0).abs() < 1e-12, "fractions must sum to 1, got {s}");
            }
        }
    }
    #[test]
    fn coalescence_detector_close_droplets() {
        let detector = CoalescenceDetector::new(0.05);
        let centers = vec![[0.0_f64, 0.0, 0.0], [0.04, 0.0, 0.0]];
        let radii = vec![0.01_f64, 0.01];
        let pairs = detector.find_pairs(&centers, &radii);
        assert_eq!(pairs.len(), 1, "one pair should coalesce");
    }
    #[test]
    fn coalescence_detector_far_droplets_no_pair() {
        let detector = CoalescenceDetector::new(0.001);
        let centers = vec![[0.0_f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let radii = vec![0.01_f64, 0.01];
        let pairs = detector.find_pairs(&centers, &radii);
        assert!(
            pairs.is_empty(),
            "well-separated droplets should not coalesce"
        );
    }
    #[test]
    fn coalescence_merged_radius_volume_conservation() {
        let r1 = 1e-3_f64;
        let r2 = 2e-3_f64;
        let r_m = CoalescenceDetector::merged_radius(r1, r2);
        let expected = (r1.powi(3) + r2.powi(3)).cbrt();
        assert!((r_m - expected).abs() < 1e-20);
    }
    #[test]
    fn coalescence_merged_center_equal_masses() {
        let c1 = [0.0_f64, 0.0, 0.0];
        let c2 = [2.0_f64, 0.0, 0.0];
        let center = CoalescenceDetector::merged_center(c1, 1.0, c2, 1.0);
        assert!((center[0] - 1.0).abs() < 1e-14, "midpoint x = 1");
    }
    #[test]
    fn cubic_spline_laplacian_zero_outside_support() {
        let lap = cubic_spline_laplacian(3.0, 1.0);
        assert_eq!(lap, 0.0, "laplacian must be 0 outside support");
    }
    #[test]
    fn cubic_spline_laplacian_finite_inside_support() {
        let lap = cubic_spline_laplacian(0.5, 1.0);
        assert!(lap.is_finite(), "laplacian must be finite inside support");
    }
}
