//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// Clamp `x` into `[lo, hi]`.
#[inline]
#[allow(dead_code)]
pub(super) fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}
/// Add a chiral off-diagonal coupling to an achiral 6×6 Voigt stiffness matrix.
///
/// Chirality couples normal stresses to shear strains (and vice-versa) via
/// an antisymmetric perturbation parameterised by `chirality` (dimensionless,
/// typically |chirality| ≪ 1 for weakly chiral materials).
///
/// Returns a new 6×6 stiffness matrix.
pub fn chiral_elastic_stiffness(c_achiral: &[[f64; 6]; 6], chirality: f64) -> [[f64; 6]; 6] {
    let mut out = *c_achiral;
    let scale = (c_achiral[0][0] + c_achiral[1][1] + c_achiral[2][2]) / 3.0;
    let delta = chirality * scale;
    #[allow(clippy::needless_range_loop)]
    for i in 0..3 {
        for j in 3..6 {
            out[i][j] += delta;
            out[j][i] += delta;
        }
    }
    out
}
/// Estimate the average mass density of an aperiodic (Penrose-tiling) metamaterial.
///
/// `phi` is the golden ratio (≈ 1.618…) governs the relative area of the two
/// Penrose tile types (fat and thin rhombi).  `unit_area` is the area of the
/// reference rhombus (m²).
///
/// Returns density in kg/m² (areal density; divide by thickness for
/// volumetric density).
pub fn penrose_tiling_density(phi: f64, unit_area: f64) -> f64 {
    let a_fat = unit_area * (72.0_f64.to_radians()).sin();
    let a_thin = unit_area * (36.0_f64.to_radians()).sin();
    (phi * a_fat + a_thin) / (phi + 1.0)
}
/// Generate a topology-optimised binary lattice with a target Poisson's ratio.
///
/// Uses a simplified SIMP-like (Solid Isotropic Material with Penalisation)
/// sensitivity iteration on a 2-D grid.  Each cell is either solid (`true`) or
/// void (`false`).  The grid is square with side length `sqrt(iterations)`,
/// clamped to a 32×32 maximum for tractability.
///
/// `target_poisson` is the desired effective Poisson's ratio (negative for auxetics).
/// `iterations` controls the number of update sweeps (minimum 1 sweep is always run).
///
/// Returns a 2-D boolean grid (outer index = row, inner = column).
pub fn topology_optimized_auxetic(target_poisson: f64, iterations: usize) -> Vec<Vec<bool>> {
    let n = ((iterations as f64).sqrt() as usize).clamp(4, 32);
    let mut grid = vec![vec![true; n]; n];
    if iterations == 0 {
        return grid;
    }
    if target_poisson < 0.0 {
        for i in (2..n - 1).step_by(3) {
            for j in (2..n - 1).step_by(3) {
                grid[i][j] = false;
            }
        }
    }
    for _iter in 0..iterations.min(200) {
        let prev = grid.clone();
        for i in 1..n - 1 {
            for j in 1..n - 1 {
                let neighbors_solid = [
                    prev[i - 1][j],
                    prev[i + 1][j],
                    prev[i][j - 1],
                    prev[i][j + 1],
                ]
                .iter()
                .filter(|&&v| v)
                .count();
                let h_solid = prev[i][j - 1] && prev[i][j + 1];
                let v_solid = prev[i - 1][j] && prev[i + 1][j];
                if target_poisson < 0.0 {
                    if neighbors_solid == 2 && (h_solid || v_solid) {
                        grid[i][j] = false;
                    }
                } else if neighbors_solid < 2 {
                    grid[i][j] = false;
                }
            }
        }
    }
    grid
}
/// Compute the Poisson's ratio of a re-entrant honeycomb from geometry.
///
/// Uses the Masters–Evans formula for the re-entrant cell:
/// ```text
/// ν = −sin(θ)·(h/l + sin(θ)) / cos²(θ)
/// ```
/// where `angle_deg` (θ) is the cell angle (degrees, negative for re-entrant),
/// `l1` = inclined rib length and `l2` = vertical rib length.
pub fn auxetic_poissons_from_geometry(angle_deg: f64, l1: f64, l2: f64) -> f64 {
    let theta = angle_deg.to_radians();
    let sin_t = theta.sin();
    let cos_t = theta.cos();
    let cos2 = cos_t * cos_t;
    if cos2.abs() < 1e-15 {
        return 0.0;
    }
    let h_over_l = l2 / l1.max(1e-30);
    -sin_t * (h_over_l + sin_t) / cos2
}
/// Compute the relative density of a regular honeycomb.
///
/// ```text
/// ρ*/ρ_s = 2·t / l
/// ```
///
/// where `t` is wall thickness (m) and `l` is cell edge length (m).
pub fn honeycomb_relative_density(t: f64, l: f64) -> f64 {
    if l < 1e-30 {
        return 1.0;
    }
    (2.0 * t / l).min(1.0)
}
/// Effective Young's modulus of a regular hexagonal honeycomb (Gibson–Ashby).
///
/// ```text
/// E* = E_s · (ρ*/ρ_s)³ / cos(θ) · (t/l)³
/// ```
/// Simplified to `E* = E_s · (rel_density)³ · f(angle)`.
pub fn honeycomb_youngs_modulus(e_s: f64, rel_density: f64, angle_deg: f64) -> f64 {
    let theta = angle_deg.to_radians();
    let cos_t = theta.cos().abs().max(1e-15);
    e_s * rel_density.powi(3) / cos_t
}
/// Linearly graded metamaterial stiffness at position `x` along a beam of `length`.
///
/// Interpolates from `e_min` at x=0 to `e_max` at x=`length`.
pub fn gradient_metamaterial_stiffness(x: f64, e_min: f64, e_max: f64, length: f64) -> f64 {
    if length < 1e-30 {
        return e_min;
    }
    let t = (x / length).clamp(0.0, 1.0);
    e_min + t * (e_max - e_min)
}
/// Acoustic bandgap centre frequency (Hz) of a mass–spring resonator.
///
/// ```text
/// f = (1 / 2π) · √(k / m)
/// ```
pub fn acoustic_bandgap_frequency(mass: f64, stiffness: f64) -> f64 {
    if mass <= 0.0 {
        return 0.0;
    }
    (stiffness / mass).sqrt() / (2.0 * PI)
}
/// Check whether a medium satisfies the negative-index condition.
///
/// Returns `true` when both `epsilon` and `mu` are strictly negative,
/// implying a negative refractive index.
pub fn negative_index_condition(epsilon: f64, mu: f64) -> bool {
    epsilon < 0.0 && mu < 0.0
}
/// Effective bulk modulus of a pentamode lattice (dilute limit).
///
/// ```text
/// K_eff = E_rod · volume_fraction / 3
/// ```
pub fn pentamode_bulk_modulus(e_rod: f64, volume_fraction: f64) -> f64 {
    e_rod * volume_fraction / 3.0
}
/// Effective thermal conductivity of a two-phase composite (Maxwell model).
///
/// ```text
/// k_eff = k1 · (k1 + k2 + v1·(k2 − k1)) / (k1 + k2 − v1·(k2 − k1))
/// ```
///
/// where `v1` is the volume fraction of phase 1.
pub fn thermal_metamaterial_flux_ratio(k1: f64, k2: f64, v1: f64) -> f64 {
    let num = k1 + k2 + v1 * (k2 - k1);
    let den = k1 + k2 - v1 * (k2 - k1);
    if den.abs() < 1e-30 {
        return k1;
    }
    k1 * num / den
}
/// Generate per-layer (E, ν) properties for a cylindrical mechanical cloak.
///
/// Linearly interpolates Young's modulus and Poisson's ratio across `n_layer`
/// concentric shells between `r_inner` and `r_outer`.  Outer layers are
/// softer (low E, high ν) and inner layers are stiffer.
///
/// Returns a `Vec` of `(E, nu)` pairs ordered from innermost to outermost layer.
pub fn mechanical_cloak_layer_properties(
    r_inner: f64,
    r_outer: f64,
    n_layer: usize,
) -> Vec<(f64, f64)> {
    if n_layer == 0 {
        return Vec::new();
    }
    let e_inner = 200.0e9f64;
    let e_outer = 1.0e9f64;
    let nu_inner = 0.10f64;
    let nu_outer = 0.49f64;
    let _r_inner = r_inner;
    let _r_outer = r_outer;
    (0..n_layer)
        .map(|i| {
            let t = if n_layer == 1 {
                0.0
            } else {
                i as f64 / (n_layer - 1) as f64
            };
            let e = e_inner + t * (e_outer - e_inner);
            let nu = nu_inner + t * (nu_outer - nu_inner);
            (e, nu)
        })
        .collect()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::metamaterials::AcousticCloakShell;
    use crate::metamaterials::AcousticMetamaterial;
    use crate::metamaterials::ActiveMetamaterial;
    use crate::metamaterials::AuxeticMaterial;
    use crate::metamaterials::DoublenegativeMedium;
    use crate::metamaterials::ElasticMetamaterial;
    use crate::metamaterials::ElectromagneticMetamaterial;
    use crate::metamaterials::PentamodeMetamaterial;
    use crate::metamaterials::PiezoResonator;
    use crate::metamaterials::ReentrantHoneycomb;
    use crate::metamaterials::SplitRingResonator;
    use crate::metamaterials::TopologicalInsulator;
    use crate::metamaterials::TopologicalInvariant;
    use crate::metamaterials::TransformationAcoustics;
    pub(super) const EPS: f64 = 1e-9;
    #[test]
    fn test_auxetic_poisson_negative() {
        let mat = AuxeticMaterial::new(-0.5, 1e9, PI / 6.0);
        assert!(
            mat.effective_poisson() < 0.0,
            "effective_poisson should be negative for auxetic material"
        );
    }
    #[test]
    fn test_auxetic_stiffness_scaling() {
        let mat = AuxeticMaterial::new(-0.3, 1e9, PI / 4.0);
        let e1 = mat.effective_stiffness(0.1);
        let e2 = mat.effective_stiffness(0.2);
        assert!(e2 > e1, "Stiffness should increase with density ratio");
    }
    #[test]
    fn test_auxetic_stiffness_bound() {
        let mat = AuxeticMaterial::new(-0.5, 1e9, 0.0);
        let e = mat.effective_stiffness(1.0);
        assert!(
            e <= mat.young_modulus + EPS,
            "Effective stiffness at density_ratio=1 should not exceed Es"
        );
    }
    #[test]
    fn test_reentrant_poisson_negative() {
        let h = ReentrantHoneycomb::new(-PI / 6.0, 1.0, 0.3, 0.1, 1e9);
        assert!(
            h.poisson_ratio_xy() < 0.0,
            "Re-entrant cell with h/l=0.3 and theta=-PI/6 should have negative Poisson ratio, got {}",
            h.poisson_ratio_xy()
        );
    }
    #[test]
    fn test_reentrant_moduli_positive() {
        let h = ReentrantHoneycomb::new(-PI / 6.0, 1.0, 0.3, 0.1, 1e9);
        assert!(h.young_modulus_x() > 0.0, "E_x must be positive");
        assert!(h.young_modulus_y() > 0.0, "E_y must be positive");
    }
    #[test]
    fn test_reentrant_shear_positive() {
        let h = ReentrantHoneycomb::new(-PI / 6.0, 1.0, 0.3, 0.1, 1e9);
        assert!(h.shear_modulus_xy() > 0.0, "G_xy must be positive");
    }
    #[test]
    fn test_reentrant_thickness_effect() {
        let h1 = ReentrantHoneycomb::new(-PI / 6.0, 1.0, 0.3, 0.05, 1e9);
        let h2 = ReentrantHoneycomb::new(-PI / 6.0, 1.0, 0.3, 0.10, 1e9);
        assert!(
            h2.young_modulus_x() > h1.young_modulus_x(),
            "Thicker ribs should give higher E_x"
        );
    }
    #[test]
    fn test_acoustic_mm_density_below_resonance() {
        let mm = AcousticMetamaterial::new(0.01, 0.001, 100.0, 1000.0, 1e6);
        let w0 = (100.0_f64 / 0.001).sqrt();
        let rho_eff = mm.effective_density(0.5 * w0);
        assert!(
            rho_eff > 0.0,
            "Effective density below resonance should be positive, got {rho_eff}"
        );
    }
    #[test]
    fn test_acoustic_mm_density_above_resonance() {
        let mm = AcousticMetamaterial::new(0.01, 0.001, 100.0, 1000.0, 1e6);
        let w0 = (100.0_f64 / 0.001).sqrt();
        let rho_eff = mm.effective_density(1.01 * w0);
        assert!(
            rho_eff < 0.0,
            "Effective density just above resonance should be negative, got {rho_eff}"
        );
    }
    #[test]
    fn test_acoustic_mm_band_gap_order() {
        let mm = AcousticMetamaterial::new(0.01, 0.001, 100.0, 1000.0, 1e6);
        let (f_lo, f_hi) = mm.band_gap_range();
        assert!(f_lo < f_hi, "f_lo={f_lo} should be < f_hi={f_hi}");
        assert!(f_lo > 0.0, "f_lo must be positive");
    }
    #[test]
    fn test_pentamode_sound_speed() {
        let pm = PentamodeMetamaterial::new(4e9, 2e9);
        let c = pm.sound_speed(1000.0);
        let expected = (4e9_f64 / 1000.0).sqrt();
        assert!(
            (c - expected).abs() < 1.0,
            "Sound speed mismatch: got {c}, expected {expected}"
        );
    }
    #[test]
    fn test_pentamode_zero_density() {
        let pm = PentamodeMetamaterial::new(4e9, 2e9);
        assert_eq!(pm.sound_speed(0.0), 0.0);
    }
    #[test]
    #[allow(clippy::needless_range_loop)]
    fn test_transform_acoustics_identity() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let cloak = TransformationAcoustics::new(identity);
        let rho_tensor = cloak.effective_density_tensor(1200.0);
        for i in 0..3 {
            assert!(
                (rho_tensor[i][i] - 1200.0).abs() < 1e-6,
                "Diagonal [{i}][{i}] = {} should be 1200",
                rho_tensor[i][i]
            );
        }
        assert!(
            (cloak.effective_bulk_modulus(1e6) - 1e6).abs() < EPS,
            "Identity jacobian: bulk modulus should equal k0"
        );
    }
    #[test]
    fn test_transform_acoustics_scaling() {
        let s = 2.0;
        let scaled = [[s, 0.0, 0.0], [0.0, s, 0.0], [0.0, 0.0, s]];
        let cloak = TransformationAcoustics::new(scaled);
        let rho_tensor = cloak.effective_density_tensor(1000.0);
        let expected_rho = 1000.0 / s;
        assert!(
            (rho_tensor[0][0] - expected_rho).abs() < 1e-6,
            "ρ_eff[0][0] = {} should be {expected_rho}",
            rho_tensor[0][0]
        );
        let expected_keff = 1e6 / (s * s * s);
        assert!(
            (cloak.effective_bulk_modulus(1e6) - expected_keff).abs() < 1e-3,
            "κ_eff mismatch"
        );
    }
    #[test]
    fn test_doublenegative_negative_index() {
        let dnm = DoublenegativeMedium::new(-2.5, -1.8);
        let n = dnm.refractive_index();
        assert!(
            n < 0.0,
            "Double-negative medium must have negative index, got {n}"
        );
    }
    #[test]
    fn test_doublenegative_positive_index() {
        let dnm = DoublenegativeMedium::new(2.25, 1.0);
        let n = dnm.refractive_index();
        assert!((n - 1.5).abs() < 1e-9, "n = sqrt(2.25*1.0) = 1.5, got {n}");
    }
    #[test]
    fn test_doublenegative_mixed_sign() {
        let dnm = DoublenegativeMedium::new(-1.0, 1.0);
        assert_eq!(dnm.refractive_index(), 0.0);
    }
    #[test]
    fn test_doublenegative_phase_velocity_negative() {
        let dnm = DoublenegativeMedium::new(-2.0, -2.0);
        let vp = dnm.phase_velocity();
        assert!(
            vp < 0.0,
            "Phase velocity should be negative for DNG medium, got {vp}"
        );
    }
    #[test]
    #[allow(clippy::needless_range_loop)]
    fn test_chiral_zero() {
        let c: [[f64; 6]; 6] = {
            let mut m = [[0.0f64; 6]; 6];
            for i in 0..6 {
                m[i][i] = 1e10;
            }
            m
        };
        let out = chiral_elastic_stiffness(&c, 0.0);
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (out[i][j] - c[i][j]).abs() < EPS,
                    "Zero chirality should not change matrix at [{i}][{j}]"
                );
            }
        }
    }
    #[test]
    #[allow(clippy::needless_range_loop)]
    fn test_chiral_nonzero() {
        let c: [[f64; 6]; 6] = {
            let mut m = [[0.0f64; 6]; 6];
            for i in 0..6 {
                m[i][i] = 1e10;
            }
            m
        };
        let out = chiral_elastic_stiffness(&c, 0.01);
        assert!(
            out[0][3].abs() > EPS,
            "Chirality should add off-diagonal coupling"
        );
    }
    #[test]
    fn test_penrose_density_positive() {
        let d = penrose_tiling_density(1.618, 1.0);
        assert!(d > 0.0, "Penrose density must be positive, got {d}");
    }
    #[test]
    fn test_penrose_density_scales_with_area() {
        let d1 = penrose_tiling_density(1.618, 1.0);
        let d2 = penrose_tiling_density(1.618, 2.0);
        assert!(
            (d2 - 2.0 * d1).abs() < 1e-12,
            "Density should scale linearly with unit_area"
        );
    }
    #[test]
    fn test_topo_grid_dimensions() {
        let grid = topology_optimized_auxetic(-0.5, 16);
        let n = grid.len();
        assert!(n >= 4, "Grid must have at least 4 rows");
        assert_eq!(grid[0].len(), n, "Grid must be square");
    }
    #[test]
    fn test_topo_auxetic_creates_voids() {
        let grid = topology_optimized_auxetic(-0.5, 200);
        let void_count: usize = grid.iter().flatten().filter(|&&v| !v).count();
        assert!(
            void_count > 0,
            "Auxetic target should create at least one void cell (iterations=200)"
        );
    }
    #[test]
    fn test_topo_zero_iterations() {
        let grid = topology_optimized_auxetic(-0.5, 0);
        let solid_count: usize = grid.iter().flatten().filter(|&&v| v).count();
        let total = grid.len() * grid[0].len();
        assert_eq!(
            solid_count, total,
            "Zero iterations should leave all cells solid"
        );
    }
    #[test]
    fn test_acoustic_mm_modulus_constant() {
        let mm = AcousticMetamaterial::new(0.01, 0.001, 100.0, 1000.0, 2e6);
        let w0 = (100.0_f64 / 0.001).sqrt();
        let m1 = mm.effective_modulus(0.5 * w0);
        let m2 = mm.effective_modulus(2.0 * w0);
        assert!(
            (m1 - m2).abs() < EPS,
            "Effective modulus should be frequency-independent in this model"
        );
    }
    #[test]
    fn test_reentrant_theta_sign() {
        let hn = ReentrantHoneycomb::new(-PI / 6.0, 1.0, 0.3, 0.1, 1e9);
        let hp = ReentrantHoneycomb::new(PI / 6.0, 1.0, 0.3, 0.1, 1e9);
        let nu_n = hn.poisson_ratio_xy();
        let nu_p = hp.poisson_ratio_xy();
        assert!(nu_n < 0.0, "nu_n should be negative: {nu_n}");
        assert!(nu_p < 0.0, "nu_p should be negative: {nu_p}");
        assert!(
            (nu_n.abs() - nu_p.abs()).abs() > 1e-9,
            "Magnitudes of ν should differ when θ flips sign: nu_n={nu_n}, nu_p={nu_p}"
        );
    }
    #[test]
    fn test_transform_acoustics_singular() {
        let singular = [[0.0f64; 3]; 3];
        let cloak = TransformationAcoustics::new(singular);
        assert_eq!(cloak.effective_bulk_modulus(1e6), 0.0);
    }
    #[test]
    fn test_auxetic_zero_angle() {
        let mat = AuxeticMaterial::new(-0.4, 1e9, 0.0);
        assert!(
            (mat.effective_poisson() - mat.poisson_ratio).abs() < EPS,
            "cell_angle=0 should give effective_poisson == poisson_ratio"
        );
    }
    #[test]
    #[allow(clippy::needless_range_loop)]
    fn test_chiral_symmetry() {
        let c: [[f64; 6]; 6] = {
            let mut m = [[0.0f64; 6]; 6];
            for i in 0..6 {
                m[i][i] = 1e10;
            }
            m
        };
        let out = chiral_elastic_stiffness(&c, 0.05);
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (out[i][j] - out[j][i]).abs() < EPS,
                    "Stiffness matrix should be symmetric at [{i}][{j}]"
                );
            }
        }
    }
    #[test]
    fn test_acoustic_mm_static_density() {
        let mm = AcousticMetamaterial::new(0.01, 0.001, 400.0, 1200.0, 1e6);
        let rho = mm.effective_density(1.0);
        assert!(
            rho.is_finite(),
            "Effective density at low omega should be finite"
        );
    }
    #[test]
    fn test_doublenegative_phase_velocity_magnitude() {
        let dnm = DoublenegativeMedium::new(-1.0, -1.0);
        let n = dnm.refractive_index();
        let vp = dnm.phase_velocity();
        let c = 2.997_924_58e8;
        assert!(
            (vp - c / n).abs() < 1.0,
            "Phase velocity magnitude mismatch: vp={vp}, c/n={}",
            c / n
        );
    }
    #[test]
    fn test_topo_positive_target_solid() {
        let grid = topology_optimized_auxetic(0.3, 25);
        let solid_count: usize = grid.iter().flatten().filter(|&&v| v).count();
        let total = grid.len() * grid[0].len();
        assert!(
            solid_count > total / 2,
            "Positive target should keep most cells solid: solid={solid_count}/{total}"
        );
    }
    #[test]
    fn test_penrose_phi_one() {
        let d = penrose_tiling_density(1.0, 1.0);
        let a_fat = (72.0_f64.to_radians()).sin();
        let a_thin = (36.0_f64.to_radians()).sin();
        let expected = (a_fat + a_thin) / 2.0;
        assert!(
            (d - expected).abs() < 1e-12,
            "Density at phi=1 should be (a_fat+a_thin)/2, got {d}"
        );
    }
    #[test]
    fn test_pentamode_speed_c11_scaling() {
        let pm1 = PentamodeMetamaterial::new(1e9, 0.5e9);
        let pm2 = PentamodeMetamaterial::new(4e9, 0.5e9);
        let c1 = pm1.sound_speed(1000.0);
        let c2 = pm2.sound_speed(1000.0);
        assert!(
            (c2 / c1 - 2.0).abs() < 1e-9,
            "Sound speed should scale as sqrt(c11): ratio={}",
            c2 / c1
        );
    }
    #[test]
    fn test_auxetic_poissons_negative_angle() {
        let nu = auxetic_poissons_from_geometry(-30.0, 1.0, 0.3);
        assert!(
            nu < 0.0,
            "Re-entrant angle should give negative Poisson's ratio, got {nu}"
        );
    }
    #[test]
    fn test_honeycomb_density_cap() {
        let rho = honeycomb_relative_density(0.5, 1.0);
        assert!(rho <= 1.0);
    }
    #[test]
    fn test_honeycomb_density_zero() {
        let rho = honeycomb_relative_density(0.0, 1.0);
        assert_eq!(rho, 0.0);
    }
    #[test]
    fn test_honeycomb_modulus_cubic() {
        let e1 = honeycomb_youngs_modulus(1e9, 0.1, 0.0);
        let e2 = honeycomb_youngs_modulus(1e9, 0.2, 0.0);
        assert!(
            (e2 / e1 - 8.0).abs() < 1e-9,
            "E should scale as rel_density³"
        );
    }
    #[test]
    fn test_gradient_at_zero() {
        let e = gradient_metamaterial_stiffness(0.0, 1e9, 10e9, 1.0);
        assert!((e - 1e9).abs() < EPS);
    }
    #[test]
    fn test_gradient_at_end() {
        let e = gradient_metamaterial_stiffness(1.0, 1e9, 10e9, 1.0);
        assert!((e - 10e9).abs() < EPS);
    }
    #[test]
    fn test_bandgap_frequency_positive() {
        let f = acoustic_bandgap_frequency(0.01, 1000.0);
        assert!(f > 0.0);
    }
    #[test]
    fn test_bandgap_frequency_formula() {
        let m = 0.01f64;
        let k = 1000.0f64;
        let f = acoustic_bandgap_frequency(m, k);
        let expected = (k / m).sqrt() / (2.0 * PI);
        assert!((f - expected).abs() < EPS);
    }
    #[test]
    fn test_negative_index_both_negative() {
        assert!(negative_index_condition(-1.0, -2.0));
    }
    #[test]
    fn test_negative_index_mixed() {
        assert!(!negative_index_condition(-1.0, 1.0));
        assert!(!negative_index_condition(1.0, -1.0));
    }
    #[test]
    fn test_pentamode_bulk_scales() {
        let k1 = pentamode_bulk_modulus(100e9, 0.1);
        let k2 = pentamode_bulk_modulus(100e9, 0.2);
        assert!((k2 - 2.0 * k1).abs() < EPS);
    }
    #[test]
    fn test_thermal_flux_v0() {
        let k = thermal_metamaterial_flux_ratio(10.0, 1.0, 0.0);
        assert!((k - 10.0).abs() < EPS);
    }
    #[test]
    fn test_cloak_empty() {
        let v = mechanical_cloak_layer_properties(0.1, 0.5, 0);
        assert!(v.is_empty());
    }
    #[test]
    fn test_cloak_layer_count() {
        let v = mechanical_cloak_layer_properties(0.1, 0.5, 5);
        assert_eq!(v.len(), 5);
    }
    #[test]
    fn test_cloak_stiffness_gradient() {
        let v = mechanical_cloak_layer_properties(0.1, 0.5, 5);
        assert!(v[0].0 > v[4].0, "Inner layer should be stiffer than outer");
    }
    fn make_elastic_meta() -> ElasticMetamaterial {
        ElasticMetamaterial::new(1e9, -0.3, 800.0, 3.0, -PI / 6.0, 200e9, 7850.0)
    }
    #[test]
    fn test_elastic_meta_longitudinal_speed_positive() {
        let m = make_elastic_meta();
        let c = m.longitudinal_wave_speed();
        assert!(c > 0.0, "c_L={c}");
    }
    #[test]
    fn test_elastic_meta_shear_speed_positive() {
        let m = make_elastic_meta();
        let cs = m.shear_wave_speed();
        assert!(cs > 0.0, "c_S={cs}");
    }
    #[test]
    fn test_elastic_meta_bulk_modulus_positive() {
        let m = make_elastic_meta();
        let k = m.bulk_modulus();
        assert!(k > 0.0, "K={k}");
    }
    #[test]
    fn test_elastic_meta_is_auxetic() {
        let m = make_elastic_meta();
        assert!(m.is_auxetic());
    }
    #[test]
    fn test_elastic_meta_amplified_density() {
        let m = make_elastic_meta();
        let rho_amp = m.amplified_density();
        assert!((rho_amp - 3.0 * 7850.0).abs() < 1.0, "rho_amp={rho_amp}");
    }
    #[test]
    fn test_elastic_meta_inertial_bandgap() {
        let m = make_elastic_meta();
        let (f0, f1) = m.inertial_bandgap(1e6);
        assert!(f0 > 0.0 && f1 > f0, "f0={f0} f1={f1}");
    }
    #[test]
    fn test_elastic_meta_compute_auxetic_poisson() {
        let mut m = make_elastic_meta();
        let nu = m.compute_auxetic_poisson(0.3);
        assert!(nu < 0.0, "nu={nu}");
    }
    fn make_srr() -> SplitRingResonator {
        SplitRingResonator::new(1e-3, 1e-4, 1e-4, 5e-3, 1e-9, 1e-12)
    }
    #[test]
    fn test_srr_resonance_positive() {
        let srr = make_srr();
        let w0 = srr.resonance_frequency();
        assert!(w0 > 0.0, "w0={w0}");
    }
    #[test]
    fn test_srr_filling_fraction_range() {
        let srr = make_srr();
        let f = srr.filling_fraction();
        assert!(f > 0.0 && f < 1.0, "f={f}");
    }
    #[test]
    fn test_srr_resonance_formula() {
        let srr = make_srr();
        let expected = 1.0 / (1e-9_f64 * 1e-12_f64).sqrt();
        let got = srr.resonance_frequency();
        assert!((got - expected).abs() < 1.0, "w0={got} expected={expected}");
    }
    fn make_em_meta() -> ElectromagneticMetamaterial {
        let srr = make_srr();
        ElectromagneticMetamaterial::new(srr, 2e11, 1e9, 1e8)
    }
    #[test]
    fn test_em_permittivity_below_plasma() {
        let m = make_em_meta();
        let eps = m.effective_permittivity(1e10);
        assert!(eps < 0.0, "eps={eps}");
    }
    #[test]
    fn test_em_permittivity_above_plasma() {
        let m = make_em_meta();
        let eps = m.effective_permittivity(1e13);
        assert!(eps > 0.0 && eps < 1.0 + 1e-3, "eps={eps}");
    }
    #[test]
    fn test_em_permeability_finite() {
        let m = make_em_meta();
        let mu = m.effective_permeability(1e11);
        assert!(mu.is_finite(), "mu={mu}");
    }
    #[test]
    fn test_em_refractive_index_finite() {
        let m = make_em_meta();
        let n = m.refractive_index(1e11);
        assert!(n.is_finite(), "n={n}");
    }
    #[test]
    fn test_em_phase_velocity_finite() {
        let m = make_em_meta();
        let vp = m.phase_velocity(1e13);
        assert!(vp.is_finite());
    }
    #[test]
    fn test_em_left_handed_bandwidth() {
        let m = make_em_meta();
        let bw = m.left_handed_bandwidth(1e10, 1e13, 1000);
        assert!(bw >= 0.0);
    }
    fn make_cloak_shell() -> AcousticCloakShell {
        AcousticCloakShell::new(0.05, 0.10, 1000.0, 2.2e9, 10)
    }
    #[test]
    fn test_cloak_shell_radial_density_in_shell() {
        let c = make_cloak_shell();
        let rho = c.radial_density(0.075);
        assert!(rho >= 0.0, "rho_r={rho}");
    }
    #[test]
    fn test_cloak_shell_radial_density_outside() {
        let c = make_cloak_shell();
        assert_eq!(c.radial_density(0.11), 0.0);
        assert_eq!(c.radial_density(0.04), 0.0);
    }
    #[test]
    fn test_cloak_shell_angular_density_positive() {
        let c = make_cloak_shell();
        let rho = c.angular_density(0.075);
        assert!(rho > 0.0, "rho_theta={rho}");
    }
    #[test]
    fn test_cloak_shell_effective_modulus() {
        let c = make_cloak_shell();
        let kappa = c.effective_modulus(0.075);
        assert!(kappa >= 0.0);
    }
    #[test]
    fn test_cloak_shell_layer_count() {
        let c = make_cloak_shell();
        let layers = c.layer_properties();
        assert_eq!(layers.len(), 10);
    }
    #[test]
    fn test_cloak_shell_layers_increasing_r() {
        let c = make_cloak_shell();
        let layers = c.layer_properties();
        for i in 1..layers.len() {
            assert!(layers[i].0 > layers[i - 1].0, "r should increase");
        }
    }
    fn make_topo_topo() -> TopologicalInsulator {
        TopologicalInsulator::new(1.0, 3.0, 1.0, 10)
    }
    fn make_topo_trivial() -> TopologicalInsulator {
        TopologicalInsulator::new(3.0, 1.0, 1.0, 10)
    }
    #[test]
    fn test_topo_chern_number_topological() {
        let t = make_topo_topo();
        assert_eq!(t.invariant.chern_number, 1);
    }
    #[test]
    fn test_topo_chern_number_trivial() {
        let t = make_topo_trivial();
        assert_eq!(t.invariant.chern_number, 0);
    }
    #[test]
    fn test_topo_band_gap_positive() {
        let t = make_topo_topo();
        assert!(t.band_gap > 0.0, "band_gap={}", t.band_gap);
    }
    #[test]
    fn test_topo_bulk_edge_correspondence_topo() {
        let t = make_topo_topo();
        assert!(t.invariant.satisfies_bulk_edge_correspondence());
    }
    #[test]
    fn test_topo_bulk_edge_correspondence_trivial() {
        let t = make_topo_trivial();
        assert!(t.invariant.satisfies_bulk_edge_correspondence());
    }
    #[test]
    fn test_topo_berry_phase_topological() {
        let t = make_topo_topo();
        assert!((t.berry_phase() - PI).abs() < 1e-10);
    }
    #[test]
    fn test_topo_berry_phase_trivial() {
        let t = make_topo_trivial();
        assert_eq!(t.berry_phase(), 0.0);
    }
    #[test]
    fn test_topo_edge_state_exists_in_topological() {
        let t = make_topo_topo();
        assert!(t.edge_state_frequency().is_some());
    }
    #[test]
    fn test_topo_edge_state_absent_trivial() {
        let t = make_topo_trivial();
        assert!(t.edge_state_frequency().is_none());
    }
    #[test]
    fn test_topo_dispersion_finite() {
        let t = make_topo_topo();
        for i in 0..=10 {
            let k = i as f64 / 10.0;
            let w = t.dispersion(k);
            assert!(w.is_finite() && w >= 0.0, "k={k} w={w}");
        }
    }
    #[test]
    fn test_topo_set_couplings() {
        let mut t = make_topo_trivial();
        t.set_couplings(1.0, 5.0);
        assert_eq!(t.invariant.chern_number, 1);
    }
    fn make_active_meta() -> ActiveMetamaterial {
        let resonators = vec![
            PiezoResonator::new(1e-3, 1e3, 100.0, 5.0, 1e-10),
            PiezoResonator::new(1e-3, 1e3, 100.0, 3.0, 1e-10),
            PiezoResonator::new(1e-3, 1e3, 100.0, -2.0, 1e-10),
        ];
        ActiveMetamaterial::new(resonators, 1e4, 1000.0, 0.5, 0.5)
    }
    #[test]
    fn test_active_meta_mean_freq_positive() {
        let m = make_active_meta();
        let f = m.mean_resonance_frequency();
        assert!(f > 0.0, "f={f}");
    }
    #[test]
    fn test_active_meta_stiffness_range() {
        let m = make_active_meta();
        let (k_min, k_max) = m.stiffness_range();
        assert!(k_max >= k_min, "k_min={k_min} k_max={k_max}");
    }
    #[test]
    fn test_active_meta_pt_symmetry_balanced() {
        let m = make_active_meta();
        assert!(m.check_pt_symmetry());
        assert!(m.pt_symmetric);
    }
    #[test]
    fn test_active_meta_pt_symmetry_unbalanced() {
        let mut m = make_active_meta();
        m.gain = 1.0;
        m.loss = 0.5;
        assert!(!m.check_pt_symmetry());
    }
    #[test]
    fn test_active_meta_exceptional_point_zero() {
        let m = make_active_meta();
        let ep = m.exceptional_point_proximity();
        assert_eq!(ep, 0.0);
    }
    #[test]
    fn test_active_meta_exceptional_point_nonzero() {
        let mut m = make_active_meta();
        m.gain = 1.0;
        m.loss = 0.5;
        let ep = m.exceptional_point_proximity();
        assert!(ep > 0.0, "ep={ep}");
    }
    #[test]
    fn test_active_meta_set_all_voltages() {
        let mut m = make_active_meta();
        m.set_all_voltages(10.0);
        for r in &m.resonators {
            assert_eq!(r.voltage, 10.0);
        }
    }
    #[test]
    fn test_active_meta_gain_loss_balance_zero() {
        let m = make_active_meta();
        assert_eq!(m.gain_loss_balance(), 0.0);
    }
    #[test]
    fn test_active_meta_effective_attenuation_zero() {
        let m = make_active_meta();
        assert_eq!(m.effective_attenuation(), 0.0);
    }
    #[test]
    fn test_piezo_resonator_elongation() {
        let r = PiezoResonator::new(1e-3, 1e3, 100.0, 5.0, 2e-10);
        let elong = r.piezo_elongation();
        assert!((elong - 1e-9).abs() < 1e-20, "elong={elong}");
    }
    #[test]
    fn test_piezo_resonator_effective_stiffness() {
        let r = PiezoResonator::new(1e-3, 1e3, 10.0, 5.0, 1e-10);
        assert!((r.effective_stiffness() - 1050.0).abs() < 1e-8);
    }
    #[test]
    fn test_topological_invariant_bulk_edge() {
        let inv = TopologicalInvariant::new(2, 0, 2);
        assert!(inv.satisfies_bulk_edge_correspondence());
    }
    #[test]
    fn test_topological_invariant_bulk_edge_fail() {
        let inv = TopologicalInvariant::new(2, 0, 1);
        assert!(!inv.satisfies_bulk_edge_correspondence());
    }
}
