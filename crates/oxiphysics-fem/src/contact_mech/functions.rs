//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{AdhesionModel, NodeToSegmentContact};

/// Complementary error function erfc(x) via rational approximation (Abramowitz & Stegun 7.1.26).
pub(super) fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    poly * (-x * x).exp()
}
/// Distance from a point to a line segment.
pub(super) fn point_segment_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let proj = NodeToSegmentContact::project_onto_segment(p, a, b);
    let d = [p[0] - proj[0], p[1] - proj[1], p[2] - proj[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::contact_mech::*;
    fn steel() -> ElasticSolid {
        ElasticSolid {
            youngs_modulus: 200e9,
            poisson_ratio: 0.3,
        }
    }
    fn hertz_equal() -> HertzContact {
        HertzContact {
            radius_a: 0.01,
            radius_b: 0.01,
            solid_a: steel(),
            solid_b: steel(),
        }
    }
    #[test]
    fn test_contact_radius_positive() {
        let h = hertz_equal();
        let a = h.contact_radius(100.0);
        assert!(
            a > 0.0,
            "contact radius must be positive for positive force"
        );
    }
    #[test]
    fn test_max_pressure_decreases_with_radius() {
        let h_small = HertzContact {
            radius_a: 0.005,
            radius_b: f64::INFINITY,
            solid_a: steel(),
            solid_b: steel(),
        };
        let h_large = HertzContact {
            radius_a: 0.05,
            radius_b: f64::INFINITY,
            solid_a: steel(),
            solid_b: steel(),
        };
        let f = 500.0;
        assert!(
            h_large.max_pressure(f) < h_small.max_pressure(f),
            "larger sphere should give lower peak pressure"
        );
    }
    #[test]
    fn test_effective_radius_equal_spheres() {
        let h = hertz_equal();
        let r_eff = h.effective_radius();
        let expected = 0.01 / 2.0;
        assert!((r_eff - expected).abs() < 1e-12);
    }
    #[test]
    fn test_pressure_distribution_at_centre() {
        let h = hertz_equal();
        let f = 100.0;
        let p0 = h.max_pressure(f);
        let p_centre = h.pressure_distribution(0.0, f);
        assert!((p_centre - p0).abs() / p0 < 1e-10, "p(0) must equal p0");
    }
    #[test]
    fn test_penalty_zero_for_positive_gap() {
        let fn_ = PenaltyContactConstraint::normal_force(1e-3, 1e6);
        assert_eq!(fn_, 0.0, "no force for open gap");
    }
    #[test]
    fn test_penalty_positive_for_negative_gap() {
        let fn_ = PenaltyContactConstraint::normal_force(-1e-4, 1e8);
        assert!(fn_ > 0.0, "penetration must produce positive contact force");
    }
    #[test]
    fn test_project_midpoint() {
        let a = [0.0, 0.0, 0.0_f64];
        let b = [2.0, 0.0, 0.0_f64];
        let p = [1.0, 1.0, 0.0_f64];
        let proj = NodeToSegmentContact::project_onto_segment(p, a, b);
        assert!(
            (proj[0] - 1.0).abs() < 1e-12,
            "x-projection should be midpoint"
        );
        assert!(proj[1].abs() < 1e-12);
        assert!(proj[2].abs() < 1e-12);
    }
    #[test]
    fn test_max_shear_depth() {
        let a = 0.01;
        let z_max = HertzianStressField::max_shear_depth(a);
        assert!((z_max - 0.48 * a).abs() < 1e-15);
    }
    #[test]
    fn test_mean_pressure_relation() {
        let h = hertz_equal();
        let f = 100.0;
        let p0 = h.max_pressure(f);
        let p_mean = h.mean_pressure(f);
        assert!((p_mean - 2.0 / 3.0 * p0).abs() / p0 < 1e-6);
    }
    #[test]
    fn test_contact_area_positive() {
        let h = hertz_equal();
        let area = h.contact_area(100.0);
        assert!(area > 0.0);
    }
    #[test]
    fn test_stored_energy_positive() {
        let h = hertz_equal();
        let u = h.stored_energy(100.0);
        assert!(u > 0.0);
    }
    #[test]
    fn test_force_from_indentation_zero() {
        let h = hertz_equal();
        assert_eq!(h.force_from_indentation(0.0), 0.0);
        assert_eq!(h.force_from_indentation(-0.01), 0.0);
    }
    #[test]
    fn test_force_indentation_roundtrip() {
        let h = hertz_equal();
        let f_in = 100.0;
        let delta = h.hertz_indentation(f_in);
        let f_out = h.force_from_indentation(delta);
        assert!(
            (f_in - f_out).abs() / f_in < 1e-6,
            "Force roundtrip: {f_in} -> delta={delta} -> {f_out}"
        );
    }
    #[test]
    fn test_pressure_profile_length() {
        let h = hertz_equal();
        let profile = h.pressure_profile(100.0, 20);
        assert_eq!(profile.len(), 20);
        assert!(profile[0].0.abs() < 1e-12);
        assert!(profile[0].1 > 0.0);
    }
    #[test]
    fn test_plane_strain_modulus() {
        let s = steel();
        let e_prime = s.plane_strain_modulus();
        let expected = 200e9 / (1.0 - 0.09);
        assert!((e_prime - expected).abs() / expected < 1e-10);
    }
    #[test]
    fn test_mindlin_stick_radius() {
        let h = hertz_equal();
        let mc = MindlinContact {
            hertz: h,
            friction: 0.5,
        };
        let fn_ = 100.0;
        let ft = 10.0;
        let c = mc.stick_radius(ft, fn_).unwrap();
        let a = mc.hertz.contact_radius(fn_);
        assert!(c > 0.0 && c < a, "Stick radius should be between 0 and a");
    }
    #[test]
    fn test_mindlin_stick_radius_gross_slip() {
        let h = hertz_equal();
        let mc = MindlinContact {
            hertz: h,
            friction: 0.5,
        };
        let fn_ = 100.0;
        let ft = 60.0;
        assert!(mc.stick_radius(ft, fn_).is_none());
    }
    #[test]
    fn test_fretting_energy_zero_amplitude() {
        let h = hertz_equal();
        let mc = MindlinContact {
            hertz: h,
            friction: 0.5,
        };
        let energy = mc.fretting_energy_per_cycle(0.0, 100.0);
        assert!(energy.abs() < 1e-20);
    }
    #[test]
    fn test_fretting_energy_positive() {
        let h = hertz_equal();
        let mc = MindlinContact {
            hertz: h,
            friction: 0.5,
        };
        let energy = mc.fretting_energy_per_cycle(10.0, 100.0);
        assert!(energy > 0.0);
    }
    #[test]
    fn test_rolling_resistance() {
        let h = hertz_equal();
        let rc = RollingContact::new(h, 0.01, 0.01);
        let fr = rc.rolling_resistance(100.0);
        assert!((fr - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_rolling_torque() {
        let h = hertz_equal();
        let rc = RollingContact::new(h, 0.01, 0.01);
        let t = rc.rolling_torque(100.0);
        assert!((t - 0.01).abs() < 1e-10);
    }
    #[test]
    fn test_creep_ratio_zero_rolling() {
        assert_eq!(RollingContact::creep_ratio(1.0, 0.0), 0.0);
    }
    #[test]
    fn test_creep_ratio_normal() {
        let cr = RollingContact::creep_ratio(0.1, 10.0);
        assert!((cr - 0.01).abs() < 1e-12);
    }
    #[test]
    fn test_carter_creep_force_direction() {
        let h = hertz_equal();
        let rc = RollingContact::new(h, 0.01, 0.01);
        let f_pos = rc.carter_creep_force(0.001, 100.0, 0.3);
        let f_neg = rc.carter_creep_force(-0.001, 100.0, 0.3);
        assert!(f_pos < 0.0, "Positive creep should give negative force");
        assert!(f_neg > 0.0, "Negative creep should give positive force");
    }
    #[test]
    fn test_contact_patch_hertzian_force() {
        let h = hertz_equal();
        let f_applied = 100.0;
        let a = h.contact_radius(f_applied);
        let p0 = h.max_pressure(f_applied);
        let patch = ContactPatch::hertzian(a, p0, 200);
        let f_computed = patch.total_force();
        assert!(
            (f_computed - f_applied).abs() / f_applied < 0.05,
            "Patch force={f_computed:.3e}, expected {f_applied:.3e}"
        );
    }
    #[test]
    fn test_contact_patch_max_pressure() {
        let patch = ContactPatch::hertzian(0.001, 1e6, 50);
        assert!((patch.max_pressure() - 1e6).abs() < 1e3);
    }
    #[test]
    fn test_fretting_contact_cycles() {
        let h = hertz_equal();
        let mc = MindlinContact {
            hertz: h,
            friction: 0.5,
        };
        let mut fc = FrettingContact::new(mc, 1e-6);
        fc.advance_cycles(1000);
        assert_eq!(fc.cycles, 1000);
        let energy = fc.total_dissipated_energy(100.0);
        assert!(energy >= 0.0);
    }
    #[test]
    fn test_wear_depth() {
        let h = FrettingContact::wear_depth(1e-4, 1e6, 1e-6, 1000, 1e9);
        assert!(h > 0.0);
    }
    #[test]
    fn test_augmented_lagrangian_no_contact() {
        let f = PenaltyContactConstraint::augmented_lagrangian_force(0.0, 0.01, 1e6);
        assert!(f.abs() < 1e-10);
    }
    #[test]
    fn test_augmented_lagrangian_contact() {
        let f = PenaltyContactConstraint::augmented_lagrangian_force(0.0, -0.01, 1e6);
        assert!(f > 0.0);
    }
    #[test]
    fn test_segment_segment_distance() {
        let a1 = [0.0, 0.0, 0.0];
        let a2 = [1.0, 0.0, 0.0];
        let b1 = [0.5, 1.0, 0.0];
        let b2 = [0.5, 2.0, 0.0];
        let d = NodeToSegmentContact::segment_segment_distance(a1, a2, b1, b2);
        assert!((d - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_axial_normal_stress_surface() {
        let p0 = 1e6;
        let a = 0.001;
        let sigma = HertzianStressField::axial_normal_stress(0.0, a, p0);
        assert!((sigma + p0).abs() < 1e-3, "sigma_z(0) = -p0");
    }
    #[test]
    fn test_axial_normal_stress_deep() {
        let p0 = 1e6;
        let a = 0.001;
        let sigma = HertzianStressField::axial_normal_stress(100.0 * a, a, p0);
        assert!(
            sigma.abs() < p0 * 1e-3,
            "Deep stress should be near zero, got {}",
            sigma
        );
    }
    #[test]
    fn test_plasticity_index() {
        let rough = RoughSurfaceContact {
            asperity_density: 1e10,
            asperity_radius: 1e-6,
            asperity_height_std: 1e-7,
        };
        let psi = rough.plasticity_index(200e9, 2e9);
        assert!(psi > 0.0);
    }
    #[test]
    fn test_max_subsurface_von_mises() {
        let vm = HertzianStressField::max_subsurface_von_mises(0.001, 1e6);
        assert!((vm - 0.31e6).abs() < 1.0);
    }
    #[test]
    fn test_line_contact_half_width() {
        let b = RollingContact::line_contact_half_width(1e4, 0.01, 200e9);
        assert!(b > 0.0);
    }
}
/// Compute adhesion energy in the JKR model.
///
/// U_adh = -π * W * a²
#[allow(dead_code)]
pub fn jkr_adhesion_energy(contact_radius: f64, work_of_adhesion: f64) -> f64 {
    -PI * work_of_adhesion * contact_radius * contact_radius
}
/// Compute contact radius as a function of load for different adhesion models.
///
/// Returns (load, contact_radius) pairs for the given load range.
#[allow(dead_code)]
pub fn contact_radius_vs_load(
    r_star: f64,
    e_star: f64,
    work_of_adhesion: f64,
    model: AdhesionModel,
    n_points: usize,
    f_min: f64,
    f_max: f64,
) -> Vec<(f64, f64)> {
    let mut result = Vec::with_capacity(n_points);
    for i in 0..n_points {
        let t = i as f64 / (n_points as f64 - 1.0).max(1.0);
        let f = f_min + t * (f_max - f_min);
        let a = match model {
            AdhesionModel::Hertz => {
                let f_pos = f.max(0.0);
                (3.0 * f_pos * r_star / (4.0 * e_star)).cbrt()
            }
            AdhesionModel::Jkr => {
                let term1 = 3.0 * PI * work_of_adhesion * r_star;
                let under_sqrt = 6.0 * PI * work_of_adhesion * r_star * f + term1 * term1;
                if under_sqrt < 0.0 {
                    0.0
                } else {
                    let a3 = (r_star / e_star) * (f + term1 + under_sqrt.sqrt());
                    if a3 <= 0.0 { 0.0 } else { a3.cbrt() }
                }
            }
            AdhesionModel::Dmt => {
                let f_eff = f + 2.0 * PI * work_of_adhesion * r_star;
                if f_eff <= 0.0 {
                    0.0
                } else {
                    (3.0 * f_eff * r_star / (4.0 * e_star)).cbrt()
                }
            }
        };
        result.push((f, a));
    }
    result
}
#[cfg(test)]
mod tests_adhesion {
    use super::*;
    use crate::contact_mech::*;
    fn steel() -> ElasticSolid {
        ElasticSolid {
            youngs_modulus: 200e9,
            poisson_ratio: 0.3,
        }
    }
    fn hertz_sphere_flat() -> HertzContact {
        HertzContact {
            radius_a: 0.01,
            radius_b: f64::INFINITY,
            solid_a: steel(),
            solid_b: steel(),
        }
    }
    #[test]
    fn test_jkr_pull_off_force_negative() {
        let jkr = JkrContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.1,
        };
        let f_po = jkr.pull_off_force();
        assert!(
            f_po < 0.0,
            "pull-off force must be negative (tensile): {f_po}"
        );
    }
    #[test]
    fn test_jkr_pull_off_magnitude() {
        let jkr = JkrContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.1,
        };
        let r_star = jkr.hertz.effective_radius();
        let expected = -1.5 * PI * 0.1 * r_star;
        let got = jkr.pull_off_force();
        assert!(
            (got - expected).abs() < 1e-10,
            "JKR pull-off: got {got}, expected {expected}"
        );
    }
    #[test]
    fn test_jkr_contact_radius_positive_load() {
        let jkr = JkrContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.1,
        };
        let a = jkr.contact_radius(100.0);
        assert!(
            a > 0.0,
            "JKR contact radius must be positive for positive load: {a}"
        );
    }
    #[test]
    fn test_jkr_contact_radius_zero_load() {
        let jkr = JkrContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.1,
        };
        let a = jkr.contact_radius(0.0);
        assert!(a > 0.0, "JKR: adhesion causes contact at zero load: {a}");
    }
    #[test]
    fn test_jkr_contact_radius_larger_than_hertz() {
        let jkr = JkrContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.5,
        };
        let f = 100.0;
        let a_jkr = jkr.contact_radius(f);
        let a_hertz = jkr.hertz.contact_radius(f);
        assert!(
            a_jkr > a_hertz,
            "JKR a ({a_jkr}) should > Hertz a ({a_hertz}) for positive W"
        );
    }
    #[test]
    fn test_jkr_zero_load_contact_area_positive() {
        let jkr = JkrContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.2,
        };
        let area = jkr.zero_load_contact_area();
        assert!(
            area > 0.0,
            "JKR contact area at zero load must be positive: {area}"
        );
    }
    #[test]
    fn test_jkr_stored_energy_finite() {
        let jkr = JkrContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.1,
        };
        let e = jkr.stored_elastic_energy(50.0);
        assert!(e.is_finite(), "JKR elastic energy must be finite: {e}");
    }
    #[test]
    fn test_dmt_pull_off_force_negative() {
        let dmt = DmtContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.1,
        };
        let f_po = dmt.pull_off_force();
        assert!(f_po < 0.0, "DMT pull-off force must be negative: {f_po}");
    }
    #[test]
    fn test_dmt_pull_off_magnitude() {
        let dmt = DmtContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.1,
        };
        let r_star = dmt.hertz.effective_radius();
        let expected = -2.0 * PI * 0.1 * r_star;
        let got = dmt.pull_off_force();
        assert!(
            (got - expected).abs() < 1e-10,
            "DMT pull-off: got {got}, expected {expected}"
        );
    }
    #[test]
    fn test_dmt_contact_radius_positive_load() {
        let dmt = DmtContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.1,
        };
        let a = dmt.contact_radius(100.0);
        assert!(a > 0.0, "DMT contact radius for positive load: {a}");
    }
    #[test]
    fn test_dmt_larger_than_hertz_same_load() {
        let dmt = DmtContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.3,
        };
        let f = 100.0;
        let a_dmt = dmt.contact_radius(f);
        let a_hertz = dmt.hertz.contact_radius(f);
        assert!(
            a_dmt > a_hertz,
            "DMT radius ({a_dmt}) should > Hertz ({a_hertz})"
        );
    }
    #[test]
    fn test_dmt_jkr_pull_off_ordering() {
        let hertz = hertz_sphere_flat();
        let w = 0.5;
        let jkr = JkrContact {
            hertz: HertzContact {
                radius_a: hertz.radius_a,
                radius_b: hertz.radius_b,
                solid_a: steel(),
                solid_b: steel(),
            },
            work_of_adhesion: w,
        };
        let dmt = DmtContact {
            hertz: HertzContact {
                radius_a: hertz.radius_a,
                radius_b: hertz.radius_b,
                solid_a: steel(),
                solid_b: steel(),
            },
            work_of_adhesion: w,
        };
        let f_jkr = jkr.pull_off_force().abs();
        let f_dmt = dmt.pull_off_force().abs();
        assert!(
            f_jkr < f_dmt,
            "JKR pull-off ({f_jkr}) should < DMT pull-off ({f_dmt})"
        );
    }
    #[test]
    fn test_maugis_parameter_positive() {
        let md = MaugisDugdaleContact {
            effective_radius: 0.01,
            reduced_modulus: 100e9,
            work_of_adhesion: 0.1,
            cohesive_stress: 100e6,
        };
        let lambda = md.maugis_parameter();
        assert!(lambda > 0.0, "Maugis lambda must be positive: {lambda}");
    }
    #[test]
    fn test_maugis_pull_off_between_dmt_jkr() {
        let md = MaugisDugdaleContact {
            effective_radius: 0.01,
            reduced_modulus: 100e9,
            work_of_adhesion: 0.1,
            cohesive_stress: 100e6,
        };
        let f_po = md.pull_off_force();
        let f_dmt = md.dmt_pull_off_force();
        let f_jkr = md.jkr_pull_off_force();
        assert!(f_po.abs() >= f_jkr.abs(), "Maugis |F_po| >= |F_jkr|");
        assert!(f_po.abs() <= f_dmt.abs(), "Maugis |F_po| <= |F_dmt|");
    }
    #[test]
    fn test_maugis_zero_force_radius_positive() {
        let md = MaugisDugdaleContact {
            effective_radius: 0.01,
            reduced_modulus: 100e9,
            work_of_adhesion: 0.1,
            cohesive_stress: 100e6,
        };
        let a0 = md.zero_force_contact_radius();
        assert!(
            a0 > 0.0,
            "Maugis zero-force contact radius must be positive: {a0}"
        );
    }
    #[test]
    fn test_jkr_adhesion_energy_negative() {
        let u = jkr_adhesion_energy(1e-4, 0.5);
        assert!(u < 0.0, "JKR adhesion energy must be negative: {u}");
    }
    #[test]
    fn test_jkr_adhesion_energy_scales_with_area() {
        let a1 = 1e-4;
        let a2 = 2e-4;
        let w = 0.5;
        let u1 = jkr_adhesion_energy(a1, w);
        let u2 = jkr_adhesion_energy(a2, w);
        assert!(
            (u2 / u1 - 4.0).abs() < 1e-9,
            "Adhesion energy scales as a²: ratio={}",
            u2 / u1
        );
    }
    #[test]
    fn test_contact_radius_vs_load_hertz_monotone() {
        let curve = contact_radius_vs_load(0.005, 100e9, 0.1, AdhesionModel::Hertz, 10, 1.0, 100.0);
        for i in 1..curve.len() {
            assert!(
                curve[i].1 >= curve[i - 1].1 - 1e-20,
                "Hertz radius should increase with load: a[{i}]={} < a[{}]={}",
                curve[i].1,
                i - 1,
                curve[i - 1].1
            );
        }
    }
    #[test]
    fn test_contact_radius_vs_load_jkr_larger_than_hertz() {
        let r_star = 0.005;
        let e_star = 100e9;
        let w = 0.5;
        let curve_hertz =
            contact_radius_vs_load(r_star, e_star, w, AdhesionModel::Hertz, 5, 10.0, 100.0);
        let curve_jkr =
            contact_radius_vs_load(r_star, e_star, w, AdhesionModel::Jkr, 5, 10.0, 100.0);
        for i in 0..5 {
            assert!(
                curve_jkr[i].1 >= curve_hertz[i].1 - 1e-20,
                "JKR radius should >= Hertz radius at same load"
            );
        }
    }
    #[test]
    fn test_contact_ld_displacements_positive() {
        let h = hertz_sphere_flat();
        let ld = ContactLoadDisplacement {
            hertz: h,
            loads: vec![10.0, 50.0, 100.0, 200.0],
        };
        for &d in ld.displacements().iter() {
            assert!(d > 0.0, "Hertz indentation must be positive: {d}");
        }
    }
    #[test]
    fn test_contact_ld_elastic_work_positive() {
        let h = hertz_sphere_flat();
        let ld = ContactLoadDisplacement {
            hertz: h,
            loads: vec![0.0, 50.0, 100.0],
        };
        let work = ld.elastic_work();
        assert!(work > 0.0, "Elastic work must be positive: {work}");
    }
    #[test]
    fn test_contact_ld_stiffnesses_positive() {
        let h = hertz_sphere_flat();
        let ld = ContactLoadDisplacement {
            hertz: h,
            loads: vec![10.0, 100.0, 1000.0],
        };
        for &k in ld.stiffnesses().iter() {
            assert!(k > 0.0, "Contact stiffness must be positive: {k}");
        }
    }
    #[test]
    fn test_contact_ld_peak_load() {
        let h = hertz_sphere_flat();
        let ld = ContactLoadDisplacement {
            hertz: h,
            loads: vec![10.0, 200.0, 50.0],
        };
        assert!((ld.peak_load() - 200.0).abs() < 1e-10);
    }
    #[test]
    fn test_extended_gw_adhesive_force_positive() {
        let gw = RoughSurfaceContact {
            asperity_density: 1e10,
            asperity_radius: 1e-6,
            asperity_height_std: 1e-7,
        };
        let egw = ExtendedGwContact {
            gw,
            hardness: 2e9,
            work_of_adhesion: 0.1,
        };
        let f_adh = egw.adhesive_force(0.0, 1e-4);
        assert!(f_adh >= 0.0, "Adhesive force must be non-negative: {f_adh}");
    }
    #[test]
    fn test_extended_gw_friction_coefficient() {
        let gw = RoughSurfaceContact {
            asperity_density: 1e10,
            asperity_radius: 1e-6,
            asperity_height_std: 1e-7,
        };
        let egw = ExtendedGwContact {
            gw,
            hardness: 3e9,
            work_of_adhesion: 0.1,
        };
        let mu = egw.friction_coefficient(1e9);
        assert!(mu > 0.0, "Friction coefficient must be positive: {mu}");
    }
    #[test]
    fn test_dmt_zero_load_area_positive() {
        let dmt = DmtContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.2,
        };
        let area = dmt.zero_load_contact_area();
        assert!(
            area > 0.0,
            "DMT zero-load contact area must be positive: {area}"
        );
    }
    #[test]
    fn test_jkr_indentation_finite() {
        let jkr = JkrContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.1,
        };
        let delta = jkr.indentation(100.0);
        assert!(delta.is_finite(), "JKR indentation must be finite: {delta}");
    }
    #[test]
    fn test_adhesion_model_enum_variants() {
        let f1 = contact_radius_vs_load(0.005, 100e9, 0.1, AdhesionModel::Hertz, 3, 10.0, 100.0);
        let f2 = contact_radius_vs_load(0.005, 100e9, 0.1, AdhesionModel::Jkr, 3, 10.0, 100.0);
        let f3 = contact_radius_vs_load(0.005, 100e9, 0.1, AdhesionModel::Dmt, 3, 10.0, 100.0);
        assert_eq!(f1.len(), 3);
        assert_eq!(f2.len(), 3);
        assert_eq!(f3.len(), 3);
    }
    #[test]
    fn test_maugis_large_lambda_jkr_like() {
        let md = MaugisDugdaleContact {
            effective_radius: 0.01,
            reduced_modulus: 100e9,
            work_of_adhesion: 0.01,
            cohesive_stress: 1e12,
        };
        let f_po = md.pull_off_force();
        let f_jkr = md.jkr_pull_off_force();
        assert!(
            (f_po - f_jkr).abs() / f_jkr.abs() < 0.01,
            "Large lambda → JKR limit: got {f_po}, expected ~{f_jkr}"
        );
    }
    #[test]
    fn test_maugis_small_lambda_dmt_like() {
        let md = MaugisDugdaleContact {
            effective_radius: 0.01,
            reduced_modulus: 100e9,
            work_of_adhesion: 0.1,
            cohesive_stress: 1.0,
        };
        let f_po = md.pull_off_force();
        let f_dmt = md.dmt_pull_off_force();
        assert!(
            (f_po - f_dmt).abs() / f_dmt.abs() < 0.01,
            "Small lambda → DMT limit: got {f_po}, expected ~{f_dmt}"
        );
    }
    #[test]
    fn test_gw_normal_force_increases_with_approach() {
        let gw = RoughSurfaceContact {
            asperity_density: 1e10,
            asperity_radius: 1e-6,
            asperity_height_std: 1e-7,
        };
        let solid = ElasticSolid {
            youngs_modulus: 200e9,
            poisson_ratio: 0.3,
        };
        let f_far = gw.normal_force(1e-6, &solid);
        let f_close = gw.normal_force(0.0, &solid);
        assert!(
            f_close >= f_far,
            "Closer separation gives more force: {f_close} vs {f_far}"
        );
    }
    #[test]
    fn test_gw_real_contact_area_decreases_with_separation() {
        let gw = RoughSurfaceContact {
            asperity_density: 1e10,
            asperity_radius: 1e-6,
            asperity_height_std: 1e-7,
        };
        let r_close = gw.real_contact_area_fraction(0.0);
        let r_far = gw.real_contact_area_fraction(1e-6);
        assert!(
            r_close >= r_far,
            "Closer approach → larger contact fraction"
        );
    }
    #[test]
    fn test_jkr_contact_radius_monotone_increasing() {
        let jkr = JkrContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.1,
        };
        let loads = [10.0_f64, 50.0, 200.0, 500.0];
        let radii: Vec<f64> = loads.iter().map(|&f| jkr.contact_radius(f)).collect();
        for i in 1..radii.len() {
            assert!(
                radii[i] >= radii[i - 1] - 1e-20,
                "JKR radius should increase with load: r[{}]={} < r[{}]={}",
                i,
                radii[i],
                i - 1,
                radii[i - 1]
            );
        }
    }
    #[test]
    fn test_dmt_contact_radius_monotone_increasing() {
        let dmt = DmtContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.1,
        };
        let loads = [10.0_f64, 50.0, 200.0, 500.0];
        let radii: Vec<f64> = loads.iter().map(|&f| dmt.contact_radius(f)).collect();
        for i in 1..radii.len() {
            assert!(
                radii[i] >= radii[i - 1] - 1e-20,
                "DMT radius should increase with load"
            );
        }
    }
    #[test]
    fn test_elastic_energy_increases_with_contact_radius() {
        let jkr = JkrContact {
            hertz: hertz_sphere_flat(),
            work_of_adhesion: 0.05,
        };
        let e_low = jkr.stored_elastic_energy(10.0);
        let e_high = jkr.stored_elastic_energy(1000.0);
        assert!(
            e_high > e_low,
            "Higher load → more stored elastic energy: {e_high} vs {e_low}"
        );
    }
}
