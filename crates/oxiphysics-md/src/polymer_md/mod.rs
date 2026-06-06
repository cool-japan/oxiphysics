// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Polymer MD simulation.
//!
//! Provides bead-spring chains, WLC, Rouse/Zimm models, block copolymer
//! microphase separation, dendrimer generation, and polymer adsorption.

mod analysis;
mod block_copolymers;
mod chain_models;
mod dynamics;
mod melt_network;
mod specialized;

pub use analysis::*;
pub use block_copolymers::*;
pub use chain_models::*;
pub use dynamics::*;
pub use melt_network::*;
pub use specialized::*;

use std::f64::consts::PI;

// --- Helper Functions ---

/// FENE spring force.
///
/// Returns F(r) = -K*r / (1 - (r/R0)²)
pub fn fene_force(r: f64, k: f64, r0: f64) -> f64 {
    let x = r / r0;
    if x >= 1.0 {
        return -1e10;
    }
    -k * r / (1.0 - x * x)
}

/// Rouse mode eigenvalue.
///
/// lambda_p = 2*k*(1-cos(p*pi/N))
pub fn rouse_eigenvalue(p: usize, n: usize, k: f64) -> f64 {
    let angle = p as f64 * PI / n as f64;
    2.0 * k * (1.0 - angle.cos())
}

/// WLC force-extension relation.
///
/// F = kT/p * (1/(4(1-r/L)²) - 1/4 + r/L)
pub fn wlc_extension(r: f64, l: f64, p: f64, kt: f64) -> f64 {
    let x = (r / l).clamp(0.0, 0.9999);
    kt / p * (1.0 / (4.0 * (1.0 - x) * (1.0 - x)) - 0.25 + x)
}

/// Flory-Huggins chi parameter for block copolymer.
///
/// Returns chi*N.
pub fn block_copolymer_chi(chi: f64, n: usize) -> f64 {
    chi * n as f64
}

// --- Vector helpers (shared across submodules) ---

pub(crate) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub(crate) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(crate) fn mag3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

pub(crate) fn norm3(a: [f64; 3]) -> [f64; 3] {
    let m = mag3(a).max(1e-30);
    [a[0] / m, a[1] / m, a[2] / m]
}

pub(crate) fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    mag3(sub3(a, b))
}

/// Compute distance between two 3D points.
pub(crate) fn dist3_vec(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub3_vec(a, b);
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

/// Subtract two 3D vectors.
pub(crate) fn sub3_vec(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Dot product of two 3D vectors.
pub(crate) fn dot3_vec(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(crate) fn sub3_arr(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub(crate) fn dot3_arr(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// ============================================================
// Tests
// ============================================================
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polymer_chain_creation() {
        let chain = PolymerChain::new(10, 1.5, 30.0, 1.0, 1.0, 10.0, 1.0, 0.001);
        assert_eq!(chain.n_beads, 10);
        assert_eq!(chain.positions.len(), 10);
    }

    #[test]
    fn test_polymer_chain_end_to_end() {
        let chain = PolymerChain::new(5, 1.5, 30.0, 1.0, 1.0, 5.0, 1.0, 0.001);
        let r = chain.end_to_end();
        assert!(r > 0.0);
    }

    #[test]
    fn test_polymer_chain_rg() {
        let chain = PolymerChain::new(10, 1.5, 30.0, 1.0, 1.0, 10.0, 1.0, 0.001);
        let rg = chain.radius_of_gyration();
        assert!(rg > 0.0);
    }

    #[test]
    fn test_polymer_chain_step() {
        let mut chain = PolymerChain::new(5, 1.5, 30.0, 0.1, 1.0, 5.0, 1.0, 0.0001);
        chain.step();
        for p in &chain.positions {
            assert!(p[0].is_finite());
        }
    }

    #[test]
    fn test_fene_force_zero() {
        let f = fene_force(0.0, 30.0, 1.5);
        assert!(f.abs() < 1e-10);
    }

    #[test]
    fn test_fene_force_increases() {
        let f1 = fene_force(0.5, 30.0, 1.5).abs();
        let f2 = fene_force(1.0, 30.0, 1.5).abs();
        assert!(f2 > f1);
    }

    #[test]
    fn test_wlc_force_positive() {
        let wlc = WormlikeChain::new(100.0, 50.0, 4.1);
        let f = wlc.force(50.0);
        assert!(f > 0.0);
    }

    #[test]
    fn test_wlc_force_increases() {
        let wlc = WormlikeChain::new(100.0, 50.0, 4.1);
        let f1 = wlc.force(50.0);
        let f2 = wlc.force(80.0);
        assert!(f2 > f1);
    }

    #[test]
    fn test_rouse_eigenvalue_p1() {
        let lam = rouse_eigenvalue(1, 100, 1.0);
        assert!(lam > 0.0);
    }

    #[test]
    fn test_rouse_model_tau() {
        let rouse = RouseModel::new(100, 1.0, 1.0, 1.0);
        let tau = rouse.tau_rouse();
        assert!(tau > 0.0);
    }

    #[test]
    fn test_rouse_mode_amplitude() {
        let rouse = RouseModel::new(100, 1.0, 1.0, 1.0);
        let amp = rouse.mode_amplitude_sq(1);
        assert!(amp > 0.0);
    }

    #[test]
    fn test_zimm_time() {
        let zimm = ZimmModel::new(100, 1.0, 1e-3, 4.1e-21);
        let tau = zimm.zimm_time();
        assert!(tau > 0.0);
    }

    #[test]
    fn test_zimm_diffusion() {
        let zimm = ZimmModel::new(100, 1.0, 1e-3, 4.1e-21);
        let d = zimm.diffusion_coefficient();
        assert!(d > 0.0);
    }

    #[test]
    fn test_block_copolymer_chin() {
        let bcp = BlockCopolymer::new(50, 50, 0.1, 1.0, 300.0);
        let chin = bcp.chi_n();
        assert!((chin - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_block_copolymer_phase_sep() {
        let bcp = BlockCopolymer::new(100, 100, 0.2, 1.0, 300.0);
        assert!(bcp.is_phase_separated()); // chi*N = 40 > 10.5
    }

    #[test]
    fn test_block_copolymer_lamellar() {
        let bcp = BlockCopolymer::new(50, 50, 0.2, 1.0, 300.0);
        let d = bcp.lamellar_spacing();
        assert!(d > 0.0);
    }

    #[test]
    fn test_polymer_crystallization_nucleation() {
        let pc = PolymerCrystallization::new(400.0, 450.0, 1e-2, 1e-3, 8000.0, 1e6, 1500.0, 230.0);
        let rate = pc.nucleation_rate();
        assert!(rate >= 0.0);
    }

    #[test]
    fn test_polymer_crystallization_lamellar() {
        let pc = PolymerCrystallization::new(400.0, 450.0, 1e-2, 1e-3, 8000.0, 1e6, 1500.0, 230.0);
        let l = pc.lamellar_thickness();
        assert!(l > 0.0);
    }

    #[test]
    fn test_polymer_melt_reptation() {
        let melt = PolymerMelt::new(200, 20, 1.0, 1.0, 1.0, 100.0);
        let tau = melt.reptation_time();
        assert!(tau > 0.0);
    }

    #[test]
    fn test_polymer_melt_viscosity() {
        let melt = PolymerMelt::new(200, 20, 1.0, 1.0, 1.0, 100.0);
        let eta = melt.viscosity();
        assert!(eta > 0.0);
    }

    #[test]
    fn test_polymer_network_modulus() {
        let net = PolymerNetwork::new(1000, 500, 1e4, 900.0, 300.0);
        let g = net.modulus();
        assert!(g > 0.0);
    }

    #[test]
    fn test_polymer_network_functionality() {
        let net = PolymerNetwork::new(1000, 500, 1e4, 900.0, 300.0);
        let f = net.functionality();
        assert!((f - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_dendrimer_end_groups() {
        let den = DendriticPolymer::new(3, 3, 2, 5, 1.0);
        assert_eq!(den.n_end_groups(), 3 * 8); // 3 * 2^3
    }

    #[test]
    fn test_dendrimer_radius() {
        let den = DendriticPolymer::new(3, 3, 2, 5, 1.0);
        let r = den.radius();
        assert!(r > 0.0);
    }

    #[test]
    fn test_polymer_adsorption_energy() {
        let ads = PolymerAdsorption::new(100, 2.0, 1.0, 1e-3);
        let e = ads.adsorption_energy();
        assert!(e < 0.0); // favorable
    }

    #[test]
    fn test_polymer_adsorption_langmuir() {
        let ads = PolymerAdsorption::new(100, 2.0, 1.0, 1e-3);
        let theta = ads.langmuir_coverage();
        assert!(theta > 0.0 && theta < 1.0);
    }

    #[test]
    fn test_wlc_energy() {
        let wlc = WormlikeChain::new(100.0, 50.0, 4.1);
        let e = wlc.energy(20.0);
        assert!(e > 0.0);
    }

    #[test]
    fn test_rouse_msd() {
        let rouse = RouseModel::new(100, 1.0, 1.0, 1.0);
        let msd = rouse.msd_center_of_mass(1.0);
        assert!(msd > 0.0);
    }

    #[test]
    fn test_block_copolymer_interfacial_width() {
        let bcp = BlockCopolymer::new(50, 50, 0.2, 1.0, 300.0);
        let w = bcp.interfacial_width();
        assert!(w > 0.0);
    }
}

// ─── Extended tests ───────────────────────────────────────────────────────────

#[cfg(test)]
mod extended_tests {
    use super::*;

    // ── WLC sampler ────────────────────────────────────────────────────────────

    #[test]
    fn wlc_sampler_mean_cos_angle() {
        let wlc = WLCSampler::new(100, 1.0, 10.0, 1.0);
        let cos_a = wlc.mean_cos_angle();
        assert!(cos_a > 0.0 && cos_a < 1.0);
    }

    #[test]
    fn wlc_sampler_tangent_correlation() {
        let wlc = WLCSampler::new(100, 1.0, 10.0, 1.0);
        let c0 = wlc.tangent_correlation(0.0);
        let c5 = wlc.tangent_correlation(5.0);
        assert!((c0 - 1.0).abs() < 1e-10);
        assert!(c5 < c0);
    }

    #[test]
    fn wlc_sampler_mean_r2_positive() {
        let wlc = WLCSampler::new(100, 1.0, 10.0, 1.0);
        assert!(wlc.mean_r2() > 0.0);
    }

    #[test]
    fn wlc_sampler_rg2_positive() {
        let wlc = WLCSampler::new(100, 1.0, 10.0, 1.0);
        assert!(wlc.radius_of_gyration_sq() > 0.0);
    }

    #[test]
    fn wlc_sampler_fit_persistence() {
        let lp = WLCSampler::fit_persistence_length(1.0, 0.9);
        assert!(lp > 0.0 && lp.is_finite());
    }

    // ── FJC ────────────────────────────────────────────────────────────────────

    #[test]
    fn fjc_contour_length() {
        let fjc = FreelyJointedChain::new(100, 1.5, 1.0);
        assert!((fjc.contour_length() - 150.0).abs() < 1e-10);
    }

    #[test]
    fn fjc_mean_r2() {
        let fjc = FreelyJointedChain::new(100, 1.5, 1.0);
        let expected = 100.0 * 1.5 * 1.5;
        assert!((fjc.mean_r2() - expected).abs() < 1e-10);
    }

    #[test]
    fn fjc_langevin_small_x() {
        let l = FreelyJointedChain::langevin(0.0001);
        assert!((l - 0.0001 / 3.0).abs() < 1e-8);
    }

    #[test]
    fn fjc_inverse_langevin_zero() {
        let l = FreelyJointedChain::inverse_langevin(0.0);
        assert!(l.abs() < 1e-10);
    }

    #[test]
    fn fjc_force_extension_positive() {
        let fjc = FreelyJointedChain::new(100, 1.0, 1.0);
        let f = fjc.force_extension(50.0);
        assert!(f > 0.0);
    }

    #[test]
    fn fjc_entropic_spring() {
        let fjc = FreelyJointedChain::new(100, 1.0, 1.0);
        let k = fjc.entropic_spring_constant();
        assert!((k - 0.03).abs() < 1e-10);
    }

    // ── FENEBond ──────────────────────────────────────────────────────────────

    #[test]
    fn fene_bond_potential_zero_at_origin() {
        let bond = FENEBond::new(30.0, 1.5, 0.96);
        assert!(bond.potential(0.0).abs() < 1e-10);
    }

    #[test]
    fn fene_bond_potential_diverges() {
        let bond = FENEBond::new(30.0, 1.5, 0.96);
        let u = bond.potential(1.5);
        assert!(u == f64::INFINITY);
    }

    #[test]
    fn fene_bond_force_positive() {
        let bond = FENEBond::new(30.0, 1.5, 0.96);
        let f = bond.force_magnitude(0.5);
        assert!(f > 0.0);
    }

    #[test]
    fn fene_wca_combined() {
        let bond = FENEBond::new(30.0, 1.5, 0.96);
        let u = bond.fene_wca_potential(1.0, 1.0, 1.0);
        assert!(u.is_finite());
    }

    // ── RouseChain ────────────────────────────────────────────────────────────

    #[test]
    fn rouse_chain_spring_force_equilibrium() {
        let chain = RouseChain::new(3, 1.0, 1.0, 1.0, 0.01);
        // Evenly spaced: force on middle bead should be 0
        let f = chain.spring_force(1);
        assert!(f.abs() < 1e-10);
    }

    #[test]
    fn rouse_chain_normal_mode_p0() {
        let chain = RouseChain::new(5, 1.0, 1.0, 1.0, 0.01);
        let x0 = chain.normal_mode(0);
        // p=0 mode is center-of-mass (should be non-zero for shifted chain)
        assert!(x0.is_finite());
    }

    #[test]
    fn rouse_chain_tau_p_decreases() {
        let chain = RouseChain::new(50, 1.0, 1.0, 1.0, 0.01);
        let tau1 = chain.tau_p(1);
        let tau2 = chain.tau_p(2);
        assert!(tau1 > tau2);
    }

    #[test]
    fn rouse_chain_deterministic_step() {
        let mut chain = RouseChain::new(5, 1.0, 1.0, 1.0, 0.01);
        chain.deterministic_step();
        for &x in &chain.positions {
            assert!(x.is_finite());
        }
    }

    // ── TubeModel ─────────────────────────────────────────────────────────────

    #[test]
    fn tube_model_z_entanglements() {
        let tube = TubeModel::new(1000, 50, 5.0, 1.0, 1.0, 1.0);
        assert!((tube.z_entanglements() - 20.0).abs() < 1e-10);
    }

    #[test]
    fn tube_model_tau_d_positive() {
        let tube = TubeModel::new(1000, 50, 5.0, 1.0, 1.0, 1.0);
        assert!(tube.tau_d() > 0.0);
    }

    #[test]
    fn tube_model_tau_d_gt_tau_rouse() {
        let tube = TubeModel::new(1000, 50, 5.0, 1.0, 1.0, 1.0);
        assert!(tube.tau_d() > tube.tau_rouse());
    }

    #[test]
    fn tube_model_diffusion_positive() {
        let tube = TubeModel::new(1000, 50, 5.0, 1.0, 1.0, 1.0);
        assert!(tube.diffusion_coefficient() > 0.0);
    }

    #[test]
    fn tube_model_msd_increases() {
        let tube = TubeModel::new(1000, 50, 5.0, 1.0, 1.0, 1.0);
        let m1 = tube.msd_regimes(0.001);
        let m2 = tube.msd_regimes(1.0);
        assert!(m2 > m1);
    }

    // ── RingPolymer ───────────────────────────────────────────────────────────

    #[test]
    fn ring_polymer_omega_positive() {
        let rp = RingPolymer::new(32, 1.0, 1.0, 1.0);
        assert!(rp.omega_p() > 0.0);
    }

    #[test]
    fn ring_polymer_centroid_zero() {
        let rp = RingPolymer::new(32, 1.0, 1.0, 1.0);
        assert!(rp.centroid().abs() < 1e-10);
    }

    #[test]
    fn ring_polymer_internal_energy_zero() {
        let rp = RingPolymer::new(32, 1.0, 1.0, 1.0);
        // All beads at same position -> zero internal energy
        assert!(rp.internal_energy().abs() < 1e-10);
    }

    // ── FloryHuggins ──────────────────────────────────────────────────────────

    #[test]
    fn flory_huggins_phi_critical_symmetric() {
        let fh = FloryHuggins::new(100, 100, 0.1);
        let pc = fh.phi_critical();
        assert!((pc - 0.5).abs() < 1e-10);
    }

    #[test]
    fn flory_huggins_free_energy_minimum() {
        let fh = FloryHuggins::new(100, 100, 0.0);
        let f50 = fh.free_energy(0.5);
        let f20 = fh.free_energy(0.2);
        assert!(f50 < f20); // mixed state favorable when chi=0
    }

    #[test]
    fn flory_huggins_spinodal_positive() {
        let fh = FloryHuggins::new(100, 100, 0.1);
        assert!(fh.chi_spinodal() > 0.0);
    }

    // ── PolymerBrush ──────────────────────────────────────────────────────────

    #[test]
    fn brush_height_positive() {
        let brush = PolymerBrush::new(100, 1.0, 1.0, 0.1, 1.0);
        assert!(brush.brush_height() > 0.0);
    }

    #[test]
    fn brush_compression_force_zero_above_height() {
        let brush = PolymerBrush::new(100, 1.0, 1.0, 0.1, 1.0);
        let h = brush.brush_height();
        assert!(brush.compression_force(h + 1.0).abs() < 1e-10);
    }

    #[test]
    fn brush_free_energy_positive() {
        let brush = PolymerBrush::new(100, 1.0, 1.0, 0.1, 1.0);
        assert!(brush.free_energy_per_chain() > 0.0);
    }

    // ── Polyelectrolyte ───────────────────────────────────────────────────────

    #[test]
    fn polyelectrolyte_manning_positive() {
        let pe = Polyelectrolyte::new(100, 0.17, 0.71, 1, -1, 0.1, 4.1);
        assert!(pe.manning_parameter() > 0.0);
    }

    #[test]
    fn polyelectrolyte_condensed_fraction_bounded() {
        let pe = Polyelectrolyte::new(100, 0.17, 0.71, 1, -1, 0.1, 4.1);
        let theta = pe.condensed_fraction();
        assert!((0.0..=1.0).contains(&theta));
    }

    #[test]
    fn polyelectrolyte_debye_length_positive() {
        let pe = Polyelectrolyte::new(100, 0.17, 0.71, 1, -1, 0.1, 4.1);
        assert!(pe.debye_length() > 0.0);
    }

    // ── GoModel ───────────────────────────────────────────────────────────────

    #[test]
    fn go_model_native_contacts() {
        let positions = vec![[0.0, 0.0, 0.0], [3.8, 0.0, 0.0], [7.6, 0.0, 0.0]];
        let go = GoModel::new(positions, 1.0, 10.0);
        assert!(go.fraction_native_contacts() >= 0.0);
    }

    #[test]
    fn go_model_rmsd_zero_for_native() {
        let positions = vec![[0.0, 0.0, 0.0], [3.8, 0.0, 0.0], [7.6, 0.0, 0.0]];
        let go = GoModel::new(positions, 1.0, 10.0);
        assert!(go.rmsd_from_native() < 1e-10);
    }

    // ── ENM ───────────────────────────────────────────────────────────────────

    #[test]
    fn enm_contacts_positive() {
        let pos = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let enm = ElasticNetworkModel::new(pos, 1.0, 8.0);
        assert!(enm.n_contacts() > 0);
    }

    #[test]
    fn enm_b_factor_positive() {
        let pos = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let enm = ElasticNetworkModel::new(pos, 1.0, 8.0);
        assert!(enm.b_factor(0) > 0.0);
    }

    // ── DebyeHuckel ───────────────────────────────────────────────────────────

    #[test]
    fn debye_huckel_repulsive() {
        let dh = DebyeHuckel::new(0.7, 1.0);
        let u = dh.potential(1.0, 1.0, 1.0);
        assert!(u > 0.0);
    }

    #[test]
    fn debye_huckel_force_positive() {
        let dh = DebyeHuckel::new(0.7, 1.0);
        let f = dh.force_magnitude(1.0, 1.0, 1.0);
        assert!(f > 0.0);
    }

    // ── SelfAvoidingWalk ──────────────────────────────────────────────────────

    #[test]
    fn saw_flory_radius_increases() {
        let s1 = SelfAvoidingWalk::new(100, 1.0);
        let s2 = SelfAvoidingWalk::new(200, 1.0);
        assert!(s2.flory_radius() > s1.flory_radius());
    }

    #[test]
    fn saw_rg_positive() {
        let s = SelfAvoidingWalk::new(100, 1.0);
        assert!(s.radius_of_gyration() > 0.0);
    }

    // ── StarPolymer ───────────────────────────────────────────────────────────

    #[test]
    fn star_polymer_radius_increases_with_f() {
        let s1 = StarPolymer::new(100, 3, 1.0, 1.0);
        let s2 = StarPolymer::new(100, 10, 1.0, 1.0);
        assert!(s2.star_radius() > s1.star_radius());
    }

    #[test]
    fn star_polymer_n_total() {
        let s = StarPolymer::new(100, 5, 1.0, 1.0);
        assert_eq!(s.n_total(), 500);
    }

    // ── CombPolymer ───────────────────────────────────────────────────────────

    #[test]
    fn comb_polymer_n_total() {
        let c = CombPolymer::new(200, 20, 10, 1.0, 10);
        assert_eq!(c.n_total(), 400);
    }

    #[test]
    fn comb_polymer_stretching_gt_one() {
        let c = CombPolymer::new(200, 20, 10, 1.0, 10);
        assert!(c.backbone_stretching() > 1.0);
    }

    // ── LangevinIntegrator ────────────────────────────────────────────────────

    #[test]
    fn langevin_kinetic_energy_zero() {
        let li = LangevinIntegrator::new(1.0, 1.0, 0.001, vec![1.0; 3]);
        let v = vec![[0.0; 3]; 3];
        assert!(li.kinetic_energy(&v) < 1e-10);
    }

    #[test]
    fn langevin_temperature_zero() {
        let li = LangevinIntegrator::new(1.0, 1.0, 0.001, vec![1.0; 3]);
        let v = vec![[0.0; 3]; 3];
        assert!(li.temperature(&v) < 1e-10);
    }
}

// ─── Tests for new polymer MD features ────────────────────────────────────────

#[cfg(test)]
mod polymer_extended_tests {
    use super::*;

    // --- KremerGrestBond ---

    #[test]
    fn test_kg_fene_potential_zero() {
        let kg = KremerGrestBond::new(30.0, 1.5, 1.0, 1.0);
        assert!((kg.fene_potential(0.0)).abs() < 1e-10);
    }

    #[test]
    fn test_kg_wca_potential_zero_outside_cutoff() {
        let kg = KremerGrestBond::new(30.0, 1.5, 1.0, 1.0);
        let r_cut = 2.0_f64.powf(1.0 / 6.0);
        assert!((kg.wca_potential(r_cut + 0.01)).abs() < 1e-10);
    }

    #[test]
    fn test_kg_total_potential_finite() {
        let kg = KremerGrestBond::new(30.0, 1.5, 1.0, 1.0);
        let u = kg.total_potential(1.0);
        assert!(u.is_finite());
    }

    #[test]
    fn test_kg_fene_force_positive() {
        let kg = KremerGrestBond::new(30.0, 1.5, 1.0, 1.0);
        let f = kg.fene_force_magnitude(0.5);
        assert!(f > 0.0);
    }

    #[test]
    fn test_kg_equilibrium_approx() {
        let kg = KremerGrestBond::new(30.0, 1.5, 1.0, 1.0);
        let r_eq = kg.equilibrium_approx();
        assert!((r_eq - 2.0_f64.powf(1.0 / 6.0)).abs() < 1e-10);
    }

    #[test]
    fn test_kg_not_over_extended() {
        let kg = KremerGrestBond::new(30.0, 1.5, 1.0, 1.0);
        assert!(!kg.is_over_extended(1.0));
    }

    #[test]
    fn test_kg_over_extended() {
        let kg = KremerGrestBond::new(30.0, 1.5, 1.0, 1.0);
        assert!(kg.is_over_extended(1.4));
    }

    // --- ShakeConstraints ---

    #[test]
    fn test_shake_no_violation_for_correct_bonds() {
        let shake = ShakeConstraints::new(3, 1.0, 1e-6, 100);
        let pos = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let max_v = shake.max_violation(&pos);
        assert!(max_v < 1e-10);
    }

    #[test]
    fn test_shake_apply_corrects_bonds() {
        let shake = ShakeConstraints::new(3, 1.0, 1e-6, 100);
        let old = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut new_pos = vec![[0.0, 0.0, 0.0], [1.1, 0.0, 0.0], [2.2, 0.0, 0.0]];
        let _iters = shake.apply(&mut new_pos, &old);
        let max_v = shake.max_violation(&new_pos);
        assert!(max_v < 0.01, "max_v = {max_v}");
    }

    // --- PolymerPressureTensor ---

    #[test]
    fn test_pressure_tensor_kinetic() {
        let pt = PolymerPressureTensor::new(100, 1000.0, 1.0);
        assert!(pt.kinetic_pressure() > 0.0);
    }

    #[test]
    fn test_pressure_tensor_isotropic() {
        let pt = PolymerPressureTensor::new(100, 1000.0, 1.0);
        let virial = [0.1, 0.0, 0.0, 0.1, 0.0, 0.1];
        let p = pt.isotropic_pressure(&virial);
        assert!(p > 0.0);
    }

    // --- EndToEndCorrelation ---

    #[test]
    fn test_e2e_correlation_self() {
        let mut e2e = EndToEndCorrelation::new(100, 1.0, 1.0, 1000);
        e2e.record([1.0, 0.0, 0.0]);
        e2e.record([1.0, 0.0, 0.0]);
        let c0 = e2e.correlation(0);
        assert!((c0 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_e2e_mean_r2() {
        let mut e2e = EndToEndCorrelation::new(100, 1.0, 1.0, 1000);
        e2e.record([3.0, 4.0, 0.0]);
        assert!((e2e.mean_r2() - 25.0).abs() < 1e-10);
    }

    // --- RadiusOfGyrationTracker ---

    #[test]
    fn test_rg_tracker_compute_rg() {
        let pos: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let rg = RadiusOfGyrationTracker::compute_rg(&pos);
        assert!(rg > 0.0);
    }

    #[test]
    fn test_rg_tracker_mean() {
        let mut tracker = RadiusOfGyrationTracker::new(10, 1.0, 0.5);
        let pos: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        tracker.record(&pos);
        tracker.record(&pos);
        assert!(tracker.mean_rg() > 0.0);
    }

    #[test]
    fn test_rg_tracker_variance_zero_identical() {
        let mut tracker = RadiusOfGyrationTracker::new(10, 1.0, 0.5);
        let pos: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        tracker.record(&pos);
        tracker.record(&pos);
        tracker.record(&pos);
        assert!(tracker.variance_rg() < 1e-10);
    }

    // --- ChainDiffusion ---

    #[test]
    fn test_chain_diffusion_rouse_positive() {
        let cd = ChainDiffusion::new(100, 1.0, 1.0, 1e-3, 1.0);
        assert!(cd.rouse_diffusion() > 0.0);
    }

    #[test]
    fn test_chain_diffusion_zimm_positive() {
        let cd = ChainDiffusion::new(100, 1.0, 1.0, 1e-3, 1.0);
        assert!(cd.zimm_diffusion() > 0.0);
    }

    #[test]
    fn test_chain_diffusion_msd_increases() {
        let cd = ChainDiffusion::new(100, 1.0, 1.0, 1e-3, 1.0);
        assert!(cd.cm_msd_rouse(2.0) > cd.cm_msd_rouse(1.0));
    }

    // --- ReptationDynamics ---

    #[test]
    fn test_reptation_tau_d0_positive() {
        let rd = ReptationDynamics::new(200, 20, 1.0, 1.0, 1.0, 0.5);
        assert!(rd.tau_d0() > 0.0);
    }

    #[test]
    fn test_reptation_clf_lt_one() {
        let rd = ReptationDynamics::new(200, 20, 1.0, 1.0, 1.0, 0.5);
        assert!(rd.clf_correction() < 1.0);
    }

    #[test]
    fn test_reptation_tau_d_clf_lt_tau_d0() {
        let rd = ReptationDynamics::new(200, 20, 1.0, 1.0, 1.0, 0.5);
        assert!(rd.tau_d_clf() < rd.tau_d0());
    }

    #[test]
    fn test_reptation_viscosity_positive() {
        let rd = ReptationDynamics::new(200, 20, 1.0, 1.0, 1.0, 0.5);
        assert!(rd.viscosity_zero_shear() > 0.0);
    }

    #[test]
    fn test_reptation_tube_survival_prob_one_at_zero() {
        let rd = ReptationDynamics::new(200, 20, 1.0, 1.0, 1.0, 0.5);
        let p = rd.tube_survival_prob(0.0, 50);
        assert!((p - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_reptation_tube_survival_decays() {
        let rd = ReptationDynamics::new(200, 20, 1.0, 1.0, 1.0, 0.5);
        let p0 = rd.tube_survival_prob(0.0, 50);
        let p1 = rd.tube_survival_prob(rd.tau_d_cr(), 50);
        assert!(p1 < p0);
    }

    // --- BondBreakingKinetics ---

    #[test]
    fn test_bell_rate_increases_with_force() {
        let bk = BondBreakingKinetics::new(1.0, 0.5, 1.0);
        let k0 = bk.k_off_bell(0.0);
        let k1 = bk.k_off_bell(5.0);
        assert!(k1 > k0);
    }

    #[test]
    fn test_survival_prob_one_at_zero() {
        let bk = BondBreakingKinetics::new(1.0, 0.5, 1.0);
        assert!((bk.survival_probability(0.0, 0.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_survival_prob_decreases() {
        let bk = BondBreakingKinetics::new(1.0, 0.5, 1.0);
        let s0 = bk.survival_probability(1.0, 0.0);
        let s1 = bk.survival_probability(1.0, 1.0);
        assert!(s1 < s0);
    }

    #[test]
    fn test_mean_rupture_force_positive() {
        let bk = BondBreakingKinetics::new(1.0, 0.5, 1.0);
        let f = bk.mean_rupture_force(100.0);
        assert!(f > 0.0);
    }

    // --- KratkyPorodChain ---

    #[test]
    fn test_kp_mean_r2_positive() {
        let kp = KratkyPorodChain::new(100.0, 10.0, 1.0, 100);
        assert!(kp.mean_r2() > 0.0);
    }

    #[test]
    fn test_kp_rg2_positive() {
        let kp = KratkyPorodChain::new(100.0, 10.0, 1.0, 100);
        assert!(kp.radius_of_gyration_sq() > 0.0);
    }

    #[test]
    fn test_kp_force_extension_increases() {
        let kp = KratkyPorodChain::new(100.0, 10.0, 1.0, 100);
        let f1 = kp.force_extension(20.0);
        let f2 = kp.force_extension(80.0);
        assert!(f2 > f1);
    }

    #[test]
    fn test_kp_tangent_correlation_decays() {
        let kp = KratkyPorodChain::new(100.0, 10.0, 1.0, 100);
        let c0 = kp.tangent_correlation(0.0);
        let c10 = kp.tangent_correlation(10.0);
        assert!((c0 - 1.0).abs() < 1e-10);
        assert!(c10 < c0);
    }

    #[test]
    fn test_kp_regime_classification() {
        let rod = KratkyPorodChain::new(1.0, 100.0, 1.0, 100);
        assert_eq!(rod.regime(), "rod");
        let coil = KratkyPorodChain::new(1000.0, 1.0, 1.0, 100);
        assert_eq!(coil.regime(), "coil");
    }

    // --- RouseModeTracker ---

    #[test]
    fn test_rouse_tracker_lambda_p_positive() {
        let rt = RouseModeTracker::new(50, 1.0, 1.0, 1.0, 0.01, 100);
        assert!(rt.lambda_p(1) > 0.0);
    }

    #[test]
    fn test_rouse_tracker_tau_p_positive() {
        let rt = RouseModeTracker::new(50, 1.0, 1.0, 1.0, 0.01, 100);
        assert!(rt.tau_p(1) > 0.0);
    }

    #[test]
    fn test_rouse_tracker_relax_all() {
        let mut rt = RouseModeTracker::new(5, 1.0, 1.0, 1.0, 0.01, 100);
        rt.amplitudes[1] = 1.0;
        rt.relax_all_modes();
        assert!(rt.amplitudes[1].abs() < 1.0); // decayed
    }

    // --- ZimmHydrodynamics ---

    #[test]
    fn test_zimm_hydro_oseen_element_zero_diagonal() {
        let zh = ZimmHydrodynamics::new(10, 1.0, 1e-3, 1.0);
        assert!(zh.oseen_element(3, 3).abs() < 1e-10);
    }

    #[test]
    fn test_zimm_hydro_diffusion_positive() {
        let zh = ZimmHydrodynamics::new(10, 1.0, 1e-3, 1.0);
        assert!(zh.diffusion_coefficient() > 0.0);
    }

    #[test]
    fn test_zimm_hydro_zimm_time_positive() {
        let zh = ZimmHydrodynamics::new(10, 1.0, 1e-3, 1.0);
        assert!(zh.zimm_time() > 0.0);
    }

    #[test]
    fn test_zimm_hydro_intrinsic_viscosity_positive() {
        let zh = ZimmHydrodynamics::new(10, 1.0, 1e-3, 1.0);
        assert!(zh.intrinsic_viscosity() > 0.0);
    }
}

// ─── Tests for new structures ─────────────────────────────────────────────────

#[cfg(test)]
mod polymer_new_struct_tests {
    use super::*;

    // ─── RadiusOfGyration tests ───────────────────────────────────────────────

    #[test]
    fn test_rg_single_bead_zero() {
        let rg = RadiusOfGyration::new(vec![[1.0, 2.0, 3.0]]);
        assert!(rg.rg_squared().abs() < 1e-12, "Single bead Rg² should be 0");
    }

    #[test]
    fn test_rg_two_beads_symmetric() {
        // Two beads at ±1 along x: COM = 0, Rg² = 1
        let rg = RadiusOfGyration::new(vec![[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let rg2 = rg.rg_squared();
        assert!((rg2 - 1.0).abs() < 1e-10, "rg2={rg2}");
    }

    #[test]
    fn test_rg_linear_chain_positive() {
        let pos: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let rg = RadiusOfGyration::new(pos);
        assert!(rg.rg() > 0.0);
    }

    #[test]
    fn test_rg_increases_with_chain_length() {
        let pos1: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let pos2: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let rg1 = RadiusOfGyration::new(pos1).rg();
        let rg2 = RadiusOfGyration::new(pos2).rg();
        assert!(rg2 > rg1, "Longer chain → larger Rg");
    }

    #[test]
    fn test_rg_principal_radii_sorted() {
        let pos: Vec<[f64; 3]> = (0..8)
            .map(|i| {
                let t = i as f64 * 0.3;
                [t, t * 0.5, t * 0.2]
            })
            .collect();
        let rg = RadiusOfGyration::new(pos);
        let lambdas = rg.principal_radii_squared();
        assert!(lambdas[0] <= lambdas[1] + 1e-12);
        assert!(lambdas[1] <= lambdas[2] + 1e-12);
    }

    #[test]
    fn test_rg_principal_radii_sum_equals_rg2() {
        let pos: Vec<[f64; 3]> = (0..6)
            .map(|i| {
                let t = i as f64;
                [t, t * 0.7, t * 0.3]
            })
            .collect();
        let rg = RadiusOfGyration::new(pos);
        let rg2 = rg.rg_squared();
        let lambdas = rg.principal_radii_squared();
        let sum = lambdas[0] + lambdas[1] + lambdas[2];
        assert!((sum - rg2).abs() < 1e-8, "sum={sum}, rg2={rg2}");
    }

    #[test]
    fn test_rg_relative_shape_anisotropy_range() {
        let pos: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let rg = RadiusOfGyration::new(pos);
        let kappa2 = rg.relative_shape_anisotropy();
        // Rod-like → κ² close to 1, sphere → κ² = 0
        assert!((0.0..=1.0 + 1e-6).contains(&kappa2), "kappa2={kappa2}");
    }

    #[test]
    fn test_rg_rod_high_anisotropy() {
        // Linear chain along x should be highly anisotropic
        let pos: Vec<[f64; 3]> = (0..20).map(|i| [i as f64, 0.0, 0.0]).collect();
        let rg = RadiusOfGyration::new(pos);
        let kappa2 = rg.relative_shape_anisotropy();
        assert!(
            kappa2 > 0.5,
            "Rod should have high anisotropy, kappa2={kappa2}"
        );
    }

    #[test]
    fn test_rg_asphericity_rod_positive() {
        let pos: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let rg = RadiusOfGyration::new(pos);
        let b = rg.asphericity();
        assert!(b > 0.0, "Rod asphericity should be positive, b={b}");
    }

    #[test]
    fn test_rg_acylindricity_nonneg() {
        let pos: Vec<[f64; 3]> = (0..6)
            .map(|i| {
                let t = i as f64;
                [t, t * 0.3, 0.0]
            })
            .collect();
        let rg = RadiusOfGyration::new(pos);
        let c = rg.acylindricity();
        assert!(c >= 0.0, "Acylindricity should be non-negative, c={c}");
    }

    #[test]
    fn test_rg_center_of_mass_correct() {
        let pos = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let rg = RadiusOfGyration::new(pos);
        let com = rg.center_of_mass();
        assert!((com[0] - 1.0).abs() < 1e-12);
        assert!(com[1].abs() < 1e-12);
        assert!(com[2].abs() < 1e-12);
    }

    #[test]
    fn test_rg_gyration_tensor_trace_equals_rg2() {
        let pos: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, (i as f64) * 0.5, 0.0]).collect();
        let rg = RadiusOfGyration::new(pos);
        let s = rg.gyration_tensor();
        let trace = s[0] + s[3] + s[5];
        let rg2 = rg.rg_squared();
        assert!((trace - rg2).abs() < 1e-10, "trace={trace}, rg2={rg2}");
    }

    // ─── PolymerEntanglement tests ────────────────────────────────────────────

    #[test]
    fn test_entanglement_ne_positive() {
        let pe = PolymerEntanglement::new(200, 10, 1.0, 0.85, 1.0);
        assert!(pe.ne > 0.0, "Ne should be positive");
    }

    #[test]
    fn test_entanglement_tube_diameter_positive() {
        let pe = PolymerEntanglement::new(200, 10, 1.0, 0.85, 1.0);
        assert!(pe.tube_diameter() > 0.0);
    }

    #[test]
    fn test_entanglement_z_entanglements() {
        let pe = PolymerEntanglement::with_ne(200, 10, 1.0, 20.0, 1.0, 0.85);
        let z = pe.z_entanglements();
        assert!(
            (z - 10.0).abs() < 1e-10,
            "Z should be N/Ne = 200/20 = 10, z={z}"
        );
    }

    #[test]
    fn test_entanglement_tube_diameter_formula() {
        let pe = PolymerEntanglement::with_ne(100, 5, 1.0, 25.0, 1.0, 0.85);
        let d_t = pe.tube_diameter();
        let expected = 1.0 * 25.0_f64.sqrt();
        assert!(
            (d_t - expected).abs() < 1e-10,
            "d_T=b*sqrt(Ne), d_T={d_t}, expected={expected}"
        );
    }

    #[test]
    fn test_entanglement_path_length_positive() {
        let pe = PolymerEntanglement::new(100, 5, 1.0, 0.85, 1.0);
        assert!(pe.primitive_path_length() > 0.0);
    }

    #[test]
    fn test_entanglement_z1_min_path_straight() {
        let positions: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let l = PolymerEntanglement::z1_min_path(&positions);
        assert!(
            (l - 4.0).abs() < 1e-10,
            "Straight chain min path = 4, l={l}"
        );
    }

    #[test]
    fn test_entanglement_z1_min_path_single() {
        let l = PolymerEntanglement::z1_min_path(&[[0.0, 0.0, 0.0]]);
        assert!((l).abs() < 1e-10);
    }

    #[test]
    fn test_entanglement_rouse_time_positive() {
        let pe = PolymerEntanglement::new(200, 10, 1.0, 0.85, 1.0);
        let tau_e = pe.rouse_time_entanglement_strand(1.0, 1.0);
        assert!(tau_e > 0.0);
    }

    #[test]
    fn test_entanglement_reptation_time_greater_rouse() {
        let pe = PolymerEntanglement::with_ne(200, 10, 1.0, 20.0, 1.0, 0.85);
        let tau_e = pe.rouse_time_entanglement_strand(1.0, 1.0);
        let tau_d = pe.reptation_time(1.0, 1.0);
        assert!(tau_d > tau_e, "Reptation time > entanglement Rouse time");
    }

    #[test]
    fn test_entanglement_path_reduction_le_one() {
        let pe = PolymerEntanglement::new(100, 5, 1.0, 0.85, 1.0);
        let ratio = pe.path_reduction_ratio();
        assert!(
            ratio <= 1.0 + 1e-10,
            "Primitive path ≤ full contour, ratio={ratio}"
        );
    }

    #[test]
    fn test_entanglement_plateau_modulus_positive() {
        let pe = PolymerEntanglement::new(200, 10, 1.0, 0.85, 1.0);
        assert!(pe.plateau_modulus > 0.0);
    }

    // ─── DendrimerModel tests ─────────────────────────────────────────────────

    #[test]
    fn test_dendrimer_terminal_groups_g0() {
        let d = DendrimerModel::new(0, 3, 2, 3, 1.0, 1.0);
        // Generation 0: f × b^0 = 3 × 1 = 3
        assert_eq!(d.n_terminal_groups(), 3);
    }

    #[test]
    fn test_dendrimer_terminal_groups_g2_pamam() {
        let d = DendrimerModel::pamam(2, 1.0, 1.0);
        // f=3, b=2, g=2: 3 × 2^2 = 12
        assert_eq!(d.n_terminal_groups(), 12);
    }

    #[test]
    fn test_dendrimer_total_monomers_positive() {
        let d = DendrimerModel::pamam(3, 1.0, 1.0);
        assert!(d.n_total_monomers() > 0);
    }

    #[test]
    fn test_dendrimer_mass_positive() {
        let d = DendrimerModel::pamam(4, 1.0, 1.0);
        assert!(d.mass() > 0.0);
    }

    #[test]
    fn test_dendrimer_rg_increases_with_generation() {
        let d1 = DendrimerModel::pamam(2, 1.0, 1.0);
        let d2 = DendrimerModel::pamam(4, 1.0, 1.0);
        assert!(d2.radius_of_gyration() > d1.radius_of_gyration());
    }

    #[test]
    fn test_dendrimer_hydrodynamic_radius_gt_rg() {
        let d = DendrimerModel::pamam(3, 1.0, 1.0);
        assert!(d.hydrodynamic_radius() >= d.radius_of_gyration());
    }

    #[test]
    fn test_dendrimer_fractal_dimension_positive() {
        let d = DendrimerModel::pamam(4, 1.0, 1.0);
        let df = d.fractal_dimension();
        assert!(df > 0.0 && df.is_finite(), "Df={df}");
    }

    #[test]
    fn test_dendrimer_intrinsic_viscosity_positive() {
        let d = DendrimerModel::pamam(3, 1.0, 1.0);
        assert!(d.intrinsic_viscosity() > 0.0);
    }

    #[test]
    fn test_dendrimer_diffusion_coefficient_positive() {
        let d = DendrimerModel::pamam(3, 1.0, 1.0);
        let dc = d.diffusion_coefficient(1e-3);
        assert!(dc > 0.0);
    }

    #[test]
    fn test_dendrimer_arm_end_to_end_distance() {
        let d = DendrimerModel::new(3, 3, 2, 4, 1.0, 1.0);
        let r_arm = d.arm_end_to_end_distance();
        let expected = 1.0 * (4.0_f64).sqrt();
        assert!((r_arm - expected).abs() < 1e-10, "r_arm={r_arm}");
    }

    #[test]
    fn test_dendrimer_critical_generation_positive() {
        let d = DendrimerModel::pamam(5, 1.0, 1.0);
        assert!(d.critical_generation() >= 1);
    }

    #[test]
    fn test_dendrimer_terminal_groups_double_with_generation() {
        // PAMAM: f=3, b=2: 3*2^g → doubling per generation
        let d1 = DendrimerModel::pamam(3, 1.0, 1.0);
        let d2 = DendrimerModel::pamam(4, 1.0, 1.0);
        let n1 = d1.n_terminal_groups();
        let n2 = d2.n_terminal_groups();
        assert_eq!(
            n2,
            2 * n1,
            "Each generation doubles terminal groups for b=2"
        );
    }
}
