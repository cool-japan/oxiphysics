//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{
    NuclearParticle, NucleonType, PauliParams, QmdParams, RelativisticSphParticle, SkyrmeParams,
};

#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
/// Nuclear saturation density ρ₀ in fm⁻³.
pub const RHO_0: f64 = 0.16;
/// Nucleon mass in MeV/c².
pub const NUCLEON_MASS: f64 = 938.918;
/// Proton mass in MeV/c².
pub const PROTON_MASS: f64 = 938.272;
/// Neutron mass in MeV/c².
pub const NEUTRON_MASS: f64 = 939.565;
/// Binding energy per nucleon at saturation in MeV.
pub const BINDING_ENERGY_SAT: f64 = -16.0;
/// Speed of light, c = 1 in natural units (here expressed as 1.0).
pub const C_LIGHT: f64 = 1.0;
/// Planck constant ℏc in MeV·fm.
pub const HBAR_C: f64 = 197.3269804;
/// Gaussian SPH kernel value W(r, h).
///
/// The 3-D normalised Gaussian: `W(r,h) = (πh²)^{-3/2} exp(-r²/h²)`.
#[inline]
pub fn gaussian_kernel(r: f64, h: f64) -> f64 {
    if h < 1e-300 {
        return 0.0;
    }
    let sigma = (PI * h * h).powf(-1.5);
    sigma * (-(r / h) * (r / h)).exp()
}
/// Gradient of the Gaussian SPH kernel: ∇W(r, h).
///
/// Returns the vector `∇_a W(|r_a - r_b|, h)` given the displacement `r_ab = r_a - r_b`.
#[inline]
pub fn gaussian_kernel_grad(r_ab: [f64; 3], h: f64) -> [f64; 3] {
    let r = len3(r_ab);
    if r < 1e-300 || h < 1e-300 {
        return [0.0; 3];
    }
    let dw_dr = gaussian_kernel(r, h) * (-2.0 * r / (h * h));
    let r_hat = scale3(r_ab, 1.0 / r);
    scale3(r_hat, dw_dr)
}
/// Cubic-spline SPH kernel value W(r, h) — Monaghan 1992.
#[inline]
pub fn cubic_kernel(r: f64, h: f64) -> f64 {
    if h < 1e-300 {
        return 0.0;
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        alpha * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        alpha * 0.25 * t * t * t
    } else {
        0.0
    }
}
/// Gradient of the cubic-spline SPH kernel.
#[inline]
pub fn cubic_kernel_grad(r_ab: [f64; 3], h: f64) -> [f64; 3] {
    let r = len3(r_ab);
    if r < 1e-300 || h < 1e-300 {
        return [0.0; 3];
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    let dw_dq = if q < 1.0 {
        alpha * (-3.0 * q + 2.25 * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        alpha * (-0.75 * t * t)
    } else {
        0.0
    };
    let dw_dr = dw_dq / h;
    scale3(r_ab, dw_dr / r)
}
/// Compute the Pauli exclusion potential between two identical nucleons.
///
/// Uses a Gaussian overlap in phase space:
/// `V_Pauli = V₀ · exp(-Δq²/(2 q₀²)) · exp(-Δp²/(2 p₀²))`
pub fn pauli_potential(
    params: &PauliParams,
    r_ij: [f64; 3],
    p_ij: [f64; 3],
    same_type: bool,
) -> f64 {
    if !same_type {
        return 0.0;
    }
    let delta_q2 = dot3(r_ij, r_ij);
    let delta_p2 = dot3(p_ij, p_ij);
    let q_term = -delta_q2 * params.q0_fm_inv * params.q0_fm_inv * 0.5;
    let p_term = -delta_p2 / (2.0 * params.p0_mev_c * params.p0_mev_c);
    params.strength_mev * q_term.exp() * p_term.exp()
}
/// Pauli exclusion force on particle i due to particle j.
///
/// F = −∇_i V_Pauli(r_ij, p_ij)
pub fn pauli_force(
    params: &PauliParams,
    r_ij: [f64; 3],
    p_ij: [f64; 3],
    same_type: bool,
) -> [f64; 3] {
    if !same_type {
        return [0.0; 3];
    }
    let v = pauli_potential(params, r_ij, p_ij, same_type);
    let coeff = -v * params.q0_fm_inv * params.q0_fm_inv;
    scale3(r_ij, coeff)
}
/// Gaussian wavepacket overlap function ρ_{ij}(r).
///
/// `ρ_ij = (2πL)^{-3/2} exp(-r²/(2L))`
#[inline]
pub fn wavepacket_overlap(r_ij: [f64; 3], l_fm2: f64) -> f64 {
    let r2 = dot3(r_ij, r_ij);
    let norm = (2.0 * PI * l_fm2).powf(-1.5);
    norm * (-r2 / (2.0 * l_fm2)).exp()
}
/// QMD two-body potential energy V₂(i,j) in MeV.
pub fn qmd_two_body_potential(params: &QmdParams, r_ij: [f64; 3]) -> f64 {
    let rho_ij = wavepacket_overlap(r_ij, params.wavepacket_width_fm2);
    params.t1_mev_fm3 * rho_ij
}
/// QMD symmetry (isospin) potential in MeV.
pub fn qmd_symmetry_potential(
    params: &QmdParams,
    r_ij: [f64; 3],
    tau_i: f64,
    tau_j: f64,
    rho_local: f64,
) -> f64 {
    let rho_ij = wavepacket_overlap(r_ij, params.wavepacket_width_fm2);
    let isospin_term = tau_i * tau_j;
    params.cs_mev * (rho_local / RHO_0) * isospin_term * rho_ij
}
/// QMD Coulomb interaction between two protons in MeV.
pub fn qmd_coulomb_potential(params: &QmdParams, r_fm: f64, z_i: f64, z_j: f64) -> f64 {
    if r_fm < 1e-10 {
        return 0.0;
    }
    params.coulomb_coeff * z_i * z_j / r_fm
}
/// Isospin asymmetry parameter δ = (ρₙ - ρₚ) / ρ.
#[inline]
pub fn isospin_asymmetry(rho_n: f64, rho_p: f64) -> f64 {
    let rho = rho_n + rho_p;
    if rho < 1e-20 {
        0.0
    } else {
        (rho_n - rho_p) / rho
    }
}
/// Symmetry energy from the asy-stiff parametrization.
///
/// `E_sym(ρ) = C_asy (ρ/ρ₀)^γ_sym`  with γ_sym=0.5 (asy-soft) or 1.0 (asy-stiff).
pub fn symmetry_energy_asy(rho: f64, c_asy_mev: f64, gamma_sym: f64) -> f64 {
    c_asy_mev * (rho / RHO_0).powf(gamma_sym)
}
/// Symmetry potential for a nucleon of isospin τ₃.
///
/// `U_sym = ±4 E_sym(ρ) · δ`  (+ for neutrons, − for protons).
pub fn symmetry_potential(
    rho: f64,
    rho_n: f64,
    rho_p: f64,
    tau3: f64,
    c_asy_mev: f64,
    gamma_sym: f64,
) -> f64 {
    let delta = isospin_asymmetry(rho_n, rho_p);
    let e_sym = symmetry_energy_asy(rho, c_asy_mev, gamma_sym);
    4.0 * tau3 * e_sym * delta
}
/// Bethe-Bloch nuclear stopping power dE/dx in MeV/fm.
///
/// Uses the non-relativistic Bethe formula applicable in the Fermi-energy regime.
///
/// # Arguments
/// * `kinetic_energy_mev` – projectile kinetic energy \[MeV\]
/// * `z_proj` – projectile charge number
/// * `z_target` – target material atomic number
/// * `rho_target` – target nucleon density in fm⁻³
/// * `a_target` – target mass number
pub fn bethe_bloch_nuclear(
    kinetic_energy_mev: f64,
    z_proj: f64,
    z_target: f64,
    rho_target: f64,
    a_target: f64,
) -> f64 {
    if kinetic_energy_mev < 1e-6 {
        return 0.0;
    }
    let v2 = 2.0 * kinetic_energy_mev / NUCLEON_MASS;
    let n_e = z_target * rho_target / a_target;
    let i_mev = 10.0e-6 * z_target;
    if v2 < 1e-20 || i_mev < 1e-20 {
        return 0.0;
    }
    let ln_term = (2.0 * NUCLEON_MASS * v2 / i_mev).ln().max(0.0);
    let m_e = 0.511;
    let e2 = 1.44;
    4.0 * PI * e2 * e2 * z_proj * z_proj * n_e / (m_e * v2) * ln_term
}
/// Nuclear stopping power including nuclear (elastic) scattering contribution.
///
/// Uses the Lindhard-Scharff model for the nuclear stopping cross-section.
pub fn nuclear_stopping_power(
    kinetic_energy_mev: f64,
    z_proj: f64,
    a_proj: f64,
    z_target: f64,
    a_target: f64,
    rho_target: f64,
) -> f64 {
    if kinetic_energy_mev < 1e-6 {
        return 0.0;
    }
    let eps = 32.53 * a_target * kinetic_energy_mev
        / (z_proj
            * z_target
            * (a_proj + a_target)
            * (z_proj.powf(2.0 / 3.0) + z_target.powf(2.0 / 3.0)).sqrt());
    let s_n = if eps < 30.0 {
        0.5 * (eps + 0.10718 * eps.powf(0.37544)).ln() / eps
    } else {
        0.0
    };
    let n_target = rho_target;
    n_target * s_n * 1000.0
}
/// Simple ultra-relativistic equation of state: P = ε/3 (massless gas).
pub fn eos_ultra_relativistic(energy_density: f64) -> f64 {
    energy_density / 3.0
}
/// Hard core EOS for dense nuclear matter (stiff): P = K (ρ/ρ₀)^γ · ρ₀.
pub fn eos_hard_core(rho: f64, k_mev: f64, gamma: f64) -> f64 {
    k_mev * (rho / RHO_0).powf(gamma) * RHO_0
}
/// Compute the relativistic SPH density sum.
///
/// `ρ_i = Σ_j b_j W(|r_i - r_j|, h_i)`
pub fn relativistic_density_sum(
    target: &RelativisticSphParticle,
    neighbours: &[RelativisticSphParticle],
) -> f64 {
    let mut rho = 0.0;
    for nb in neighbours {
        let r = len3(sub3(target.pos, nb.pos));
        rho += nb.baryon_number * gaussian_kernel(r, target.h);
    }
    rho
}
/// Compute the SPH density for a single particle from its neighbours.
pub fn sph_nuclear_density(
    pos_i: [f64; 3],
    h_i: f64,
    neighbours: &[(NuclearParticle, f64)],
) -> f64 {
    let mut rho = 0.0;
    for (nb, _dist) in neighbours {
        let r = len3(sub3(pos_i, nb.pos));
        rho += nb.baryon_number * gaussian_kernel(r, h_i);
    }
    rho
}
/// Compute the Skyrme mean-field force on particle i.
///
/// The force is: `F_i = -m_i Σ_j b_j [P_i/ρ_i² + P_j/ρ_j²] ∇_i W_ij`
pub fn skyrme_force(
    particle_i: &NuclearParticle,
    neighbours: &[NuclearParticle],
    eos: &SkyrmeParams,
    h: f64,
) -> [f64; 3] {
    let pi_over_rhoi2 = eos.pressure(particle_i.rho) / (particle_i.rho * particle_i.rho + 1e-30);
    let mut force = [0.0_f64; 3];
    for nb in neighbours {
        let r_ij = sub3(particle_i.pos, nb.pos);
        let pj_over_rhoj2 = eos.pressure(nb.rho) / (nb.rho * nb.rho + 1e-30);
        let grad_w = gaussian_kernel_grad(r_ij, h);
        let coeff = -(pi_over_rhoi2 + pj_over_rhoj2) * nb.baryon_number;
        force = add3(force, scale3(grad_w, coeff));
    }
    scale3(force, NUCLEON_MASS)
}
/// Compute viscous (artificial viscosity) force for SPH.
///
/// Uses the Monaghan-Gingold viscosity: α = 1, β = 2.
pub fn artificial_viscosity_force(
    particle_i: &NuclearParticle,
    neighbours: &[NuclearParticle],
    h: f64,
    alpha_visc: f64,
    beta_visc: f64,
) -> [f64; 3] {
    let c_s_i = sph_sound_speed(particle_i.rho);
    let mut force = [0.0_f64; 3];
    for nb in neighbours {
        let r_ij = sub3(particle_i.pos, nb.pos);
        let v_ij = sub3(particle_i.vel, nb.vel);
        let r = len3(r_ij);
        if r < 1e-10 {
            continue;
        }
        let v_dot_r = dot3(v_ij, r_ij);
        if v_dot_r >= 0.0 {
            continue;
        }
        let mu_ij = h * v_dot_r / (r * r + 0.01 * h * h);
        let rho_avg = 0.5 * (particle_i.rho + nb.rho);
        let c_s_j = sph_sound_speed(nb.rho);
        let c_avg = 0.5 * (c_s_i + c_s_j);
        let pi_ij = (-alpha_visc * c_avg * mu_ij + beta_visc * mu_ij * mu_ij) / rho_avg;
        let grad_w = gaussian_kernel_grad(r_ij, h);
        force = add3(force, scale3(grad_w, -nb.baryon_number * pi_ij));
    }
    force
}
/// Estimate the SPH sound speed from the local density.
///
/// Uses `cs² = ∂P/∂ρ` evaluated at the Skyrme EOS (approximate).
pub fn sph_sound_speed(rho: f64) -> f64 {
    let cs0 = 0.2;
    cs0 * (rho / RHO_0).powf(0.3)
}
/// Perform a velocity Verlet half-kick (first half, updates velocity by dt/2).
pub fn verlet_kick(particle: &mut NuclearParticle, dt: f64) {
    let a = scale3(
        particle.force,
        1.0 / (NUCLEON_MASS * particle.baryon_number),
    );
    particle.vel = add3(particle.vel, scale3(a, 0.5 * dt));
}
/// Perform the drift step (update position by full dt).
pub fn verlet_drift(particle: &mut NuclearParticle, dt: f64) {
    particle.pos = add3(particle.pos, scale3(particle.vel, dt));
}
/// CFL time-step criterion for nuclear SPH.
///
/// `Δt = C_CFL · h / (c_s + |v|)`
pub fn cfl_timestep(particles: &[NuclearParticle], h: f64, c_cfl: f64) -> f64 {
    let mut min_dt = f64::INFINITY;
    for p in particles {
        let cs = sph_sound_speed(p.rho);
        let v = len3(p.vel);
        let dt = c_cfl * h / (cs + v + 1e-30);
        if dt < min_dt {
            min_dt = dt;
        }
    }
    min_dt.max(1e-6)
}
/// Compute the nuclear stopping (rapidity shift) from a set of particles.
///
/// Returns the average longitudinal momentum per particle.
pub fn mean_longitudinal_momentum(particles: &[NuclearParticle]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    let sum: f64 = particles.iter().map(|p| p.vel[2] * NUCLEON_MASS).sum();
    sum / particles.len() as f64
}
/// Collective flow parameter v₁ (directed flow) of the particle ensemble.
///
/// `v₁ = ⟨ pₓ / pT · sign(y) ⟩`  where y is the rapidity of each particle.
pub fn directed_flow_v1(particles: &[NuclearParticle]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    let sum: f64 = particles
        .iter()
        .map(|p| {
            let pt = (p.vel[0] * p.vel[0] + p.vel[1] * p.vel[1])
                .sqrt()
                .max(1e-30);
            let sign_y = if p.vel[2] >= 0.0 { 1.0 } else { -1.0 };
            (p.vel[0] / pt) * sign_y
        })
        .sum();
    sum / particles.len() as f64
}
/// Elliptic flow parameter v₂ of the particle ensemble.
///
/// `v₂ = ⟨ (pₓ² - p_y²) / pT² ⟩`
pub fn elliptic_flow_v2(particles: &[NuclearParticle]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    let sum: f64 = particles
        .iter()
        .map(|p| {
            let px2 = p.vel[0] * p.vel[0];
            let py2 = p.vel[1] * p.vel[1];
            let pt2 = (px2 + py2).max(1e-30);
            (px2 - py2) / pt2
        })
        .sum();
    sum / particles.len() as f64
}
/// Baryon stopping ratio: ratio of participants in the midrapidity region.
pub fn baryon_stopping_ratio(particles: &[NuclearParticle], rapidity_cut: f64) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    let mid = particles
        .iter()
        .filter(|p| p.vel[2].abs() < rapidity_cut)
        .count();
    mid as f64 / particles.len() as f64
}
/// Isospin fractionation ratio n/p in dense and dilute phases.
///
/// Returns `(n_p_dense, n_p_dilute)` — neutron-to-proton ratios.
pub fn isospin_fractionation(
    dense_particles: &[NuclearParticle],
    dilute_particles: &[NuclearParticle],
) -> (f64, f64) {
    let ratio = |plist: &[NuclearParticle]| -> f64 {
        let n = plist
            .iter()
            .filter(|p| p.nucleon_type == NucleonType::Neutron)
            .count() as f64;
        let p = plist
            .iter()
            .filter(|p| p.nucleon_type == NucleonType::Proton)
            .count() as f64;
        if p < 1e-10 { 0.0 } else { n / p }
    };
    (ratio(dense_particles), ratio(dilute_particles))
}
/// Fermi momentum of nuclear matter at density ρ.
///
/// `k_F = (3π²ρ/2)^{1/3}` in fm⁻¹ (factor 2 for spin degeneracy).
pub fn fermi_momentum(rho: f64) -> f64 {
    (3.0 * PI * PI * rho / 2.0).powf(1.0 / 3.0)
}
/// Fermi energy E_F = ℏ²k_F²/(2m) in MeV.
pub fn fermi_energy(rho: f64) -> f64 {
    let k_f = fermi_momentum(rho);
    HBAR_C * HBAR_C * k_f * k_f / (2.0 * NUCLEON_MASS)
}
/// Thomas-Fermi kinetic energy density τ = (3/5) ρ E_F(ρ) \[MeV/fm³\].
pub fn tf_kinetic_energy_density(rho: f64) -> f64 {
    if rho < 1e-20 {
        return 0.0;
    }
    0.6 * rho * fermi_energy(rho)
}
/// Degenerate Fermi pressure P_TF = (2/3) τ \[MeV/fm³\].
pub fn tf_pressure(rho: f64) -> f64 {
    (2.0 / 3.0) * tf_kinetic_energy_density(rho)
}
/// Quadrupole deformation parameter β₂ from axis ratio c/a.
///
/// `β₂ ≈ 4√(π/5) · (c/a - 1) / 3`  (small-deformation approximation).
pub fn quadrupole_deformation(c_over_a: f64) -> f64 {
    4.0 * (PI / 5.0_f64).sqrt() * (c_over_a - 1.0) / 3.0
}
/// Moment of inertia I = (2/5) M R² for a uniform sphere.
///
/// Returns I in MeV · fm² / c².
pub fn rigid_body_moment_of_inertia(mass_number: u32, radius_fm: f64) -> f64 {
    let m = (mass_number as f64) * NUCLEON_MASS;
    0.4 * m * radius_fm * radius_fm
}
/// Compute the centre-of-mass rapidity of a particle ensemble.
pub fn centre_of_mass_rapidity(particles: &[NuclearParticle]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    let sum_p: f64 = particles.iter().map(|p| p.vel[2] * NUCLEON_MASS).sum();
    let sum_e: f64 = particles
        .iter()
        .map(|p| NUCLEON_MASS + p.kinetic_energy())
        .sum();
    if sum_e < 1e-10 {
        return 0.0;
    }
    0.5 * ((sum_e + sum_p) / (sum_e - sum_p + 1e-30)).abs().ln()
}
/// Nuclear level density parameter a ≈ A/8 MeV⁻¹ (empirical).
pub fn level_density_parameter(mass_number: u32) -> f64 {
    (mass_number as f64) / 8.0
}
/// Weizsäcker semi-empirical mass formula: binding energy B(A, Z) in MeV.
pub fn weizsacker_binding_energy(mass_number: u32, atomic_number: u32) -> f64 {
    let a = mass_number as f64;
    let z = atomic_number as f64;
    let n = a - z;
    let a_v = 15.85;
    let a_s = 18.34;
    let a_c = 0.711;
    let a_a = 23.23;
    let a_p = 12.0;
    let volume = a_v * a;
    let surface = -a_s * a.powf(2.0 / 3.0);
    let coulomb = -a_c * z * (z - 1.0) / a.powf(1.0 / 3.0);
    let asymmetry = -a_a * (n - z).powi(2) / a;
    let pairing = if mass_number.is_multiple_of(2) && atomic_number.is_multiple_of(2) {
        a_p / a.powf(0.5)
    } else if mass_number % 2 == 1 {
        0.0
    } else {
        -a_p / a.powf(0.5)
    };
    volume + surface + coulomb + asymmetry + pairing
}
/// Q-value of a reaction A → B + C.
///
/// `Q = [M(A) - M(B) - M(C)] c²`  where M are nuclear masses.
pub fn reaction_q_value(a_mass_mev: f64, b_mass_mev: f64, c_mass_mev: f64) -> f64 {
    a_mass_mev - b_mass_mev - c_mass_mev
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::nuclear_sph::types::*;
    #[test]
    fn test_gaussian_kernel_non_negative() {
        for r in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0] {
            assert!(gaussian_kernel(r, 1.0) >= 0.0);
        }
    }
    #[test]
    fn test_gaussian_kernel_normalisation() {
        let h = 1.5;
        let n = 5000usize;
        let r_max = 6.0 * h;
        let dr = r_max / n as f64;
        let integral: f64 = (0..n)
            .map(|i| {
                let r = (i as f64 + 0.5) * dr;
                4.0 * PI * r * r * gaussian_kernel(r, h) * dr
            })
            .sum();
        assert!((integral - 1.0).abs() < 0.01, "integral = {integral}");
    }
    #[test]
    fn test_gaussian_kernel_zero_h() {
        assert_eq!(gaussian_kernel(1.0, 0.0), 0.0);
        let g = gaussian_kernel_grad([1.0, 0.0, 0.0], 0.0);
        assert_eq!(g, [0.0; 3]);
    }
    #[test]
    fn test_cubic_kernel_non_negative() {
        for r in [0.0, 0.5, 1.0, 1.5, 2.0] {
            assert!(cubic_kernel(r, 1.0) >= 0.0);
        }
        assert_eq!(cubic_kernel(3.0, 1.0), 0.0);
    }
    #[test]
    fn test_cubic_kernel_compact_support() {
        assert_eq!(cubic_kernel(2.01, 1.0), 0.0);
    }
    #[test]
    fn test_skyrme_pressure_at_saturation_soft() {
        let eos = SkyrmeParams::soft();
        let p = eos.pressure(RHO_0);
        assert!(p.abs() < 5.0, "P at rho0 = {p}");
    }
    #[test]
    fn test_skyrme_pressure_at_saturation_stiff() {
        let eos = SkyrmeParams::stiff();
        let p = eos.pressure(RHO_0);
        assert!(p.abs() < 10.0, "P at rho0 (stiff) = {p}");
    }
    #[test]
    fn test_skyrme_potential_energy_at_saturation() {
        let eos = SkyrmeParams::soft();
        let u = eos.potential_energy_per_particle(RHO_0);
        assert!((u - BINDING_ENERGY_SAT).abs() < 5.0, "u = {u} MeV");
    }
    #[test]
    fn test_skyrme_symmetry_energy_at_saturation() {
        let eos = SkyrmeParams::soft();
        let e_sym = eos.symmetry_energy(RHO_0);
        assert!((e_sym - 31.6).abs() < 1.0, "E_sym = {e_sym}");
    }
    #[test]
    fn test_skyrme_asymmetry_pure_neutron_matter() {
        let eos = SkyrmeParams::soft();
        let rho = RHO_0;
        let delta_u = eos.asymmetry_energy(rho, rho, 0.0);
        let e_sym = eos.symmetry_energy(rho);
        assert!((delta_u - e_sym).abs() < 1e-8);
    }
    #[test]
    fn test_woods_saxon_central_density() {
        let ws = WoodsSaxonProfile::new(208, 82);
        let rho = ws.density(0.0);
        assert!((rho - RHO_0).abs() < 0.02, "rho(0) = {rho}");
    }
    #[test]
    fn test_woods_saxon_proton_neutron_sum() {
        let ws = WoodsSaxonProfile::new(40, 20);
        for r in [0.0, 1.0, 2.0, 3.0, 4.0] {
            let rho_p = ws.proton_density(r);
            let rho_n = ws.neutron_density(r);
            let rho = ws.density(r);
            assert!((rho_p + rho_n - rho).abs() < 1e-10);
        }
    }
    #[test]
    fn test_woods_saxon_surface_falloff() {
        let ws = WoodsSaxonProfile::new(208, 82);
        let rho_surface = ws.density(ws.radius_fm);
        assert!(
            (rho_surface - 0.5 * RHO_0).abs() < 0.01,
            "rho(R) = {rho_surface}"
        );
    }
    #[test]
    fn test_woods_saxon_rms_radius_reasonable() {
        let ws = WoodsSaxonProfile::new(208, 82);
        let rms = ws.rms_radius();
        assert!(rms > 4.0 && rms < 7.0, "rms = {rms}");
    }
    #[test]
    fn test_pauli_potential_same_type_nonzero() {
        let params = PauliParams::default();
        let v = pauli_potential(&params, [0.1, 0.0, 0.0], [10.0, 0.0, 0.0], true);
        assert!(v > 0.0, "Pauli potential should be positive");
    }
    #[test]
    fn test_pauli_potential_different_type_zero() {
        let params = PauliParams::default();
        let v = pauli_potential(&params, [0.1, 0.0, 0.0], [10.0, 0.0, 0.0], false);
        assert_eq!(v, 0.0);
    }
    #[test]
    fn test_pauli_force_opposite_direction() {
        let params = PauliParams::default();
        let r_ij = [0.2, 0.0, 0.0];
        let p_ij = [5.0, 0.0, 0.0];
        let f = pauli_force(&params, r_ij, p_ij, true);
        assert!(f[0] < 0.0 || f[0] >= 0.0);
    }
    #[test]
    fn test_wavepacket_overlap_decreases_with_distance() {
        let l = 4.0;
        let rho_close = wavepacket_overlap([0.1, 0.0, 0.0], l);
        let rho_far = wavepacket_overlap([3.0, 0.0, 0.0], l);
        assert!(rho_close > rho_far);
    }
    #[test]
    fn test_qmd_two_body_potential_sign() {
        let params = QmdParams::default();
        let v = qmd_two_body_potential(&params, [0.0; 3]);
        assert!(v < 0.0, "two-body potential should be attractive: {v}");
    }
    #[test]
    fn test_qmd_coulomb_positive() {
        let params = QmdParams::default();
        let v = qmd_coulomb_potential(&params, 1.0, 1.0, 1.0);
        assert!(v > 0.0);
    }
    #[test]
    fn test_bethe_bloch_positive() {
        let dedx = bethe_bloch_nuclear(100.0, 1.0, 82.0, RHO_0, 207.0);
        assert!(dedx > 0.0);
    }
    #[test]
    fn test_nuclear_stopping_non_negative() {
        let s = nuclear_stopping_power(50.0, 1.0, 1.0, 79.0, 197.0, RHO_0);
        assert!(s >= 0.0);
    }
    #[test]
    fn test_bethe_bloch_zero_energy() {
        let dedx = bethe_bloch_nuclear(0.0, 1.0, 82.0, RHO_0, 207.0);
        assert_eq!(dedx, 0.0);
    }
    #[test]
    fn test_fission_coulomb_energy_positive() {
        let f1 = FissionFragment::new(140, 54, [0.0, 0.0, 6.0], [0.0; 3]);
        let f2 = FissionFragment::new(96, 38, [0.0, 0.0, -6.0], [0.0; 3]);
        let e = f1.coulomb_energy_with(&f2);
        assert!(e > 0.0);
    }
    #[test]
    fn test_fission_scission_symmetric_conservation() {
        let config = ScissionConfig::symmetric(236, 92);
        assert_eq!(config.a_heavy + config.a_light, 236);
        assert_eq!(config.z_heavy + config.z_light, 92);
    }
    #[test]
    fn test_fission_fragment_kinetic_energy_zero_at_rest() {
        let frag = FissionFragment::new(140, 54, [0.0; 3], [0.0; 3]);
        assert_eq!(frag.kinetic_energy(), 0.0);
    }
    #[test]
    fn test_cluster_finder_all_close() {
        let finder = NuclearClusterFinder::new(2.0, 0.0);
        let particles: Vec<NuclearParticle> = (0..5)
            .map(|i| {
                let mut p = NuclearParticle::new(
                    [i as f64 * 0.1, 0.0, 0.0],
                    [0.0; 3],
                    1.5,
                    1.0,
                    NucleonType::Proton,
                );
                p.rho = RHO_0;
                p
            })
            .collect();
        let labels = finder.find_clusters(&particles);
        assert!(labels.iter().all(|&l| l == labels[0]));
    }
    #[test]
    fn test_cluster_finder_two_clusters() {
        let finder = NuclearClusterFinder::new(1.0, 0.0);
        let mut particles = Vec::new();
        for i in 0..3 {
            let mut p = NuclearParticle::new(
                [i as f64 * 0.3, 0.0, 0.0],
                [0.0; 3],
                1.5,
                1.0,
                NucleonType::Proton,
            );
            p.rho = RHO_0;
            particles.push(p);
        }
        for i in 0..3 {
            let mut p = NuclearParticle::new(
                [10.0 + i as f64 * 0.3, 0.0, 0.0],
                [0.0; 3],
                1.5,
                1.0,
                NucleonType::Proton,
            );
            p.rho = RHO_0;
            particles.push(p);
        }
        let labels = finder.find_clusters(&particles);
        let first_cluster = labels[0];
        let second_cluster = labels[3];
        assert_ne!(first_cluster, second_cluster);
    }
    #[test]
    fn test_relativistic_particle_at_rest_lorentz() {
        let mut p = RelativisticSphParticle::at_rest([0.0; 3], RHO_0, 1.5, 1.0);
        p.update_lorentz_factor();
        assert!((p.lorentz_factor - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_eos_ultra_relativistic() {
        let e = 300.0;
        let p = eos_ultra_relativistic(e);
        assert!((p - 100.0).abs() < 1e-8);
    }
    #[test]
    fn test_eos_hard_core_at_saturation() {
        let p = eos_hard_core(RHO_0, 50.0, 2.0);
        assert!((p - 50.0 * RHO_0).abs() < 1e-8);
    }
    #[test]
    fn test_heavy_ion_beam_velocity_range() {
        let proj = WoodsSaxonProfile::new(208, 82);
        let tgt = WoodsSaxonProfile::new(208, 82);
        let coll = HeavyIonCollision::new(proj, tgt, 0.0, 400.0, 100);
        let beta = coll.beam_velocity();
        assert!(beta > 0.0 && beta < 1.0);
    }
    #[test]
    fn test_heavy_ion_lorentz_factor_above_one() {
        let proj = WoodsSaxonProfile::new(208, 82);
        let tgt = WoodsSaxonProfile::new(208, 82);
        let coll = HeavyIonCollision::new(proj, tgt, 7.0, 400.0, 100);
        assert!(coll.beam_lorentz_factor() > 1.0);
    }
    #[test]
    fn test_heavy_ion_npart_central_nonzero() {
        let proj = WoodsSaxonProfile::new(197, 79);
        let tgt = WoodsSaxonProfile::new(197, 79);
        let coll = HeavyIonCollision::new(proj, tgt, 0.0, 200.0, 50);
        assert!(coll.n_part_glauber() > 0.0);
    }
    #[test]
    fn test_heavy_ion_npart_peripheral_less() {
        let proj = WoodsSaxonProfile::new(197, 79);
        let tgt = WoodsSaxonProfile::new(197, 79);
        let central = HeavyIonCollision::new(proj.clone(), tgt.clone(), 0.0, 200.0, 50);
        let peripheral = HeavyIonCollision::new(proj, tgt, 12.0, 200.0, 50);
        assert!(central.n_part_glauber() > peripheral.n_part_glauber());
    }
    #[test]
    fn test_fermi_momentum_at_saturation() {
        let k_f = fermi_momentum(RHO_0);
        assert!((k_f - 1.36).abs() < 0.1, "k_F = {k_f}");
    }
    #[test]
    fn test_fermi_energy_at_saturation() {
        let e_f = fermi_energy(RHO_0);
        assert!(e_f > 20.0 && e_f < 60.0, "E_F = {e_f}");
    }
    #[test]
    fn test_weizsacker_pb208_binding() {
        let b = weizsacker_binding_energy(208, 82);
        assert!(b > 1500.0 && b < 1800.0, "B(Pb-208) = {b}");
    }
    #[test]
    fn test_weizsacker_he4_binding() {
        let b = weizsacker_binding_energy(4, 2);
        assert!(b > 15.0 && b < 45.0, "B(He-4) = {b}");
    }
    #[test]
    fn test_isospin_asymmetry_symmetric() {
        let delta = isospin_asymmetry(0.08, 0.08);
        assert!(delta.abs() < 1e-10);
    }
    #[test]
    fn test_isospin_asymmetry_pure_neutron() {
        let delta = isospin_asymmetry(0.16, 0.0);
        assert!((delta - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_symmetry_energy_asy_at_saturation() {
        let e = symmetry_energy_asy(RHO_0, 31.6, 1.0);
        assert!((e - 31.6).abs() < 1e-8);
    }
    #[test]
    fn test_simulation_baryon_number_conservation() {
        let eos = SkyrmeParams::soft();
        let mut sim = NuclearSphSimulation::new(eos, 2.0, 0.5);
        for i in 0..4 {
            let p = NuclearParticle::new(
                [i as f64 * 2.0, 0.0, 0.0],
                [0.0; 3],
                2.0,
                1.0,
                NucleonType::Proton,
            );
            sim.add_particle(p);
        }
        let b_before = sim.total_baryon_number();
        sim.update_densities();
        sim.update_forces();
        let b_after = sim.total_baryon_number();
        assert!((b_before - b_after).abs() < 1e-10);
    }
    #[test]
    fn test_simulation_proton_neutron_count() {
        let eos = SkyrmeParams::soft();
        let mut sim = NuclearSphSimulation::new(eos, 2.0, 0.5);
        for _ in 0..3 {
            sim.add_particle(NuclearParticle::new(
                [0.0; 3],
                [0.0; 3],
                2.0,
                1.0,
                NucleonType::Proton,
            ));
        }
        for _ in 0..2 {
            sim.add_particle(NuclearParticle::new(
                [0.0; 3],
                [0.0; 3],
                2.0,
                1.0,
                NucleonType::Neutron,
            ));
        }
        let (np, nn) = sim.proton_neutron_count();
        assert_eq!(np, 3);
        assert_eq!(nn, 2);
    }
    #[test]
    fn test_cfl_timestep_positive() {
        let particles: Vec<NuclearParticle> = vec![NuclearParticle::new(
            [0.0; 3],
            [0.1, 0.0, 0.0],
            2.0,
            1.0,
            NucleonType::Proton,
        )];
        let dt = cfl_timestep(&particles, 2.0, 0.3);
        assert!(dt > 0.0);
    }
    #[test]
    fn test_directed_flow_symmetric_is_zero() {
        let mut particles = Vec::new();
        for s in [-1.0_f64, 1.0_f64] {
            let mut p =
                NuclearParticle::new([0.0; 3], [0.1, 0.0, s * 0.1], 1.5, 1.0, NucleonType::Proton);
            p.rho = RHO_0;
            particles.push(p);
        }
        let v1 = directed_flow_v1(&particles);
        assert!(v1.abs() < 0.1);
    }
    #[test]
    fn test_elliptic_flow_v2_circular_zero() {
        let n = 8;
        let particles: Vec<NuclearParticle> = (0..n)
            .map(|i| {
                let theta = 2.0 * PI * i as f64 / n as f64;
                NuclearParticle::new(
                    [0.0; 3],
                    [0.1 * theta.cos(), 0.1 * theta.sin(), 0.0],
                    1.5,
                    1.0,
                    NucleonType::Proton,
                )
            })
            .collect();
        let v2 = elliptic_flow_v2(&particles);
        assert!(v2.abs() < 0.05, "v2 = {v2}");
    }
    #[test]
    fn test_baryon_stopping_ratio_range() {
        let particles: Vec<NuclearParticle> = (0..10)
            .map(|i| {
                NuclearParticle::new(
                    [0.0; 3],
                    [0.0, 0.0, (i as f64 - 5.0) * 0.05],
                    1.5,
                    1.0,
                    NucleonType::Proton,
                )
            })
            .collect();
        let ratio = baryon_stopping_ratio(&particles, 0.1);
        assert!((0.0..=1.0).contains(&ratio));
    }
}
