//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{CrankNicolsonResult, ThermalMesh1D};

/// Solve a tridiagonal system using the Thomas algorithm.
pub fn thomas_algorithm(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = b.len();
    assert_eq!(a.len(), n);
    assert_eq!(c.len(), n);
    assert_eq!(d.len(), n);
    let mut b_mod = b.to_vec();
    let mut d_mod = d.to_vec();
    let mut x = vec![0.0_f64; n];
    for i in 1..n {
        let m = a[i] / b_mod[i - 1];
        b_mod[i] -= m * c[i - 1];
        d_mod[i] -= m * d_mod[i - 1];
    }
    x[n - 1] = d_mod[n - 1] / b_mod[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (d_mod[i] - c[i] * x[i + 1]) / b_mod[i];
    }
    x
}
/// Solve a transient heat conduction problem using the theta method.
///
/// (C/dt + theta*K) T_{n+1} = (C/dt - (1-theta)*K) T_n + theta*Q_{n+1} + (1-theta)*Q_n
///
/// theta = 0: forward Euler (explicit)
/// theta = 0.5: Crank-Nicolson
/// theta = 1: backward Euler (implicit)
///
/// Returns the temperature history at each time step.
pub fn transient_theta_method(
    mesh: &mut ThermalMesh1D,
    dt: f64,
    n_steps: usize,
    theta: f64,
    heat_sources: &[f64],
) -> Vec<Vec<f64>> {
    assert!((0.0..=1.0).contains(&theta));
    let n = mesh.n_nodes;
    let k = mesh.assemble_conductance_matrix();
    let c = mesh.assemble_capacitance_vector();
    let mut history = Vec::with_capacity(n_steps + 1);
    history.push(mesh.temperatures.clone());
    for _step in 0..n_steps {
        let mut a = vec![0.0_f64; n];
        let mut b = vec![0.0_f64; n];
        let mut cv = vec![0.0_f64; n];
        let mut rhs = vec![0.0_f64; n];
        for i in 0..n {
            b[i] = c[i] / dt + theta * k[i][i];
            if i > 0 {
                a[i] = theta * k[i][i - 1];
            }
            if i < n - 1 {
                cv[i] = theta * k[i][i + 1];
            }
            let kt_i: f64 = k[i]
                .iter()
                .zip(mesh.temperatures.iter())
                .map(|(&kij, &tj)| kij * tj)
                .sum();
            rhs[i] = c[i] / dt * mesh.temperatures[i] - (1.0 - theta) * kt_i + heat_sources[i];
        }
        for &(node, temp) in &mesh.bc {
            a[node] = 0.0;
            b[node] = 1.0;
            cv[node] = 0.0;
            rhs[node] = temp;
        }
        mesh.temperatures = thomas_algorithm(&a, &b, &cv, &rhs);
        history.push(mesh.temperatures.clone());
    }
    history
}
/// Advance the phase-change front position using Stefan's condition.
///
/// `ds/dt = q_flux / (rho * L)`  (heat balance at the front)
///
/// Returns the new front position after one explicit time step `dt`.
///
/// # Arguments
/// * `s`         - current front position \[m\]
/// * `q_flux`    - heat flux at the front \[W/m²\]
/// * `rho`       - density \[kg/m³\]
/// * `latent_heat` - latent heat of fusion L \[J/kg\]
/// * `dt`        - time step \[s\]
pub fn stefan_front_position(s: f64, q_flux: f64, rho: f64, latent_heat: f64, dt: f64) -> f64 {
    if latent_heat.abs() < 1e-30 || rho.abs() < 1e-30 {
        return s;
    }
    s + q_flux / (rho * latent_heat) * dt
}
/// Neumann exact solution for the one-phase Stefan problem.
///
/// A semi-infinite region, initially at the melt temperature, has its boundary
/// held at `T_melt + delta_t` (`delta_t > 0`); the melt front then advances as
/// `s(t) = 2 * lambda * sqrt(alpha * t)`, where `alpha` is the thermal
/// diffusivity of the melt.  The dimensionless front coefficient `lambda` is the
/// positive root of the transcendental equation
///
/// `lambda * exp(lambda²) * erf(lambda) = St / sqrt(pi)`,
///
/// with the Stefan number `St = cp * delta_t / latent_heat`.  `lambda` depends
/// on the material only through `St` (the diffusivity and the elapsed time enter
/// the front position, not `lambda`).
///
/// # Arguments
/// * `cp`          - specific heat of the melt \[J/(kg K)\]
/// * `delta_t`     - superheat `T_boundary - T_melt` \[K\]
/// * `latent_heat` - latent heat of fusion `L` \[J/kg\]
///
/// # Returns
/// `lambda >= 0`.  `0.0` when `St == 0` (no superheat, so no front motion).
/// `NaN` when an input is not finite, when `latent_heat <= 0`, or when `St < 0`
/// (a sub-melt boundary temperature has no one-phase melting solution).
pub fn stefan_neumann_lambda(cp: f64, delta_t: f64, latent_heat: f64) -> f64 {
    if !(cp.is_finite() && delta_t.is_finite() && latent_heat.is_finite()) || latent_heat <= 0.0 {
        return f64::NAN;
    }
    let stefan_number = cp * delta_t / latent_heat;
    if stefan_number.is_nan() || stefan_number < 0.0 {
        return f64::NAN;
    }
    if stefan_number == 0.0 {
        return 0.0;
    }
    let sqrt_pi = std::f64::consts::PI.sqrt();
    let target = stefan_number / sqrt_pi;
    // g(lambda) = lambda e^{lambda^2} erf(lambda) - target is increasing and
    // convex on lambda > 0 with g(0) < 0, so it has exactly one root.
    let residual = |l: f64| l * (l * l).exp() * erf_series(l) - target;
    let slope = |l: f64| {
        let e = (l * l).exp();
        e * erf_series(l) * (1.0 + 2.0 * l * l) + 2.0 * l / sqrt_pi
    };
    // Bracket the root.  For St <= 1/2, lambda = sqrt(2 St) <= 1 already gives
    // g >= 0 (g(l) >= 2/sqrt(pi) l^2 on [0, 1]) and lies within a factor of two
    // of the small-St root sqrt(St / 2); e^{lambda^2} overflows beyond lambda ~ 26.
    const LAMBDA_MAX: f64 = 26.0;
    let mut hi = (2.0 * stefan_number).sqrt().min(1.0);
    while residual(hi) < 0.0 {
        if hi >= LAMBDA_MAX {
            return f64::NAN;
        }
        hi = (hi * 2.0).min(LAMBDA_MAX);
    }
    let mut lo = 0.0_f64;
    // Newton started on the right of a convex increasing function converges
    // monotonically; bisection keeps every iterate inside the bracket anyway.
    let mut lambda = hi;
    for _ in 0..200 {
        let g = residual(lambda);
        if g == 0.0 {
            break;
        }
        if g < 0.0 {
            lo = lambda;
        } else {
            hi = lambda;
        }
        let mut next = lambda - g / slope(lambda);
        if !next.is_finite() || next <= lo || next >= hi {
            next = 0.5 * (lo + hi);
        }
        let converged = (next - lambda).abs() <= 1e-14 * next.abs();
        lambda = next;
        if converged {
            break;
        }
    }
    lambda
}
/// Melt-front position of the one-phase Stefan problem at time `t`,
/// `s(t) = 2 * lambda * sqrt(alpha * t)`, with the Neumann coefficient
/// `lambda = stefan_neumann_lambda(cp, delta_t, latent_heat)`.
///
/// # Arguments
/// * `cp`          - specific heat of the melt \[J/(kg K)\]
/// * `delta_t`     - superheat `T_boundary - T_melt` \[K\]
/// * `latent_heat` - latent heat of fusion `L` \[J/kg\]
/// * `alpha`       - thermal diffusivity of the melt \[m²/s\]
/// * `t`           - elapsed time \[s\]
///
/// Returns `NaN` when [`stefan_neumann_lambda`] does, or when `alpha` or `t` is
/// negative or not finite.
pub fn stefan_neumann_front_position(
    cp: f64,
    delta_t: f64,
    latent_heat: f64,
    alpha: f64,
    t: f64,
) -> f64 {
    if !(alpha.is_finite() && t.is_finite()) || alpha < 0.0 || t < 0.0 {
        return f64::NAN;
    }
    2.0 * stefan_neumann_lambda(cp, delta_t, latent_heat) * (alpha * t).sqrt()
}
/// Error function to double precision.
///
/// Uses the non-alternating Maclaurin series
/// `erf(x) = 2/sqrt(pi) * exp(-x²) * sum_n 2^n x^(2n+1) / (2n+1)!!`, whose terms
/// are all positive (no cancellation), and saturates to `±1` beyond `|x| = 6`
/// where `1 - |erf(x)| < 2e-17`.
pub(super) fn erf_series(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let ax = x.abs();
    if ax >= 6.0 {
        return x.signum();
    }
    let x2 = ax * ax;
    let mut term = ax;
    let mut sum = ax;
    for n in 1..400 {
        term *= 2.0 * x2 / (2.0 * n as f64 + 1.0);
        sum += term;
        if term <= sum * f64::EPSILON * 1e-2 {
            break;
        }
    }
    let value = 2.0 / std::f64::consts::PI.sqrt() * (-x2).exp() * sum;
    if x < 0.0 { -value } else { value }
}
/// Latent heat source term for the enthalpy method.
///
/// Returns `rho * L * df_l/dt` at the current temperature node.
/// Near the melting point, releases latent heat over a small temperature range.
///
/// `f_l` = liquid fraction = `max(0, min(1, (T - T_solidus) / width))`
pub fn latent_heat_source(temperature: f64, t_melt: f64, rho: f64, latent_heat: f64) -> f64 {
    let width = 0.5;
    let t_solidus = t_melt - width;
    let t_liquidus = t_melt + width;
    if temperature < t_solidus || temperature > t_liquidus {
        return 0.0;
    }
    rho * latent_heat / (2.0 * width)
}
/// Effective specific heat with latent heat contribution near melting.
///
/// `Cp_eff = Cp + rho * L * df_l/dT`
///
/// where `df_l/dT = 1 / (2 * width)` within the mushy zone.
pub fn effective_specific_heat_melting(
    temperature: f64,
    cp: f64,
    latent_heat: f64,
    width: f64,
    t_melt: f64,
) -> f64 {
    let t_solidus = t_melt - width;
    let t_liquidus = t_melt + width;
    if (t_solidus..=t_liquidus).contains(&temperature) {
        cp + latent_heat / (2.0 * width)
    } else {
        cp
    }
}
/// Advance the solidification front by one time step using Stefan's condition.
///
/// `ds/dt = k * dT/dx|_front / (rho * L)`
///
/// # Arguments
/// * `s`       - current front position \[m\]
/// * `q_grad`  - temperature gradient `dT/dx` at front \[K/m\]
/// * `k`       - conductivity \[W/(m K)\]
/// * `rho`     - density \[kg/m³\]
/// * `lat`     - latent heat \[J/kg\]
/// * `nu`      - kinematic viscosity (unused, for API symmetry)
/// * `dt`      - time step \[s\]
pub fn solidification_front_advance(
    s: f64,
    q_grad: f64,
    k: f64,
    rho: f64,
    lat: f64,
    _nu: f64,
    dt: f64,
) -> f64 {
    if lat.abs() < 1e-30 || rho.abs() < 1e-30 {
        return s;
    }
    s + k * q_grad / (rho * lat) * dt
}
/// Liquid fraction at a given temperature (for mushy-zone tracking).
///
/// `f_l = clamp((T - T_solidus) / (T_liquidus - T_solidus), 0, 1)`
///
/// `width` = half-width of mushy zone \[K\].
pub fn solidification_fraction(temperature: f64, t_solidus: f64, width: f64) -> f64 {
    let t_liquidus = t_solidus + 2.0 * width;
    ((temperature - t_solidus) / (t_liquidus - t_solidus)).clamp(0.0, 1.0)
}
/// Estimate front position from temperature profile using linear interpolation.
///
/// Finds the node where the temperature crosses `t_melt`.
pub fn find_solidification_front(temperatures: &[f64], positions: &[f64], t_melt: f64) -> f64 {
    let n = temperatures.len();
    if n < 2 {
        return positions[0];
    }
    for i in 0..n - 1 {
        let t0 = temperatures[i];
        let t1 = temperatures[i + 1];
        if (t0 - t_melt) * (t1 - t_melt) <= 0.0 {
            let alpha = (t_melt - t0) / (t1 - t0 + 1e-30);
            return positions[i] + alpha * (positions[i + 1] - positions[i]);
        }
    }
    positions[n - 1]
}
/// View factor between two equal, parallel, directly-opposed rectangles of side `a`
/// separated by distance `h` (simplified formula for infinite parallel plates):
///
/// `F12 = [sqrt(4/H^2 + 1) - 1] / H` where `H = h/a`.
///
/// For the limiting case H → 0 (large plates), `F12 → 1`.
pub fn view_factor_parallel_plates_equal(side: f64, separation: f64) -> f64 {
    if side < 1e-30 {
        return 0.0;
    }
    let h_ratio = separation / side;
    if h_ratio < 1e-10 {
        return 1.0;
    }
    ((h_ratio * h_ratio + 1.0).sqrt() - h_ratio).clamp(0.0, 1.0)
}
/// Radiation heat transfer coefficient (linearized Stefan-Boltzmann).
///
/// `h_r = epsilon * sigma * (T1^2 + T2^2) * (T1 + T2)`
pub fn radiation_heat_transfer_coefficient(t1: f64, t2: f64, emissivity: f64, sigma: f64) -> f64 {
    emissivity * sigma * (t1 * t1 + t2 * t2) * (t1 + t2)
}
/// Gray-body heat exchange between two parallel surfaces.
///
/// `q = sigma * (T1^4 - T2^4) / (1/eps1 + 1/eps2 - 1) * A`
pub fn gray_body_heat_exchange(
    t1: f64,
    t2: f64,
    eps1: f64,
    eps2: f64,
    area: f64,
    sigma: f64,
) -> f64 {
    let denom = 1.0 / eps1 + 1.0 / eps2 - 1.0;
    if denom.abs() < 1e-30 {
        return 0.0;
    }
    sigma * (t1.powi(4) - t2.powi(4)) / denom * area
}
/// Convective thermal resistance: `R_conv = 1 / (h * A)` \[K/W\].
pub fn convective_resistance(h: f64, area: f64) -> f64 {
    if (h * area).abs() < 1e-30 {
        return f64::INFINITY;
    }
    1.0 / (h * area)
}
/// Conductive thermal resistance: `R_cond = L / (k * A)` \[K/W\].
pub fn conductive_resistance(k: f64, area: f64, length: f64) -> f64 {
    if (k * area).abs() < 1e-30 {
        return f64::INFINITY;
    }
    length / (k * area)
}
/// Total thermal resistance in series: `R_total = sum(R_i)`.
pub fn series_thermal_resistance(resistances: &[f64]) -> f64 {
    resistances.iter().sum()
}
/// Total thermal resistance in parallel: `1/R_total = sum(1/R_i)`.
pub fn parallel_thermal_resistance(resistances: &[f64]) -> f64 {
    let inv_sum: f64 = resistances
        .iter()
        .map(|&r| if r.abs() > 1e-30 { 1.0 / r } else { 0.0 })
        .sum();
    if inv_sum.abs() < 1e-30 {
        return f64::INFINITY;
    }
    1.0 / inv_sum
}
/// Overall heat transfer coefficient U for a flat wall with fluid on both sides.
///
/// `U = 1 / (1/h1 + L/k + 1/h2)`
///
/// # Arguments
/// * `h1`  - convective coefficient on side 1 \[W/(m²K)\]
/// * `h2`  - convective coefficient on side 2 \[W/(m²K)\]
/// * `area` - area \[m²\] (for normalisation; result is per unit area here)
/// * `k`   - wall conductivity \[W/(mK)\]
/// * `l`   - wall thickness \[m\]
pub fn overall_heat_transfer_coefficient(h1: f64, h2: f64, area: f64, k: f64, l: f64) -> f64 {
    let _ = area;
    let r_total = 1.0 / h1 + l / k + 1.0 / h2;
    if r_total.abs() < 1e-30 {
        return 0.0;
    }
    1.0 / r_total
}
/// Heat exchanger effectiveness using the NTU method (counter-flow, Cr < 1).
///
/// `eps = (1 - exp(-NTU*(1-Cr))) / (1 - Cr * exp(-NTU*(1-Cr)))`
///
/// For `Cr = 1`: `eps = NTU / (1 + NTU)`.
pub fn ntu_effectiveness(ntu: f64, cr: f64) -> f64 {
    if (cr - 1.0).abs() < 1e-10 {
        ntu / (1.0 + ntu)
    } else {
        let exp_term = (-ntu * (1.0 - cr)).exp();
        (1.0 - exp_term) / (1.0 - cr * exp_term)
    }
}
/// Log Mean Temperature Difference (LMTD) for a heat exchanger.
///
/// `LMTD = (ΔT1 - ΔT2) / ln(ΔT1/ΔT2)`
///
/// For counter-flow: `ΔT1 = T_h,in - T_c,out`, `ΔT2 = T_h,out - T_c,in`.
pub fn log_mean_temperature_difference(
    t_h_in: f64,
    t_c_in: f64,
    t_h_out: f64,
    t_c_out: f64,
) -> f64 {
    let dt1 = t_h_in - t_c_out;
    let dt2 = t_h_out - t_c_in;
    if (dt1 - dt2).abs() < 1e-8 {
        return (dt1 + dt2) / 2.0;
    }
    if dt1 <= 0.0 || dt2 <= 0.0 {
        return (dt1.abs() + dt2.abs()) / 2.0;
    }
    (dt1 - dt2) / (dt1 / dt2).ln()
}
/// Compute the Biot number for a body with convection.
///
/// Bi = h * L_c / k
///
/// If Bi << 0.1, the lumped capacitance method is valid.
///
/// # Arguments
/// * `h`        - convection coefficient (W/(m^2 K))
/// * `l_c`      - characteristic length = V / A_s (m)
/// * `k`        - thermal conductivity (W/(m K))
pub fn biot_number(h: f64, l_c: f64, k: f64) -> f64 {
    if k.abs() < 1e-30 {
        return f64::INFINITY;
    }
    h * l_c / k
}
/// Check whether the lumped capacitance assumption is valid.
///
/// Returns `true` if Bi < 0.1 (standard criterion).
pub fn lumped_capacitance_valid(h: f64, l_c: f64, k: f64) -> bool {
    biot_number(h, l_c, k) < 0.1
}
/// Lumped capacitance transient solution: T(t) = T_inf + (T_0 - T_inf) * exp(-t / tau)
///
/// where tau = rho * c_p * V / (h * A_s) = m * c_p / (h * A_s)
///
/// # Arguments
/// * `t`        - time (s)
/// * `t_0`      - initial temperature (K)
/// * `t_inf`    - ambient temperature (K)
/// * `tau`      - thermal time constant (s)
pub fn lumped_capacitance_temperature(t: f64, t_0: f64, t_inf: f64, tau: f64) -> f64 {
    if tau.abs() < 1e-30 {
        return t_inf;
    }
    t_inf + (t_0 - t_inf) * (-t / tau).exp()
}
/// Thermal time constant for lumped capacitance: tau = rho * cp * V / (h * As).
pub fn thermal_time_constant(rho: f64, cp: f64, volume: f64, h: f64, area_surface: f64) -> f64 {
    let denom = h * area_surface;
    if denom.abs() < 1e-30 {
        return f64::INFINITY;
    }
    rho * cp * volume / denom
}
/// Thermal strain vector (Voigt notation) for isotropic material.
///
/// epsilon_th = alpha * (T - T_ref) * \[1, 1, 1, 0, 0, 0\]^T
///
/// # Arguments
/// * `alpha`  - coefficient of thermal expansion (1/K)
/// * `t`      - current temperature (K)
/// * `t_ref`  - reference (stress-free) temperature (K)
pub fn thermal_strain_isotropic(alpha: f64, t: f64, t_ref: f64) -> [f64; 6] {
    let eps_th = alpha * (t - t_ref);
    [eps_th, eps_th, eps_th, 0.0, 0.0, 0.0]
}
/// Thermal stress for isotropic material: sigma_th = -C * epsilon_th.
///
/// For isotropic material (plane-stress simplification for diagonal terms):
/// sigma_th = -E * alpha * (T - T_ref) / (1 - 2*nu)  (hydrostatic)
///
/// Returns the equivalent thermal pressure (negative = compressive under heating).
pub fn thermal_stress_isotropic(e_modulus: f64, nu: f64, alpha: f64, t: f64, t_ref: f64) -> f64 {
    let factor = e_modulus * alpha * (t - t_ref) / (1.0 - 2.0 * nu);
    -factor
}
/// Thermal load vector for a 1-D rod element due to temperature change.
///
/// f_th = E * A * alpha * delta_T * \[-1, 1\]
///
/// # Arguments
/// * `e_modulus` - Young's modulus (Pa)
/// * `area`      - cross-section area (m^2)
/// * `alpha`     - CTE (1/K)
/// * `delta_t`   - temperature change (K)
pub fn thermal_load_vector_1d(e_modulus: f64, area: f64, alpha: f64, delta_t: f64) -> [f64; 2] {
    let f = e_modulus * area * alpha * delta_t;
    [-f, f]
}
/// Solve the transient heat equation C*dT/dt + K*T = f using Crank-Nicolson (theta=0.5).
///
/// System:
/// (C/dt + 0.5*K) T_{n+1} = (C/dt - 0.5*K) T_n + f
///
/// # Arguments
/// * `k_global`     - global conductance matrix n×n (dense)
/// * `c_lumped`     - lumped capacitance vector of length n
/// * `t_init`       - initial temperature vector
/// * `heat_sources` - constant heat source vector
/// * `dirichlet`    - Dirichlet BCs: (node, temperature)
/// * `dt`           - time step (s)
/// * `n_steps`      - number of time steps
pub fn crank_nicolson_transient(
    k_global: &[Vec<f64>],
    c_lumped: &[f64],
    t_init: &[f64],
    heat_sources: &[f64],
    dirichlet: &[(usize, f64)],
    dt: f64,
    n_steps: usize,
) -> CrankNicolsonResult {
    let n = t_init.len();
    assert_eq!(k_global.len(), n);
    assert_eq!(c_lumped.len(), n);
    assert_eq!(heat_sources.len(), n);
    let mut temperatures = t_init.to_vec();
    let mut history = Vec::with_capacity(n_steps + 1);
    history.push(temperatures.clone());
    for _step in 0..n_steps {
        let mut a_vec = vec![0.0_f64; n];
        let mut b_vec = vec![0.0_f64; n];
        let mut c_vec = vec![0.0_f64; n];
        let mut rhs = vec![0.0_f64; n];
        for i in 0..n {
            b_vec[i] = c_lumped[i] / dt + 0.5 * k_global[i][i];
            if i > 0 {
                a_vec[i] = 0.5 * k_global[i][i - 1];
            }
            if i < n - 1 {
                c_vec[i] = 0.5 * k_global[i][i + 1];
            }
            let mut kt_i = 0.0;
            for j in 0..n {
                kt_i += k_global[i][j] * temperatures[j];
            }
            rhs[i] = c_lumped[i] / dt * temperatures[i] - 0.5 * kt_i + heat_sources[i];
        }
        for &(node, temp) in dirichlet {
            a_vec[node] = 0.0;
            b_vec[node] = 1.0;
            c_vec[node] = 0.0;
            rhs[node] = temp;
        }
        temperatures = thomas_algorithm(&a_vec, &b_vec, &c_vec, &rhs);
        history.push(temperatures.clone());
    }
    CrankNicolsonResult {
        final_temperatures: temperatures,
        history,
    }
}
/// Enthalpy-method update for phase change in a single node.
///
/// Converts the nodal temperature to specific enthalpy, advances it by one
/// explicit step with the net heat input, and converts back to temperature.
///
/// The specific enthalpy (J/kg) is
///
/// `e(T) = cp * T + L * f_l(T)`,
///
/// with the liquid fraction `f_l(T)` rising linearly from 0 at the solidus
/// `T_melt - mush_width` to 1 at the liquidus `T_melt + mush_width`.  One step
/// adds the heat `q_net * dt` to the node mass `rho * volume`:
///
/// `e_new = e(T) + q_net * dt / (rho * volume)`.
///
/// `e(T)` is continuous, strictly increasing and piecewise linear (kinks at the
/// solidus and the liquidus), so `T_new = e^{-1}(e_new)` is evaluated exactly
/// on the branch that contains `e_new`; there is no iteration that could stall
/// or cycle across the kinks.  For `mush_width <= 0` (sharp interface) `e(T)`
/// jumps by `L` at `T_melt` and the temperature stays at `T_melt` while `e_new`
/// lies on that plateau.
///
/// # Arguments
/// * `temperature`  - current nodal temperature (K)
/// * `q_net`        - net heat supply rate (W)
/// * `rho`          - density (kg/m^3)
/// * `cp`           - specific heat (J/(kg K))
/// * `latent_heat`  - latent heat L (J/kg)
/// * `t_melt`       - melting temperature (K)
/// * `mush_width`   - mushy zone half-width (K)
/// * `dt`           - time step (s)
/// * `volume`       - node volume (m^3)
pub fn enthalpy_update(
    temperature: f64,
    q_net: f64,
    rho: f64,
    cp: f64,
    latent_heat: f64,
    t_melt: f64,
    mush_width: f64,
    dt: f64,
    volume: f64,
) -> f64 {
    let half_width = mush_width.max(0.0);
    let t_sol = t_melt - half_width;
    let t_liq = t_melt + half_width;
    let cp_safe = cp.max(1e-30);
    let liquid_fraction = |t: f64| -> f64 {
        if half_width > 0.0 {
            ((t - t_sol) / (t_liq - t_sol)).clamp(0.0, 1.0)
        } else if t > t_melt {
            1.0
        } else {
            0.0
        }
    };
    let e = cp_safe * temperature + latent_heat * liquid_fraction(temperature);
    let mass = (rho * volume).max(1e-30);
    let e_new = e + q_net * dt / mass;
    let e_solidus = cp_safe * t_sol;
    let e_liquidus = cp_safe * t_liq + latent_heat;
    if e_new <= e_solidus {
        e_new / cp_safe
    } else if e_new >= e_liquidus {
        (e_new - latent_heat) / cp_safe
    } else if half_width > 0.0 {
        let mush_slope = cp_safe + latent_heat / (t_liq - t_sol);
        t_sol + (e_new - e_solidus) / mush_slope
    } else {
        t_melt
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermal::*;
    #[test]
    fn test_diffusivity() {
        let elem = ThermalElement::new(1.0, 1.0, 1.0, 1.0);
        assert!((elem.diffusivity() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_steady_state_linear_profile() {
        let n = 5;
        let length = 1.0_f64;
        let mut mesh = ThermalMesh1D::new_uniform(n, length, 1.0, 1.0, 1.0);
        mesh.set_temperature_bc(0, 0.0);
        mesh.set_temperature_bc(n - 1, 100.0);
        let q = vec![0.0_f64; n];
        mesh.steady_state(&q);
        for (i, &t) in mesh.temperatures.iter().enumerate() {
            let expected = 100.0 * i as f64 / (n - 1) as f64;
            assert!(
                (t - expected).abs() < 1e-9,
                "node {i}: got {t}, expected {expected}"
            );
        }
    }
    #[test]
    fn test_explicit_uniform_stays_uniform() {
        let n = 5;
        let mut mesh = ThermalMesh1D::new_uniform(n, 1.0, 1.0, 1.0, 1.0);
        for t in mesh.temperatures.iter_mut() {
            *t = 50.0;
        }
        let q = vec![0.0_f64; n];
        mesh.step_explicit(0.01, &q);
        for (i, &t) in mesh.temperatures.iter().enumerate() {
            assert!(
                (t - 50.0).abs() < 1e-12,
                "node {i}: temperature changed to {t}"
            );
        }
    }
    #[test]
    fn test_thomas_algorithm_3x3() {
        let a = [0.0, -1.0, -1.0];
        let b = [2.0, 2.0, 2.0];
        let c = [-1.0, -1.0, 0.0];
        let d = [1.0, 0.0, 1.0];
        let x = thomas_algorithm(&a, &b, &c, &d);
        for (i, &xi) in x.iter().enumerate() {
            assert!((xi - 1.0).abs() < 1e-12, "x[{i}] = {xi}, expected 1.0");
        }
    }
    #[test]
    fn test_heat_flux_linear_profile() {
        let n = 4;
        let length = 1.0_f64;
        let k = 2.0_f64;
        let mut mesh = ThermalMesh1D::new_uniform(n, length, k, 1.0, 1.0);
        mesh.set_temperature_bc(0, 0.0);
        mesh.set_temperature_bc(n - 1, 30.0);
        let q = vec![0.0_f64; n];
        mesh.steady_state(&q);
        let expected_flux = -k * 30.0 / length;
        for e in 0..mesh.elements.len() {
            let flux = mesh.heat_flux_at(e);
            assert!(
                (flux - expected_flux).abs() < 1e-9,
                "element {e}: flux = {flux}, expected {expected_flux}"
            );
        }
    }
    #[test]
    fn test_backward_euler_converges_to_steady_state() {
        let n = 5;
        let mut mesh = ThermalMesh1D::new_uniform(n, 1.0, 1.0, 1.0, 1.0);
        mesh.set_temperature_bc(0, 0.0);
        mesh.set_temperature_bc(n - 1, 100.0);
        let q = vec![0.0_f64; n];
        for _ in 0..1000 {
            mesh.step_implicit(0.1, &q);
        }
        for (i, &t) in mesh.temperatures.iter().enumerate() {
            let expected = 100.0 * i as f64 / (n - 1) as f64;
            assert!(
                (t - expected).abs() < 0.1,
                "node {i}: got {t}, expected {expected}"
            );
        }
    }
    #[test]
    fn test_convection_bc_heat_flux() {
        let conv = ConvectionBC::new(0, 10.0, 0.01, 300.0);
        let q = conv.heat_flux(400.0);
        assert!((q - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_convection_bc_at_ambient() {
        let conv = ConvectionBC::new(0, 10.0, 0.01, 300.0);
        let q = conv.heat_flux(300.0);
        assert!(q.abs() < 1e-12);
    }
    #[test]
    fn test_convection_bc_stiffness_and_load() {
        let conv = ConvectionBC::new(0, 10.0, 0.01, 300.0);
        assert!((conv.stiffness_contribution() - 0.1).abs() < 1e-12);
        assert!((conv.load_contribution() - 30.0).abs() < 1e-12);
    }
    #[test]
    fn test_radiation_bc_heat_flux() {
        let rad = RadiationBC::new(0, 0.8, 0.01, 300.0);
        let q = rad.heat_flux(500.0);
        assert!(q > 0.0, "radiation should emit heat when T > T_surr");
        let q2 = rad.heat_flux(200.0);
        assert!(q2 < 0.0, "radiation should absorb heat when T < T_surr");
    }
    #[test]
    fn test_radiation_bc_zero_at_equilibrium() {
        let rad = RadiationBC::new(0, 0.8, 0.01, 300.0);
        let q = rad.heat_flux(300.0);
        assert!(q.abs() < 1e-10, "no radiation flux at equilibrium");
    }
    #[test]
    fn test_radiation_linearized_coefficient() {
        let rad = RadiationBC::new(0, 1.0, 1.0, 300.0);
        let h = rad.linearized_coefficient(300.0);
        let expected = 4.0 * 5.670374419e-8 * 300.0_f64.powi(3);
        assert!((h - expected).abs() / expected < 1e-10);
    }
    #[test]
    fn test_thermal_contact_resistance_heat_flux() {
        let tcr = ThermalContactResistance::new(0, 1, 100.0, 0.01);
        let q = tcr.heat_flux(500.0, 300.0);
        assert!((q - 200.0).abs() < 1e-10);
    }
    #[test]
    fn test_thermal_contact_resistance_inverse() {
        let tcr = ThermalContactResistance::new(0, 1, 100.0, 0.01);
        let g = tcr.conductance();
        let r = tcr.resistance();
        assert!((g * r - 1.0).abs() < 1e-12, "G * R should be 1");
    }
    #[test]
    fn test_steady_state_with_convection() {
        let n = 5;
        let mut mesh = ThermalMesh1D::new_uniform(n, 1.0, 1.0, 1.0, 1.0);
        mesh.set_temperature_bc(0, 100.0);
        mesh.add_convection_bc(n - 1, 10.0, 1.0, 20.0);
        let q = vec![0.0_f64; n];
        mesh.steady_state(&q);
        assert!((mesh.temperatures[0] - 100.0).abs() < 1e-9);
        assert!(
            mesh.temperatures[n - 1] < 100.0,
            "convection end should be cooler"
        );
        assert!(mesh.temperatures[n - 1] > 20.0, "should not reach ambient");
    }
    #[test]
    fn test_steady_state_with_contact_resistance() {
        let n = 3;
        let mut mesh = ThermalMesh1D::new_uniform(n, 1.0, 1.0, 1.0, 1.0);
        mesh.set_temperature_bc(0, 100.0);
        mesh.set_temperature_bc(n - 1, 0.0);
        mesh.add_contact_resistance(1, 2, 0.5, 1.0);
        let q = vec![0.0_f64; n];
        mesh.steady_state(&q);
        assert!(
            mesh.temperatures[1] > 0.0 && mesh.temperatures[1] < 100.0,
            "T[1] = {}",
            mesh.temperatures[1]
        );
    }
    #[test]
    fn test_total_heat_content() {
        let n = 3;
        let mut mesh = ThermalMesh1D::new_uniform(n, 1.0, 1.0, 1.0, 1.0);
        for t in mesh.temperatures.iter_mut() {
            *t = 100.0;
        }
        let q_total = mesh.total_heat_content();
        assert!(
            q_total > 0.0,
            "heat content should be positive for positive temperatures"
        );
    }
    #[test]
    fn test_max_min_temperature() {
        let n = 5;
        let mut mesh = ThermalMesh1D::new_uniform(n, 1.0, 1.0, 1.0, 1.0);
        mesh.set_temperature_bc(0, 0.0);
        mesh.set_temperature_bc(n - 1, 100.0);
        let q = vec![0.0; n];
        mesh.steady_state(&q);
        assert!((mesh.min_temperature() - 0.0).abs() < 1e-9);
        assert!((mesh.max_temperature() - 100.0).abs() < 1e-9);
    }
    #[test]
    fn test_cfl_limit() {
        let mesh = ThermalMesh1D::new_uniform(5, 1.0, 1.0, 1.0, 1.0);
        let dt_crit = mesh.cfl_limit();
        assert!(dt_crit > 0.0);
        assert!(dt_crit.is_finite());
        assert!((dt_crit - 0.03125).abs() < 1e-10);
    }
    #[test]
    fn test_theta_method_backward_euler() {
        let n = 5;
        let mut mesh = ThermalMesh1D::new_uniform(n, 1.0, 1.0, 1.0, 1.0);
        mesh.set_temperature_bc(0, 0.0);
        mesh.set_temperature_bc(n - 1, 100.0);
        let q = vec![0.0; n];
        let history = transient_theta_method(&mut mesh, 0.1, 100, 1.0, &q);
        assert_eq!(history.len(), 101);
        for (i, &t) in mesh.temperatures.iter().enumerate() {
            let expected = 100.0 * i as f64 / (n - 1) as f64;
            assert!(
                (t - expected).abs() < 0.5,
                "node {i}: got {t}, expected {expected}"
            );
        }
    }
    #[test]
    fn test_theta_method_crank_nicolson() {
        let n = 5;
        let mut mesh = ThermalMesh1D::new_uniform(n, 1.0, 1.0, 1.0, 1.0);
        mesh.set_temperature_bc(0, 0.0);
        mesh.set_temperature_bc(n - 1, 100.0);
        let q = vec![0.0; n];
        let history = transient_theta_method(&mut mesh, 0.01, 500, 0.5, &q);
        assert_eq!(history.len(), 501);
        for (i, &t) in mesh.temperatures.iter().enumerate() {
            let expected = 100.0 * i as f64 / (n - 1) as f64;
            assert!(
                (t - expected).abs() < 1.0,
                "node {i}: got {t}, expected {expected}"
            );
        }
    }
    #[test]
    fn test_theta_method_preserves_bc() {
        let n = 5;
        let mut mesh = ThermalMesh1D::new_uniform(n, 1.0, 1.0, 1.0, 1.0);
        mesh.set_temperature_bc(0, 50.0);
        mesh.set_temperature_bc(n - 1, 200.0);
        let q = vec![0.0; n];
        let _history = transient_theta_method(&mut mesh, 0.1, 10, 1.0, &q);
        assert!(
            (mesh.temperatures[0] - 50.0).abs() < 1e-9,
            "BC at node 0 should be preserved"
        );
        assert!(
            (mesh.temperatures[n - 1] - 200.0).abs() < 1e-9,
            "BC at node n-1 should be preserved"
        );
    }
    #[test]
    fn test_stefan_front_moves_forward() {
        let front0 = 0.0;
        let front1 = stefan_front_position(front0, 100.0, 300.0, 0.1, 1.0);
        assert!(
            front1 > front0,
            "Stefan front should move: {front0} → {front1}"
        );
    }
    #[test]
    fn test_stefan_front_zero_heat_flux() {
        let front0 = 0.5;
        let front1 = stefan_front_position(front0, 0.0, 300.0, 0.1, 1.0);
        assert!(
            (front1 - front0).abs() < 1e-12,
            "Zero flux: front should not move"
        );
    }
    #[test]
    fn test_stefan_neumann_solution_positive() {
        let lambda = stefan_neumann_lambda(1.0, 50.0, 100.0);
        assert!(lambda > 0.0, "Neumann lambda should be positive: {lambda}");
    }
    #[test]
    fn test_stefan_neumann_lambda_st_0_1_is_0_2200() {
        // Water-like melt: cp = 1000 J/(kg K), superheat 33.4 K, L = 334 kJ/kg
        // gives St = cp dT / L = 0.1 and lambda = 0.2200163 (mpmath root of
        // lambda exp(lambda^2) erf(lambda) = St / sqrt(pi)).
        let lambda = stefan_neumann_lambda(1000.0, 33.4, 334_000.0);
        assert!(
            (lambda - 0.220_016_272_743).abs() < 1e-9,
            "lambda = {lambda}"
        );
        assert!((lambda - 0.2200).abs() < 1e-4);
    }
    #[test]
    fn test_stefan_neumann_lambda_matches_reference_roots() {
        // (St, lambda) pairs from a 30-digit mpmath root solve.
        let reference = [
            (0.01, 0.070_593_276_560_429_3),
            (0.05, 0.156_820_922_334_547_3),
            (0.2, 0.306_423_905_361_211_2),
            (0.5, 0.464_785_920_646_244_4),
            (1.0, 0.620_062_633_313_595_5),
            (2.0, 0.800_601_362_805_608_3),
            (5.0, 1.059_687_014_281_902_4),
        ];
        for (st, expected) in reference {
            // cp * dT / L = St with cp = 2000, L = 100 kJ/kg
            let delta_t = st * 100_000.0 / 2000.0;
            let lambda = stefan_neumann_lambda(2000.0, delta_t, 100_000.0);
            assert!(
                (lambda - expected).abs() < 1e-11,
                "St = {st}: lambda = {lambda}, expected {expected}"
            );
        }
    }
    #[test]
    fn test_stefan_neumann_lambda_satisfies_transcendental_equation() {
        let sqrt_pi = std::f64::consts::PI.sqrt();
        for st in [1e-6, 1e-3, 0.03, 0.7, 3.0, 20.0, 400.0] {
            let lambda = stefan_neumann_lambda(1.0, st, 1.0);
            let lhs = lambda * (lambda * lambda).exp() * erf_series(lambda);
            let rhs = st / sqrt_pi;
            assert!(
                (lhs - rhs).abs() <= 1e-12 * rhs,
                "St = {st}: lambda = {lambda}, lhs = {lhs}, rhs = {rhs}"
            );
        }
    }
    #[test]
    fn test_stefan_neumann_lambda_uses_cp_delta_t_over_latent_heat() {
        // Only the Stefan number cp dT / L matters; the old code evaluated
        // cp / L and ignored the temperature difference entirely.
        let base = stefan_neumann_lambda(4000.0, 10.0, 200_000.0);
        let double_cp = stefan_neumann_lambda(8000.0, 10.0, 200_000.0);
        let double_dt = stefan_neumann_lambda(4000.0, 20.0, 200_000.0);
        let half_l = stefan_neumann_lambda(4000.0, 10.0, 100_000.0);
        assert!((double_cp - double_dt).abs() < 1e-14);
        assert!((double_cp - half_l).abs() < 1e-14);
        assert!(double_dt > base);
        assert!(
            (stefan_neumann_lambda(4000.0, 10.0, 200_000.0) - stefan_neumann_lambda(1.0, 0.2, 1.0))
                .abs()
                < 1e-14
        );
    }
    #[test]
    fn test_stefan_neumann_lambda_small_stefan_number_limit() {
        // lambda -> sqrt(St / 2) as St -> 0.
        let st = 1e-8;
        let lambda = stefan_neumann_lambda(1.0, st, 1.0);
        assert!(
            (lambda / (0.5 * st).sqrt() - 1.0).abs() < 1e-6,
            "lambda = {lambda}"
        );
    }
    #[test]
    fn test_stefan_neumann_lambda_is_increasing_in_superheat() {
        let mut previous = 0.0;
        for k in 1..=40 {
            let lambda = stefan_neumann_lambda(1000.0, 5.0 * k as f64, 300_000.0);
            assert!(lambda > previous, "k = {k}: {lambda} <= {previous}");
            previous = lambda;
        }
    }
    #[test]
    fn test_stefan_neumann_lambda_zero_superheat_is_zero() {
        assert_eq!(stefan_neumann_lambda(1000.0, 0.0, 334_000.0), 0.0);
        assert_eq!(stefan_neumann_lambda(0.0, 25.0, 334_000.0), 0.0);
    }
    #[test]
    fn test_stefan_neumann_front_position_follows_two_lambda_sqrt_alpha_t() {
        // St = 0.1 (lambda = 0.2200163), alpha = 1e-7 m^2/s, one hour.
        let s = stefan_neumann_front_position(1000.0, 33.4, 334_000.0, 1.0e-7, 3600.0);
        let expected = 2.0 * 0.220_016_272_743 * (1.0e-7_f64 * 3600.0).sqrt();
        assert!((s - expected).abs() < 1e-11, "s = {s}, expected {expected}");
        // The front starts at the boundary and grows like sqrt(t).
        assert_eq!(
            stefan_neumann_front_position(1000.0, 33.4, 334_000.0, 1.0e-7, 0.0),
            0.0
        );
        let s4 = stefan_neumann_front_position(1000.0, 33.4, 334_000.0, 1.0e-7, 4.0 * 3600.0);
        assert!((s4 / s - 2.0).abs() < 1e-12);
        // No superheat, no motion; invalid inputs propagate NaN.
        assert_eq!(
            stefan_neumann_front_position(1000.0, 0.0, 334_000.0, 1.0e-7, 3600.0),
            0.0
        );
        assert!(stefan_neumann_front_position(1000.0, 33.4, 0.0, 1.0e-7, 3600.0).is_nan());
        assert!(stefan_neumann_front_position(1000.0, 33.4, 334_000.0, -1.0, 3600.0).is_nan());
        assert!(stefan_neumann_front_position(1000.0, 33.4, 334_000.0, 1.0e-7, -1.0).is_nan());
        assert!(stefan_neumann_front_position(1000.0, 33.4, 334_000.0, f64::NAN, 1.0).is_nan());
    }
    #[test]
    fn test_stefan_neumann_lambda_rejects_invalid_input() {
        assert!(stefan_neumann_lambda(1000.0, 10.0, 0.0).is_nan());
        assert!(stefan_neumann_lambda(1000.0, 10.0, -334_000.0).is_nan());
        assert!(stefan_neumann_lambda(1000.0, -10.0, 334_000.0).is_nan());
        assert!(stefan_neumann_lambda(f64::NAN, 10.0, 334_000.0).is_nan());
        assert!(stefan_neumann_lambda(1000.0, f64::INFINITY, 334_000.0).is_nan());
        assert!(stefan_neumann_lambda(1000.0, 10.0, f64::INFINITY).is_nan());
    }
    #[test]
    fn test_erf_series_matches_known_values() {
        // High-precision reference values of erf.
        let reference = [
            (0.0, 0.0),
            (0.1, 0.112_462_916_018_284_9),
            (0.5, 0.520_499_877_813_046_5),
            (1.0, 0.842_700_792_949_714_9),
            (2.0, 0.995_322_265_018_952_7),
            (3.0, 0.999_977_909_503_001_4),
        ];
        for (x, expected) in reference {
            assert!((erf_series(x) - expected).abs() < 1e-15, "erf({x})");
            assert!((erf_series(-x) + expected).abs() < 1e-15, "erf(-{x})");
        }
        assert_eq!(erf_series(7.0), 1.0);
        assert_eq!(erf_series(-7.0), -1.0);
        assert!(erf_series(f64::NAN).is_nan());
    }
    #[test]
    fn test_latent_heat_source_at_melting_point() {
        let qs = latent_heat_source(0.0, 273.15, 1.0, 334000.0);
        assert!(qs.is_finite(), "Latent heat source = {qs}");
    }
    #[test]
    fn test_latent_heat_source_away_from_melting() {
        let qs = latent_heat_source(100.0, 273.15, 1.0, 334000.0);
        assert!(qs.abs() < 1e-20, "No latent heat far from melting: {qs}");
    }
    #[test]
    fn test_effective_specific_heat_at_melting() {
        let cp_eff = effective_specific_heat_melting(273.15, 1000.0, 334000.0, 0.1, 273.15);
        assert!(
            cp_eff > 1000.0,
            "Effective Cp should be boosted at melting: {cp_eff}"
        );
    }
    #[test]
    fn test_solidification_front_advances() {
        let s0 = 0.01;
        let s1 = solidification_front_advance(s0, 1.0, 80.0, 1000.0, 334000.0, 1e-3, 0.001);
        assert!(s1 > s0, "Solidification front should advance: {s0} → {s1}");
    }
    #[test]
    fn test_solidification_fraction_range() {
        let fl = solidification_fraction(280.0, 273.15, 3.0);
        assert!(
            (0.0..=1.0).contains(&fl),
            "Fraction should be in [0,1]: {fl}"
        );
    }
    #[test]
    fn test_solidification_fraction_below_solidus() {
        let fl = solidification_fraction(250.0, 273.15, 3.0);
        assert!(
            fl.abs() < 1e-12,
            "Below solidus: fully solid (fl=0), fl = {fl}"
        );
    }
    #[test]
    fn test_solidification_fraction_above_liquidus() {
        let fl = solidification_fraction(290.0, 273.15, 0.01);
        assert!(
            (fl - 1.0).abs() < 1e-10,
            "Above liquidus: fully liquid (fl=1), fl = {fl}"
        );
    }
    #[test]
    fn test_view_factor_parallel_plates_between_zero_and_one() {
        let vf = view_factor_parallel_plates_equal(0.5, 1.0);
        assert!(
            vf > 0.0 && vf <= 1.0,
            "View factor should be in (0,1]: {vf}"
        );
    }
    #[test]
    fn test_view_factor_large_plates_approaches_one() {
        let vf = view_factor_parallel_plates_equal(1000.0, 1.0);
        assert!((vf - 1.0).abs() < 0.01, "Large plate ratio → VF ≈ 1: {vf}");
    }
    #[test]
    fn test_radiation_heat_transfer_coefficient_positive() {
        let h_r = radiation_heat_transfer_coefficient(500.0, 400.0, 0.8, 5.670374419e-8);
        assert!(
            h_r > 0.0,
            "Radiation heat transfer coefficient should be positive: {h_r}"
        );
    }
    #[test]
    fn test_gray_body_heat_exchange_positive() {
        let q = gray_body_heat_exchange(600.0, 400.0, 0.8, 0.9, 1.0, 5.670374419e-8);
        assert!(q > 0.0, "Heat should flow from hot to cold: {q}");
    }
    #[test]
    fn test_gray_body_heat_exchange_zero_at_equilibrium() {
        let q = gray_body_heat_exchange(400.0, 400.0, 0.8, 0.9, 1.0, 5.670374419e-8);
        assert!(q.abs() < 1e-10, "No net heat exchange at equilibrium: {q}");
    }
    #[test]
    fn test_convective_resistance_positive() {
        let r = convective_resistance(50.0, 0.1);
        assert!(r > 0.0, "Convective resistance should be positive: {r}");
    }
    #[test]
    fn test_conductive_resistance_positive() {
        let r = conductive_resistance(1.0, 0.01, 1.0);
        assert!(r > 0.0, "Conductive resistance should be positive: {r}");
    }
    #[test]
    fn test_total_thermal_resistance_sum() {
        let r_conv = convective_resistance(50.0, 0.1);
        let r_cond = conductive_resistance(1.0, 0.01, 1.0);
        let r_total = series_thermal_resistance(&[r_conv, r_cond]);
        assert!(
            (r_total - r_conv - r_cond).abs() < 1e-12,
            "Series: {r_total}"
        );
    }
    #[test]
    fn test_overall_heat_transfer_coefficient() {
        let u = overall_heat_transfer_coefficient(50.0, 50.0, 1.0, 0.01, 1.0);
        assert!(u > 0.0, "Overall U should be positive: {u}");
    }
    #[test]
    fn test_ntu_method_effectiveness() {
        let eps = ntu_effectiveness(2.0, 0.5);
        assert!(
            eps > 0.0 && eps < 1.0,
            "Effectiveness should be in (0,1): {eps}"
        );
    }
    #[test]
    fn test_log_mean_temperature_difference() {
        let lmtd = log_mean_temperature_difference(80.0, 20.0, 60.0, 10.0);
        assert!(lmtd > 0.0, "LMTD should be positive: {lmtd}");
    }
    #[test]
    fn test_log_mean_temperature_equal_diffs_gives_that_value() {
        let lmtd = log_mean_temperature_difference(100.0, 30.0, 80.0, 50.0);
        assert!(
            (lmtd - 50.0).abs() < 1e-6,
            "Equal ΔT LMTD = {lmtd}, expected 50"
        );
    }
}
