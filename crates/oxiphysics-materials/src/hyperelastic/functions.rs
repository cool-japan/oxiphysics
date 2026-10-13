//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{MooneyRivlin, Ogden};

/// Helper: compute invariants from a 3x3 deformation gradient F.
pub(super) fn deformation_invariants(f: &[[f64; 3]; 3]) -> (f64, f64, f64) {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = (0..3).map(|k| f[k][i] * f[k][j]).sum();
        }
    }
    let i1 = c[0][0] + c[1][1] + c[2][2];
    let c2_trace: f64 = (0..3)
        .flat_map(|i| (0..3).map(move |k| c[i][k] * c[k][i]))
        .sum();
    let i2 = 0.5 * (i1 * i1 - c2_trace);
    let det_f = f[0][0] * (f[1][1] * f[2][2] - f[1][2] * f[2][1])
        - f[0][1] * (f[1][0] * f[2][2] - f[1][2] * f[2][0])
        + f[0][2] * (f[1][0] * f[2][1] - f[1][1] * f[2][0]);
    let i3 = det_f * det_f;
    (i1, i2, i3)
}
/// Helper: determinant of a 3x3 matrix.
pub(super) fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
/// Helper: inverse of a 3x3 matrix.
pub(super) fn inv3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let d = det3(m);
    if d.abs() < 1e-30 {
        return [[0.0; 3]; 3];
    }
    let inv_d = 1.0 / d;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_d,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_d,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_d,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_d,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_d,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_d,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_d,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_d,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_d,
        ],
    ]
}
/// Helper: transpose of a 3x3 matrix.
pub(super) fn transpose3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}
/// Compute the principal stretches λ₁ ≥ λ₂ ≥ λ₃ ≥ 0 of a deformation gradient F.
///
/// The principal stretches are the singular values of F, i.e. the square roots
/// of the eigenvalues of the right Cauchy-Green tensor C = FᵀF.  They are
/// obtained with a one-sided (Hestenes) Jacobi iteration that orthogonalises
/// the columns of F by plane rotations and takes the column norms.  Working on
/// F instead of C avoids squaring the condition number and, unlike the
/// closed-form cubic solution, keeps full accuracy for repeated stretches
/// (uniaxial, equibiaxial and hydrostatic states).
pub(super) fn principal_stretches(f: &[[f64; 3]; 3]) -> (f64, f64, f64) {
    const MAX_SWEEPS: usize = 40;
    const COLUMN_PAIRS: [(usize, usize); 3] = [(0, 1), (0, 2), (1, 2)];
    let mut a = *f;
    for _ in 0..MAX_SWEEPS {
        let mut rotated = false;
        for &(p, q) in COLUMN_PAIRS.iter() {
            let alpha: f64 = a.iter().map(|row| row[p] * row[p]).sum();
            let beta: f64 = a.iter().map(|row| row[q] * row[q]).sum();
            let gamma: f64 = a.iter().map(|row| row[p] * row[q]).sum();
            if gamma == 0.0 || gamma.abs() <= f64::EPSILON * (alpha * beta).sqrt() {
                continue;
            }
            rotated = true;
            // Rotation angle that annihilates the inner product of columns p and q.
            let zeta = (beta - alpha) / (2.0 * gamma);
            let t = zeta.signum() / (zeta.abs() + (1.0 + zeta * zeta).sqrt());
            let c = 1.0 / (1.0 + t * t).sqrt();
            let s = c * t;
            for row in a.iter_mut() {
                let (x_p, x_q) = (row[p], row[q]);
                row[p] = c * x_p - s * x_q;
                row[q] = s * x_p + c * x_q;
            }
        }
        if !rotated {
            break;
        }
    }
    let mut stretches = [0.0_f64; 3];
    for (j, stretch) in stretches.iter_mut().enumerate() {
        *stretch = a.iter().map(|row| row[j] * row[j]).sum::<f64>().sqrt();
    }
    stretches.sort_by(|x, y| y.total_cmp(x));
    (stretches[0], stretches[1], stretches[2])
}
/// Predict uniaxial Cauchy stress for a Neo-Hookean material (special case of Ogden α=2).
///
/// For incompressible uniaxial tension: λ2 = λ3 = 1/√λ1.
/// σ = μ*(λ1² - λ1⁻¹)
pub fn neo_hookean_uniaxial_stress(mu: f64, lambda: f64) -> f64 {
    mu * (lambda * lambda - 1.0 / lambda)
}
/// Predict equibiaxial Cauchy stress for a Neo-Hookean material.
///
/// Incompressible biaxial: λ1=λ2=λ, λ3=1/λ².
/// σ = μ*(λ² - λ⁻⁴)
pub fn neo_hookean_biaxial_stress(mu: f64, lambda: f64) -> f64 {
    mu * (lambda * lambda - lambda.powi(-4))
}
/// Uniaxial Mooney-Rivlin Cauchy stress.
///
/// σ = 2*(C10 + C01/λ)*(λ² - 1/λ)
pub fn mooney_rivlin_uniaxial_stress(c10: f64, c01: f64, lambda: f64) -> f64 {
    2.0 * (c10 + c01 / lambda) * (lambda * lambda - 1.0 / lambda)
}
/// Biaxial Mooney-Rivlin Cauchy stress.
///
/// σ_biaxial = 2*(C10 + λ²·C01)*(λ² - λ⁻⁴)
pub fn mooney_rivlin_biaxial_stress(c10: f64, c01: f64, lambda: f64) -> f64 {
    2.0 * (c10 + lambda * lambda * c01) * (lambda * lambda - lambda.powi(-4))
}
/// Cauchy stress from strain energy W using finite difference dW/dλ.
///
/// For uniaxial incompressible: σ = λ * dW/dλ (first Piola-Kirchhoff / lambda approx).
pub fn cauchy_stress_uniaxial_fd<F: Fn(f64) -> f64>(w_fn: F, lambda: f64) -> f64 {
    let h = 1e-6;
    let dw = (w_fn(lambda + h) - w_fn(lambda - h)) / (2.0 * h);
    lambda * dw
}
/// Check incompressibility: returns `|J - 1|` where J = det(F).
pub fn incompressibility_error(f: &[[f64; 3]; 3]) -> f64 {
    (det3(f) - 1.0).abs()
}
/// Build a uniaxial-tension deformation gradient for stretch λ along X axis.
///
/// Assumes incompressibility: λ2 = λ3 = 1/√λ.
pub fn uniaxial_deformation_gradient(lambda: f64) -> [[f64; 3]; 3] {
    let inv_sqrt_lambda = 1.0 / lambda.sqrt();
    [
        [lambda, 0.0, 0.0],
        [0.0, inv_sqrt_lambda, 0.0],
        [0.0, 0.0, inv_sqrt_lambda],
    ]
}
/// Build an equibiaxial-tension deformation gradient.
///
/// λ1=λ2=λ, λ3=1/λ² (incompressible).
pub fn biaxial_deformation_gradient(lambda: f64) -> [[f64; 3]; 3] {
    let lambda3 = 1.0 / (lambda * lambda);
    [[lambda, 0.0, 0.0], [0.0, lambda, 0.0], [0.0, 0.0, lambda3]]
}
/// Cavitation criterion for incompressible hyperelastic materials.
///
/// Under spherically symmetric (triaxial) tension the material can undergo
/// sudden cavitation when the applied pressure `p_cav` is reached.
///
/// For a Neo-Hookean material (shear modulus μ):
///
///   p_cav = 5μ/2   (Ball 1982)
///
/// For a general Ogden material the result is computed numerically.
pub fn neo_hookean_cavitation_pressure(shear_modulus: f64) -> f64 {
    2.5 * shear_modulus
}
/// Estimate cavitation pressure for a Mooney-Rivlin material.
///
/// p_cav ≈ 5/2 (C10 + C01)   (leading order, same structure as Neo-Hookean).
pub fn mooney_rivlin_cavitation_pressure(c10: f64, c01: f64) -> f64 {
    2.5 * (c10 + c01)
}
/// Check the Drucker stability condition for a Neo-Hookean material
/// under uniaxial stretch λ.
///
/// The material is Drucker stable if dσ/dε > 0 everywhere in the
/// loading path.  For the incompressible Neo-Hookean:
///
///   dσ/dλ = μ(2λ + λ⁻²) > 0  for all λ > 0.
pub fn neo_hookean_drucker_stable(shear_modulus: f64, lambda: f64) -> bool {
    let d_sigma_d_lambda = shear_modulus * (2.0 * lambda + lambda.powi(-2));
    d_sigma_d_lambda > 0.0
}
/// Check the Baker-Ericksen inequality for principal Cauchy stresses.
///
/// The inequality requires: (σᵢ - σⱼ)(λᵢ - λⱼ) ≥ 0 for i ≠ j.
pub fn baker_ericksen_satisfied(stresses: &[f64; 3], stretches: &[f64; 3]) -> bool {
    for i in 0..3 {
        for j in 0..3 {
            if i != j {
                let ds = stresses[i] - stresses[j];
                let dl = stretches[i] - stretches[j];
                if ds * dl < -1e-12 {
                    return false;
                }
            }
        }
    }
    true
}
/// Compute the uniaxial Cauchy stress for an incompressible Ogden material.
///
/// For incompressible uniaxial tension (λ₂ = λ₃ = 1/√λ):
///
///   σ = Σ_p μ_p (λ^α_p - λ^(-α_p/2))
pub fn ogden_uniaxial_stress_incompressible(ogden: &Ogden, lambda: f64) -> f64 {
    ogden
        .mu
        .iter()
        .zip(ogden.alpha.iter())
        .map(|(&mu_p, &alpha_p)| mu_p * (lambda.powf(alpha_p) - lambda.powf(-alpha_p / 2.0)))
        .sum()
}
/// Compute the equibiaxial Cauchy stress for an incompressible Ogden material.
///
/// λ₁ = λ₂ = λ, λ₃ = λ⁻²:
///
///   σ = Σ_p μ_p (λ^α_p - λ^(-2 α_p))
pub fn ogden_biaxial_stress_incompressible(ogden: &Ogden, lambda: f64) -> f64 {
    ogden
        .mu
        .iter()
        .zip(ogden.alpha.iter())
        .map(|(&mu_p, &alpha_p)| mu_p * (lambda.powf(alpha_p) - lambda.powf(-2.0 * alpha_p)))
        .sum()
}
/// Compute numerical tangent stiffness for a hyperelastic model (volumetric).
///
/// Uses finite differences on the strain energy density to estimate
/// the tangent bulk modulus K = -V * dP/dV at a given stretch λ_vol.
///
/// For isotropic uniform dilation λ₁=λ₂=λ₃=λ, J = λ³.
pub fn tangent_bulk_modulus_fd<F: Fn(f64) -> f64>(w_fn: F, lambda_vol: f64) -> f64 {
    let h = lambda_vol * 1e-5;
    let lp = lambda_vol + h;
    let lm = lambda_vol - h;
    let dwdlam = (w_fn(lp) - w_fn(lm)) / (2.0 * h);
    let j = lambda_vol * lambda_vol * lambda_vol;
    let p_hi = -(w_fn(lp + h) - w_fn(lp - h)) / (2.0 * h) / (3.0 * lp * lp);
    let p_lo = -(w_fn(lm + h) - w_fn(lm - h)) / (2.0 * h) / (3.0 * lm * lm);
    let dj = 3.0 * lambda_vol * lambda_vol * 2.0 * h;
    let _ = dwdlam;
    -j * (p_hi - p_lo) / dj
}
/// Compute numerical tangent shear modulus for a hyperelastic model.
///
/// Uses finite differences dσ/dγ around zero shear.
/// The model is queried via a simple shear deformation gradient.
pub fn tangent_shear_modulus_fd<F: Fn(&[[f64; 3]; 3]) -> f64>(w_fn: F, gamma: f64) -> f64 {
    let h = 1e-6;
    let f_plus = [[1.0, gamma + h, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let f_minus = [[1.0, gamma - h, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    (w_fn(&f_plus) - w_fn(&f_minus)) / (2.0 * h)
}
/// Incompressible Cauchy stress from strain energy via FD (triaxial).
///
/// Returns the 3 principal Cauchy stresses (σ₁, σ₂, σ₃) for principal
/// stretches (λ₁, λ₂, λ₃) using finite differences on W.
pub fn principal_cauchy_stress_incompressible<F: Fn(f64, f64, f64) -> f64>(
    w_fn: &F,
    lambda1: f64,
    lambda2: f64,
    lambda3: f64,
    lagrange_multiplier: f64,
) -> [f64; 3] {
    let h = 1e-6;
    let dw1 =
        (w_fn(lambda1 + h, lambda2, lambda3) - w_fn(lambda1 - h, lambda2, lambda3)) / (2.0 * h);
    let dw2 =
        (w_fn(lambda1, lambda2 + h, lambda3) - w_fn(lambda1, lambda2 - h, lambda3)) / (2.0 * h);
    let dw3 =
        (w_fn(lambda1, lambda2, lambda3 + h) - w_fn(lambda1, lambda2, lambda3 - h)) / (2.0 * h);
    let j = lambda1 * lambda2 * lambda3;
    [
        lambda1 * dw1 / j - lagrange_multiplier,
        lambda2 * dw2 / j - lagrange_multiplier,
        lambda3 * dw3 / j - lagrange_multiplier,
    ]
}
/// Solve for the Lagrange multiplier for incompressible uniaxial loading.
///
/// For uniaxial tension (σ₂ = σ₃ = 0):
/// p = λ₂ * ∂W/∂λ₂ / J
pub fn lagrange_multiplier_uniaxial<F: Fn(f64, f64, f64) -> f64>(w_fn: &F, lambda1: f64) -> f64 {
    let lambda2 = 1.0 / lambda1.sqrt();
    let lambda3 = lambda2;
    let j = lambda1 * lambda2 * lambda3;
    let h = 1e-6;
    let dw2 =
        (w_fn(lambda1, lambda2 + h, lambda3) - w_fn(lambda1, lambda2 - h, lambda3)) / (2.0 * h);
    lambda2 * dw2 / j
}
/// Full 4th-order constitutive tangent of the Neo-Hookean material in
/// Voigt notation (6×6 matrix, flat row-major).
///
/// The spatial (Eulerian) tangent moduli are:
///   c_ijkl = μ/J * (δ_ik δ_jl + δ_il δ_jk) - 2μ/J * (b̄_ijkl - ...)
///            + K * (2J-1) * δ_ij δ_kl - K * (J-1) * (δ_ij δ_kl - ...)
///
/// In the simplified linearized form for the isotropic Neo-Hookean:
///   C_ep (Voigt 6×6) = λ_eff * (I⊗I) + 2 μ_eff * I_sym
///
/// where λ_eff = K - 2G/3,  μ_eff = G = μ/J (push-forward),
/// and J = det(F).
///
/// For the reference configuration (material tangent):
///   λ = K - 2/3 * μ,  using current J.
pub fn neo_hookean_constitutive_tangent(
    shear_modulus: f64,
    bulk_modulus: f64,
    f: &[[f64; 3]; 3],
) -> [f64; 36] {
    let j = det3(f);
    let j_safe = if j.abs() < 1e-30 { 1.0 } else { j };
    let mu_eff = shear_modulus / j_safe;
    let k_eff = bulk_modulus;
    let lam = k_eff - 2.0 / 3.0 * mu_eff;
    let mut c = [0.0_f64; 36];
    for i in 0..3 {
        for j_idx in 0..3 {
            c[i * 6 + j_idx] = lam;
        }
        c[i * 6 + i] += 2.0 * mu_eff;
    }
    c[3 * 6 + 3] = mu_eff;
    c[4 * 6 + 4] = mu_eff;
    c[5 * 6 + 5] = mu_eff;
    c
}
/// Compute the Hessian of the Mooney-Rivlin strain energy W with respect to
/// the Green-Lagrange strain tensor components (Voigt notation).
///
/// The Hessian W,CC gives the incremental stress response — this is the
/// material tangent modulus in the reference configuration.
///
/// For the Mooney-Rivlin model W = C10*(Ī₁-3) + C01*(Ī₂-3) + D1*(J-1)²:
///   d²W/dC² is approximated via central finite differences with step h.
///
/// Returns the 6×6 symmetric Hessian in Voigt notation (flat row-major).
pub fn mooney_rivlin_strain_energy_hessian(
    mr: &MooneyRivlin,
    c_voigt: &[f64; 6],
    h: f64,
) -> [f64; 36] {
    let voigt_to_c3 = |v: &[f64; 6]| -> [[f64; 3]; 3] {
        [[v[0], v[5], v[4]], [v[5], v[1], v[3]], [v[4], v[3], v[2]]]
    };
    let w_from_c = |v: &[f64; 6]| -> f64 {
        let c3 = voigt_to_c3(v);
        let i1 = c3[0][0] + c3[1][1] + c3[2][2];
        let i2 = c3[0][0] * c3[1][1] + c3[1][1] * c3[2][2] + c3[0][0] * c3[2][2]
            - c3[0][1] * c3[0][1]
            - c3[1][2] * c3[1][2]
            - c3[0][2] * c3[0][2];
        let i3 = c3[0][0] * (c3[1][1] * c3[2][2] - c3[1][2] * c3[2][1])
            - c3[0][1] * (c3[1][0] * c3[2][2] - c3[1][2] * c3[2][0])
            + c3[0][2] * (c3[1][0] * c3[2][1] - c3[1][1] * c3[2][0]);
        let j = i3.max(0.0).sqrt();
        let j_23 = if j > 1e-30 { j.powf(-2.0 / 3.0) } else { 1.0 };
        let i1_bar = j_23 * i1;
        let i2_bar = j_23 * j_23 * i2;
        let d1 = mr.bulk_modulus / 2.0;
        mr.c10 * (i1_bar - 3.0) + mr.c01 * (i2_bar - 3.0) + d1 * (j - 1.0).powi(2)
    };
    let w0 = w_from_c(c_voigt);
    let mut hess = [0.0_f64; 36];
    for i in 0..6 {
        for j_idx in 0..6 {
            let d2w = if i == j_idx {
                let mut cpi = *c_voigt;
                cpi[i] += h;
                let mut cmi = *c_voigt;
                cmi[i] -= h;
                (w_from_c(&cpi) - 2.0 * w0 + w_from_c(&cmi)) / (h * h)
            } else {
                let mut cpp = *c_voigt;
                cpp[i] += h;
                cpp[j_idx] += h;
                let mut cmm = *c_voigt;
                cmm[i] -= h;
                cmm[j_idx] -= h;
                let mut cpm = *c_voigt;
                cpm[i] += h;
                cpm[j_idx] -= h;
                let mut cmp = *c_voigt;
                cmp[i] -= h;
                cmp[j_idx] += h;
                (w_from_c(&cpp) - w_from_c(&cpm) - w_from_c(&cmp) + w_from_c(&cmm)) / (4.0 * h * h)
            };
            hess[i * 6 + j_idx] = d2w;
        }
    }
    hess
}
/// Compute principal Cauchy stresses from principal stretches for an Ogden material.
///
/// For the strain energy
///   W = Σ_p (μ_p/α_p) (λ̄₁^α_p + λ̄₂^α_p + λ̄₃^α_p − 3) + K/2 (J − 1)²,
/// with the isochoric stretches λ̄_i = J^(−1/3) λ_i and J = λ₁λ₂λ₃, the
/// principal Cauchy stresses σ_i = λ_i (∂W/∂λ_i) / J are
///
///   σ_i = (1/J) Σ_p μ_p (λ̄_i^α_p − ⅓ Σ_j λ̄_j^α_p) + K (J − 1).
///
/// The mean subtracted in each term is the deviatoric projection of the
/// isochoric response, so the stress vanishes in the reference configuration
/// λ = (1, 1, 1) and the isochoric part is traceless in the sense
/// Σ_i (σ_i − K(J − 1)) = 0.
///
/// For an incompressible state (J = 1) the stress differences are
/// σ_i − σ_j = Σ_p μ_p (λ_i^α_p − λ_j^α_p); in particular a uniaxial test
/// (λ, λ^(−1/2), λ^(−1/2)) gives σ₁ − σ₂ = Σ_p μ_p (λ^α_p − λ^(−α_p/2)).
///
/// # Arguments
/// * `lambda1, lambda2, lambda3` — principal stretches
///
/// # Returns
/// `[sigma1, sigma2, sigma3]` — principal Cauchy stresses (Pa)
pub fn ogden_principal_stresses(
    ogden: &Ogden,
    lambda1: f64,
    lambda2: f64,
    lambda3: f64,
) -> [f64; 3] {
    let j = lambda1 * lambda2 * lambda3;
    let j13 = if j > 1e-30 { j.powf(-1.0 / 3.0) } else { 1.0 };
    let l = [j13 * lambda1, j13 * lambda2, j13 * lambda3];
    let mut sigma_iso = [0.0_f64; 3];
    for (mu_p, alpha_p) in ogden.mu.iter().zip(ogden.alpha.iter()) {
        let powers = [
            l[0].powf(*alpha_p),
            l[1].powf(*alpha_p),
            l[2].powf(*alpha_p),
        ];
        let mean_power = (powers[0] + powers[1] + powers[2]) / 3.0;
        for (sigma_k, power_k) in sigma_iso.iter_mut().zip(powers.iter()) {
            *sigma_k += mu_p / j * (power_k - mean_power);
        }
    }
    let p_vol = ogden.bulk_modulus * (j - 1.0);
    [
        sigma_iso[0] + p_vol,
        sigma_iso[1] + p_vol,
        sigma_iso[2] + p_vol,
    ]
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::DruckerPrager;
    use crate::J2Plasticity;
    use crate::JwlEos;
    use crate::MieGruneisenEos;
    use crate::hyperelastic::ArrudaBoyce;
    use crate::hyperelastic::BilinearCohesiveZone;
    use crate::hyperelastic::BlatzKo;
    use crate::hyperelastic::Fung;
    use crate::hyperelastic::Gent;
    use crate::hyperelastic::Hencky;
    use crate::hyperelastic::HolzapfelGasserOgden;
    use crate::hyperelastic::HolzapfelOgden;
    use crate::hyperelastic::Varga;
    use crate::hyperelastic::Yeoh;
    #[test]
    fn test_mooney_rivlin_identity_zero_stress() {
        let mr = MooneyRivlin::new(0.5e6, 0.1e6, 1.0e9);
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let p = mr.first_piola_kirchhoff_stress(&identity);
        for (i, row) in p.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    val.abs() < 1.0,
                    "P[{i}][{j}] = {val} should be ~0 at identity"
                );
            }
        }
    }
    #[test]
    fn test_mooney_rivlin_shear_modulus() {
        let mr = MooneyRivlin::new(0.5e6, 0.1e6, 1.0e9);
        assert!((mr.shear_modulus() - 1.2e6).abs() < 1.0);
    }
    #[test]
    fn test_ogden_identity_zero_energy() {
        let og = Ogden::one_term(1.0e6, 2.0, 1.0e9);
        let energy = og.strain_energy_from_stretches(1.0, 1.0, 1.0);
        assert!(
            energy.abs() < 1.0e-6,
            "Energy at identity should be ~0, got {energy}"
        );
    }
    #[test]
    fn test_ogden_shear_modulus() {
        let og = Ogden::one_term(2.0e6, 2.0, 1.0e9);
        assert!((og.shear_modulus() - 2.0e6).abs() < 1.0);
    }
    #[test]
    fn test_j2_von_mises_uniaxial() {
        let stress = [100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let vm = J2Plasticity::von_mises_stress(&stress);
        assert!(
            (vm - 100.0e6).abs() < 1.0e3,
            "VM should be 100 MPa, got {vm}"
        );
    }
    #[test]
    fn test_j2_elastic_no_yield() {
        let j2 = J2Plasticity::new(200.0e9, 0.3, 250.0e6, 1.0e9);
        let trial = [100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let (stress, dp) = j2.return_mapping(&trial, 0.0);
        assert!((dp).abs() < 1e-10, "Should be elastic");
        assert!((stress[0] - 100.0e6).abs() < 1.0);
    }
    #[test]
    fn test_j2_plastic_return() {
        let j2 = J2Plasticity::new(200.0e9, 0.3, 250.0e6, 1.0e9);
        let trial = [500.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let (_stress, dp) = j2.return_mapping(&trial, 0.0);
        assert!(dp > 0.0, "Should yield");
    }
    #[test]
    fn test_drucker_prager_yield() {
        let dp = DruckerPrager::new(30.0_f64.to_radians(), 1.0e6, 30.0e9, 0.2);
        let hydrostatic = [-100.0e6, -100.0e6, -100.0e6, 0.0, 0.0, 0.0];
        let f = dp.yield_function(&hydrostatic);
        assert!(
            f < 0.0,
            "Hydrostatic compression should be below yield, F={f}"
        );
    }
    #[test]
    fn test_jwl_eos() {
        let jwl = JwlEos::new(3.712e11, 3.231e9, 4.15, 0.95, 0.30, 1630.0);
        let p = jwl.pressure(1630.0, 0.0);
        assert!(p.is_finite(), "JWL pressure should be finite");
    }
    #[test]
    fn test_mie_gruneisen_reference() {
        let mg = MieGruneisenEos::new(2700.0, 5386.0, 1.339, 2.0);
        let p = mg.pressure(2700.0, 0.0);
        assert!(p.abs() < 1.0, "Pressure at reference should be ~0, got {p}");
    }
    #[test]
    fn test_principal_stretches_identity() {
        let f = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (l1, l2, l3) = principal_stretches(&f);
        assert!((l1 - 1.0).abs() < 1e-6);
        assert!((l2 - 1.0).abs() < 1e-6);
        assert!((l3 - 1.0).abs() < 1e-6);
    }
    fn matmul3(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let mut c = [[0.0; 3]; 3];
        for (i, row) in c.iter_mut().enumerate() {
            for (j, entry) in row.iter_mut().enumerate() {
                *entry = (0..3).map(|k| a[i][k] * b[k][j]).sum();
            }
        }
        c
    }
    /// Rotation about a (not necessarily unit) axis by `angle` radians (Rodrigues).
    fn rotation(axis: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
        let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        let (x, y, z) = (axis[0] / norm, axis[1] / norm, axis[2] / norm);
        let (s, c) = angle.sin_cos();
        let t = 1.0 - c;
        [
            [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
            [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
            [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
        ]
    }
    fn assert_stretches(got: (f64, f64, f64), want: [f64; 3], tol: f64) {
        let got = [got.0, got.1, got.2];
        for (g, w) in got.iter().zip(want.iter()) {
            assert!(
                (g - w).abs() < tol,
                "principal stretches {got:?} differ from {want:?}"
            );
        }
    }
    #[test]
    fn test_principal_stretches_diagonal_gradient() {
        let f = [[1.2, 0.0, 0.0], [0.0, 0.9, 0.0], [0.0, 0.0, 1.1]];
        // Returned in descending order.
        assert_stretches(principal_stretches(&f), [1.2, 1.1, 0.9], 1e-14);
    }
    #[test]
    fn test_principal_stretches_uniaxial_gradient() {
        let f = uniaxial_deformation_gradient(1.5);
        let lateral = 1.0 / 1.5_f64.sqrt();
        let (l1, l2, l3) = principal_stretches(&f);
        assert_stretches((l1, l2, l3), [1.5, lateral, lateral], 1e-14);
        assert!((l1 * l2 * l3 - 1.0).abs() < 1e-14, "J must be 1");
    }
    #[test]
    fn test_principal_stretches_equal_stretches_are_not_averaged() {
        // The former closed-form fallback returned I1/3 = 1.0176 for all three
        // principal stretches of F = 1.2 I.
        let f = [[1.2, 0.0, 0.0], [0.0, 1.2, 0.0], [0.0, 0.0, 1.2]];
        assert_stretches(principal_stretches(&f), [1.2, 1.2, 1.2], 1e-14);
    }
    #[test]
    fn test_principal_stretches_simple_shear() {
        let gamma = 0.3_f64;
        let f = [[1.0, gamma, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (l1, l2, l3) = principal_stretches(&f);
        assert!(
            (l1 * l2 * l3 - 1.0).abs() < 1e-14,
            "simple shear is isochoric"
        );
        // tr(FᵀF) = 3 + γ².
        let sum_sq = l1 * l1 + l2 * l2 + l3 * l3;
        assert!(
            (sum_sq - (3.0 + gamma * gamma)).abs() < 1e-14,
            "Σλ² = {sum_sq}"
        );
        // Eigenvalues of the in-plane block [[1, γ], [γ, 1 + γ²]] of C, plus the free direction.
        let mean = 1.0 + 0.5 * gamma * gamma;
        let radius = (0.25 * gamma.powi(4) + gamma * gamma).sqrt();
        assert_stretches(
            (l1, l2, l3),
            [(mean + radius).sqrt(), 1.0, (mean - radius).sqrt()],
            1e-14,
        );
    }
    #[test]
    fn test_principal_stretches_invariant_under_rotations() {
        let d = [[1.2, 0.0, 0.0], [0.0, 0.9, 0.0], [0.0, 0.0, 1.1]];
        let r1 = rotation([1.0, 2.0, 3.0], 0.7);
        let r2 = rotation([-2.0, 0.5, 1.0], 1.9);
        let expected = [1.2, 1.1, 0.9];
        // F = D R: the columns of F are not orthogonal, so FᵀF is a full matrix.
        assert_stretches(principal_stretches(&matmul3(&d, &r1)), expected, 1e-13);
        // F = R1 D R2 (general polar decomposition).
        let general = matmul3(&r1, &matmul3(&d, &r2));
        assert_stretches(principal_stretches(&general), expected, 1e-13);
        // Wide range of stretches.
        let wide = [[10.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.1]];
        let wide_f = matmul3(&r2, &matmul3(&wide, &r1));
        assert_stretches(principal_stretches(&wide_f), [10.0, 1.0, 0.1], 1e-12);
    }
    #[test]
    fn test_principal_stretches_match_cauchy_green_invariants() {
        let f = [[1.1, 0.2, -0.1], [0.05, 0.95, 0.3], [0.0, -0.2, 1.05]];
        let (l1, l2, l3) = principal_stretches(&f);
        assert!(l1 >= l2 && l2 >= l3 && l3 > 0.0, "{l1} {l2} {l3}");
        let (i1, i2, i3) = deformation_invariants(&f);
        let (s1, s2, s3) = (l1 * l1, l2 * l2, l3 * l3);
        assert!((s1 + s2 + s3 - i1).abs() < 1e-13);
        assert!((s1 * s2 + s2 * s3 + s1 * s3 - i2).abs() < 1e-13);
        assert!((s1 * s2 * s3 - i3).abs() < 1e-13);
        assert!((l1 * l2 * l3 - det3(&f).abs()).abs() < 1e-13);
    }
    #[test]
    fn test_principal_stretches_rank_deficient_gradient() {
        let f = [[2.0, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.0]];
        assert_stretches(principal_stretches(&f), [2.0, 0.5, 0.0], 1e-15);
    }
    #[test]
    fn test_isochoric_energies_from_uniaxial_gradient_are_nonzero() {
        // Regression: the old eigenvalue fallback returned (I1/3, I1/3, I1/3) for repeated
        // stretches, which carries no isochoric strain energy in any isotropic model.
        let lambda = 1.5_f64;
        let lateral = 1.0 / lambda.sqrt();
        let f = uniaxial_deformation_gradient(lambda);
        let yeoh = Yeoh::new(0.5e6, 0.0, 0.0, 1.0e9);
        let i1_bar = lambda * lambda + 2.0 / lambda;
        let expected = 0.5e6 * (i1_bar - 3.0);
        let w = yeoh.strain_energy_density(&f);
        assert!(
            (w - expected).abs() / expected < 1e-9,
            "Yeoh W = {w}, expected {expected}"
        );
        let gent = Gent::new(1.0e6, 50.0, 1.0e9);
        let w_gent = gent.strain_energy_from_stretches(lambda, lateral, lateral);
        assert!(w_gent > 0.0);
        assert!((gent.strain_energy_density(&f) - w_gent).abs() / w_gent < 1e-9);
        let ab = ArrudaBoyce::new(1.0e6, 7.0, 1.0e9);
        let w_ab = ab.strain_energy_from_stretches(lambda, lateral, lateral);
        assert!(w_ab > 0.0);
        assert!((ab.strain_energy_density(&f) - w_ab).abs() / w_ab < 1e-9);
    }
    #[test]
    fn test_arruda_boyce_identity_zero_energy() {
        let ab = ArrudaBoyce::new(1.0e6, 7.0, 1.0e9);
        let w = ab.strain_energy_from_stretches(1.0, 1.0, 1.0);
        assert!(
            w.abs() < 1.0e-3,
            "AB energy at identity should be ~0, got {w}"
        );
    }
    #[test]
    fn test_arruda_boyce_positive_for_tension() {
        let ab = ArrudaBoyce::new(1.0e6, 7.0, 1.0e9);
        let w = ab.strain_energy_from_stretches(1.5, 1.0 / 1.5_f64.sqrt(), 1.0 / 1.5_f64.sqrt());
        assert!(w > 0.0, "AB energy should be positive under tension");
    }
    #[test]
    fn test_arruda_boyce_shear_modulus() {
        let ab = ArrudaBoyce::new(3.0e6, 5.0, 1.0e9);
        assert!((ab.shear_modulus() - 3.0e6).abs() < 1.0);
    }
    #[test]
    fn test_arruda_boyce_from_deformation_gradient() {
        let ab = ArrudaBoyce::new(1.0e6, 7.0, 1.0e9);
        let f = uniaxial_deformation_gradient(1.2);
        let w = ab.strain_energy_density(&f);
        assert!(w > 0.0, "AB energy from F should be positive");
    }
    #[test]
    fn test_gent_identity_zero_energy() {
        let g = Gent::new(1.0e6, 50.0, 1.0e9);
        let w = g.strain_energy_from_stretches(1.0, 1.0, 1.0);
        assert!(
            w.abs() < 1.0e-3,
            "Gent energy at identity should be ~0, got {w}"
        );
    }
    #[test]
    fn test_gent_energy_increases_with_stretch() {
        let g = Gent::new(1.0e6, 50.0, 1.0e9);
        let lambda1 = 1.2_f64;
        let inv1 = 1.0 / lambda1.sqrt();
        let w1 = g.strain_energy_from_stretches(lambda1, inv1, inv1);
        let lambda2 = 1.5_f64;
        let inv2 = 1.0 / lambda2.sqrt();
        let w2 = g.strain_energy_from_stretches(lambda2, inv2, inv2);
        assert!(w2 > w1, "Gent energy should increase with stretch");
    }
    #[test]
    fn test_gent_shear_modulus() {
        let g = Gent::new(2.0e6, 100.0, 1.0e9);
        assert!((g.shear_modulus() - 2.0e6).abs() < 1.0);
    }
    #[test]
    fn test_yeoh_identity_zero_energy() {
        let y = Yeoh::new(0.5e6, -0.01e6, 0.001e6, 1.0e9);
        let w = y.strain_energy_from_stretches(1.0, 1.0, 1.0);
        assert!(w.abs() < 1.0e-3, "Yeoh at identity should be ~0, got {w}");
    }
    #[test]
    fn test_yeoh_shear_modulus() {
        let y = Yeoh::new(0.5e6, -0.01e6, 0.001e6, 1.0e9);
        assert!((y.shear_modulus() - 1.0e6).abs() < 1.0);
    }
    #[test]
    fn test_yeoh_energy_from_f() {
        let y = Yeoh::new(0.5e6, 0.0, 0.0, 1.0e9);
        let f = uniaxial_deformation_gradient(1.5);
        let w = y.strain_energy_density(&f);
        assert!(w > 0.0, "Yeoh energy should be positive");
    }
    #[test]
    fn test_neo_hookean_uniaxial_stress_positive() {
        let sigma = neo_hookean_uniaxial_stress(1.0e6, 1.5);
        assert!(sigma > 0.0, "uniaxial stress should be positive: {sigma}");
    }
    #[test]
    fn test_neo_hookean_biaxial_stress_positive() {
        let sigma = neo_hookean_biaxial_stress(1.0e6, 1.3);
        assert!(sigma > 0.0, "biaxial stress should be positive: {sigma}");
    }
    #[test]
    fn test_neo_hookean_zero_stress_at_identity() {
        let sigma = neo_hookean_uniaxial_stress(1.0e6, 1.0);
        assert!(sigma.abs() < 1e-6, "stress at λ=1 should be ~0: {sigma}");
    }
    #[test]
    fn test_mooney_rivlin_uniaxial_stress() {
        let sigma = mooney_rivlin_uniaxial_stress(0.5e6, 0.1e6, 1.5);
        assert!(
            sigma > 0.0,
            "MR uniaxial stress should be positive: {sigma}"
        );
    }
    #[test]
    fn test_mooney_rivlin_biaxial_stress() {
        let sigma = mooney_rivlin_biaxial_stress(0.5e6, 0.1e6, 1.3);
        assert!(sigma > 0.0, "MR biaxial stress should be positive: {sigma}");
    }
    #[test]
    fn test_cauchy_stress_fd_neo_hookean() {
        let mu = 1.0e6;
        let lambda = 1.4_f64;
        let inv = 1.0 / lambda.sqrt();
        let w_fn = |lam: f64| {
            let i1 = lam * lam + 2.0 / lam;
            mu / 2.0 * (i1 - 3.0)
        };
        let sigma_fd = cauchy_stress_uniaxial_fd(w_fn, lambda);
        let sigma_exact = neo_hookean_uniaxial_stress(mu, lambda);
        let _ = inv;
        assert!(
            (sigma_fd - sigma_exact).abs() / sigma_exact.abs() < 0.01,
            "FD stress {sigma_fd} vs exact {sigma_exact}"
        );
    }
    #[test]
    fn test_uniaxial_gradient_incompressible() {
        let f = uniaxial_deformation_gradient(1.5);
        let err = incompressibility_error(&f);
        assert!(err < 1e-12, "J should be 1 for incompressible F, err={err}");
    }
    #[test]
    fn test_biaxial_gradient_incompressible() {
        let f = biaxial_deformation_gradient(1.3);
        let err = incompressibility_error(&f);
        assert!(err < 1e-12, "J should be 1 for biaxial F, err={err}");
    }
    #[test]
    fn test_incompressibility_error_identity() {
        let f = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!(incompressibility_error(&f) < 1e-15);
    }
    #[test]
    fn test_ogden_rubber_default_positive_energy() {
        let og = Ogden::rubber_default();
        let f = uniaxial_deformation_gradient(1.2);
        let w = og.strain_energy_density(&f);
        assert!(
            w > 0.0,
            "rubber default Ogden energy should be positive: {w}"
        );
    }
    #[test]
    fn test_hencky_identity_zero_energy() {
        let h = Hencky::from_young_poisson(200.0e9, 0.3);
        let w = h.strain_energy_from_stretches(1.0, 1.0, 1.0);
        assert!(w.abs() < 1e-6, "Hencky at identity should be ~0, got {w}");
    }
    #[test]
    fn test_hencky_young_poisson_roundtrip() {
        let e_in = 200.0e9;
        let nu_in = 0.3;
        let h = Hencky::from_young_poisson(e_in, nu_in);
        let e_out = h.young_modulus();
        let nu_out = h.poisson_ratio();
        assert!((e_out - e_in).abs() / e_in < 1e-10, "E={e_out}");
        assert!((nu_out - nu_in).abs() < 1e-10, "ν={nu_out}");
    }
    #[test]
    fn test_hencky_bulk_modulus_positive() {
        let h = Hencky::from_young_poisson(200.0e9, 0.3);
        assert!(h.bulk_modulus() > 0.0);
    }
    #[test]
    fn test_hencky_energy_positive_under_tension() {
        let h = Hencky::from_young_poisson(200.0e9, 0.3);
        let w = h.strain_energy_from_stretches(1.1, 0.99, 0.99);
        assert!(
            w > 0.0,
            "Hencky energy should be positive under tension, got {w}"
        );
    }
    #[test]
    fn test_hencky_principal_stress_identity() {
        let h = Hencky::from_young_poisson(200.0e9, 0.3);
        let sigma = h.principal_cauchy_stress(1.0, 1.0, 1.0);
        for s in sigma {
            assert!(s.abs() < 1e-3, "stress at identity should be ~0: {s}");
        }
    }
    #[test]
    fn test_hencky_energy_from_gradient() {
        let h = Hencky::from_young_poisson(200.0e9, 0.3);
        let f = uniaxial_deformation_gradient(1.05);
        let w = h.strain_energy_density(&f);
        assert!(w > 0.0, "Hencky energy from F should be positive: {w}");
    }
    #[test]
    fn test_blatz_ko_identity_zero_energy() {
        let bk = BlatzKo::foam_default(1.0e5);
        let w = bk.strain_energy_from_stretches(1.0, 1.0, 1.0);
        assert!(w.abs() < 1.0, "Blatz-Ko at identity should be ~0, got {w}");
    }
    #[test]
    fn test_blatz_ko_young_modulus() {
        let bk = BlatzKo::new(1.0e5, 0.5, 0.25);
        let e = bk.young_modulus();
        assert!((e - 2.5e5).abs() < 1.0, "E={e}");
    }
    #[test]
    fn test_blatz_ko_positive_energy_under_compression() {
        let bk = BlatzKo::foam_default(1.0e5);
        let w = bk.strain_energy_from_stretches(0.8, 0.8, 0.8);
        assert!(w >= 0.0, "Blatz-Ko energy under compression: {w}");
    }
    #[test]
    fn test_varga_identity_zero_energy() {
        let v = Varga::new(1.0e6, 1.0e9);
        let w = v.strain_energy_from_stretches(1.0, 1.0, 1.0);
        assert!(w.abs() < 1.0, "Varga at identity: {w}");
    }
    #[test]
    fn test_varga_uniaxial_stress_zero_at_identity() {
        let v = Varga::new(1.0e6, 1.0e9);
        let sigma = v.uniaxial_stress_incompressible(1.0);
        assert!(sigma.abs() < 1e-6, "σ(λ=1)={sigma}");
    }
    #[test]
    fn test_varga_uniaxial_stress_positive_for_tension() {
        let v = Varga::new(1.0e6, 1.0e9);
        let sigma = v.uniaxial_stress_incompressible(1.5);
        assert!(sigma > 0.0, "σ(λ=1.5)={sigma}");
    }
    #[test]
    fn test_varga_energy_positive_under_tension() {
        let v = Varga::new(1.0e6, 1.0e9);
        let w = v.strain_energy_from_stretches(1.5, 1.0 / 1.5_f64.sqrt(), 1.0 / 1.5_f64.sqrt());
        assert!(w > 0.0, "Varga energy under tension: {w}");
    }
    #[test]
    fn test_holzapfel_ogden_identity_positive_energy() {
        let ho = HolzapfelOgden::myocardium_default();
        let f_id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let w = ho.strain_energy_density(&f_id);
        assert!(w.abs() < 1.0, "W at identity should be ~0, got {w}");
    }
    #[test]
    fn test_holzapfel_ogden_fibre_active_tension() {
        let ho = HolzapfelOgden::myocardium_default();
        let f = uniaxial_deformation_gradient(1.2);
        assert!(
            ho.fibre_active(&f, 1),
            "fibre 1 should be active in tension"
        );
    }
    #[test]
    fn test_holzapfel_ogden_fibre_inactive_compression() {
        let ho = HolzapfelOgden::myocardium_default();
        let f = [[0.8, 0.0, 0.0], [0.0, 1.118, 0.0], [0.0, 0.0, 1.118]];
        assert!(
            !ho.fibre_active(&f, 1),
            "fibre 1 should be inactive in compression"
        );
    }
    #[test]
    fn test_holzapfel_ogden_energy_increases_with_stretch() {
        let ho = HolzapfelOgden::myocardium_default();
        let f1 = uniaxial_deformation_gradient(1.1);
        let f2 = uniaxial_deformation_gradient(1.3);
        let w1 = ho.strain_energy_density(&f1);
        let w2 = ho.strain_energy_density(&f2);
        assert!(w2 > w1, "W should increase with stretch: {w1} < {w2}");
    }
    #[test]
    fn test_neo_hookean_drucker_stable_for_positive_lambda() {
        let mu = 1.0e6;
        for &lam in &[0.5, 1.0, 2.0, 5.0] {
            assert!(neo_hookean_drucker_stable(mu, lam), "not stable at λ={lam}");
        }
    }
    #[test]
    fn test_neo_hookean_cavitation_pressure() {
        let mu = 1.0e6;
        let p_cav = neo_hookean_cavitation_pressure(mu);
        assert!((p_cav - 2.5e6).abs() < 1.0, "p_cav={p_cav}");
    }
    #[test]
    fn test_mooney_rivlin_cavitation_pressure() {
        let p = mooney_rivlin_cavitation_pressure(0.5e6, 0.1e6);
        assert!((p - 1.5e6).abs() < 1.0, "p_cav={p}");
    }
    #[test]
    fn test_baker_ericksen_isotropic_satisfied() {
        let stresses = [1.0e6, 1.0e6, 1.0e6];
        let stretches = [1.1, 1.1, 1.1];
        assert!(baker_ericksen_satisfied(&stresses, &stretches));
    }
    #[test]
    fn test_baker_ericksen_ordered_satisfied() {
        let stresses = [3.0e6, 2.0e6, 1.0e6];
        let stretches = [1.3, 1.2, 1.1];
        assert!(baker_ericksen_satisfied(&stresses, &stretches));
    }
    #[test]
    fn test_baker_ericksen_violated() {
        let stresses = [3.0e6, 2.0e6, 1.0e6];
        let stretches = [1.1, 1.2, 1.3];
        assert!(!baker_ericksen_satisfied(&stresses, &stretches));
    }
    #[test]
    fn test_j2_consistent_tangent_less_than_elastic() {
        let j2 = J2Plasticity::new(200.0e9, 0.3, 250.0e6, 2.0e9);
        let c_ep = j2.consistent_tangent_uniaxial();
        assert!(c_ep < j2.young_modulus, "C_ep should be < E");
        assert!(c_ep > 0.0, "C_ep should be positive");
    }
    #[test]
    fn test_j2_back_stress_increment() {
        let j2 = J2Plasticity::new(200.0e9, 0.3, 250.0e6, 2.0e9);
        let dep = [0.001, -0.0005, -0.0005, 0.0, 0.0, 0.0];
        let dalpha = j2.back_stress_increment(10.0e9, &dep);
        assert!(dalpha[0].abs() > 0.0, "back stress xx non-zero");
    }
    #[test]
    fn test_j2_equivalent_plastic_strain_increment() {
        let dep = [0.0, 0.0, 0.0, 0.001, 0.0, 0.0];
        let dep_bar = J2Plasticity::equivalent_plastic_strain_increment(&dep);
        assert!(dep_bar > 0.0, "Δε̄_p should be positive: {dep_bar}");
    }
    #[test]
    fn test_j2_plane_stress_elastic() {
        let j2 = J2Plasticity::new(200.0e9, 0.3, 250.0e6, 1.0e9);
        let ([sxx, syy, sxy], dp) = j2.return_mapping_plane_stress(100.0e6, 0.0, 0.0, 0.0);
        assert!((dp).abs() < 1e-15, "should be elastic");
        assert!((sxx - 100.0e6).abs() < 1.0);
        let _ = (syy, sxy);
    }
    #[test]
    fn test_j2_plane_stress_plastic() {
        let j2 = J2Plasticity::new(200.0e9, 0.3, 250.0e6, 1.0e9);
        let (_, dp) = j2.return_mapping_plane_stress(400.0e6, 0.0, 0.0, 0.0);
        assert!(dp > 0.0, "should yield");
    }
    fn plane_stress_von_mises(s: [f64; 3]) -> f64 {
        J2Plasticity::von_mises_stress(&[s[0], s[1], 0.0, s[2], 0.0, 0.0])
    }
    /// C P σ for the plane-stress elasticity matrix C and projection P (engineering shear).
    fn plane_stress_c_p_sigma(young: f64, nu: f64, s: [f64; 3]) -> [f64; 3] {
        let p = [
            (2.0 * s[0] - s[1]) / 3.0,
            (2.0 * s[1] - s[0]) / 3.0,
            2.0 * s[2],
        ];
        let c = young / (1.0 - nu * nu);
        [
            c * (p[0] + nu * p[1]),
            c * (nu * p[0] + p[1]),
            c * (1.0 - nu) / 2.0 * p[2],
        ]
    }
    fn assert_on_hardened_yield_surface(j2: &J2Plasticity, s: [f64; 3], ep_bar: f64, dp: f64) {
        let vm = plane_stress_von_mises(s);
        let sigma_y = j2.current_yield_stress(ep_bar + dp);
        assert!(
            (vm - sigma_y).abs() / sigma_y < 1e-9,
            "von Mises {vm} must equal the hardened yield stress {sigma_y}"
        );
    }
    #[test]
    fn test_j2_plane_stress_uniaxial_trial_returns_to_hardened_yield_surface() {
        let j2 = J2Plasticity::new(200.0e9, 0.3, 250.0e6, 1.0e9);
        let (s, dp) = j2
            .try_return_mapping_plane_stress(400.0e6, 0.0, 0.0, 0.0)
            .unwrap();
        assert!(dp > 0.0);
        assert_on_hardened_yield_surface(&j2, s, 0.0, dp);
        // The closest-point return of a uniaxial trial leaves a non-zero transverse stress.
        assert!(s[1].abs() > 1.0e6, "σ_yy = {}", s[1]);
        assert!(s[0] < 400.0e6 && s[0] > 250.0e6, "σ_xx = {}", s[0]);
        assert_eq!(s[2], 0.0);
        // Reference values from an independent bracketing (Brent) solution of the same
        // closest-point problem.
        assert!(
            (s[0] - 2.608432203e8).abs() <= 2e-9 * 2.608432203e8,
            "σ_xx = {}",
            s[0]
        );
        assert!(
            (s[1] - 2.162695755e7).abs() <= 2e-9 * 2.162695755e7,
            "σ_yy = {}",
            s[1]
        );
        assert!(
            (dp - 7.302646442e-4).abs() <= 2e-9 * 7.302646442e-4,
            "Δε̄_p = {dp}"
        );
    }
    #[test]
    fn test_j2_plane_stress_mixed_trial_with_prior_hardening() {
        let j2 = J2Plasticity::new(210.0e9, 0.29, 300.0e6, 2.0e9);
        let ep_bar = 0.002;
        let (s, dp) = j2
            .try_return_mapping_plane_stress(350.0e6, -200.0e6, 120.0e6, ep_bar)
            .unwrap();
        assert!(dp > 0.0);
        assert_on_hardened_yield_surface(&j2, s, ep_bar, dp);
    }
    #[test]
    fn test_j2_plane_stress_pure_shear_and_equibiaxial_trials() {
        let j2 = J2Plasticity::new(200.0e9, 0.3, 250.0e6, 1.0e9);
        // Pure shear: only the d2 branch of the update acts.
        let (s, dp) = j2
            .try_return_mapping_plane_stress(0.0, 0.0, 300.0e6, 0.0)
            .unwrap();
        assert_on_hardened_yield_surface(&j2, s, 0.0, dp);
        assert!(s[0].abs() < 1e-3 && s[1].abs() < 1e-3, "{s:?}");
        assert!(s[2] > 0.0 && s[2] < 300.0e6);
        // Equibiaxial: only the d1 branch acts and the stress stays equibiaxial.
        let (s, dp) = j2
            .try_return_mapping_plane_stress(350.0e6, 350.0e6, 0.0, 0.0)
            .unwrap();
        assert_on_hardened_yield_surface(&j2, s, 0.0, dp);
        assert!((s[0] - s[1]).abs() < 1e-3 * s[0].abs());
        assert!(s[2].abs() < 1e-6);
    }
    #[test]
    fn test_j2_plane_stress_satisfies_the_backward_euler_flow_rule() {
        // σ_trial − σ = Δγ C P σ with Δγ = 3 Δε̄_p / (2 σ_vm): the update is the closest-point
        // projection, not merely a point on the yield surface.
        let (young, nu) = (200.0e9, 0.3);
        let j2 = J2Plasticity::new(young, nu, 250.0e6, 1.0e9);
        let trial = [400.0e6, -150.0e6, 90.0e6];
        let (s, dp) = j2
            .try_return_mapping_plane_stress(trial[0], trial[1], trial[2], 0.001)
            .unwrap();
        let delta_gamma = 1.5 * dp / plane_stress_von_mises(s);
        let flow = plane_stress_c_p_sigma(young, nu, s);
        for ((trial_i, s_i), flow_i) in trial.iter().zip(s.iter()).zip(flow.iter()) {
            let expected = trial_i - delta_gamma * flow_i;
            assert!(
                (s_i - expected).abs() <= 1e-8 * trial_i.abs().max(1.0e6),
                "σ = {s_i} but σ_trial − Δγ C P σ = {expected}"
            );
        }
    }
    #[test]
    fn test_j2_plane_stress_perfect_plasticity_and_monotonic_in_trial() {
        let j2 = J2Plasticity::new(200.0e9, 0.3, 250.0e6, 0.0);
        let mut previous = 0.0;
        for &sxx in &[300.0e6, 400.0e6, 600.0e6, 1.0e9] {
            let (s, dp) = j2
                .try_return_mapping_plane_stress(sxx, 0.5 * sxx, 0.0, 0.0)
                .unwrap();
            assert!(
                (plane_stress_von_mises(s) - 250.0e6).abs() / 250.0e6 < 1e-9,
                "perfectly plastic von Mises must stay at the initial yield stress"
            );
            assert!(dp > previous, "Δε̄_p must grow with the overstress");
            previous = dp;
        }
    }
    #[test]
    fn test_j2_plane_stress_elastic_trial_is_unchanged() {
        let j2 = J2Plasticity::new(200.0e9, 0.3, 250.0e6, 1.0e9);
        let (s, dp) = j2
            .try_return_mapping_plane_stress(100.0e6, -50.0e6, 40.0e6, 0.0)
            .unwrap();
        assert_eq!(s, [100.0e6, -50.0e6, 40.0e6]);
        assert_eq!(dp, 0.0);
    }
    #[test]
    fn test_j2_plane_stress_invalid_input_is_reported() {
        use crate::Error;
        let bad_modulus = J2Plasticity::new(-1.0, 0.3, 250.0e6, 1.0e9);
        assert!(matches!(
            bad_modulus.try_return_mapping_plane_stress(400.0e6, 0.0, 0.0, 0.0),
            Err(Error::ParameterOutOfRange { .. })
        ));
        let bad_nu = J2Plasticity::new(200.0e9, 1.0, 250.0e6, 1.0e9);
        assert!(matches!(
            bad_nu.try_return_mapping_plane_stress(400.0e6, 0.0, 0.0, 0.0),
            Err(Error::ParameterOutOfRange { .. })
        ));
        let j2 = J2Plasticity::new(200.0e9, 0.3, 250.0e6, 1.0e9);
        assert!(matches!(
            j2.try_return_mapping_plane_stress(f64::NAN, 0.0, 0.0, 0.0),
            Err(Error::InadmissibleState(_))
        ));
        // The infallible wrapper reports invalid input as NaN instead of a bogus state.
        let (s, dp) = bad_modulus.return_mapping_plane_stress(400.0e6, 0.0, 0.0, 0.0);
        assert!(s.iter().all(|v| v.is_nan()) && dp.is_nan());
        // Valid input: wrapper and fallible variant agree.
        let (s_wrapped, dp_wrapped) = j2.return_mapping_plane_stress(400.0e6, 0.0, 0.0, 0.0);
        let (s_try, dp_try) = j2
            .try_return_mapping_plane_stress(400.0e6, 0.0, 0.0, 0.0)
            .unwrap();
        assert_eq!(s_wrapped, s_try);
        assert_eq!(dp_wrapped, dp_try);
    }
    #[test]
    fn test_ogden_uniaxial_stress_positive_tension() {
        let og = Ogden::one_term(1.0e6, 2.0, 1.0e9);
        let s = ogden_uniaxial_stress_incompressible(&og, 1.5);
        assert!(s > 0.0, "uniaxial stress should be positive: {s}");
    }
    #[test]
    fn test_ogden_biaxial_stress_positive_tension() {
        let og = Ogden::one_term(1.0e6, 2.0, 1.0e9);
        let s = ogden_biaxial_stress_incompressible(&og, 1.3);
        assert!(s > 0.0, "biaxial stress should be positive: {s}");
    }
    #[test]
    fn test_cohesive_peak_traction() {
        let cz = BilinearCohesiveZone::new(1.0e9, 1.0e-4, 5.0e-4);
        assert!(
            (cz.peak_traction() - 1.0e5).abs() < 1.0,
            "T_max={}",
            cz.peak_traction()
        );
    }
    #[test]
    fn test_cohesive_fracture_toughness() {
        let cz = BilinearCohesiveZone::new(1.0e9, 1.0e-4, 5.0e-4);
        let gc = cz.fracture_toughness();
        assert!((gc - 25.0).abs() < 0.1, "G_c={gc}");
    }
    #[test]
    fn test_cohesive_traction_at_zero() {
        let cz = BilinearCohesiveZone::new(1.0e9, 1.0e-4, 5.0e-4);
        assert!(cz.traction(0.0).abs() < 1e-10);
    }
    #[test]
    fn test_cohesive_traction_linear_region() {
        let cz = BilinearCohesiveZone::new(1.0e9, 1.0e-4, 5.0e-4);
        let t = cz.traction(5.0e-5);
        assert!((t - 5.0e4).abs() < 1.0, "traction in linear region: {t}");
    }
    #[test]
    fn test_cohesive_traction_at_failure() {
        let cz = BilinearCohesiveZone::new(1.0e9, 1.0e-4, 5.0e-4);
        let t = cz.traction(5.0e-4 + 1e-10);
        assert!(t == 0.0, "traction after failure: {t}");
    }
    #[test]
    fn test_cohesive_damage_zero_at_initiation() {
        let cz = BilinearCohesiveZone::new(1.0e9, 1.0e-4, 5.0e-4);
        assert!((cz.damage_variable(1.0e-4)).abs() < 1e-10);
    }
    #[test]
    fn test_cohesive_damage_one_at_failure() {
        let cz = BilinearCohesiveZone::new(1.0e9, 1.0e-4, 5.0e-4);
        assert!((cz.damage_variable(5.0e-4) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_fung_identity_zero_energy() {
        let f = Fung::aorta();
        let w = f.strain_energy_from_stretches(1.0, 1.0, 1.0);
        assert!(
            w.abs() < 1.0e-3,
            "Fung energy at identity should be ~0, got {w}"
        );
    }
    #[test]
    fn test_fung_energy_positive_for_tension() {
        let f = Fung::aorta();
        let lambda = 1.3_f64;
        let inv = 1.0 / lambda.sqrt();
        let w = f.strain_energy_from_stretches(lambda, inv, inv);
        assert!(w > 0.0, "Fung energy should be positive under tension: {w}");
    }
    #[test]
    fn test_fung_energy_increases_with_stretch() {
        let fung = Fung::myocardium();
        let inv1 = 1.0_f64 / 1.2_f64.sqrt();
        let inv2 = 1.0_f64 / 1.5_f64.sqrt();
        let w1 = fung.strain_energy_from_stretches(1.2, inv1, inv1);
        let w2 = fung.strain_energy_from_stretches(1.5, inv2, inv2);
        assert!(
            w2 > w1,
            "Fung energy should increase with stretch: {w1} < {w2}"
        );
    }
    #[test]
    fn test_fung_cauchy_stress_positive_for_tension() {
        let f = Fung::aorta();
        let sigma = f.uniaxial_cauchy_stress(1.3);
        assert!(
            sigma > 0.0,
            "Fung uniaxial Cauchy stress should be positive: {sigma}"
        );
    }
    #[test]
    fn test_fung_cauchy_stress_zero_at_identity() {
        let f = Fung::aorta();
        let sigma = f.uniaxial_cauchy_stress(1.0);
        assert!(
            sigma.abs() < 1.0,
            "Fung stress at λ=1 should be ~0: {sigma}"
        );
    }
    #[test]
    fn test_fung_small_strain_shear_modulus() {
        let f = Fung::new(100.0, 2.0, 0.0, 1.0e6);
        assert!((f.small_strain_shear_modulus() - 200.0).abs() < 1e-10);
    }
    #[test]
    fn test_fung_energy_from_gradient() {
        let f = Fung::myocardium();
        let grad = uniaxial_deformation_gradient(1.2);
        let w = f.strain_energy_density(&grad);
        assert!(w > 0.0, "Fung energy from gradient should be positive: {w}");
    }
    #[test]
    fn test_hgo_identity_near_zero_energy() {
        let hgo = HolzapfelGasserOgden::aorta_adventitia();
        let f_id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let w = hgo.strain_energy_density(&f_id);
        assert!(w.abs() < 1.0, "HGO at identity should be ~0, got {w}");
    }
    #[test]
    fn test_hgo_energy_positive_for_tension_along_fiber() {
        let hgo = HolzapfelGasserOgden::aorta_adventitia();
        let f_tens = uniaxial_deformation_gradient(1.3);
        let w = hgo.strain_energy_density(&f_tens);
        assert!(
            w > 0.0,
            "HGO energy should be positive under fiber tension: {w}"
        );
    }
    #[test]
    fn test_hgo_energy_increases_with_stretch() {
        let hgo = HolzapfelGasserOgden::aorta_media();
        let f1 = uniaxial_deformation_gradient(1.1);
        let f2 = uniaxial_deformation_gradient(1.3);
        let w1 = hgo.strain_energy_density(&f1);
        let w2 = hgo.strain_energy_density(&f2);
        assert!(
            w2 > w1,
            "HGO energy should increase with stretch: {w1} < {w2}"
        );
    }
    #[test]
    fn test_hgo_fibers_active_under_tension() {
        let hgo = HolzapfelGasserOgden::aorta_adventitia();
        let f_tens = uniaxial_deformation_gradient(1.3);
        let (active1, _active2) = hgo.fibers_active(&f_tens);
        assert!(
            active1,
            "Fiber 1 should be active under X-direction tension"
        );
    }
    #[test]
    fn test_hgo_volumetric_penalty_at_compression() {
        let hgo = HolzapfelGasserOgden::aorta_adventitia();
        let lam = 0.9_f64;
        let f_comp = [[lam, 0.0, 0.0], [0.0, lam, 0.0], [0.0, 0.0, lam]];
        let w = hgo.strain_energy_density(&f_comp);
        assert!(
            w > 0.0,
            "HGO volumetric penalty should give positive energy under compression"
        );
    }
    #[test]
    fn test_tangent_bulk_modulus_positive() {
        let mr = MooneyRivlin::new(0.5e6, 0.1e6, 1.0e9);
        let k = tangent_bulk_modulus_fd(
            |lam| {
                let f = biaxial_deformation_gradient(lam);
                mr.strain_energy_density(&f)
            },
            1.0,
        );
        assert!(k.is_finite(), "bulk modulus should be finite: {k}");
    }
    #[test]
    fn test_tangent_shear_modulus_fd_neohookean() {
        let mu = 1.0e6_f64;
        let tau_at_0 = tangent_shear_modulus_fd(
            |f| {
                let (i1, _i2, _i3) = deformation_invariants(f);
                mu / 2.0 * (i1 - 3.0)
            },
            0.0,
        );
        assert!(
            tau_at_0.abs() < 1.0,
            "τ at γ=0 for Neo-Hookean should be ~0: {tau_at_0}"
        );
    }
    #[test]
    fn test_principal_cauchy_stress_incompressible_uniaxial() {
        let mu = 1.0e6_f64;
        let lambda = 1.5_f64;
        let inv = 1.0 / lambda.sqrt();
        let p = lagrange_multiplier_uniaxial(
            &|l1, l2, l3| {
                let i1 = l1 * l1 + l2 * l2 + l3 * l3;
                mu / 2.0 * (i1 - 3.0)
            },
            lambda,
        );
        let sigma = principal_cauchy_stress_incompressible(
            &|l1, l2, l3| {
                let i1 = l1 * l1 + l2 * l2 + l3 * l3;
                mu / 2.0 * (i1 - 3.0)
            },
            lambda,
            inv,
            inv,
            p,
        );
        assert!(sigma[1].abs() < 1.0, "σ₂ should be ~0: {}", sigma[1]);
        assert!(sigma[2].abs() < 1.0, "σ₃ should be ~0: {}", sigma[2]);
        let sigma_exact = neo_hookean_uniaxial_stress(mu, lambda);
        assert!(
            (sigma[0] - sigma_exact).abs() / sigma_exact.abs() < 0.01,
            "σ₁ = {} vs exact {}",
            sigma[0],
            sigma_exact
        );
    }
    /// At identity F, tangent equals isotropic elastic stiffness.
    #[test]
    fn test_neo_hookean_tangent_identity_deformation() {
        let mu = 1.0e6_f64;
        let k = 3.0e6_f64;
        let f_id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let c = neo_hookean_constitutive_tangent(mu, k, &f_id);
        let c11_expected = k + 4.0 / 3.0 * mu;
        assert!(
            (c[0] - c11_expected).abs() / c11_expected < 1e-8,
            "C11 at identity: {} vs {}",
            c[0],
            c11_expected
        );
    }
    /// Tangent is symmetric for isotropic Neo-Hookean.
    #[test]
    fn test_neo_hookean_tangent_symmetric() {
        let f = [[1.1, 0.0, 0.0], [0.0, 0.9, 0.0], [0.0, 0.0, 1.0 / 0.99]];
        let c = neo_hookean_constitutive_tangent(1.0e6, 3.0e6, &f);
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (c[i * 6 + j] - c[j * 6 + i]).abs() < 1e-6,
                    "Tangent not symmetric at [{i},{j}]: {} vs {}",
                    c[i * 6 + j],
                    c[j * 6 + i]
                );
            }
        }
    }
    /// Shear tangent components (C44, C55, C66) are positive.
    #[test]
    fn test_neo_hookean_tangent_shear_positive() {
        let f_id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let c = neo_hookean_constitutive_tangent(1.0e6, 3.0e6, &f_id);
        assert!(
            c[3 * 6 + 3] > 0.0,
            "C44 should be positive: {}",
            c[3 * 6 + 3]
        );
        assert!(
            c[4 * 6 + 4] > 0.0,
            "C55 should be positive: {}",
            c[4 * 6 + 4]
        );
        assert!(
            c[5 * 6 + 5] > 0.0,
            "C66 should be positive: {}",
            c[5 * 6 + 5]
        );
    }
    /// Hessian at identity state is positive definite (diagonal > 0).
    #[test]
    fn test_mooney_rivlin_hessian_diagonal_positive() {
        let mr = MooneyRivlin::new(0.5e6, 0.1e6, 2.0e9);
        let c_voigt = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0_f64];
        let hess = mooney_rivlin_strain_energy_hessian(&mr, &c_voigt, 1e-5);
        for i in 0..3 {
            assert!(
                hess[i * 6 + i] > 0.0,
                "Hessian diagonal[{i}] should be positive: {}",
                hess[i * 6 + i]
            );
        }
    }
    /// Hessian is approximately symmetric.
    #[test]
    fn test_mooney_rivlin_hessian_symmetric() {
        let mr = MooneyRivlin::new(0.5e6, 0.1e6, 2.0e9);
        let c_voigt = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0_f64];
        let hess = mooney_rivlin_strain_energy_hessian(&mr, &c_voigt, 1e-5);
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (hess[i * 6 + j] - hess[j * 6 + i]).abs() < 1.0,
                    "Hessian not symmetric at [{i},{j}]: {} vs {}",
                    hess[i * 6 + j],
                    hess[j * 6 + i]
                );
            }
        }
    }
    /// At identity (no deformation) every principal stress vanishes: the Ogden stress
    /// carries the deviatoric projection of the isochoric response.
    #[test]
    fn test_ogden_principal_stresses_identity() {
        let mu = 1.0e6;
        let ogden = Ogden::one_term(mu, 2.0, 1.0e9);
        let sigma = ogden_principal_stresses(&ogden, 1.0, 1.0, 1.0);
        for (i, s) in sigma.iter().enumerate() {
            assert!(
                s.abs() < 1.0e-6 * mu,
                "σ_{} = {s} should vanish at the identity",
                i + 1
            );
        }
        let rubber = Ogden::rubber_default();
        for s in ogden_principal_stresses(&rubber, 1.0, 1.0, 1.0) {
            assert!(s.abs() < 1.0e-6 * mu, "rubber σ = {s} at the identity");
        }
    }
    /// Uniaxial stress state (λ, λ^-1/2, λ^-1/2), J = 1:
    /// σ₁ − σ₂ = Σ μ_p (λ^α_p − λ^(−α_p/2)).
    #[test]
    fn test_ogden_principal_stresses_uniaxial_matches_incompressible_formula() {
        let ogden = Ogden::rubber_default();
        for &lambda in &[0.7_f64, 1.1, 1.4, 2.0, 3.0] {
            let lateral = 1.0 / lambda.sqrt();
            let sigma = ogden_principal_stresses(&ogden, lambda, lateral, lateral);
            let expected = ogden_uniaxial_stress_incompressible(&ogden, lambda);
            assert!(expected != 0.0);
            assert!(
                ((sigma[0] - sigma[1]) - expected).abs() <= 1e-9 * expected.abs(),
                "λ = {lambda}: σ₁ − σ₂ = {}, expected {expected}",
                sigma[0] - sigma[1]
            );
            assert!(
                (sigma[1] - sigma[2]).abs() <= 1e-9 * expected.abs(),
                "λ = {lambda}: the two lateral stresses must coincide"
            );
        }
    }
    /// The Neo-Hookean special case (α = 2) reproduces μ (λ² − 1/λ).
    #[test]
    fn test_ogden_principal_stresses_neo_hookean_limit() {
        let mu = 0.8e6;
        let ogden = Ogden::one_term(mu, 2.0, 1.0e10);
        for &lambda in &[1.2_f64, 1.8] {
            let lateral = 1.0 / lambda.sqrt();
            let sigma = ogden_principal_stresses(&ogden, lambda, lateral, lateral);
            let expected = neo_hookean_uniaxial_stress(mu, lambda);
            assert!(((sigma[0] - sigma[1]) - expected).abs() <= 1e-9 * expected.abs());
        }
    }
    /// Σ σ_i = 3 K (J − 1): the isochoric Ogden stress is deviatoric.
    #[test]
    fn test_ogden_principal_stresses_isochoric_part_is_traceless() {
        let k = 2.0e6;
        let ogden = Ogden::new(vec![6.3e5, 1.2e3, -1.0e4], vec![1.3, 5.0, -2.0], k);
        let (l1, l2, l3) = (1.3, 0.9, 1.1);
        let volumetric = 3.0 * k * (l1 * l2 * l3 - 1.0);
        let sigma = ogden_principal_stresses(&ogden, l1, l2, l3);
        let trace: f64 = sigma.iter().sum();
        assert!(
            (trace - volumetric).abs() <= 1e-9 * volumetric.abs(),
            "tr σ = {trace}, 3K(J − 1) = {volumetric}"
        );
    }
    /// σ_i = λ_i (∂W/∂λ_i) / J for the compressible Ogden energy.
    #[test]
    fn test_ogden_principal_stresses_match_energy_derivative() {
        let ogden = Ogden::new(vec![6.3e5, 1.2e3, -1.0e4], vec![1.3, 5.0, -2.0], 2.0e6);
        let lambda = [1.3_f64, 0.9, 1.1];
        let j = lambda[0] * lambda[1] * lambda[2];
        let sigma = ogden_principal_stresses(&ogden, lambda[0], lambda[1], lambda[2]);
        for (i, &sigma_i) in sigma.iter().enumerate() {
            let h = 1.0e-6 * lambda[i];
            let mut plus = lambda;
            plus[i] += h;
            let mut minus = lambda;
            minus[i] -= h;
            let dw = (ogden.strain_energy_from_stretches(plus[0], plus[1], plus[2])
                - ogden.strain_energy_from_stretches(minus[0], minus[1], minus[2]))
                / (2.0 * h);
            let expected = lambda[i] * dw / j;
            assert!(
                (sigma_i - expected).abs() <= 1e-6 * expected.abs() + 1.0e-3,
                "σ_{} = {sigma_i}, λ ∂W/∂λ / J = {expected}",
                i + 1
            );
        }
    }
    /// Under uniaxial tension (λ₁ > 1), σ₁ > σ₂ for Ogden.
    #[test]
    fn test_ogden_principal_stresses_uniaxial_tension() {
        let ogden = Ogden::one_term(1.0e6, 2.0, 1.0e10);
        let lambda = 1.5_f64;
        let inv = 1.0 / lambda.sqrt();
        let sigma = ogden_principal_stresses(&ogden, lambda, inv, inv);
        assert!(
            sigma[0] > sigma[1],
            "σ₁ ({}) should exceed σ₂ ({}) under uniaxial tension",
            sigma[0],
            sigma[1]
        );
    }
    /// Principal stresses increase with stretch for Ogden.
    #[test]
    fn test_ogden_principal_stresses_increase_with_stretch() {
        let ogden = Ogden::one_term(1.0e6, 2.0, 1.0e10);
        let inv1 = 1.0 / 1.2_f64.sqrt();
        let inv2 = 1.0 / 1.5_f64.sqrt();
        let s1 = ogden_principal_stresses(&ogden, 1.2, inv1, inv1);
        let s2 = ogden_principal_stresses(&ogden, 1.5, inv2, inv2);
        assert!(
            s2[0] > s1[0],
            "σ₁ should increase with stretch: {} at 1.2 vs {} at 1.5",
            s1[0],
            s2[0]
        );
    }
    /// Rubber default: principal stresses are finite and well-ordered under biaxial.
    #[test]
    fn test_ogden_rubber_default_biaxial() {
        let ogden = Ogden::rubber_default();
        let lambda = 1.3_f64;
        let lambda3 = 1.0 / (lambda * lambda);
        let sigma = ogden_principal_stresses(&ogden, lambda, lambda, lambda3);
        assert!(sigma[0].is_finite(), "σ₁ should be finite: {}", sigma[0]);
        assert!(sigma[1].is_finite(), "σ₂ should be finite: {}", sigma[1]);
        assert!(
            (sigma[0] - sigma[1]).abs() < 1.0,
            "σ₁ and σ₂ should be equal for equibiaxial: {} vs {}",
            sigma[0],
            sigma[1]
        );
    }
}
