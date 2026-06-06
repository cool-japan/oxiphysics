//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{BucklingProblem, SlendernessClass};

/// Euler column buckling critical load.
///
/// P_cr = π² E I / (k L)²
///
/// * `E` – Young's modulus
/// * `I` – second moment of area
/// * `L` – column length
/// * `k_factor` – effective length factor (1.0 = pinned-pinned, 0.5 = fixed-fixed, etc.)
pub fn column_buckling_load(e: f64, i: f64, l: f64, k_factor: f64) -> f64 {
    PI * PI * e * i / (k_factor * l).powi(2)
}
/// Plate buckling coefficient for simply supported rectangular plate.
///
/// k = (m b/a + n² a/(m b))²
///
/// * `a` – plate length (loading direction)
/// * `b` – plate width
/// * `m` – half-waves in loading direction
/// * `n` – half-waves in transverse direction
pub fn plate_buckling_coefficient(a: f64, b: f64, m: u32, n: u32) -> f64 {
    let mf = m as f64;
    let nf = n as f64;
    let term = mf * b / a + nf * nf * a / (mf * b);
    term * term
}
/// Critical buckling stress for a plate.
///
/// σ_cr = k π² E / (12 (1 − ν²)) (t/b)²
///
/// * `e`  – Young's modulus
/// * `nu` – Poisson's ratio
/// * `t`  – plate thickness
/// * `b`  – plate width (shorter dimension)
/// * `k`  – plate buckling coefficient
pub fn critical_stress(e: f64, nu: f64, t: f64, b: f64, k: f64) -> f64 {
    k * PI * PI * e / (12.0 * (1.0 - nu * nu)) * (t / b).powi(2)
}
/// Koiter imperfection sensitivity.
///
/// λ_imperfect ≈ λ_perfect × (1 − c √|δ|)  where c = 0.5
///
/// * `load_factor`            – perfect-structure buckling load factor
/// * `imperfection_amplitude` – normalized imperfection amplitude δ
pub fn imperfection_sensitivity(load_factor: f64, imperfection_amplitude: f64) -> f64 {
    pub(super) const C: f64 = 0.5;
    load_factor * (1.0 - C * imperfection_amplitude.abs().sqrt())
}
/// Assemble beam K and Kg for a pinned-pinned Euler-Bernoulli column and return
/// the lowest buckling load multiplier.
///
/// * `n_elements` – number of beam elements
/// * `E`          – Young's modulus
/// * `I`          – second moment of area
/// * `L`          – total column length
/// * `P`          – reference axial load (Kg is assembled for this load)
pub fn column_fem_buckling(n_elements: usize, e: f64, i_mom: f64, l: f64, p: f64) -> f64 {
    let le = l / n_elements as f64;
    let n_nodes = n_elements + 1;
    let total_dof = 2 * n_nodes;
    let ei = e * i_mom;
    let mut k_full = vec![0.0f64; total_dof * total_dof];
    let mut kg_full = vec![0.0f64; total_dof * total_dof];
    for elem in 0..n_elements {
        let n0 = elem;
        let n1 = elem + 1;
        let dofs = [2 * n0, 2 * n0 + 1, 2 * n1, 2 * n1 + 1];
        let c = ei / (le * le * le);
        let ke: [[f64; 4]; 4] = [
            [12.0 * c, 6.0 * c * le, -12.0 * c, 6.0 * c * le],
            [
                6.0 * c * le,
                4.0 * c * le * le,
                -6.0 * c * le,
                2.0 * c * le * le,
            ],
            [-12.0 * c, -6.0 * c * le, 12.0 * c, -6.0 * c * le],
            [
                6.0 * c * le,
                2.0 * c * le * le,
                -6.0 * c * le,
                4.0 * c * le * le,
            ],
        ];
        let cg = p / (30.0 * le);
        let kge: [[f64; 4]; 4] = [
            [36.0 * cg, 3.0 * cg * le, -36.0 * cg, 3.0 * cg * le],
            [
                3.0 * cg * le,
                4.0 * cg * le * le,
                -3.0 * cg * le,
                -cg * le * le,
            ],
            [-36.0 * cg, -3.0 * cg * le, 36.0 * cg, -3.0 * cg * le],
            [
                3.0 * cg * le,
                -cg * le * le,
                -3.0 * cg * le,
                4.0 * cg * le * le,
            ],
        ];
        for a in 0..4 {
            for b in 0..4 {
                k_full[dofs[a] * total_dof + dofs[b]] += ke[a][b];
                kg_full[dofs[a] * total_dof + dofs[b]] += kge[a][b];
            }
        }
    }
    let fixed: Vec<usize> = vec![0, 2 * n_elements];
    let free: Vec<usize> = (0..total_dof).filter(|d| !fixed.contains(d)).collect();
    let nf = free.len();
    let mut k_red = vec![0.0f64; nf * nf];
    let mut kg_red = vec![0.0f64; nf * nf];
    for (ia, &ga) in free.iter().enumerate() {
        for (ib, &gb) in free.iter().enumerate() {
            k_red[ia * nf + ib] = k_full[ga * total_dof + gb];
            kg_red[ia * nf + ib] = kg_full[ga * total_dof + gb];
        }
    }
    let mut prob = BucklingProblem::new(nf);
    for ia in 0..nf {
        for ib in 0..nf {
            prob.set_k(ia, ib, k_red[ia * nf + ib]);
            prob.set_kg(ia, ib, kg_red[ia * nf + ib]);
        }
    }
    let factors = prob.solve_buckling_load_factors(1);
    factors.into_iter().next().unwrap_or(f64::NAN)
}
pub(super) fn gaussian_solve(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut mat = vec![0.0f64; n * (n + 1)];
    for i in 0..n {
        for j in 0..n {
            mat[i * (n + 1) + j] = a[i * n + j];
        }
        mat[i * (n + 1) + n] = b[i];
    }
    for col in 0..n {
        let mut max_row = col;
        let mut max_val = mat[col * (n + 1) + col].abs();
        for row in (col + 1)..n {
            let v = mat[row * (n + 1) + col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_val < 1e-300 {
            return None;
        }
        if max_row != col {
            for k in 0..=(n) {
                mat.swap(col * (n + 1) + k, max_row * (n + 1) + k);
            }
        }
        let pivot = mat[col * (n + 1) + col];
        for row in (col + 1)..n {
            let factor = mat[row * (n + 1) + col] / pivot;
            for k in col..=(n) {
                let val = mat[col * (n + 1) + k] * factor;
                mat[row * (n + 1) + k] -= val;
            }
        }
    }
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let mut s = mat[i * (n + 1) + n];
        for j in (i + 1)..n {
            s -= mat[i * (n + 1) + j] * x[j];
        }
        x[i] = s / mat[i * (n + 1) + i];
    }
    Some(x)
}
/// Assemble a simple truss geometric stiffness matrix for `n_elements`
/// bar elements of length `le` carrying axial load `p`.
///
/// Returns an (n_nodes x n_nodes) dense matrix stored row-major.
pub fn truss_geometric_stiffness(n_elements: usize, le: f64, p: f64) -> Vec<f64> {
    let n_nodes = n_elements + 1;
    let n = n_nodes;
    let mut kg = vec![0.0f64; n * n];
    let cg = p / le;
    for e in 0..n_elements {
        let i = e;
        let j = e + 1;
        kg[i * n + i] += cg;
        kg[i * n + j] -= cg;
        kg[j * n + i] -= cg;
        kg[j * n + j] += cg;
    }
    kg
}
/// Assemble a simple truss elastic stiffness matrix for `n_elements`
/// bar elements of length `le` with axial stiffness `ea`.
pub fn truss_elastic_stiffness(n_elements: usize, le: f64, ea: f64) -> Vec<f64> {
    let n_nodes = n_elements + 1;
    let n = n_nodes;
    let mut k = vec![0.0f64; n * n];
    let c = ea / le;
    for e in 0..n_elements {
        let i = e;
        let j = e + 1;
        k[i * n + i] += c;
        k[i * n + j] -= c;
        k[j * n + i] -= c;
        k[j * n + j] += c;
    }
    k
}
/// Single step of the arc-length (Riks) method for 1D post-buckling analysis.
///
/// Given a load-displacement curve approximated as a polynomial, trace the
/// equilibrium path beyond the critical point.
///
/// * `k_tangent` - tangent stiffness (scalar for 1D)
/// * `load` - current load level
/// * `disp` - current displacement
/// * `ds` - arc-length increment
///
/// Returns `(new_disp, new_load)`.
pub fn arc_length_step_1d(k_tangent: f64, load: f64, disp: f64, ds: f64) -> (f64, f64) {
    if k_tangent.abs() < 1e-15 {
        return (disp + ds, load);
    }
    let d_disp = ds / (1.0 + 1.0 / (k_tangent * k_tangent)).sqrt();
    let d_load = d_disp / k_tangent;
    (disp + d_disp, load + d_load)
}
/// Trace a post-buckling path using the arc-length method.
///
/// * `k_fn` - function that returns tangent stiffness given displacement
/// * `initial_load` - starting load
/// * `initial_disp` - starting displacement
/// * `ds` - arc-length step size
/// * `n_steps` - number of steps to trace
///
/// Returns a vector of `(displacement, load)` pairs.
pub fn trace_post_buckling_path(
    k_fn: impl Fn(f64) -> f64,
    initial_load: f64,
    initial_disp: f64,
    ds: f64,
    n_steps: usize,
) -> Vec<(f64, f64)> {
    let mut path = Vec::with_capacity(n_steps + 1);
    let mut disp = initial_disp;
    let mut load = initial_load;
    path.push((disp, load));
    for _ in 0..n_steps {
        let kt = k_fn(disp);
        let (d, l) = arc_length_step_1d(kt, load, disp, ds);
        disp = d;
        load = l;
        path.push((disp, load));
    }
    path
}
/// Compute imperfect buckling loads for a range of imperfection amplitudes.
///
/// Uses the Koiter formula: lambda_imp = lambda_perf * (1 - c * sqrt(|delta|))
///
/// Returns a vector of `(amplitude, reduced_load_factor)` pairs.
pub fn imperfection_sensitivity_curve(
    perfect_load_factor: f64,
    amplitudes: &[f64],
    c: f64,
) -> Vec<(f64, f64)> {
    amplitudes
        .iter()
        .map(|&delta| {
            let reduced = perfect_load_factor * (1.0 - c * delta.abs().sqrt());
            (delta, reduced.max(0.0))
        })
        .collect()
}
/// Multi-mode imperfection sensitivity (combined imperfection of first two modes).
///
/// lambda_imp = lambda_perf * (1 - c1*sqrt(|d1|) - c2*sqrt(|d2|))
pub fn multimode_imperfection_sensitivity(
    perfect_load_factor: f64,
    d1: f64,
    c1: f64,
    d2: f64,
    c2: f64,
) -> f64 {
    let reduced = perfect_load_factor * (1.0 - c1 * d1.abs().sqrt() - c2 * d2.abs().sqrt());
    reduced.max(0.0)
}
/// Detect snap-through by monitoring the tangent stiffness along a load path.
///
/// Returns the index in the path where the tangent stiffness first becomes
/// negative (indicating snap-through), or `None` if no snap-through occurs.
///
/// * `stiffness_values` - tangent stiffness at each load step
pub fn detect_snap_through(stiffness_values: &[f64]) -> Option<usize> {
    for (i, &k) in stiffness_values.iter().enumerate() {
        if k < 0.0 {
            return Some(i);
        }
    }
    None
}
/// Analyze load-displacement data for snap-through behavior.
///
/// Returns `Some((snap_load, snap_disp))` at the point where the load
/// starts decreasing, or `None` if the load is monotonically increasing.
pub fn snap_through_from_path(path: &[(f64, f64)]) -> Option<(f64, f64)> {
    if path.len() < 2 {
        return None;
    }
    for i in 1..path.len() {
        let (d_prev, l_prev) = path[i - 1];
        let (_d_curr, l_curr) = path[i];
        if l_curr < l_prev {
            return Some((l_prev, d_prev));
        }
    }
    None
}
/// Compute the limit point (maximum load) from a load-displacement path.
pub fn find_limit_point(path: &[(f64, f64)]) -> Option<(f64, f64)> {
    path.iter()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .copied()
}
/// Assemble the full geometric stiffness matrix for a cantilever beam
/// under tip load.
///
/// Uses Euler-Bernoulli beam theory with consistent geometric stiffness.
/// Returns a dense `(2*n_nodes) x (2*n_nodes)` matrix.
pub fn beam_geometric_stiffness(n_elements: usize, le: f64, p: f64) -> Vec<f64> {
    let n_nodes = n_elements + 1;
    let total_dof = 2 * n_nodes;
    let mut kg = vec![0.0f64; total_dof * total_dof];
    let cg = p / (30.0 * le);
    for elem in 0..n_elements {
        let n0 = elem;
        let n1 = elem + 1;
        let dofs = [2 * n0, 2 * n0 + 1, 2 * n1, 2 * n1 + 1];
        let kge: [[f64; 4]; 4] = [
            [36.0 * cg, 3.0 * cg * le, -36.0 * cg, 3.0 * cg * le],
            [
                3.0 * cg * le,
                4.0 * cg * le * le,
                -3.0 * cg * le,
                -cg * le * le,
            ],
            [-36.0 * cg, -3.0 * cg * le, 36.0 * cg, -3.0 * cg * le],
            [
                3.0 * cg * le,
                -cg * le * le,
                -3.0 * cg * le,
                4.0 * cg * le * le,
            ],
        ];
        for a in 0..4 {
            for b in 0..4 {
                kg[dofs[a] * total_dof + dofs[b]] += kge[a][b];
            }
        }
    }
    kg
}
/// Compute the critical load for a fixed-free (cantilever) column.
///
/// Uses FEM with Euler-Bernoulli elements and fixed-free boundary conditions.
/// The analytical solution is P_cr = pi^2*EI/(4*L^2).
pub fn cantilever_fem_buckling(n_elements: usize, e: f64, i_mom: f64, l: f64, p: f64) -> f64 {
    let le = l / n_elements as f64;
    let n_nodes = n_elements + 1;
    let total_dof = 2 * n_nodes;
    let ei = e * i_mom;
    let mut k_full = vec![0.0f64; total_dof * total_dof];
    let mut kg_full = vec![0.0f64; total_dof * total_dof];
    for elem in 0..n_elements {
        let n0 = elem;
        let n1 = elem + 1;
        let dofs = [2 * n0, 2 * n0 + 1, 2 * n1, 2 * n1 + 1];
        let c = ei / (le * le * le);
        let ke: [[f64; 4]; 4] = [
            [12.0 * c, 6.0 * c * le, -12.0 * c, 6.0 * c * le],
            [
                6.0 * c * le,
                4.0 * c * le * le,
                -6.0 * c * le,
                2.0 * c * le * le,
            ],
            [-12.0 * c, -6.0 * c * le, 12.0 * c, -6.0 * c * le],
            [
                6.0 * c * le,
                2.0 * c * le * le,
                -6.0 * c * le,
                4.0 * c * le * le,
            ],
        ];
        let cg = p / (30.0 * le);
        let kge: [[f64; 4]; 4] = [
            [36.0 * cg, 3.0 * cg * le, -36.0 * cg, 3.0 * cg * le],
            [
                3.0 * cg * le,
                4.0 * cg * le * le,
                -3.0 * cg * le,
                -cg * le * le,
            ],
            [-36.0 * cg, -3.0 * cg * le, 36.0 * cg, -3.0 * cg * le],
            [
                3.0 * cg * le,
                -cg * le * le,
                -3.0 * cg * le,
                4.0 * cg * le * le,
            ],
        ];
        for a in 0..4 {
            for b in 0..4 {
                k_full[dofs[a] * total_dof + dofs[b]] += ke[a][b];
                kg_full[dofs[a] * total_dof + dofs[b]] += kge[a][b];
            }
        }
    }
    let fixed: Vec<usize> = vec![0, 1];
    let free: Vec<usize> = (0..total_dof).filter(|d| !fixed.contains(d)).collect();
    let nf = free.len();
    let mut k_red = vec![0.0f64; nf * nf];
    let mut kg_red = vec![0.0f64; nf * nf];
    for (ia, &ga) in free.iter().enumerate() {
        for (ib, &gb) in free.iter().enumerate() {
            k_red[ia * nf + ib] = k_full[ga * total_dof + gb];
            kg_red[ia * nf + ib] = kg_full[ga * total_dof + gb];
        }
    }
    let mut prob = BucklingProblem::new(nf);
    for ia in 0..nf {
        for ib in 0..nf {
            prob.set_k(ia, ib, k_red[ia * nf + ib]);
            prob.set_kg(ia, ib, kg_red[ia * nf + ib]);
        }
    }
    let factors = prob.solve_buckling_load_factors(1);
    factors.into_iter().next().unwrap_or(f64::NAN)
}
/// Consistent geometric stiffness matrix for an Euler-Bernoulli beam element.
///
/// DOF order: \[w_i, θ_i, w_j, θ_j\]
///
/// * `N` – axial compressive force (positive = compression)
/// * `L` – element length
///
/// Returns a 4×4 matrix (row-major).
pub fn geometric_stiffness_beam(n: f64, l: f64) -> [[f64; 4]; 4] {
    let c = n / (30.0 * l);
    [
        [36.0 * c, 3.0 * c * l, -36.0 * c, 3.0 * c * l],
        [3.0 * c * l, 4.0 * c * l * l, -3.0 * c * l, -c * l * l],
        [-36.0 * c, -3.0 * c * l, 36.0 * c, -3.0 * c * l],
        [3.0 * c * l, -c * l * l, -3.0 * c * l, 4.0 * c * l * l],
    ]
}
/// Johnson parabolic formula for short-column buckling.
///
/// σ_cr = σ_y − (σ_y² / (4 π² E)) * λ²
///
/// Valid when λ < λ_c = π * sqrt(2 E / σ_y).
///
/// * `sigma_y`    – yield stress
/// * `e`          – Young's modulus
/// * `slenderness` – slenderness ratio λ = (K L) / r
pub fn johnson_parabolic_formula(sigma_y: f64, e: f64, slenderness: f64) -> f64 {
    (sigma_y - (sigma_y * sigma_y / (4.0 * PI * PI * e)) * slenderness * slenderness).max(0.0)
}
/// Engesser tangent-modulus theory critical stress.
///
/// σ_cr = π² E_t / λ²
///
/// * `_sigma`     – current stress (kept for API; used via E_t externally)
/// * `_e`         – elastic modulus (kept for API completeness)
/// * `e_t`        – tangent modulus at the current stress
/// * `slenderness` – slenderness ratio λ
pub fn tangent_modulus_theory(_sigma: f64, _e: f64, e_t: f64, slenderness: f64) -> f64 {
    if slenderness.abs() < 1e-15 {
        return f64::INFINITY;
    }
    PI * PI * e_t / (slenderness * slenderness)
}
/// Timoshenko beam buckling – accounts for transverse shear deformation.
///
/// The critical load for a simply supported column:
///
/// P_cr_T = P_E / (1 + P_E / (κ A G))
///
/// where P_E = π² E I / L² (Euler load) and κ A G is the shear rigidity.
///
/// * `e`  – Young's modulus
/// * `i`  – second moment of area
/// * `l`  – effective length
/// * `a`  – cross-sectional area
/// * `g`  – shear modulus
/// * `kappa` – shear correction factor (≈ 5/6 for rectangular cross-section)
pub fn timoshenko_critical_load(e: f64, i: f64, l: f64, a: f64, g: f64, kappa: f64) -> f64 {
    let p_euler = PI * PI * e * i / (l * l);
    let shear_stiffness = kappa * a * g;
    if shear_stiffness.abs() < 1e-30 {
        return p_euler;
    }
    p_euler / (1.0 + p_euler / shear_stiffness)
}
/// Ratio of Timoshenko to Euler critical load.
///
/// Returns P_cr_T / P_E = 1 / (1 + P_E / (κ A G)).
/// Indicates how much shear deformation reduces the critical load.
pub fn timoshenko_to_euler_ratio(e: f64, i: f64, l: f64, a: f64, g: f64, kappa: f64) -> f64 {
    let p_euler = PI * PI * e * i / (l * l);
    let shear_stiffness = kappa * a * g;
    if shear_stiffness.abs() < 1e-30 {
        return 1.0;
    }
    1.0 / (1.0 + p_euler / shear_stiffness)
}
/// Timoshenko beam buckling for a cantilever (fixed-free) column.
///
/// Uses effective length L_eff = 2L for the Euler term.
pub fn timoshenko_cantilever_buckling(e: f64, i: f64, l: f64, a: f64, g: f64, kappa: f64) -> f64 {
    timoshenko_critical_load(e, i, 2.0 * l, a, g, kappa)
}
/// Timoshenko beam element stiffness matrix (4×4).
///
/// DOF: \[w_i, θ_i, w_j, θ_j\]
///
/// Includes shear deformation for thick beams.
///
/// * `e`     – Young's modulus
/// * `i`     – second moment of area
/// * `g`     – shear modulus
/// * `a`     – cross-sectional area
/// * `kappa` – shear correction factor
/// * `l`     – element length
pub fn timoshenko_beam_stiffness(
    e: f64,
    i: f64,
    g: f64,
    a: f64,
    kappa: f64,
    l: f64,
) -> [[f64; 4]; 4] {
    let ei = e * i;
    let kag = kappa * a * g;
    let phi = if kag.abs() > 1e-30 {
        12.0 * ei / (kag * l * l)
    } else {
        0.0
    };
    let d = 1.0 + phi;
    let k11 = 12.0 * ei / (l * l * l * d);
    let k12 = 6.0 * ei / (l * l * d);
    let k22 = (4.0 + phi) * ei / (l * d);
    let k24 = (2.0 - phi) * ei / (l * d);
    [
        [k11, k12, -k11, k12],
        [k12, k22, -k12, k24],
        [-k11, -k12, k11, -k12],
        [k12, k24, -k12, k22],
    ]
}
/// Slenderness ratio below which Timoshenko theory significantly differs from Euler.
///
/// Returns the slenderness λ = L/r above which P_T/P_E > 0.99 (shear effect < 1%).
pub fn timoshenko_slenderness_threshold(e: f64, g: f64, kappa: f64, r_gyration: f64) -> f64 {
    let ratio = 0.01 / 0.99;
    if kappa * g * ratio < 1e-30 {
        return 0.0;
    }
    (PI * PI * e / (kappa * g * ratio)).sqrt() * r_gyration
}
/// Critical temperature rise for thermal buckling of a simply-supported plate.
///
/// ΔT_cr = k · π² / (α · b²) · (t/b)²  × D/(E h α)
///
/// Simplified: ΔT_cr = k · π² · t² / (12(1−ν²) α L²)
///
/// * `alpha` – thermal expansion coefficient (1/K)
/// * `nu`    – Poisson's ratio
/// * `t`     – plate thickness
/// * `l`     – plate dimension (shorter side)
/// * `k_plate` – plate buckling coefficient (use k=4 for simply-supported)
pub fn thermal_buckling_temperature(alpha: f64, nu: f64, t: f64, l: f64, k_plate: f64) -> f64 {
    k_plate * PI * PI * t * t / (12.0 * (1.0 - nu * nu) * alpha * l * l)
}
/// Thermal buckling load for a column with fixed ends.
///
/// P_cr_T = α · E · A · ΔT / (1 + α · E · A / P_cr_Euler)
///
/// Simplified: since the thermal force N_T = E A α ΔT acts as a compressive load,
/// buckling occurs when N_T ≥ P_cr_Euler.
///
/// Returns the critical temperature rise ΔT_cr.
pub fn column_thermal_buckling_delta_t(
    e: f64,
    i: f64,
    l: f64,
    a: f64,
    alpha: f64,
    k_factor: f64,
) -> f64 {
    let p_cr = PI * PI * e * i / (k_factor * l).powi(2);
    p_cr / (e * a * alpha)
}
/// Critical temperature rise for constrained thermal expansion of a beam.
///
/// When both ends are fully fixed and the beam cannot expand, the thermal
/// compressive force is N_T = E·A·α·ΔT, and buckling occurs when N_T ≥ P_cr.
pub fn constrained_beam_thermal_buckling(e: f64, i: f64, l: f64, alpha: f64) -> f64 {
    column_thermal_buckling_delta_t(e, i, l, 1.0, alpha, 0.5)
}
/// Interaction formula for combined axial load and bending moment.
///
/// Beam-column interaction according to the linear interaction formula:
///
/// P/P_cr + M/M_cr ≤ 1
///
/// Returns the utilization ratio (1.0 = critical, < 1.0 = safe, > 1.0 = failed).
pub fn combined_loading_utilization(p: f64, p_cr: f64, m: f64, m_cr: f64) -> f64 {
    p / p_cr.max(1e-30) + m / m_cr.max(1e-30)
}
/// Modified interaction formula (Dutheil) for eccentric compression:
///
/// P/P_cr + (M + P·e) / M_cr ≤ 1
///
/// where e is the eccentricity.
pub fn eccentric_compression_utilization(p: f64, p_cr: f64, m0: f64, m_cr: f64, e: f64) -> f64 {
    let m_total = m0 + p * e;
    combined_loading_utilization(p, p_cr, m_total, m_cr)
}
/// Axial load amplification factor for beam-columns (moment amplification).
///
/// Used to compute second-order moments from first-order analysis:
///
/// M_II = M_I / (1 − P/P_cr)
///
/// Returns `None` if P ≥ P_cr (unstable).
pub fn moment_amplification_factor(p: f64, p_cr: f64) -> Option<f64> {
    if p >= p_cr {
        return None;
    }
    Some(1.0 / (1.0 - p / p_cr))
}
/// Multi-mode buckling load considering interaction between two buckling modes.
///
/// Uses Dunkerley's method (lower bound):
///
/// 1/P_cr_combined = 1/P_cr1 + 1/P_cr2
pub fn dunkerley_combined_buckling(p_cr1: f64, p_cr2: f64) -> f64 {
    if p_cr1.abs() < 1e-30 || p_cr2.abs() < 1e-30 {
        return 0.0;
    }
    1.0 / (1.0 / p_cr1 + 1.0 / p_cr2)
}
/// Southwell plot method: estimates critical load from pre-buckling measurements.
///
/// The Southwell plot linearizes the load-displacement response:
///   w / P = w / P_cr + e_0 / P_cr
///
/// Given two measurement pairs (P_i, w_i), estimates P_cr as:
///
///   P_cr = (w2 − w1) / (w2/P2 − w1/P1)
///
/// Returns `None` if the denominator is too small.
pub fn southwell_critical_load(p1: f64, w1: f64, p2: f64, w2: f64) -> Option<f64> {
    let denom = w2 / p2 - w1 / p1;
    if denom.abs() < 1e-30 {
        return None;
    }
    Some((w2 - w1) / denom)
}
/// Reserve factor: ratio of critical load to applied load.
///
/// RF > 1.0 means the structure is safe against buckling.
pub fn reserve_factor(p_cr: f64, p_applied: f64) -> f64 {
    if p_applied.abs() < 1e-30 {
        return f64::INFINITY;
    }
    p_cr / p_applied
}
/// Check if a column satisfies the buckling safety requirement.
///
/// Returns `true` if P_cr / P_applied ≥ safety_factor.
pub fn buckling_safety_check(p_cr: f64, p_applied: f64, safety_factor: f64) -> bool {
    reserve_factor(p_cr, p_applied) >= safety_factor
}
/// Classify a column by slenderness.
///
/// * `slenderness` – slenderness ratio λ = K·L / r
/// * `sigma_y`     – yield stress
/// * `e`           – Young's modulus
pub fn classify_slenderness(slenderness: f64, sigma_y: f64, e: f64) -> SlendernessClass {
    let lambda_e = PI * (2.0 * e / sigma_y).sqrt();
    if slenderness <= 0.0 {
        SlendernessClass::Short
    } else if slenderness < lambda_e {
        SlendernessClass::Medium
    } else {
        SlendernessClass::Long
    }
}
