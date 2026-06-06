//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use crate::colormap::{Colormap, map_scalar};
use crate::primitives::{Color, LinePrimitive};
use oxiphysics_core::math::Vec3;
use std::f64::consts::PI;

/// A trajectory line following the direction of a principal stress.
#[derive(Debug, Clone)]
pub struct StressTrajectory {
    /// Sequence of 3D points along the trajectory.
    pub points: Vec<[f64; 3]>,
    /// Which principal stress (0=max, 1=mid, 2=min) this follows.
    pub principal_index: usize,
}
/// Parameters for a superquadric tensor glyph.
#[derive(Debug, Clone)]
pub struct SuperquadricGlyph {
    /// Centre position.
    pub position: [f64; 3],
    /// Semi-axes (eigenvalue magnitudes).
    pub semi_axes: [f64; 3],
    /// Principal directions (each row is a unit eigenvector).
    pub axes: [[f64; 3]; 3],
    /// Superquadric exponents (alpha, beta) controlling shape.
    pub exponents: [f64; 2],
    /// Color.
    pub color: crate::primitives::Color,
}
impl SuperquadricGlyph {
    /// Create from a stress tensor at a given position.
    pub fn from_stress(
        stress: &StressTensor,
        position: [f64; 3],
        scale: f64,
        colormap: Colormap,
    ) -> Self {
        let principals = stress.principal_stresses();
        let axes = principal_eigenvectors(stress.voigt);
        let vm = stress.von_mises();
        let max_abs = principals.iter().map(|&s| s.abs()).fold(0.0f64, f64::max);
        let semi_axes = principals.map(|s| (s.abs() * scale).max(1e-6));
        let lin = (principals[0] - principals[1]).abs() / (max_abs + 1e-30);
        let pla = (principals[1] - principals[2]).abs() / (max_abs + 1e-30);
        let alpha = 1.0 - lin;
        let beta = 1.0 - pla;
        Self {
            position,
            semi_axes,
            axes,
            exponents: [alpha.clamp(0.1, 2.0), beta.clamp(0.1, 2.0)],
            color: map_scalar(vm, 0.0, max_abs.max(1e-30), colormap),
        }
    }
    /// Evaluate the superquadric implicit function at a local point.
    /// Value < 1 is inside, > 1 is outside.
    pub fn implicit(&self, local_point: [f64; 3]) -> f64 {
        let [ax, ay, az] = self.semi_axes;
        let [alpha, beta] = self.exponents;
        let px = (local_point[0] / ax.max(1e-30)).abs().powf(2.0 / beta);
        let py = (local_point[1] / ay.max(1e-30)).abs().powf(2.0 / beta);
        let pz_base = (px + py).powf(beta / alpha);
        let pz = (local_point[2] / az.max(1e-30)).abs().powf(2.0 / alpha);
        pz_base + pz
    }
}
/// A tensor glyph representing a stress/strain state at a point.
///
/// The glyph is an ellipsoid whose axes are aligned with the principal
/// stress directions and scaled by the principal stress magnitudes.
#[derive(Debug, Clone)]
pub struct TensorGlyph {
    /// Position of the glyph center.
    pub position: [f64; 3],
    /// Principal stress magnitudes (eigenvalues).
    pub magnitudes: [f64; 3],
    /// Principal directions (eigenvectors) as columns.
    /// Each entry is a 3-component direction vector.
    pub directions: [[f64; 3]; 3],
    /// Color assigned to this glyph.
    pub color: Color,
}
/// Which stress invariant to use for coloring.
#[derive(Debug, Clone, Copy)]
pub enum StressInvariant {
    /// First invariant I1 = trace(sigma).
    I1,
    /// Second deviatoric invariant J2.
    J2,
    /// Third deviatoric invariant J3.
    J3,
    /// Von Mises equivalent stress.
    VonMises,
    /// Hydrostatic stress.
    Hydrostatic,
    /// Stress triaxiality.
    Triaxiality,
    /// Lode angle.
    LodeAngle,
}
/// A stress path records the evolution of stress state over time.
///
/// Typically plotted in p-q space (hydrostatic stress vs deviatoric stress)
/// or in principal stress space.
#[derive(Debug, Clone)]
pub struct StressPath {
    /// Time values.
    pub times: Vec<f64>,
    /// Stress states at each time.
    pub states: Vec<StressTensor>,
}
impl StressPath {
    /// Create a new empty stress path.
    pub fn new() -> Self {
        Self {
            times: Vec::new(),
            states: Vec::new(),
        }
    }
    /// Record a stress state at the given time.
    pub fn record(&mut self, time: f64, stress: StressTensor) {
        self.times.push(time);
        self.states.push(stress);
    }
    /// Number of recorded states.
    pub fn len(&self) -> usize {
        self.states.len()
    }
    /// Whether the path is empty.
    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }
    /// Extract p-q data for plotting (hydrostatic stress vs von Mises stress).
    ///
    /// Returns pairs (p, q) for each recorded state.
    pub fn pq_data(&self) -> Vec<(f64, f64)> {
        self.states
            .iter()
            .map(|s| (s.hydrostatic(), s.von_mises()))
            .collect()
    }
    /// Extract triaxiality over time.
    pub fn triaxiality_history(&self) -> Vec<(f64, f64)> {
        self.times
            .iter()
            .zip(self.states.iter())
            .map(|(&t, s)| (t, s.triaxiality()))
            .collect()
    }
    /// Extract Lode angle over time.
    pub fn lode_angle_history(&self) -> Vec<(f64, f64)> {
        self.times
            .iter()
            .zip(self.states.iter())
            .map(|(&t, s)| (t, s.lode_angle()))
            .collect()
    }
    /// Maximum von Mises stress encountered along the path.
    pub fn max_von_mises(&self) -> f64 {
        self.states
            .iter()
            .map(|s| s.von_mises())
            .fold(0.0_f64, f64::max)
    }
    /// Generate line primitives for stress path visualization in 3D.
    ///
    /// Maps each stress state to a 3D point using principal stresses as coordinates.
    pub fn to_principal_space_lines(&self, scale: f64) -> Vec<LinePrimitive> {
        if self.states.len() < 2 {
            return Vec::new();
        }
        let mut lines = Vec::with_capacity(self.states.len() - 1);
        for i in 0..self.states.len() - 1 {
            let p0 = self.states[i].principal_stresses();
            let p1 = self.states[i + 1].principal_stresses();
            lines.push(LinePrimitive {
                start: Vec3::new(p0[0] * scale, p0[1] * scale, p0[2] * scale),
                end: Vec3::new(p1[0] * scale, p1[1] * scale, p1[2] * scale),
                color: Color::white(),
            });
        }
        lines
    }
}
/// A field of stress tensors defined at integration points with mesh connectivity.
#[derive(Debug, Clone)]
pub struct StressField {
    /// Stress tensors at each point/node.
    pub stresses: Vec<StressTensor>,
    /// Element connectivity: each element is a list of node indices.
    pub connectivity: Vec<Vec<usize>>,
}
impl StressField {
    /// Create a new stress field.
    pub fn new(stresses: Vec<StressTensor>, connectivity: Vec<Vec<usize>>) -> Self {
        Self {
            stresses,
            connectivity,
        }
    }
    /// Compute the average stress across all nodes.
    pub fn average_stress(&self) -> StressTensor {
        if self.stresses.is_empty() {
            return StressTensor::zero();
        }
        let mut avg = [0.0_f64; 6];
        for s in &self.stresses {
            for (a, &sv) in avg.iter_mut().zip(s.voigt.iter()) {
                *a += sv;
            }
        }
        let n = self.stresses.len() as f64;
        for v in avg.iter_mut() {
            *v /= n;
        }
        StressTensor::new(avg)
    }
    /// Maximum von Mises stress in the field.
    pub fn max_von_mises(&self) -> f64 {
        self.stresses
            .iter()
            .map(|s| s.von_mises())
            .fold(0.0_f64, f64::max)
    }
}
/// Data for Mohr's circle visualization.
///
/// For a 3D stress state there are three Mohr's circles defined by pairs
/// of principal stresses.
#[derive(Debug, Clone, Copy)]
pub struct MohrCircleData {
    /// Principal stresses in descending order \[s1, s2, s3\].
    pub principals: [f64; 3],
    /// Centers of the three circles: (s1+s3)/2, (s1+s2)/2, (s2+s3)/2.
    pub centers: [f64; 3],
    /// Radii of the three circles: (s1-s3)/2, (s1-s2)/2, (s2-s3)/2.
    pub radii: [f64; 3],
}
impl MohrCircleData {
    /// Compute Mohr's circle data from a stress tensor in Voigt notation.
    pub fn from_voigt(s: [f64; 6]) -> Self {
        let principals = principal_stresses(s);
        let s1 = principals[0];
        let s2 = principals[1];
        let s3 = principals[2];
        let centers = [(s1 + s3) / 2.0, (s1 + s2) / 2.0, (s2 + s3) / 2.0];
        let radii = [(s1 - s3) / 2.0, (s1 - s2) / 2.0, (s2 - s3) / 2.0];
        Self {
            principals,
            centers,
            radii,
        }
    }
    /// Maximum shear stress (radius of the outer circle).
    pub fn max_shear(&self) -> f64 {
        self.radii[0]
    }
    /// Mean normal stress (center of the outer circle).
    pub fn mean_normal(&self) -> f64 {
        self.centers[0]
    }
    /// Generate points on the outer Mohr's circle for plotting.
    ///
    /// Returns `n_points` pairs of (sigma, tau) on the circle.
    pub fn outer_circle_points(&self, n_points: usize) -> Vec<(f64, f64)> {
        let center = self.centers[0];
        let radius = self.radii[0];
        let mut points = Vec::with_capacity(n_points);
        for i in 0..n_points {
            let theta = 2.0 * PI * (i as f64) / (n_points as f64);
            let sigma = center + radius * theta.cos();
            let tau = radius * theta.sin();
            points.push((sigma, tau));
        }
        points
    }
}
/// Symmetric stress tensor in Voigt notation: \[s11, s22, s33, s12, s23, s13\].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StressTensor {
    /// Voigt components \[s11, s22, s33, s12, s23, s13\].
    pub voigt: [f64; 6],
}
impl StressTensor {
    /// Create a new stress tensor from Voigt components.
    pub fn new(voigt: [f64; 6]) -> Self {
        Self { voigt }
    }
    /// Zero stress tensor.
    pub fn zero() -> Self {
        Self { voigt: [0.0; 6] }
    }
    /// Create a uniaxial stress state.
    pub fn uniaxial(sigma: f64) -> Self {
        Self {
            voigt: [sigma, 0.0, 0.0, 0.0, 0.0, 0.0],
        }
    }
    /// Create a hydrostatic (isotropic) stress state.
    pub fn hydrostatic_state(pressure: f64) -> Self {
        Self {
            voigt: [pressure, pressure, pressure, 0.0, 0.0, 0.0],
        }
    }
    /// Create a pure shear stress state.
    pub fn pure_shear(tau: f64) -> Self {
        Self {
            voigt: [0.0, 0.0, 0.0, tau, 0.0, 0.0],
        }
    }
    /// Compute the three principal stresses (eigenvalues) via the Cardano formula.
    pub fn principal_stresses(&self) -> [f64; 3] {
        principal_stresses(self.voigt)
    }
    /// Von Mises equivalent stress.
    pub fn von_mises(&self) -> f64 {
        von_mises_stress(self.voigt)
    }
    /// Hydrostatic (mean normal) stress.
    pub fn hydrostatic(&self) -> f64 {
        hydrostatic_stress(self.voigt)
    }
    /// Deviatoric stress in Voigt notation.
    pub fn deviatoric(&self) -> [f64; 6] {
        deviatoric_stress(self.voigt)
    }
    /// Stress triaxiality: eta = sigma_hydro / sigma_vm.
    pub fn triaxiality(&self) -> f64 {
        stress_triaxiality(self.voigt)
    }
    /// Lode angle in \[0, pi/3\].
    pub fn lode_angle(&self) -> f64 {
        lode_angle(self.voigt)
    }
    /// First stress invariant I1 = trace(sigma).
    pub fn i1(&self) -> f64 {
        self.voigt[0] + self.voigt[1] + self.voigt[2]
    }
    /// Second stress invariant I2.
    pub fn i2(&self) -> f64 {
        let s = self.voigt;
        s[0] * s[1] + s[1] * s[2] + s[2] * s[0] - s[3] * s[3] - s[4] * s[4] - s[5] * s[5]
    }
    /// Third stress invariant I3 = det(sigma).
    pub fn i3(&self) -> f64 {
        let s = self.voigt;
        let mat = [[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]];
        det3(mat)
    }
    /// Frobenius norm of the stress tensor.
    pub fn frobenius_norm(&self) -> f64 {
        let s = self.voigt;
        (s[0] * s[0] + s[1] * s[1] + s[2] * s[2] + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]))
            .sqrt()
    }
    /// Octahedral shear stress.
    ///
    /// `tau_oct = sqrt(2) / 3 * sigma_vm`
    pub fn octahedral_shear(&self) -> f64 {
        2.0_f64.sqrt() / 3.0 * self.von_mises()
    }
}
/// Identifies which yield criterion to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum YieldCriterion {
    /// von Mises (J2) yield surface.
    VonMises,
    /// Tresca (maximum shear stress) yield surface.
    Tresca,
    /// Mohr-Coulomb pressure-dependent yield surface.
    MohrCoulomb,
}
/// Plastic strain state at a material point.
#[derive(Debug, Clone, Copy)]
pub struct PlasticStrain {
    /// Accumulated equivalent plastic strain (dimensionless).
    pub equivalent_plastic_strain: f64,
    /// Plastic strain in Voigt notation.
    pub voigt: [f64; 6],
    /// Is this material point yielded?
    pub yielded: bool,
}
impl PlasticStrain {
    /// Zero plastic strain (elastic state).
    pub fn zero() -> Self {
        Self {
            equivalent_plastic_strain: 0.0,
            voigt: [0.0; 6],
            yielded: false,
        }
    }
    /// Create with given equivalent plastic strain.
    pub fn new(eps_eq: f64, voigt: [f64; 6]) -> Self {
        Self {
            equivalent_plastic_strain: eps_eq,
            voigt,
            yielded: eps_eq > 0.0,
        }
    }
    /// Increment equivalent plastic strain.
    pub fn accumulate(&mut self, d_eps: f64) {
        self.equivalent_plastic_strain += d_eps.abs();
        self.yielded = self.equivalent_plastic_strain > 0.0;
    }
}
/// Hydrodynamic stress state (pressure + viscous deviatoric).
#[derive(Debug, Clone, Copy)]
pub struct HydroStress {
    /// Thermodynamic pressure p (Pa, positive = compression).
    pub pressure: f64,
    /// Viscous stress in Voigt notation.
    pub viscous: [f64; 6],
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
}
impl HydroStress {
    /// Create a hydrostatic (no-viscosity) state.
    pub fn hydrostatic(pressure: f64) -> Self {
        Self {
            pressure,
            viscous: [0.0; 6],
            viscosity: 0.0,
        }
    }
    /// Total Cauchy stress tensor in Voigt notation.
    pub fn total_stress(&self) -> StressTensor {
        let p = self.pressure;
        StressTensor::new([
            -p + self.viscous[0],
            -p + self.viscous[1],
            -p + self.viscous[2],
            self.viscous[3],
            self.viscous[4],
            self.viscous[5],
        ])
    }
    /// Dynamic pressure contribution: 0.5 * rho * v^2 (stagnation).
    pub fn dynamic_pressure(density: f64, speed: f64) -> f64 {
        0.5 * density * speed * speed
    }
    /// Von Mises stress of the viscous part.
    pub fn viscous_von_mises(&self) -> f64 {
        von_mises_stress(self.viscous)
    }
}
/// A contact force event at a surface point.
#[derive(Debug, Clone)]
pub struct ContactForce {
    /// Contact point in world space.
    pub position: [f64; 3],
    /// Contact normal (outward from the object).
    pub normal: [f64; 3],
    /// Normal force magnitude (N).
    pub normal_force: f64,
    /// Tangential (friction) force magnitude (N).
    pub tangential_force: f64,
    /// Body index A.
    pub body_a: usize,
    /// Body index B.
    pub body_b: usize,
}
impl ContactForce {
    /// Create a new contact force event.
    pub fn new(
        position: [f64; 3],
        normal: [f64; 3],
        normal_force: f64,
        tangential_force: f64,
        body_a: usize,
        body_b: usize,
    ) -> Self {
        Self {
            position,
            normal,
            normal_force,
            tangential_force,
            body_a,
            body_b,
        }
    }
    /// Total force magnitude.
    pub fn total_force(&self) -> f64 {
        (self.normal_force * self.normal_force + self.tangential_force * self.tangential_force)
            .sqrt()
    }
    /// Coefficient of friction (tangential/normal). Returns 0 if normal is zero.
    pub fn friction_ratio(&self) -> f64 {
        if self.normal_force < 1e-30 {
            return 0.0;
        }
        self.tangential_force / self.normal_force
    }
}
/// Damage variable D in \[0, 1\]: D=0 undamaged, D=1 fully failed.
#[derive(Debug, Clone)]
pub struct DamageField {
    /// Per-node damage variables.
    pub damage: Vec<f64>,
}
impl DamageField {
    /// Create a new damage field of `n` nodes (all undamaged initially).
    pub fn new(n: usize) -> Self {
        Self {
            damage: vec![0.0; n],
        }
    }
    /// Update damage from a Palmgren-Miner cumulative damage rule.
    ///
    /// `applied_cycles[i]` is the number of cycles applied at stress level i.
    /// `nf[i]` is the cycles-to-failure at that stress level.
    /// Increments `self.damage[node]` by Σ(n_i / N_f_i).
    pub fn accumulate_miner(&mut self, node: usize, applied_cycles: &[f64], nf: &[f64]) {
        if node >= self.damage.len() {
            return;
        }
        let inc: f64 = applied_cycles
            .iter()
            .zip(nf.iter())
            .map(|(&n, &nfi)| if nfi > 0.0 { n / nfi } else { 0.0 })
            .sum();
        self.damage[node] = (self.damage[node] + inc).min(1.0);
    }
    /// Return the maximum damage variable in the field.
    pub fn max_damage(&self) -> f64 {
        self.damage.iter().cloned().fold(0.0_f64, f64::max)
    }
    /// Colours based on damage value: blue=undamaged, red=fully damaged.
    pub fn to_colors(&self, colormap: crate::colormap::Colormap) -> Vec<Color> {
        self.damage
            .iter()
            .map(|&d| map_scalar(d, 0.0, 1.0, colormap))
            .collect()
    }
}
/// Strain tensor in Voigt notation: \[e11, e22, e33, e12, e23, e13\].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StrainTensor {
    /// Voigt components \[e11, e22, e33, e12, e23, e13\].
    pub voigt: [f64; 6],
}
impl StrainTensor {
    /// Create a new strain tensor.
    pub fn new(voigt: [f64; 6]) -> Self {
        Self { voigt }
    }
    /// Zero strain.
    pub fn zero() -> Self {
        Self { voigt: [0.0; 6] }
    }
    /// Volumetric strain: ev = e11 + e22 + e33.
    pub fn volumetric(&self) -> f64 {
        self.voigt[0] + self.voigt[1] + self.voigt[2]
    }
    /// Equivalent (von Mises) strain.
    pub fn von_mises_strain(&self) -> f64 {
        let e = self.voigt;
        let term1 = (e[0] - e[1]).powi(2) + (e[1] - e[2]).powi(2) + (e[2] - e[0]).powi(2);
        let term2 = 6.0 * (e[3].powi(2) + e[4].powi(2) + e[5].powi(2));
        ((term1 + term2) / 3.0).sqrt() * std::f64::consts::SQRT_2 / std::f64::consts::SQRT_2
    }
}
