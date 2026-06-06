//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// Auxetic material with negative Poisson's ratio.
///
/// These materials expand laterally when stretched, arising from re-entrant
/// or chiral microstructures.
#[derive(Debug, Clone)]
pub struct AuxeticMaterial {
    /// Poisson's ratio (negative for auxetics, typically −1…0).
    pub poisson_ratio: f64,
    /// Young's modulus of the solid constituent (Pa).
    pub young_modulus: f64,
    /// Characteristic cell angle (radians); drives the re-entrant geometry.
    pub cell_angle: f64,
}
impl AuxeticMaterial {
    /// Create a new [`AuxeticMaterial`].
    pub fn new(poisson_ratio: f64, young_modulus: f64, cell_angle: f64) -> Self {
        Self {
            poisson_ratio,
            young_modulus,
            cell_angle,
        }
    }
    /// Effective Poisson's ratio, modified by the cell angle.
    ///
    /// The sign and magnitude are modulated by `sin(cell_angle) * cos(cell_angle)`
    /// to capture the re-entrant geometry effect.
    pub fn effective_poisson(&self) -> f64 {
        let theta = self.cell_angle;
        let projection = (theta.sin() * theta.cos()).clamp(-1.0, 1.0);
        self.poisson_ratio * (1.0 + projection.abs())
    }
    /// Effective stiffness using the Gibson–Ashby scaling for re-entrant honeycomb.
    ///
    /// `density_ratio` is the relative density ρ*/ρ_s (0 < density_ratio ≤ 1).
    /// Returns the effective Young's modulus (Pa).
    pub fn effective_stiffness(&self, density_ratio: f64) -> f64 {
        let exponent = 2.0 + self.effective_poisson().abs().min(1.0);
        let c1 = 0.5;
        self.young_modulus * c1 * density_ratio.powf(exponent)
    }
}
/// Topological insulator (mechanical or acoustic).
///
/// Models the Su-Schrieffer-Heeger (SSH) chain or Chern insulator in the
/// context of phononic/mechanical metamaterials.  Computes the band structure,
/// band gap, and Chern number.
#[derive(Debug, Clone)]
pub struct TopologicalInsulator {
    /// Intra-cell coupling stiffness k₁ (N/m or normalised).
    pub k_intra: f64,
    /// Inter-cell coupling stiffness k₂ (N/m or normalised).
    pub k_inter: f64,
    /// On-site mass m (kg or normalised).
    pub mass: f64,
    /// Number of unit cells.
    pub n_cells: usize,
    /// Topological invariant.
    pub invariant: TopologicalInvariant,
    /// Computed band gap (normalised frequency units).
    pub band_gap: f64,
}
impl TopologicalInsulator {
    /// Create a new [`TopologicalInsulator`] SSH chain.
    ///
    /// # Arguments
    /// * `k_intra` — intra-cell spring constant
    /// * `k_inter` — inter-cell spring constant
    /// * `mass`    — site mass
    /// * `n_cells` — number of unit cells
    ///
    /// The Chern number is 1 when `k_inter > k_intra` (topological phase)
    /// and 0 otherwise (trivial phase).
    pub fn new(k_intra: f64, k_inter: f64, mass: f64, n_cells: usize) -> Self {
        let (chern, z2, n_edge) = if k_inter > k_intra {
            (1, 1u8, 1usize)
        } else {
            (0, 0u8, 0usize)
        };
        let invariant = TopologicalInvariant::new(chern, z2, n_edge);
        let gap = Self::compute_band_gap_static(k_intra, k_inter, mass);
        Self {
            k_intra,
            k_inter,
            mass,
            n_cells,
            invariant,
            band_gap: gap,
        }
    }
    /// Band gap (normalised frequency): Δω = 2|k₂ − k₁| / m.
    fn compute_band_gap_static(k_intra: f64, k_inter: f64, mass: f64) -> f64 {
        if mass <= 0.0 {
            return 0.0;
        }
        2.0 * (k_inter - k_intra).abs() / mass
    }
    /// Recompute and update the band gap.
    pub fn update_band_gap(&mut self) -> f64 {
        self.band_gap = Self::compute_band_gap_static(self.k_intra, self.k_inter, self.mass);
        self.band_gap
    }
    /// Dispersion relation ω(k) at Bloch wave vector `k` (in units of π/a).
    ///
    /// SSH model: ω²(k) = (k₁² + k₂² + 2 k₁ k₂ cos(k π)) / m²
    pub fn dispersion(&self, k_norm: f64) -> f64 {
        let k1 = self.k_intra;
        let k2 = self.k_inter;
        let m = self.mass;
        if m <= 0.0 {
            return 0.0;
        }
        let ka = k_norm * PI;
        let omega_sq = (k1 * k1 + k2 * k2 + 2.0 * k1 * k2 * ka.cos()) / (m * m);
        omega_sq.sqrt()
    }
    /// Edge state frequency (appears in the band gap for topological phase).
    ///
    /// In the SSH model the edge state sits at ω_edge = sqrt(k₁ · k₂) / m.
    pub fn edge_state_frequency(&self) -> Option<f64> {
        if self.invariant.chern_number == 0 {
            return None;
        }
        let freq = (self.k_intra * self.k_inter).sqrt() / self.mass;
        Some(freq)
    }
    /// Berry phase for the lower band (0 or π).
    ///
    /// In the SSH model: Berry phase = 0 for trivial, π for topological.
    pub fn berry_phase(&self) -> f64 {
        if self.k_inter > self.k_intra { PI } else { 0.0 }
    }
    /// Change coupling constants and update topology.
    pub fn set_couplings(&mut self, k_intra: f64, k_inter: f64) {
        self.k_intra = k_intra;
        self.k_inter = k_inter;
        let (chern, z2, n_edge) = if k_inter > k_intra {
            (1, 1u8, 1usize)
        } else {
            (0, 0u8, 0usize)
        };
        self.invariant = TopologicalInvariant::new(chern, z2, n_edge);
        self.update_band_gap();
    }
}
/// Transformation-acoustics cloak based on a coordinate transformation.
///
/// Given the Jacobian **J** of the spatial mapping, computes the anisotropic
/// effective density tensor and effective bulk modulus required by the cloak.
#[derive(Debug, Clone)]
pub struct TransformationAcoustics {
    /// 3×3 Jacobian matrix of the coordinate transformation.
    pub jacobian: [[f64; 3]; 3],
}
impl TransformationAcoustics {
    /// Create a new [`TransformationAcoustics`] cloak.
    pub fn new(jacobian: [[f64; 3]; 3]) -> Self {
        Self { jacobian }
    }
    /// Determinant of the 3×3 Jacobian.
    fn det_j(&self) -> f64 {
        let j = &self.jacobian;
        j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1])
            - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
            + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0])
    }
    /// Effective density tensor ρ_ij = ρ₀ · (J·J^T)_ij / det(J).
    ///
    /// `rho0` is the background density (kg/m³).
    pub fn effective_density_tensor(&self, rho0: f64) -> [[f64; 3]; 3] {
        let j = &self.jacobian;
        let det = self.det_j();
        let mut jjt = [[0.0f64; 3]; 3];
        for (i, jjt_row) in jjt.iter_mut().enumerate() {
            for (k, jjt_ik) in jjt_row.iter_mut().enumerate() {
                for (l, &j_kl) in j[k].iter().enumerate() {
                    *jjt_ik += j[i][l] * j_kl;
                }
            }
        }
        let scale = if det.abs() < 1e-30 { 0.0 } else { rho0 / det };
        let mut out = [[0.0f64; 3]; 3];
        for (out_row, jjt_row) in out.iter_mut().zip(jjt.iter()) {
            for (o, &jv) in out_row.iter_mut().zip(jjt_row.iter()) {
                *o = scale * jv;
            }
        }
        out
    }
    /// Effective bulk modulus κ_eff = κ₀ / det(J).
    ///
    /// `k0` is the background bulk modulus (Pa).
    pub fn effective_bulk_modulus(&self, k0: f64) -> f64 {
        let det = self.det_j();
        if det.abs() < 1e-30 { 0.0 } else { k0 / det }
    }
}
/// Elastic metamaterial with extreme or negative effective elastic properties.
///
/// Models include:
/// - Pentamode lattice (near-fluid solid, five near-zero eigenvalues)
/// - Inertial amplification (large effective mass through lever geometry)
/// - Negative Poisson's ratio (auxetic, re-entrant rib geometry)
#[derive(Debug, Clone)]
pub struct ElasticMetamaterial {
    /// Effective Young's modulus (Pa). Can be very small (pentamode) or large.
    pub effective_modulus: f64,
    /// Effective Poisson's ratio. Negative for auxetic unit cells.
    pub effective_poisson: f64,
    /// Effective mass density (kg/m³). Can be amplified via lever mechanism.
    pub effective_density: f64,
    /// Inertial amplification factor (ratio of effective to host mass density).
    pub inertial_amplification: f64,
    /// Rib angle for the auxetic re-entrant geometry (radians).
    pub rib_angle: f64,
    /// Host material Young's modulus E_s (Pa).
    pub host_modulus: f64,
    /// Host material density (kg/m³).
    pub host_density: f64,
}
impl ElasticMetamaterial {
    /// Create a new [`ElasticMetamaterial`].
    pub fn new(
        effective_modulus: f64,
        effective_poisson: f64,
        effective_density: f64,
        inertial_amplification: f64,
        rib_angle: f64,
        host_modulus: f64,
        host_density: f64,
    ) -> Self {
        Self {
            effective_modulus,
            effective_poisson,
            effective_density,
            inertial_amplification,
            rib_angle,
            host_modulus,
            host_density,
        }
    }
    /// Longitudinal wave speed (m/s) in the effective medium.
    ///
    /// ```text
    /// c_L = sqrt(E_eff / ρ_eff)
    /// ```
    pub fn longitudinal_wave_speed(&self) -> f64 {
        if self.effective_density <= 0.0 {
            return 0.0;
        }
        (self.effective_modulus / self.effective_density).sqrt()
    }
    /// Shear wave speed (m/s) in the effective medium.
    ///
    /// ```text
    /// c_S = sqrt(G_eff / ρ_eff)
    /// ```
    ///
    /// where `G_eff = E_eff / (2 (1 + ν_eff))`.
    pub fn shear_wave_speed(&self) -> f64 {
        if self.effective_density <= 0.0 {
            return 0.0;
        }
        let denom = 2.0 * (1.0 + self.effective_poisson);
        if denom <= 0.0 {
            return 0.0;
        }
        let g_eff = self.effective_modulus / denom;
        (g_eff / self.effective_density).sqrt()
    }
    /// Effective bulk modulus K_eff (Pa).
    ///
    /// ```text
    /// K = E / (3 (1 − 2ν))
    /// ```
    pub fn bulk_modulus(&self) -> f64 {
        let denom = 3.0 * (1.0 - 2.0 * self.effective_poisson);
        if denom.abs() < 1e-15 {
            return f64::INFINITY;
        }
        self.effective_modulus / denom
    }
    /// Inertially amplified effective density (kg/m³).
    ///
    /// In an inertial amplification unit cell, the effective density is
    /// multiplied by the amplification factor.
    pub fn amplified_density(&self) -> f64 {
        self.host_density * self.inertial_amplification
    }
    /// Band gap frequency range (Hz) due to inertial amplification.
    ///
    /// The band gap opens near the local resonance frequency of the amplified
    /// mass connected to a spring of stiffness `k` (N/m).
    ///
    /// Returns `(f_lower, f_upper)` in Hz.
    pub fn inertial_bandgap(&self, k: f64) -> (f64, f64) {
        let m_eff = self.amplified_density();
        if m_eff <= 0.0 || k <= 0.0 {
            return (0.0, 0.0);
        }
        let omega0 = (k / m_eff).sqrt();
        let f0 = omega0 / (2.0 * PI);
        let f1 = f0 * (1.0 + 1.0 / self.inertial_amplification).sqrt();
        (f0, f1)
    }
    /// Check whether this metamaterial has negative Poisson's ratio.
    pub fn is_auxetic(&self) -> bool {
        self.effective_poisson < 0.0
    }
    /// Re-entrant honeycomb Poisson's ratio from rib angle (Masters–Evans).
    ///
    /// Updates and returns `self.effective_poisson` based on `self.rib_angle`.
    pub fn compute_auxetic_poisson(&mut self, h_over_l: f64) -> f64 {
        let theta = self.rib_angle;
        let sin_t = theta.sin();
        let cos2 = theta.cos().powi(2);
        if cos2.abs() < 1e-15 {
            self.effective_poisson = 0.0;
            return 0.0;
        }
        let nu = -sin_t * (h_over_l + sin_t) / cos2;
        self.effective_poisson = nu;
        nu
    }
}
/// Piezo-actuated resonator in an active metamaterial unit cell.
#[derive(Debug, Clone)]
pub struct PiezoResonator {
    /// Resonator mass (kg).
    pub mass: f64,
    /// Passive spring stiffness (N/m).
    pub k_passive: f64,
    /// Piezo-induced stiffness change ΔK (N/m; positive or negative).
    pub dk_piezo: f64,
    /// Control voltage (V).
    pub voltage: f64,
    /// Piezo coupling coefficient d₃₃ (m/V).
    pub d33: f64,
}
impl PiezoResonator {
    /// Create a new [`PiezoResonator`].
    pub fn new(mass: f64, k_passive: f64, dk_piezo: f64, voltage: f64, d33: f64) -> Self {
        Self {
            mass,
            k_passive,
            dk_piezo,
            voltage,
            d33,
        }
    }
    /// Effective stiffness at current voltage.
    pub fn effective_stiffness(&self) -> f64 {
        self.k_passive + self.dk_piezo * self.voltage
    }
    /// Resonance frequency (Hz) at current voltage.
    pub fn resonance_frequency(&self) -> f64 {
        let k_eff = self.effective_stiffness();
        if self.mass <= 0.0 || k_eff <= 0.0 {
            return 0.0;
        }
        (k_eff / self.mass).sqrt() / (2.0 * PI)
    }
    /// Elongation produced by piezo at `voltage` (m).
    pub fn piezo_elongation(&self) -> f64 {
        self.d33 * self.voltage
    }
}
/// Active metamaterial with piezo-actuated resonators and programmable stiffness.
///
/// Supports parity-time (PT) symmetry via balanced gain and loss, enabling
/// exceptional-point physics.
#[derive(Debug, Clone)]
pub struct ActiveMetamaterial {
    /// Array of piezo-actuated resonators (unit cells).
    pub resonators: Vec<PiezoResonator>,
    /// Host medium stiffness (N/m).
    pub k_host: f64,
    /// Host medium density (kg/m³).
    pub host_density: f64,
    /// Gain coefficient of the active elements (Neper/m).
    pub gain: f64,
    /// Loss coefficient of the passive elements (Neper/m).
    pub loss: f64,
    /// Whether PT symmetry is enforced (`gain == loss`).
    pub pt_symmetric: bool,
}
impl ActiveMetamaterial {
    /// Create a new [`ActiveMetamaterial`].
    pub fn new(
        resonators: Vec<PiezoResonator>,
        k_host: f64,
        host_density: f64,
        gain: f64,
        loss: f64,
    ) -> Self {
        let pt_symmetric = (gain - loss).abs() < 1e-12;
        Self {
            resonators,
            k_host,
            host_density,
            gain,
            loss,
            pt_symmetric,
        }
    }
    /// Mean resonance frequency across all resonators (Hz).
    pub fn mean_resonance_frequency(&self) -> f64 {
        if self.resonators.is_empty() {
            return 0.0;
        }
        self.resonators
            .iter()
            .map(|r| r.resonance_frequency())
            .sum::<f64>()
            / self.resonators.len() as f64
    }
    /// Bandwidth of programmable stiffness (min to max k_eff) across all resonators (N/m).
    pub fn stiffness_range(&self) -> (f64, f64) {
        if self.resonators.is_empty() {
            return (0.0, 0.0);
        }
        let k_min = self
            .resonators
            .iter()
            .map(|r| r.effective_stiffness())
            .fold(f64::INFINITY, f64::min);
        let k_max = self
            .resonators
            .iter()
            .map(|r| r.effective_stiffness())
            .fold(f64::NEG_INFINITY, f64::max);
        (k_min, k_max)
    }
    /// PT symmetry check: returns true when |gain − loss| < tolerance.
    pub fn check_pt_symmetry(&self) -> bool {
        (self.gain - self.loss).abs() < 1e-10
    }
    /// Exceptional point proximity: how close gain is to loss (normalised).
    ///
    /// Returns 0 at the exceptional point (gain = loss), increases away from it.
    pub fn exceptional_point_proximity(&self) -> f64 {
        let avg = 0.5 * (self.gain + self.loss);
        if avg == 0.0 {
            return 0.0;
        }
        (self.gain - self.loss).abs() / avg
    }
    /// Set all resonator voltages to `v` (programmatic stiffness adjustment).
    pub fn set_all_voltages(&mut self, v: f64) {
        for r in self.resonators.iter_mut() {
            r.voltage = v;
        }
    }
    /// Net gain-loss balance for PT symmetry: (gain − loss).
    ///
    /// Zero means the system is exactly PT symmetric.
    pub fn gain_loss_balance(&self) -> f64 {
        self.gain - self.loss
    }
    /// Estimate the effective wave attenuation (Neper/m) in the host.
    ///
    /// For a balanced PT system (`gain = loss`) the wave propagates without
    /// net attenuation in the linear regime.
    pub fn effective_attenuation(&self) -> f64 {
        self.loss - self.gain
    }
}
/// Double-negative electromagnetic/acoustic medium (negative index material).
///
/// Both effective permittivity ε and permeability μ (or their acoustic
/// analogues) are simultaneously negative, leading to a negative refractive
/// index.
#[derive(Debug, Clone)]
pub struct DoublenegativeMedium {
    /// Effective permittivity / density analogue (can be negative).
    pub epsilon_eff: f64,
    /// Effective permeability / compressibility analogue (can be negative).
    pub mu_eff: f64,
}
impl DoublenegativeMedium {
    /// Create a new [`DoublenegativeMedium`].
    pub fn new(epsilon_eff: f64, mu_eff: f64) -> Self {
        Self {
            epsilon_eff,
            mu_eff,
        }
    }
    /// Refractive index n = −sqrt(ε·μ) (negative when both parameters negative).
    pub fn refractive_index(&self) -> f64 {
        let product = self.epsilon_eff * self.mu_eff;
        if product < 0.0 {
            0.0
        } else {
            let sign = if self.epsilon_eff < 0.0 && self.mu_eff < 0.0 {
                -1.0
            } else {
                1.0
            };
            sign * product.sqrt()
        }
    }
    /// Phase velocity v_p = c / n  (m/s, c = speed of light or reference speed).
    ///
    /// Uses c = 299_792_458 m/s.  Returns 0 if the refractive index is zero.
    pub fn phase_velocity(&self) -> f64 {
        let n = self.refractive_index();
        if n.abs() < 1e-30 {
            0.0
        } else {
            2.997_924_58e8 / n
        }
    }
}
/// Topological invariant descriptor for a band structure.
#[derive(Debug, Clone)]
pub struct TopologicalInvariant {
    /// Chern number Z (integer topological invariant).
    pub chern_number: i32,
    /// Z₂ topological index (0 = trivial, 1 = topological).
    pub z2_index: u8,
    /// Number of protected edge states per edge.
    pub n_edge_states: usize,
}
impl TopologicalInvariant {
    /// Create a new topological invariant.
    pub fn new(chern_number: i32, z2_index: u8, n_edge_states: usize) -> Self {
        Self {
            chern_number,
            z2_index,
            n_edge_states,
        }
    }
    /// Bulk-edge correspondence: the number of edge states equals |Chern number|.
    pub fn satisfies_bulk_edge_correspondence(&self) -> bool {
        self.n_edge_states == self.chern_number.unsigned_abs() as usize
    }
}
/// Locally resonant acoustic metamaterial (single-resonator model).
///
/// Consists of a periodic array of mass-spring resonators embedded in a host
/// medium.  Near the resonance frequency, the effective density becomes
/// negative, opening a sub-wavelength band gap.
#[derive(Debug, Clone)]
pub struct AcousticMetamaterial {
    /// Lattice constant / unit cell size (m).
    pub unit_cell_size: f64,
    /// Mass of each internal resonator (kg).
    pub resonator_mass: f64,
    /// Stiffness of each resonator spring (N/m).
    pub resonator_stiffness: f64,
    /// Density of the host medium (kg/m³).
    pub host_density: f64,
    /// Bulk modulus of the host medium (Pa).
    pub host_modulus: f64,
}
impl AcousticMetamaterial {
    /// Create a new [`AcousticMetamaterial`].
    pub fn new(
        unit_cell_size: f64,
        resonator_mass: f64,
        resonator_stiffness: f64,
        host_density: f64,
        host_modulus: f64,
    ) -> Self {
        Self {
            unit_cell_size,
            resonator_mass,
            resonator_stiffness,
            host_density,
            host_modulus,
        }
    }
    /// Angular resonance frequency ω₀ = sqrt(k/m) in rad/s.
    #[inline]
    fn omega0(&self) -> f64 {
        (self.resonator_stiffness / self.resonator_mass).sqrt()
    }
    /// Effective density at angular frequency `omega` (rad/s).
    ///
    /// Uses the single-resonator formula:
    /// ρ_eff = ρ_host · \[1 − ω₀²/(ω² − ω₀²)\]
    ///
    /// Near resonance (ω → ω₀) the effective density diverges and can become
    /// negative just above ω₀.
    pub fn effective_density(&self, omega: f64) -> f64 {
        let w0 = self.omega0();
        let denom = omega * omega - w0 * w0;
        if denom.abs() < 1e-30 {
            f64::INFINITY
        } else {
            self.host_density * (1.0 - w0 * w0 / denom)
        }
    }
    /// Effective bulk modulus at angular frequency `omega` (rad/s).
    ///
    /// For this model the bulk modulus is not frequency-dispersive (host only).
    pub fn effective_modulus(&self, _omega: f64) -> f64 {
        self.host_modulus
    }
    /// Band-gap frequency range in Hz: (f_lower, f_upper).
    ///
    /// Returns the interval where ρ_eff < 0 (sub-wavelength band gap).
    /// f_lower = ω₀ / (2π), f_upper is estimated from the mass ratio.
    pub fn band_gap_range(&self) -> (f64, f64) {
        let w0 = self.omega0();
        let f0 = w0 / (2.0 * PI);
        let cell_vol = self.unit_cell_size.powi(3);
        let mass_ratio = self.resonator_mass / (self.host_density * cell_vol);
        let f1 = f0 * (1.0 + mass_ratio).sqrt();
        (f0, f1)
    }
}
/// Split-ring resonator (SRR) parameters for electromagnetic metamaterials.
#[derive(Debug, Clone)]
pub struct SplitRingResonator {
    /// Ring radius (m).
    pub radius: f64,
    /// Gap width (m).
    pub gap_width: f64,
    /// Ring wire thickness (m).
    pub wire_thickness: f64,
    /// Unit cell size (m).
    pub unit_cell: f64,
    /// Ring inductance (H).
    pub inductance: f64,
    /// Gap capacitance (F).
    pub capacitance: f64,
}
impl SplitRingResonator {
    /// Create a new SRR.
    pub fn new(
        radius: f64,
        gap_width: f64,
        wire_thickness: f64,
        unit_cell: f64,
        inductance: f64,
        capacitance: f64,
    ) -> Self {
        Self {
            radius,
            gap_width,
            wire_thickness,
            unit_cell,
            inductance,
            capacitance,
        }
    }
    /// LC resonance angular frequency ω_0 (rad/s).
    ///
    /// ```text
    /// ω₀ = 1 / √(L · C)
    /// ```
    pub fn resonance_frequency(&self) -> f64 {
        if self.inductance <= 0.0 || self.capacitance <= 0.0 {
            return 0.0;
        }
        1.0 / (self.inductance * self.capacitance).sqrt()
    }
    /// Filling fraction F = π r² / a².
    pub fn filling_fraction(&self) -> f64 {
        let a = self.unit_cell;
        if a <= 0.0 {
            return 0.0;
        }
        PI * self.radius * self.radius / (a * a)
    }
}
/// Re-entrant honeycomb unit-cell geometry.
///
/// Parameters follow the notation of Masters & Evans (1996):
/// `theta` – cell angle (rad, negative for re-entrant), `l` – inclined rib
/// length, `h` – vertical rib length, `t` – rib thickness, `es` – solid
/// Young's modulus.
#[derive(Debug, Clone)]
pub struct ReentrantHoneycomb {
    /// Cell angle in radians (negative → auxetic behaviour).
    pub theta: f64,
    /// Length of the inclined ribs (m).
    pub l: f64,
    /// Length of the vertical ribs (m).
    pub h: f64,
    /// Rib thickness (m).
    pub t: f64,
    /// Young's modulus of the solid material (Pa).
    pub es: f64,
}
impl ReentrantHoneycomb {
    /// Create a new [`ReentrantHoneycomb`] unit cell.
    pub fn new(theta: f64, l: f64, h: f64, t: f64, es: f64) -> Self {
        Self { theta, l, h, t, es }
    }
    /// In-plane Poisson's ratio ν_xy (can be negative for re-entrant cells).
    ///
    /// Uses the Masters–Evans closed-form solution based on flexure of the
    /// inclined ribs.
    pub fn poisson_ratio_xy(&self) -> f64 {
        let sin_t = self.theta.sin();
        let cos_t = self.theta.cos();
        let h_over_l = self.h / self.l;
        let numerator = sin_t * (h_over_l + sin_t);
        let denominator = cos_t * cos_t;
        if denominator.abs() < 1e-15 {
            0.0
        } else {
            -numerator / denominator
        }
    }
    /// Normalised Young's modulus E_x (Pa) in the x-direction.
    ///
    /// Returns the absolute value of the Masters–Evans expression to ensure
    /// physically meaningful (positive) stiffness regardless of the sign of
    /// the cell angle.
    pub fn young_modulus_x(&self) -> f64 {
        let cos_t = self.theta.cos();
        let sin_t = self.theta.sin();
        let h_over_l = self.h / self.l;
        let t_over_l = self.t / self.l;
        let denom = (h_over_l + sin_t) * sin_t * sin_t;
        if denom.abs() < 1e-30 {
            0.0
        } else {
            (t_over_l.powi(3) * cos_t / denom * self.es).abs()
        }
    }
    /// Normalised Young's modulus E_y (Pa) in the y-direction.
    ///
    /// Returns the absolute value to ensure positive stiffness.
    pub fn young_modulus_y(&self) -> f64 {
        let cos_t = self.theta.cos();
        let sin_t = self.theta.sin();
        let h_over_l = self.h / self.l;
        let t_over_l = self.t / self.l;
        let denom = cos_t.powi(3);
        if denom.abs() < 1e-30 {
            0.0
        } else {
            (t_over_l.powi(3) * (h_over_l + sin_t) / denom * self.es).abs()
        }
    }
    /// In-plane shear modulus G_xy (Pa).
    ///
    /// Returns the absolute value to ensure positive stiffness.
    pub fn shear_modulus_xy(&self) -> f64 {
        let cos_t = self.theta.cos();
        let sin_t = self.theta.sin();
        let h_over_l = self.h / self.l;
        let t_over_l = self.t / self.l;
        let num = t_over_l.powi(3) * (h_over_l + sin_t);
        let denom = h_over_l * h_over_l * (1.0 + 2.0 * h_over_l) * cos_t;
        if denom.abs() < 1e-30 {
            0.0
        } else {
            (num / denom * self.es).abs()
        }
    }
}
/// Pentamode metamaterial — a solid with five near-zero eigenvalues in its
/// stiffness tensor, behaving like an acoustic fluid.
///
/// Described by Norris (2008) and used for acoustic cloaking.  The two
/// independent stiffness coefficients are c11 and c12 in Voigt notation.
#[derive(Debug, Clone)]
pub struct PentamodeMetamaterial {
    /// Stiffness coefficient C₁₁ (Pa).
    pub c11: f64,
    /// Stiffness coefficient C₁₂ (Pa).
    pub c12: f64,
}
impl PentamodeMetamaterial {
    /// Create a new [`PentamodeMetamaterial`].
    pub fn new(c11: f64, c12: f64) -> Self {
        Self { c11, c12 }
    }
    /// Compressional wave speed (m/s) at a given density `density` (kg/m³).
    ///
    /// For a pentamode medium c_p = sqrt(C₁₁ / ρ).
    pub fn sound_speed(&self, density: f64) -> f64 {
        if density <= 0.0 {
            0.0
        } else {
            (self.c11 / density).sqrt()
        }
    }
}
/// Electromagnetic metamaterial based on split-ring resonators and wire arrays.
///
/// Supports computation of effective permittivity ε(ω) and permeability μ(ω),
/// and identification of left-handed frequency bands.
#[derive(Debug, Clone)]
pub struct ElectromagneticMetamaterial {
    /// Embedded SRR array.
    pub srr: SplitRingResonator,
    /// Plasma frequency of the effective wire medium (rad/s).
    pub plasma_frequency: f64,
    /// Damping/loss rate for the permittivity (rad/s).
    pub gamma_e: f64,
    /// Damping/loss rate for the permeability (rad/s).
    pub gamma_m: f64,
}
impl ElectromagneticMetamaterial {
    /// Create a new [`ElectromagneticMetamaterial`].
    pub fn new(srr: SplitRingResonator, plasma_frequency: f64, gamma_e: f64, gamma_m: f64) -> Self {
        Self {
            srr,
            plasma_frequency,
            gamma_e,
            gamma_m,
        }
    }
    /// Effective permittivity at angular frequency ω (Drude model, real part).
    ///
    /// ```text
    /// ε(ω) = 1 − ω_p² / (ω² + i·γ_e·ω) ≈ 1 − ω_p² / (ω² + γ_e²)  [real part]
    /// ```
    pub fn effective_permittivity(&self, omega: f64) -> f64 {
        if omega == 0.0 {
            return f64::NEG_INFINITY;
        }
        let wp = self.plasma_frequency;
        let ge = self.gamma_e;
        1.0 - wp * wp / (omega * omega + ge * ge)
    }
    /// Effective permeability at angular frequency ω (Lorentz resonator, real part).
    ///
    /// ```text
    /// μ(ω) = 1 − F·ω² / (ω² − ω₀² + i·γ_m·ω)  [real part]
    /// ```
    pub fn effective_permeability(&self, omega: f64) -> f64 {
        let w0 = self.srr.resonance_frequency();
        let f = self.srr.filling_fraction();
        let gm = self.gamma_m;
        let denom = (omega * omega - w0 * w0).powi(2) + (gm * omega).powi(2);
        if denom < 1e-60 {
            return 1.0;
        }
        let real_num = omega * omega * (omega * omega - w0 * w0);
        1.0 - f * real_num / denom
    }
    /// Effective refractive index n = sqrt(ε · μ).
    ///
    /// Returns a negative value when both ε and μ are negative (left-handed medium).
    pub fn refractive_index(&self, omega: f64) -> f64 {
        let eps = self.effective_permittivity(omega);
        let mu = self.effective_permeability(omega);
        let product = eps * mu;
        if product < 0.0 {
            return 0.0;
        }
        let sign = if eps < 0.0 && mu < 0.0 { -1.0 } else { 1.0 };
        sign * product.sqrt()
    }
    /// Check whether the material is left-handed (negative index) at frequency `omega`.
    pub fn is_left_handed(&self, omega: f64) -> bool {
        self.effective_permittivity(omega) < 0.0 && self.effective_permeability(omega) < 0.0
    }
    /// Phase velocity v_p = c / n (m/s), where c = 2.997924e8 m/s.
    pub fn phase_velocity(&self, omega: f64) -> f64 {
        let n = self.refractive_index(omega);
        if n.abs() < 1e-30 {
            return 0.0;
        }
        2.997_924_58e8 / n
    }
    /// Bandwidth of the left-handed frequency window (rad/s).
    ///
    /// Scans a range \[omega_min, omega_max\] with `n_pts` points and returns the
    /// width of the largest contiguous band where both ε < 0 and μ < 0.
    pub fn left_handed_bandwidth(&self, omega_min: f64, omega_max: f64, n_pts: usize) -> f64 {
        if n_pts < 2 || omega_max <= omega_min {
            return 0.0;
        }
        let dw = (omega_max - omega_min) / (n_pts - 1) as f64;
        let mut max_width = 0.0f64;
        let mut current_width = 0.0f64;
        for i in 0..n_pts {
            let w = omega_min + i as f64 * dw;
            if self.is_left_handed(w) {
                current_width += dw;
            } else {
                max_width = max_width.max(current_width);
                current_width = 0.0;
            }
        }
        max_width.max(current_width)
    }
}
/// Acoustic cloaking shell designed via coordinate transformation.
///
/// A spherical cloak maps the region `r_1 < r < r_2` in virtual space to the
/// annulus `r_1 < r < r_2` in physical space, while pushing all wave energy
/// around the inner core.
#[derive(Debug, Clone)]
pub struct AcousticCloakShell {
    /// Inner radius of the cloak shell (m).
    pub r_inner: f64,
    /// Outer radius of the cloak shell (m).
    pub r_outer: f64,
    /// Background medium density ρ₀ (kg/m³).
    pub background_density: f64,
    /// Background medium bulk modulus κ₀ (Pa).
    pub background_modulus: f64,
    /// Number of discrete shell layers.
    pub n_layers: usize,
}
impl AcousticCloakShell {
    /// Create a new [`AcousticCloakShell`].
    pub fn new(
        r_inner: f64,
        r_outer: f64,
        background_density: f64,
        background_modulus: f64,
        n_layers: usize,
    ) -> Self {
        Self {
            r_inner,
            r_outer,
            background_density,
            background_modulus,
            n_layers,
        }
    }
    /// Effective radial density at radius `r` within the cloak shell.
    ///
    /// From linear coordinate transformation:
    ///
    /// ```text
    /// ρ_r(r) = ρ₀ · (r − r₁)² / (r · (r₂ − r₁) / r₂) · (r₂ / (r₂ − r₁))²
    /// ```
    ///
    /// Simplified closed-form (Cummer & Schurig, 2007):
    ///
    /// ```text
    /// ρ_r = ρ₀ · (r₂/(r₂−r₁))² · ((r−r₁)/r)²
    /// ```
    pub fn radial_density(&self, r: f64) -> f64 {
        let r1 = self.r_inner;
        let r2 = self.r_outer;
        if r <= r1 || r >= r2 || (r2 - r1).abs() < 1e-15 {
            return 0.0;
        }
        let factor = r2 / (r2 - r1);
        self.background_density * factor * factor * ((r - r1) / r).powi(2)
    }
    /// Effective angular density at radius `r`.
    ///
    /// ```text
    /// ρ_θ = ρ₀ · (r₂ / (r₂ − r₁))²
    /// ```
    pub fn angular_density(&self, _r: f64) -> f64 {
        let r1 = self.r_inner;
        let r2 = self.r_outer;
        if (r2 - r1).abs() < 1e-15 {
            return 0.0;
        }
        let factor = r2 / (r2 - r1);
        self.background_density * factor * factor
    }
    /// Effective bulk modulus at radius `r`.
    ///
    /// ```text
    /// κ(r) = κ₀ · (r₂/(r₂−r₁))² · ((r−r₁)/r)²
    /// ```
    pub fn effective_modulus(&self, r: f64) -> f64 {
        let r1 = self.r_inner;
        let r2 = self.r_outer;
        if r <= r1 || r >= r2 || (r2 - r1).abs() < 1e-15 {
            return 0.0;
        }
        let factor = r2 / (r2 - r1);
        self.background_modulus * factor * factor * ((r - r1) / r).powi(2)
    }
    /// Generate layer-by-layer material properties for the discretized shell.
    ///
    /// Returns a `Vec` of `(r_center, rho_r, rho_theta, kappa)` for each layer.
    pub fn layer_properties(&self) -> Vec<(f64, f64, f64, f64)> {
        if self.n_layers == 0 {
            return Vec::new();
        }
        let r1 = self.r_inner;
        let r2 = self.r_outer;
        let dr = (r2 - r1) / self.n_layers as f64;
        (0..self.n_layers)
            .map(|i| {
                let r = r1 + (i as f64 + 0.5) * dr;
                let rho_r = self.radial_density(r);
                let rho_t = self.angular_density(r);
                let kappa = self.effective_modulus(r);
                (r, rho_r, rho_t, kappa)
            })
            .collect()
    }
}
