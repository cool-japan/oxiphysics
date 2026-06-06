// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Additive manufacturing FEM module.
//!
//! Simulates layer-by-layer deposition processes such as SLM, EBM, and DED:
//!
//! - [`LayerActivation`]: element birth/death for incremental layer deposition
//! - [`ThermalHistoryFem`]: Goldak double-ellipsoid heat source and melt-pool tracking
//! - [`ResidualStressFem`]: thermal stress accumulation, plasticity, distortion prediction
//! - [`PowderBedFusion`]: SLM/EBM scan strategy, energy density, melting threshold
//! - [`DirectedEnergyDeposition`]: DED clad layer bonding and dilution ratio
//! - [`SupportStructure`]: overhang detection, support topology and removal stress
//! - [`MicrostructureEvolution`]: Hillert grain growth, JMAK phase transformation, texture
//! - [`PorosityModel`]: lack-of-fusion, keyhole, and gas porosity; effective properties
//! - [`DistortionCompensation`]: inverse distortion pre-deformation and scaling factors
//! - [`BuildProcessOptimization`]: scan speed/power optimisation for residual-stress minimisation

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Reference temperature for residual stress (K) — ambient.
const T_AMBIENT: f64 = 298.15;

// ---------------------------------------------------------------------------
// Vector helpers
// ---------------------------------------------------------------------------

/// Dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean norm of a 3-vector.
#[inline]
pub fn norm3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

/// Normalise a 3-vector (returns zero vector if norm < eps).
#[inline]
pub fn normalise3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1.0e-300 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}

/// Element-wise addition of two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Element-wise subtraction of two 3-vectors.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
pub fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

// ---------------------------------------------------------------------------
// LayerActivation
// ---------------------------------------------------------------------------

/// Element activation state used in layer-by-layer deposition simulation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ElementState {
    /// Element is not yet deposited (dormant — contributes no stiffness).
    Inactive,
    /// Element is being deposited in the current layer (just activated).
    JustActivated,
    /// Element has been deposited and has fully solidified.
    Active,
}

/// Manages element birth/death scheduling for layer-by-layer AM simulation.
///
/// Each element is assigned a layer index. Elements are activated when the
/// current simulation layer reaches their layer index.
#[derive(Debug, Clone)]
pub struct LayerActivation {
    /// Total number of elements in the mesh.
    pub n_elements: usize,
    /// Layer index for each element (0-based).
    pub element_layer: Vec<usize>,
    /// Current activation state for each element.
    pub element_state: Vec<ElementState>,
    /// Current build layer index.
    pub current_layer: usize,
    /// Total number of build layers.
    pub n_layers: usize,
    /// Layer thickness (m).
    pub layer_thickness: f64,
    /// Deposition time per layer (s).
    pub time_per_layer: f64,
}

impl LayerActivation {
    /// Create a new `LayerActivation` manager.
    ///
    /// # Arguments
    /// * `n_elements`      – total mesh element count
    /// * `element_layer`   – layer index per element
    /// * `n_layers`        – total layer count
    /// * `layer_thickness` – physical layer thickness (m)
    /// * `time_per_layer`  – deposition time per layer (s)
    pub fn new(
        n_elements: usize,
        element_layer: Vec<usize>,
        n_layers: usize,
        layer_thickness: f64,
        time_per_layer: f64,
    ) -> Self {
        assert_eq!(element_layer.len(), n_elements);
        let element_state = vec![ElementState::Inactive; n_elements];
        Self {
            n_elements,
            element_layer,
            element_state,
            current_layer: 0,
            n_layers,
            layer_thickness,
            time_per_layer,
        }
    }

    /// Advance to the next build layer, activating all elements in that layer.
    ///
    /// Returns the indices of newly activated elements.
    pub fn advance_layer(&mut self) -> Vec<usize> {
        let mut newly_activated = Vec::new();
        // Mark previously just-activated elements as fully active.
        for s in self.element_state.iter_mut() {
            if *s == ElementState::JustActivated {
                *s = ElementState::Active;
            }
        }
        // Activate elements belonging to the new layer.
        for (idx, &layer) in self.element_layer.iter().enumerate() {
            if layer == self.current_layer && self.element_state[idx] == ElementState::Inactive {
                self.element_state[idx] = ElementState::JustActivated;
                newly_activated.push(idx);
            }
        }
        if self.current_layer < self.n_layers.saturating_sub(1) {
            self.current_layer += 1;
        }
        newly_activated
    }

    /// Return the fraction of active elements (0.0 – 1.0).
    pub fn active_fraction(&self) -> f64 {
        let active = self
            .element_state
            .iter()
            .filter(|&&s| s == ElementState::Active || s == ElementState::JustActivated)
            .count();
        active as f64 / self.n_elements as f64
    }

    /// Return the current build height (m).
    pub fn build_height(&self) -> f64 {
        self.current_layer as f64 * self.layer_thickness
    }

    /// Return the elapsed simulation time (s).
    pub fn elapsed_time(&self) -> f64 {
        self.current_layer as f64 * self.time_per_layer
    }

    /// Check whether a given element index is currently active.
    pub fn is_active(&self, element_idx: usize) -> bool {
        matches!(
            self.element_state[element_idx],
            ElementState::Active | ElementState::JustActivated
        )
    }
}

// ---------------------------------------------------------------------------
// ThermalHistoryFem
// ---------------------------------------------------------------------------

/// Parameters for the Goldak double-ellipsoid moving heat source.
///
/// The heat flux distribution is:
///
/// Q(x,y,z,t) = (6√3 * f * P) / (π√π * a*b*c) * exp(-3((x-vt)²/a² + y²/b² + z²/c²))
///
/// where front (f=1) and rear (r=2) ellipsoid halves use lengths `c_front` and `c_rear`.
#[derive(Debug, Clone)]
pub struct GoldakSource {
    /// Laser/beam power (W).
    pub power: f64,
    /// Absorption efficiency (0–1).
    pub efficiency: f64,
    /// Ellipsoid semi-axis in x (scan direction, m).
    pub a: f64,
    /// Ellipsoid semi-axis in y (transverse, m).
    pub b: f64,
    /// Front ellipsoid semi-axis in z (depth, m).
    pub c_front: f64,
    /// Rear ellipsoid semi-axis in z (depth, m).
    pub c_rear: f64,
    /// Fraction of heat deposited in front ellipsoid.
    pub f_front: f64,
    /// Fraction of heat deposited in rear ellipsoid.
    pub f_rear: f64,
}

impl GoldakSource {
    /// Create a Goldak double-ellipsoid source with symmetric fractions (f_f=0.6, f_r=1.4).
    pub fn new(power: f64, efficiency: f64, a: f64, b: f64, c_front: f64, c_rear: f64) -> Self {
        Self {
            power,
            efficiency,
            a,
            b,
            c_front,
            c_rear,
            f_front: 0.6,
            f_rear: 1.4,
        }
    }

    /// Evaluate volumetric heat flux (W m⁻³) at point `p` when heat source is at `source_pos`.
    ///
    /// Uses the front or rear ellipsoid based on the sign of (p\[0\] - source_pos\[0\]).
    pub fn heat_flux(&self, p: [f64; 3], source_pos: [f64; 3]) -> f64 {
        let effective_power = self.power * self.efficiency;
        let dx = p[0] - source_pos[0];
        let dy = p[1] - source_pos[1];
        let dz = p[2] - source_pos[2];

        let (c, f) = if dx >= 0.0 {
            (self.c_front, self.f_front)
        } else {
            (self.c_rear, self.f_rear)
        };

        let prefactor =
            6.0 * 3.0_f64.sqrt() * f * effective_power / (PI * PI.sqrt() * self.a * self.b * c);
        let exponent =
            -3.0 * (dx * dx / (self.a * self.a) + dy * dy / (self.b * self.b) + dz * dz / (c * c));
        prefactor * exponent.exp()
    }
}

/// Thermal history FEM state for tracking melt pool and temperature evolution.
#[derive(Debug, Clone)]
pub struct ThermalHistoryFem {
    /// Number of nodes.
    pub n_nodes: usize,
    /// Node positions (m).
    pub positions: Vec<[f64; 3]>,
    /// Current nodal temperatures (K).
    pub temperatures: Vec<f64>,
    /// Peak temperature recorded at each node (K).
    pub peak_temperatures: Vec<f64>,
    /// Goldak heat source parameters.
    pub source: GoldakSource,
    /// Liquidus temperature (K).
    pub t_liquidus: f64,
    /// Solidus temperature (K).
    pub t_solidus: f64,
    /// Thermal conductivity (W m⁻¹ K⁻¹).
    pub conductivity: f64,
    /// Specific heat capacity (J kg⁻¹ K⁻¹).
    pub specific_heat: f64,
    /// Density (kg m⁻³).
    pub density: f64,
}

impl ThermalHistoryFem {
    /// Construct a new `ThermalHistoryFem` instance.
    pub fn new(
        positions: Vec<[f64; 3]>,
        source: GoldakSource,
        t_liquidus: f64,
        t_solidus: f64,
        conductivity: f64,
        specific_heat: f64,
        density: f64,
    ) -> Self {
        let n = positions.len();
        Self {
            n_nodes: n,
            positions,
            temperatures: vec![T_AMBIENT; n],
            peak_temperatures: vec![T_AMBIENT; n],
            source,
            t_liquidus,
            t_solidus,
            conductivity,
            specific_heat,
            density,
        }
    }

    /// Apply a single explicit thermal step with a moving heat source.
    ///
    /// Uses a simple lumped-capacity approach: dT/dt = Q/(ρ·c_p) – h·(T - T_amb)/ρ/c_p/V_node
    /// where h is the convective coefficient and V_node is the lumped nodal volume.
    pub fn step(&mut self, source_pos: [f64; 3], dt: f64, h_conv: f64, node_volume: f64) {
        let rho_cp = self.density * self.specific_heat;
        for i in 0..self.n_nodes {
            let q = self.source.heat_flux(self.positions[i], source_pos);
            let t = self.temperatures[i];
            let dt_node = (q / rho_cp - h_conv * (t - T_AMBIENT) / rho_cp / node_volume) * dt;
            self.temperatures[i] = (t + dt_node).max(T_AMBIENT);
            if self.temperatures[i] > self.peak_temperatures[i] {
                self.peak_temperatures[i] = self.temperatures[i];
            }
        }
    }

    /// Return indices of nodes currently in the melt pool (T ≥ T_liquidus).
    pub fn melt_pool_nodes(&self) -> Vec<usize> {
        self.temperatures
            .iter()
            .enumerate()
            .filter(|&(_, &t)| t >= self.t_liquidus)
            .map(|(i, _)| i)
            .collect()
    }

    /// Return fraction of nodes above solidus (partially melted + melt pool).
    pub fn mushy_zone_fraction(&self) -> f64 {
        let mushy = self
            .temperatures
            .iter()
            .filter(|&&t| t >= self.t_solidus)
            .count();
        mushy as f64 / self.n_nodes as f64
    }

    /// Thermal diffusivity α = k / (ρ·c_p) (m² s⁻¹).
    pub fn thermal_diffusivity(&self) -> f64 {
        self.conductivity / (self.density * self.specific_heat)
    }

    /// Compute temperature-dependent conductivity using linear interpolation.
    ///
    /// k(T) = k_ref * (1 + β*(T - T_ref))
    pub fn temperature_dependent_conductivity(&self, temperature: f64, beta: f64) -> f64 {
        let t_ref = 298.15;
        self.conductivity * (1.0 + beta * (temperature - t_ref))
    }
}

// ---------------------------------------------------------------------------
// ResidualStressFem
// ---------------------------------------------------------------------------

/// Accumulated thermal stress state at a material point.
#[derive(Debug, Clone)]
pub struct StressState {
    /// Stress tensor in Voigt notation \[σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz\] (Pa).
    pub stress: [f64; 6],
    /// Plastic strain in Voigt notation.
    pub plastic_strain: [f64; 6],
    /// Accumulated equivalent plastic strain.
    pub eq_plastic_strain: f64,
}

impl StressState {
    /// Create a zero stress state.
    pub fn zero() -> Self {
        Self {
            stress: [0.0; 6],
            plastic_strain: [0.0; 6],
            eq_plastic_strain: 0.0,
        }
    }

    /// Von Mises equivalent stress (Pa).
    pub fn von_mises(&self) -> f64 {
        let s = &self.stress;
        let dev_xx = s[0] - (s[0] + s[1] + s[2]) / 3.0;
        let dev_yy = s[1] - (s[0] + s[1] + s[2]) / 3.0;
        let dev_zz = s[2] - (s[0] + s[1] + s[2]) / 3.0;
        (0.5 * (dev_xx * dev_xx
            + dev_yy * dev_yy
            + dev_zz * dev_zz
            + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5])))
            .sqrt()
            * 3.0_f64.sqrt()
    }

    /// Hydrostatic pressure (Pa) — positive in tension.
    pub fn pressure(&self) -> f64 {
        (self.stress[0] + self.stress[1] + self.stress[2]) / 3.0
    }
}

/// FEM model for residual stress accumulation during additive manufacturing.
#[derive(Debug, Clone)]
pub struct ResidualStressFem {
    /// Number of integration points.
    pub n_points: usize,
    /// Stress states at each integration point.
    pub states: Vec<StressState>,
    /// Young's modulus (Pa).
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Yield stress at reference temperature (Pa).
    pub yield_stress_ref: f64,
    /// Thermal expansion coefficient (K⁻¹).
    pub thermal_expansion: f64,
    /// Isotropic hardening modulus (Pa).
    pub hardening_modulus: f64,
}

impl ResidualStressFem {
    /// Create a new residual stress FEM with given material properties.
    pub fn new(
        n_points: usize,
        young_modulus: f64,
        poisson_ratio: f64,
        yield_stress_ref: f64,
        thermal_expansion: f64,
        hardening_modulus: f64,
    ) -> Self {
        Self {
            n_points,
            states: vec![StressState::zero(); n_points],
            young_modulus,
            poisson_ratio,
            yield_stress_ref,
            thermal_expansion,
            hardening_modulus,
        }
    }

    /// Lame's first parameter λ (Pa).
    pub fn lame_lambda(&self) -> f64 {
        self.young_modulus * self.poisson_ratio
            / ((1.0 + self.poisson_ratio) * (1.0 - 2.0 * self.poisson_ratio))
    }

    /// Shear modulus G (Pa).
    pub fn shear_modulus(&self) -> f64 {
        self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Yield stress at temperature T using linear thermal softening.
    ///
    /// σ_y(T) = σ_y0 * max(0, 1 - (T - T_ref) / (T_melt - T_ref))
    pub fn yield_stress_at_temperature(&self, temperature: f64, t_melt: f64) -> f64 {
        let t_ref = T_AMBIENT;
        let factor = 1.0 - (temperature - t_ref) / (t_melt - t_ref);
        self.yield_stress_ref * factor.max(0.0)
    }

    /// Perform a radial-return plasticity correction at point `idx`.
    ///
    /// Given a thermal strain increment `d_eps_th`, update the stress state
    /// and return the plastic strain increment magnitude.
    pub fn radial_return(
        &mut self,
        point_idx: usize,
        d_eps_th: f64,
        temperature: f64,
        t_melt: f64,
    ) -> f64 {
        let g = self.shear_modulus();
        let lam = self.lame_lambda();
        let sigma_y = self.yield_stress_at_temperature(temperature, t_melt)
            + self.hardening_modulus * self.states[point_idx].eq_plastic_strain;

        // Isotropic thermal eigenstrain: Δε_th * I
        let d_sig_th = (3.0 * lam + 2.0 * g) * self.thermal_expansion * d_eps_th;

        // Update hydrostatic stress
        self.states[point_idx].stress[0] -= d_sig_th;
        self.states[point_idx].stress[1] -= d_sig_th;
        self.states[point_idx].stress[2] -= d_sig_th;

        // Trial von Mises
        let vm_trial = self.states[point_idx].von_mises();

        if vm_trial <= sigma_y {
            return 0.0;
        }

        // Plastic correction factor
        let d_gamma = (vm_trial - sigma_y) / (3.0 * g + self.hardening_modulus);
        let scale = 1.0 - 3.0 * g * d_gamma / vm_trial;

        // Scale deviatoric stress
        let p = self.states[point_idx].pressure();
        for k in 0..6 {
            let dev = if k < 3 {
                self.states[point_idx].stress[k] - p
            } else {
                self.states[point_idx].stress[k]
            };
            let corrected = dev * scale;
            if k < 3 {
                self.states[point_idx].stress[k] = corrected + p;
            } else {
                self.states[point_idx].stress[k] = corrected;
            }
        }
        self.states[point_idx].eq_plastic_strain += d_gamma;
        d_gamma
    }

    /// Maximum von Mises stress across all integration points (Pa).
    pub fn max_von_mises(&self) -> f64 {
        self.states
            .iter()
            .map(|s| s.von_mises())
            .fold(0.0_f64, f64::max)
    }

    /// Mean residual stress (hydrostatic) across all points (Pa).
    pub fn mean_residual_pressure(&self) -> f64 {
        let sum: f64 = self.states.iter().map(|s| s.pressure()).sum();
        sum / self.n_points as f64
    }
}

// ---------------------------------------------------------------------------
// PowderBedFusion
// ---------------------------------------------------------------------------

/// Scan strategy type for powder bed fusion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ScanStrategy {
    /// Unidirectional scan (all tracks in +x direction).
    Unidirectional,
    /// Bidirectional (alternating +x/−x) scan.
    Bidirectional,
    /// Island/checkerboard scan pattern.
    Island,
    /// Spiral scan from outside to centre.
    Spiral,
    /// Alternating 90° rotation between layers.
    Rotating90,
}

/// Powder bed fusion (SLM / EBM) process simulator.
#[derive(Debug, Clone)]
pub struct PowderBedFusion {
    /// Laser/beam power (W).
    pub power: f64,
    /// Scan speed (m s⁻¹).
    pub scan_speed: f64,
    /// Hatch spacing (m).
    pub hatch_spacing: f64,
    /// Layer thickness (m).
    pub layer_thickness: f64,
    /// Powder absorptivity (0–1).
    pub absorptivity: f64,
    /// Powder bed bulk density (kg m⁻³).
    pub powder_density: f64,
    /// Melting point of powder material (K).
    pub melting_point: f64,
    /// Chosen scan strategy.
    pub scan_strategy: ScanStrategy,
}

impl PowderBedFusion {
    /// Create a new `PowderBedFusion` process.
    pub fn new(
        power: f64,
        scan_speed: f64,
        hatch_spacing: f64,
        layer_thickness: f64,
        absorptivity: f64,
        powder_density: f64,
        melting_point: f64,
        scan_strategy: ScanStrategy,
    ) -> Self {
        Self {
            power,
            scan_speed,
            hatch_spacing,
            layer_thickness,
            absorptivity,
            powder_density,
            melting_point,
            scan_strategy,
        }
    }

    /// Volumetric energy density E_v (J m⁻³).
    ///
    /// E_v = P / (v · h · t)
    pub fn energy_density(&self) -> f64 {
        self.power / (self.scan_speed * self.hatch_spacing * self.layer_thickness)
    }

    /// Linear energy density E_l (J m⁻¹).
    pub fn linear_energy_density(&self) -> f64 {
        self.power / self.scan_speed
    }

    /// Check whether the volumetric energy density exceeds the melting threshold.
    ///
    /// Uses an approximate threshold:
    /// E_th = ρ · c_p · (T_melt - T_amb) · π / (4 · η)
    pub fn is_above_melting_threshold(&self, specific_heat: f64, efficiency: f64) -> bool {
        let e_th = self.powder_density * specific_heat * (self.melting_point - T_AMBIENT) * PI
            / (4.0 * efficiency);
        self.energy_density() > e_th
    }

    /// Estimate melt pool depth using the Rosenthal point-source solution (m).
    ///
    /// r_melt ≈ sqrt(2 * η * P / (π * e * k * v * (T_melt - T0)))
    pub fn melt_pool_depth(&self, conductivity: f64, efficiency: f64) -> f64 {
        let numerator = 2.0 * efficiency * self.power;
        let denominator = PI
            * std::f64::consts::E
            * conductivity
            * self.scan_speed
            * (self.melting_point - T_AMBIENT);
        (numerator / denominator).sqrt()
    }

    /// Number of scan tracks per layer for a given part width (m).
    pub fn tracks_per_layer(&self, part_width: f64) -> usize {
        (part_width / self.hatch_spacing).ceil() as usize
    }

    /// Total scan length per layer for a given part width and depth (m).
    pub fn total_scan_length_per_layer(&self, width: f64, depth: f64) -> f64 {
        let n_tracks = self.tracks_per_layer(width);
        n_tracks as f64 * depth
    }
}

// ---------------------------------------------------------------------------
// DirectedEnergyDeposition
// ---------------------------------------------------------------------------

/// Directed energy deposition (DED / LENS / LMD) process model.
#[derive(Debug, Clone)]
pub struct DirectedEnergyDeposition {
    /// Laser power (W).
    pub power: f64,
    /// Travel speed (m s⁻¹).
    pub travel_speed: f64,
    /// Powder feed rate (kg s⁻¹).
    pub powder_feed_rate: f64,
    /// Powder catchment efficiency (0–1).
    pub catchment_efficiency: f64,
    /// Substrate absorptivity (0–1).
    pub substrate_absorptivity: f64,
    /// Clad layer height (m).
    pub clad_height: f64,
    /// Clad layer width (m).
    pub clad_width: f64,
}

impl DirectedEnergyDeposition {
    /// Create a new DED process model.
    pub fn new(
        power: f64,
        travel_speed: f64,
        powder_feed_rate: f64,
        catchment_efficiency: f64,
        substrate_absorptivity: f64,
        clad_height: f64,
        clad_width: f64,
    ) -> Self {
        Self {
            power,
            travel_speed,
            powder_feed_rate,
            catchment_efficiency,
            substrate_absorptivity,
            clad_height,
            clad_width,
        }
    }

    /// Effective deposited power on the substrate (W).
    pub fn effective_power(&self) -> f64 {
        self.power * self.substrate_absorptivity
    }

    /// Deposition rate (kg s⁻¹) accounting for catchment efficiency.
    pub fn deposition_rate(&self) -> f64 {
        self.powder_feed_rate * self.catchment_efficiency
    }

    /// Dilution ratio D = (melted substrate volume) / (total melt pool volume).
    ///
    /// Uses an empirical model: D ≈ η·P / (v·A_clad·ρ·H_f + η·P)
    /// where A_clad = height × width and H_f is the fusion enthalpy.
    pub fn dilution_ratio(&self, density: f64, fusion_enthalpy: f64) -> f64 {
        let ep = self.effective_power();
        let a_clad = self.clad_height * self.clad_width;
        let denom = ep + self.travel_speed * a_clad * density * fusion_enthalpy;
        ep / denom
    }

    /// Clad cross-sectional area (m²).
    pub fn clad_area(&self) -> f64 {
        // Approximate as semi-ellipse
        PI * self.clad_height * self.clad_width / 4.0
    }

    /// Bonding strength factor (0–1) based on dilution ratio.
    ///
    /// Optimal bonding is achieved at dilution D ≈ 0.1–0.3; outside this range strength degrades.
    pub fn bonding_factor(&self, density: f64, fusion_enthalpy: f64) -> f64 {
        let d = self.dilution_ratio(density, fusion_enthalpy);
        if d < 0.05 {
            d / 0.05
        } else if d <= 0.35 {
            1.0
        } else {
            (1.0 - (d - 0.35) / 0.65).max(0.0)
        }
    }
}

// ---------------------------------------------------------------------------
// SupportStructure
// ---------------------------------------------------------------------------

/// Support structure topology model for overhang regions in AM.
#[derive(Debug, Clone)]
pub struct SupportStructure {
    /// Critical overhang angle below which support is required (radians from horizontal).
    pub critical_angle: f64,
    /// Support material fraction (porosity complement, 0–1).
    pub material_fraction: f64,
    /// Support element stiffness scaling factor (0–1).
    pub stiffness_scale: f64,
    /// Number of support elements.
    pub n_support_elements: usize,
    /// Stress in support elements (Pa).
    pub support_stress: Vec<f64>,
}

impl SupportStructure {
    /// Create a new `SupportStructure` model.
    pub fn new(
        critical_angle_deg: f64,
        material_fraction: f64,
        stiffness_scale: f64,
        n_support_elements: usize,
    ) -> Self {
        Self {
            critical_angle: critical_angle_deg.to_radians(),
            material_fraction,
            stiffness_scale,
            n_support_elements,
            support_stress: vec![0.0; n_support_elements],
        }
    }

    /// Determine whether a surface facet with normal `facet_normal` requires support.
    ///
    /// A facet needs support if the angle between its normal and the build direction
    /// (assumed to be +z) is less than the critical angle.
    pub fn needs_support(&self, facet_normal: [f64; 3]) -> bool {
        let build_dir = [0.0_f64, 0.0, 1.0];
        let cos_angle = dot3(facet_normal, build_dir).abs();
        let angle = cos_angle.acos();
        angle < self.critical_angle
    }

    /// Compute the support volume fraction needed for a given overhang angle (radians).
    ///
    /// Uses a linear model: V_sup = V_0 * max(0, 1 - θ/θ_c)
    pub fn required_volume_fraction(&self, overhang_angle_rad: f64) -> f64 {
        let ratio = overhang_angle_rad / self.critical_angle;
        (1.0 - ratio).max(0.0) * self.material_fraction
    }

    /// Estimate the stress in support element `idx` after part removal (Pa).
    pub fn removal_stress(&self, element_idx: usize, part_weight: f64, support_area: f64) -> f64 {
        let base_stress = part_weight / support_area;
        self.support_stress[element_idx] + base_stress * self.stiffness_scale
    }

    /// Maximum support stress across all elements (Pa).
    pub fn max_support_stress(&self) -> f64 {
        self.support_stress.iter().cloned().fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// MicrostructureEvolution
// ---------------------------------------------------------------------------

/// JMAK (Johnson-Mehl-Avrami-Kolmogorov) phase transformation model.
#[derive(Debug, Clone)]
pub struct JmakTransformation {
    /// Avrami exponent n.
    pub n_avrami: f64,
    /// Rate constant k (s⁻ⁿ).
    pub k_rate: f64,
}

impl JmakTransformation {
    /// Create a JMAK model with given exponent and rate constant.
    pub fn new(n_avrami: f64, k_rate: f64) -> Self {
        Self { n_avrami, k_rate }
    }

    /// Transformed fraction X at time t (s).
    ///
    /// X(t) = 1 − exp(−k · tⁿ)
    pub fn transformed_fraction(&self, time: f64) -> f64 {
        let exponent = -self.k_rate * time.powf(self.n_avrami);
        1.0 - exponent.exp()
    }

    /// Incubation time (s) for a given initial transformed fraction X_0.
    pub fn incubation_time(&self, x0: f64) -> f64 {
        let ln_term = -(1.0 - x0).ln();
        (ln_term / self.k_rate).powf(1.0 / self.n_avrami)
    }
}

/// Hillert grain growth model: dR/dt = M * (1/R_c - 1/R)
#[derive(Debug, Clone)]
pub struct HillertGrainGrowth {
    /// Grain boundary mobility pre-factor (m² s⁻¹).
    pub mobility_prefactor: f64,
    /// Grain boundary energy (J m⁻²).
    pub boundary_energy: f64,
    /// Activation energy for grain boundary migration (J mol⁻¹).
    pub activation_energy: f64,
    /// Universal gas constant (J mol⁻¹ K⁻¹).
    pub gas_constant: f64,
}

impl HillertGrainGrowth {
    /// Create a new grain growth model.
    pub fn new(mobility_prefactor: f64, boundary_energy: f64, activation_energy: f64) -> Self {
        Self {
            mobility_prefactor,
            boundary_energy,
            activation_energy,
            gas_constant: 8.314,
        }
    }

    /// Grain boundary mobility at temperature T (m² s⁻¹).
    pub fn mobility(&self, temperature: f64) -> f64 {
        self.mobility_prefactor
            * (-self.activation_energy / (self.gas_constant * temperature)).exp()
    }

    /// Mean grain radius after time dt at temperature T (m).
    ///
    /// Normal grain growth: R² - R0² = 2·M·γ·t
    pub fn grain_radius_after(&self, r0: f64, dt: f64, temperature: f64) -> f64 {
        let m = self.mobility(temperature);
        let r_sq = r0 * r0 + 2.0 * m * self.boundary_energy * dt;
        r_sq.sqrt()
    }

    /// Critical grain size (m) — grains above this size grow, below shrink.
    ///
    /// R_c = 2 * mean_radius (Hillert mean-field approximation)
    pub fn critical_radius(&self, mean_radius: f64) -> f64 {
        2.0 * mean_radius
    }
}

/// Microstructure evolution model combining grain growth and phase transformation.
#[derive(Debug, Clone)]
pub struct MicrostructureEvolution {
    /// Current mean grain radius (m).
    pub grain_radius: f64,
    /// Current transformed (solid) fraction.
    pub transformed_fraction: f64,
    /// Grain growth model.
    pub grain_growth: HillertGrainGrowth,
    /// Phase transformation model.
    pub jmak: JmakTransformation,
    /// Texture index (1 = random, higher = stronger texture).
    pub texture_index: f64,
}

impl MicrostructureEvolution {
    /// Create a new microstructure evolution model.
    pub fn new(
        initial_grain_radius: f64,
        grain_growth: HillertGrainGrowth,
        jmak: JmakTransformation,
    ) -> Self {
        Self {
            grain_radius: initial_grain_radius,
            transformed_fraction: 0.0,
            grain_growth,
            jmak,
            texture_index: 1.0,
        }
    }

    /// Advance the microstructure by time step `dt` at `temperature`.
    pub fn step(&mut self, dt: f64, temperature: f64, cumulative_time: f64) {
        // Update grain radius
        self.grain_radius =
            self.grain_growth
                .grain_radius_after(self.grain_radius, dt, temperature);
        // Update transformed fraction
        self.transformed_fraction = self.jmak.transformed_fraction(cumulative_time + dt);
        // Texture index increases with solidification direction preference
        self.texture_index = 1.0 + self.transformed_fraction * 2.5;
    }

    /// Hall-Petch yield strength contribution (Pa).
    ///
    /// σ_HP = σ_0 + k_HP / sqrt(d)  where d = 2 * R
    pub fn hall_petch_strength(&self, sigma0: f64, k_hp: f64) -> f64 {
        let diameter = 2.0 * self.grain_radius;
        sigma0 + k_hp / diameter.sqrt()
    }
}

// ---------------------------------------------------------------------------
// PorosityModel
// ---------------------------------------------------------------------------

/// Type of porosity defect in AM parts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PorosityType {
    /// Lack-of-fusion porosity from insufficient energy input.
    LackOfFusion,
    /// Keyhole porosity from excessive energy causing deep keyhole collapse.
    Keyhole,
    /// Gas porosity from trapped gas (from powder or atmosphere).
    Gas,
}

/// Porosity model tracking defect fraction and effective mechanical properties.
#[derive(Debug, Clone)]
pub struct PorosityModel {
    /// Total porosity fraction (0–1).
    pub porosity: f64,
    /// Lack-of-fusion porosity fraction.
    pub lof_fraction: f64,
    /// Keyhole porosity fraction.
    pub keyhole_fraction: f64,
    /// Gas porosity fraction.
    pub gas_fraction: f64,
    /// Dense material Young's modulus (Pa).
    pub dense_young_modulus: f64,
    /// Dense material yield stress (Pa).
    pub dense_yield_stress: f64,
}

impl PorosityModel {
    /// Create a porosity model for a given energy density and process parameters.
    pub fn new(dense_young_modulus: f64, dense_yield_stress: f64) -> Self {
        Self {
            porosity: 0.0,
            lof_fraction: 0.0,
            keyhole_fraction: 0.0,
            gas_fraction: 0.0,
            dense_young_modulus,
            dense_yield_stress,
        }
    }

    /// Update porosity fractions from process parameters.
    ///
    /// Uses empirical thresholds:
    /// - LoF dominant below E_v < E_threshold_low
    /// - Keyhole dominant above E_v > E_threshold_high
    /// - Gas porosity is approximately constant
    pub fn update_from_energy_density(
        &mut self,
        energy_density: f64,
        e_low: f64,
        e_high: f64,
        gas_base: f64,
    ) {
        self.gas_fraction = gas_base;
        if energy_density < e_low {
            // Lack-of-fusion regime
            let deficit = (e_low - energy_density) / e_low;
            self.lof_fraction = 0.15 * deficit * deficit;
            self.keyhole_fraction = 0.0;
        } else if energy_density > e_high {
            // Keyhole regime
            let excess = (energy_density - e_high) / e_high;
            self.keyhole_fraction = 0.08 * excess;
            self.lof_fraction = 0.0;
        } else {
            self.lof_fraction = 0.0;
            self.keyhole_fraction = 0.0;
        }
        self.porosity = self.lof_fraction + self.keyhole_fraction + self.gas_fraction;
    }

    /// Effective Young's modulus using the Ramakrishnan-Arunachalam model.
    ///
    /// E_eff = E_0 * (1 - p)^b where b ≈ 2.0 for spherical pores.
    pub fn effective_young_modulus(&self) -> f64 {
        self.dense_young_modulus * (1.0 - self.porosity).powf(2.0)
    }

    /// Effective yield stress using the empirical power-law model.
    ///
    /// σ_y_eff = σ_y0 * (1 - p)^(3/2)
    pub fn effective_yield_stress(&self) -> f64 {
        self.dense_yield_stress * (1.0 - self.porosity).powf(1.5)
    }

    /// Relative density (1 - porosity).
    pub fn relative_density(&self) -> f64 {
        1.0 - self.porosity
    }
}

// ---------------------------------------------------------------------------
// DistortionCompensation
// ---------------------------------------------------------------------------

/// Distortion compensation model for AM pre-deformation.
#[derive(Debug, Clone)]
pub struct DistortionCompensation {
    /// Number of nodes in the part mesh.
    pub n_nodes: usize,
    /// Nominal (design) node positions (m).
    pub nominal_positions: Vec<[f64; 3]>,
    /// Predicted distorted positions after AM (m).
    pub distorted_positions: Vec<[f64; 3]>,
    /// Pre-compensated positions (feed-forward correction) (m).
    pub compensated_positions: Vec<[f64; 3]>,
    /// Global scaling factors \[sx, sy, sz\] applied before compensation.
    pub scaling_factors: [f64; 3],
    /// Compensation fraction (0 = no compensation, 1 = full inverse).
    pub compensation_fraction: f64,
}

impl DistortionCompensation {
    /// Create a new distortion compensation model.
    pub fn new(nominal_positions: Vec<[f64; 3]>, compensation_fraction: f64) -> Self {
        let n = nominal_positions.len();
        let distorted_positions = nominal_positions.clone();
        let compensated_positions = nominal_positions.clone();
        Self {
            n_nodes: n,
            nominal_positions,
            distorted_positions,
            compensated_positions,
            scaling_factors: [1.0; 3],
            compensation_fraction,
        }
    }

    /// Set predicted distorted positions from simulation output.
    pub fn set_distorted_positions(&mut self, positions: Vec<[f64; 3]>) {
        assert_eq!(positions.len(), self.n_nodes);
        self.distorted_positions = positions;
        self.update_compensated();
    }

    /// Recompute pre-compensated positions from current distorted positions.
    fn update_compensated(&mut self) {
        for i in 0..self.n_nodes {
            let nom = self.nominal_positions[i];
            let dist = self.distorted_positions[i];
            let delta = sub3(dist, nom);
            let correction = scale3(delta, -self.compensation_fraction);
            self.compensated_positions[i] = add3(nom, correction);
        }
    }

    /// Apply scaling factors to the compensated positions.
    pub fn apply_scaling(&mut self) {
        let sf = self.scaling_factors;
        for pos in self.compensated_positions.iter_mut() {
            pos[0] *= sf[0];
            pos[1] *= sf[1];
            pos[2] *= sf[2];
        }
    }

    /// Root-mean-square distortion across all nodes (m).
    pub fn rms_distortion(&self) -> f64 {
        let sum_sq: f64 = self
            .distorted_positions
            .iter()
            .zip(self.nominal_positions.iter())
            .map(|(&d, &n)| {
                let delta = sub3(d, n);
                dot3(delta, delta)
            })
            .sum();
        (sum_sq / self.n_nodes as f64).sqrt()
    }

    /// Maximum distortion magnitude across all nodes (m).
    pub fn max_distortion(&self) -> f64 {
        self.distorted_positions
            .iter()
            .zip(self.nominal_positions.iter())
            .map(|(&d, &n)| norm3(sub3(d, n)))
            .fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// BuildProcessOptimization
// ---------------------------------------------------------------------------

/// Cost function weights for build process optimisation.
#[derive(Debug, Clone)]
pub struct OptimizationWeights {
    /// Weight for residual stress minimisation.
    pub stress_weight: f64,
    /// Weight for distortion minimisation.
    pub distortion_weight: f64,
    /// Weight for build time minimisation.
    pub time_weight: f64,
    /// Weight for porosity minimisation.
    pub porosity_weight: f64,
}

impl OptimizationWeights {
    /// Create equal-weight objective.
    pub fn uniform() -> Self {
        Self {
            stress_weight: 1.0,
            distortion_weight: 1.0,
            time_weight: 1.0,
            porosity_weight: 1.0,
        }
    }

    /// Normalise weights so they sum to 1.
    pub fn normalise(&self) -> Self {
        let total =
            self.stress_weight + self.distortion_weight + self.time_weight + self.porosity_weight;
        Self {
            stress_weight: self.stress_weight / total,
            distortion_weight: self.distortion_weight / total,
            time_weight: self.time_weight / total,
            porosity_weight: self.porosity_weight / total,
        }
    }
}

/// Build process parameter set for a single design point.
#[derive(Debug, Clone)]
pub struct ProcessParameters {
    /// Laser power (W).
    pub power: f64,
    /// Scan speed (m s⁻¹).
    pub scan_speed: f64,
    /// Hatch spacing (m).
    pub hatch_spacing: f64,
    /// Layer thickness (m).
    pub layer_thickness: f64,
}

impl ProcessParameters {
    /// Create a new process parameter set.
    pub fn new(power: f64, scan_speed: f64, hatch_spacing: f64, layer_thickness: f64) -> Self {
        Self {
            power,
            scan_speed,
            hatch_spacing,
            layer_thickness,
        }
    }

    /// Volumetric energy density (J m⁻³).
    pub fn energy_density(&self) -> f64 {
        self.power / (self.scan_speed * self.hatch_spacing * self.layer_thickness)
    }
}

/// Build process optimisation driver using a simple grid-search approach.
#[derive(Debug, Clone)]
pub struct BuildProcessOptimization {
    /// Optimisation objective weights.
    pub weights: OptimizationWeights,
    /// Candidate parameter sets.
    pub candidates: Vec<ProcessParameters>,
    /// Objective values for each candidate.
    pub objective_values: Vec<f64>,
}

impl BuildProcessOptimization {
    /// Create optimiser with given weights and candidate set.
    pub fn new(weights: OptimizationWeights, candidates: Vec<ProcessParameters>) -> Self {
        let n = candidates.len();
        Self {
            weights,
            candidates,
            objective_values: vec![f64::INFINITY; n],
        }
    }

    /// Evaluate the weighted objective function for candidate `idx`.
    ///
    /// Uses simplified proxy models:
    /// - Stress ∝ E_v / E_opt
    /// - Distortion ∝ sqrt(E_v)
    /// - Time ∝ 1 / (v * h * t)
    /// - Porosity: zero in a window, else proportional to distance from E_opt
    pub fn evaluate_candidate(&mut self, idx: usize, e_opt: f64) {
        let p = &self.candidates[idx];
        let ev = p.energy_density();
        let w = self.weights.normalise();

        let stress_term = (ev / e_opt - 1.0).abs();
        let distortion_term = (ev / e_opt).sqrt();
        let time_term = 1.0 / (p.scan_speed * p.hatch_spacing * p.layer_thickness) * 1.0e-9;
        let porosity_term = if (ev - e_opt).abs() < 0.2 * e_opt {
            0.0
        } else {
            (ev - e_opt).abs() / e_opt
        };

        self.objective_values[idx] = w.stress_weight * stress_term
            + w.distortion_weight * distortion_term
            + w.time_weight * time_term
            + w.porosity_weight * porosity_term;
    }

    /// Evaluate all candidates and return the index of the best one.
    pub fn optimise(&mut self, e_opt: f64) -> usize {
        for i in 0..self.candidates.len() {
            self.evaluate_candidate(i, e_opt);
        }
        self.objective_values
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0)
    }

    /// Return the best process parameters found after calling `optimise`.
    pub fn best_parameters(&self) -> Option<&ProcessParameters> {
        let best_idx = self
            .objective_values
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)?;
        Some(&self.candidates[best_idx])
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- LayerActivation tests ---

    #[test]
    fn test_layer_activation_initial_state() {
        let layers = vec![0, 0, 1, 1, 2, 2];
        let la = LayerActivation::new(6, layers, 3, 5.0e-5, 10.0);
        assert_eq!(la.n_elements, 6);
        assert!(
            la.element_state
                .iter()
                .all(|&s| s == ElementState::Inactive)
        );
    }

    #[test]
    fn test_layer_activation_advance_layer_zero() {
        let layers = vec![0, 0, 1, 1, 2, 2];
        let mut la = LayerActivation::new(6, layers, 3, 5.0e-5, 10.0);
        let activated = la.advance_layer();
        assert_eq!(activated.len(), 2);
        assert!(activated.contains(&0));
        assert!(activated.contains(&1));
    }

    #[test]
    fn test_layer_activation_active_fraction_increases() {
        let layers = vec![0, 1, 2, 3];
        let mut la = LayerActivation::new(4, layers, 4, 1.0e-4, 5.0);
        la.advance_layer();
        let f1 = la.active_fraction();
        la.advance_layer();
        let f2 = la.active_fraction();
        assert!(f2 >= f1);
    }

    #[test]
    fn test_layer_activation_build_height() {
        let layers = vec![0, 1, 2];
        let mut la = LayerActivation::new(3, layers, 3, 1.0e-4, 5.0);
        la.advance_layer();
        la.advance_layer();
        let h = la.build_height();
        assert!((h - 2.0e-4).abs() < 1.0e-15, "h = {:.6e}", h);
    }

    #[test]
    fn test_layer_activation_elapsed_time() {
        let layers = vec![0, 1];
        let mut la = LayerActivation::new(2, layers, 2, 1.0e-4, 5.0);
        la.advance_layer();
        assert!((la.elapsed_time() - 5.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_layer_activation_is_active() {
        let layers = vec![0, 1];
        let mut la = LayerActivation::new(2, layers, 2, 1.0e-4, 5.0);
        la.advance_layer();
        assert!(la.is_active(0));
        assert!(!la.is_active(1));
    }

    // --- GoldakSource / ThermalHistoryFem tests ---

    #[test]
    fn test_goldak_source_zero_at_distance() {
        let src = GoldakSource::new(200.0, 0.35, 1.5e-4, 1.0e-4, 1.0e-4, 2.0e-4);
        let q = src.heat_flux([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert!(q < 1.0, "Heat flux far from source should be negligible");
    }

    #[test]
    fn test_goldak_source_max_at_source_pos() {
        let src = GoldakSource::new(200.0, 0.35, 1.5e-4, 1.0e-4, 1.0e-4, 2.0e-4);
        let q_center = src.heat_flux([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let q_away = src.heat_flux([5.0e-4, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert!(q_center > q_away);
    }

    #[test]
    fn test_thermal_history_fem_diffusivity() {
        let positions = vec![[0.0; 3]; 4];
        let src = GoldakSource::new(200.0, 0.35, 1.5e-4, 1.0e-4, 1.0e-4, 2.0e-4);
        let thf = ThermalHistoryFem::new(positions, src, 1923.0, 1723.0, 12.0, 620.0, 8000.0);
        let alpha = thf.thermal_diffusivity();
        let expected = 12.0 / (8000.0 * 620.0);
        assert!((alpha - expected).abs() < 1.0e-20);
    }

    #[test]
    fn test_thermal_history_fem_initial_temperatures() {
        let positions = vec![[0.0; 3], [0.01, 0.0, 0.0]];
        let src = GoldakSource::new(200.0, 0.35, 1.5e-4, 1.0e-4, 1.0e-4, 2.0e-4);
        let thf = ThermalHistoryFem::new(positions, src, 1923.0, 1723.0, 12.0, 620.0, 8000.0);
        assert!((thf.temperatures[0] - T_AMBIENT).abs() < 1.0e-10);
    }

    #[test]
    fn test_thermal_history_fem_melt_pool_nodes_empty() {
        let positions = vec![[0.0; 3]; 10];
        let src = GoldakSource::new(1.0, 0.01, 1.5e-4, 1.0e-4, 1.0e-4, 2.0e-4);
        let thf = ThermalHistoryFem::new(positions, src, 1923.0, 1723.0, 12.0, 620.0, 8000.0);
        assert!(thf.melt_pool_nodes().is_empty());
    }

    #[test]
    fn test_thermal_history_fem_step_increases_temperature() {
        let positions = vec![[0.0; 3]];
        let src = GoldakSource::new(2000.0, 0.8, 1.5e-4, 1.0e-4, 1.0e-4, 2.0e-4);
        let mut thf = ThermalHistoryFem::new(positions, src, 1923.0, 1723.0, 12.0, 620.0, 8000.0);
        let t_init = thf.temperatures[0];
        thf.step([0.0; 3], 1.0e-4, 10.0, 1.0e-8);
        assert!(thf.temperatures[0] > t_init);
    }

    // --- ResidualStressFem tests ---

    #[test]
    fn test_residual_stress_lame_lambda() {
        let fem = ResidualStressFem::new(4, 200.0e9, 0.3, 450.0e6, 11.7e-6, 5.0e9);
        let lam = fem.lame_lambda();
        assert!(lam > 0.0);
    }

    #[test]
    fn test_residual_stress_shear_modulus() {
        let fem = ResidualStressFem::new(4, 200.0e9, 0.3, 450.0e6, 11.7e-6, 5.0e9);
        let g = fem.shear_modulus();
        let expected = 200.0e9 / (2.0 * 1.3);
        assert!((g - expected).abs() < 1.0);
    }

    #[test]
    fn test_residual_stress_yield_at_ambient() {
        let fem = ResidualStressFem::new(4, 200.0e9, 0.3, 450.0e6, 11.7e-6, 5.0e9);
        let sy = fem.yield_stress_at_temperature(T_AMBIENT, 1800.0);
        assert!((sy - 450.0e6).abs() < 1.0);
    }

    #[test]
    fn test_residual_stress_yield_at_melt() {
        let fem = ResidualStressFem::new(4, 200.0e9, 0.3, 450.0e6, 11.7e-6, 5.0e9);
        let sy = fem.yield_stress_at_temperature(1800.0, 1800.0);
        assert!(
            sy < 1.0,
            "yield stress should be near zero at melt: {:.4e}",
            sy
        );
    }

    #[test]
    fn test_von_mises_uniaxial() {
        let mut s = StressState::zero();
        s.stress[0] = 100.0e6;
        let vm = s.von_mises();
        assert!((vm - 100.0e6).abs() < 1.0);
    }

    #[test]
    fn test_von_mises_zero() {
        let s = StressState::zero();
        assert!(s.von_mises() < 1.0e-6);
    }

    // --- PowderBedFusion tests ---

    #[test]
    fn test_pbf_energy_density() {
        let pbf = PowderBedFusion::new(
            200.0,
            0.8,
            1.0e-4,
            5.0e-5,
            0.35,
            4000.0,
            1923.0,
            ScanStrategy::Bidirectional,
        );
        let ev = pbf.energy_density();
        let expected = 200.0 / (0.8 * 1.0e-4 * 5.0e-5);
        assert!((ev - expected).abs() < 1.0);
    }

    #[test]
    fn test_pbf_linear_energy_density() {
        let pbf = PowderBedFusion::new(
            200.0,
            0.8,
            1.0e-4,
            5.0e-5,
            0.35,
            4000.0,
            1923.0,
            ScanStrategy::Unidirectional,
        );
        assert!((pbf.linear_energy_density() - 250.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_pbf_tracks_per_layer() {
        let pbf = PowderBedFusion::new(
            200.0,
            0.8,
            1.0e-4,
            5.0e-5,
            0.35,
            4000.0,
            1923.0,
            ScanStrategy::Island,
        );
        let tracks = pbf.tracks_per_layer(1.0e-3);
        assert_eq!(tracks, 10);
    }

    #[test]
    fn test_pbf_melt_pool_depth_positive() {
        let pbf = PowderBedFusion::new(
            200.0,
            0.8,
            1.0e-4,
            5.0e-5,
            0.35,
            4000.0,
            1923.0,
            ScanStrategy::Bidirectional,
        );
        let depth = pbf.melt_pool_depth(12.0, 0.35);
        assert!(depth > 0.0);
    }

    // --- DirectedEnergyDeposition tests ---

    #[test]
    fn test_ded_effective_power() {
        let ded = DirectedEnergyDeposition::new(2000.0, 0.01, 5.0e-4, 0.7, 0.5, 1.0e-3, 3.0e-3);
        assert!((ded.effective_power() - 1000.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_ded_deposition_rate() {
        let ded = DirectedEnergyDeposition::new(2000.0, 0.01, 5.0e-4, 0.7, 0.5, 1.0e-3, 3.0e-3);
        assert!((ded.deposition_rate() - 3.5e-4).abs() < 1.0e-15);
    }

    #[test]
    fn test_ded_dilution_ratio_bounds() {
        let ded = DirectedEnergyDeposition::new(2000.0, 0.01, 5.0e-4, 0.7, 0.5, 1.0e-3, 3.0e-3);
        let d = ded.dilution_ratio(8000.0, 2.5e5);
        assert!(d > 0.0 && d < 1.0);
    }

    // --- SupportStructure tests ---

    #[test]
    fn test_support_overhang_detection_horizontal() {
        let sup = SupportStructure::new(45.0, 0.3, 0.5, 10);
        // A horizontal face (normal pointing down) should need support
        assert!(sup.needs_support([0.0, 0.0, -1.0]));
    }

    #[test]
    fn test_support_overhang_detection_vertical() {
        let sup = SupportStructure::new(45.0, 0.3, 0.5, 10);
        // A vertical face (normal in x) should NOT need support
        assert!(!sup.needs_support([1.0, 0.0, 0.0]));
    }

    // --- MicrostructureEvolution tests ---

    #[test]
    fn test_jmak_zero_at_t_zero() {
        let jmak = JmakTransformation::new(2.5, 0.01);
        assert!((jmak.transformed_fraction(0.0) - 0.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_jmak_approaches_unity() {
        let jmak = JmakTransformation::new(2.5, 0.01);
        let x = jmak.transformed_fraction(1000.0);
        assert!(x > 0.99);
    }

    #[test]
    fn test_hillert_grain_growth_increases() {
        let hgg = HillertGrainGrowth::new(1.0e-9, 0.5, 100.0e3);
        let r0 = 10.0e-6;
        let r1 = hgg.grain_radius_after(r0, 1.0, 1500.0);
        assert!(r1 >= r0);
    }

    #[test]
    fn test_hall_petch_strength_decreases_with_grain_size() {
        let hgg = HillertGrainGrowth::new(1.0e-9, 0.5, 100.0e3);
        let jmak = JmakTransformation::new(2.5, 0.01);
        let mut me = MicrostructureEvolution::new(1.0e-6, hgg.clone(), jmak.clone());
        let hp1 = me.hall_petch_strength(100.0e6, 500.0e3);
        me.grain_radius = 10.0e-6;
        let hp2 = me.hall_petch_strength(100.0e6, 500.0e3);
        assert!(hp1 > hp2);
    }

    // --- PorosityModel tests ---

    #[test]
    fn test_porosity_model_relative_density() {
        let mut pm = PorosityModel::new(200.0e9, 450.0e6);
        pm.update_from_energy_density(60.0e6, 50.0e6, 150.0e6, 0.001);
        assert!((pm.relative_density() + pm.porosity - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_porosity_model_effective_modulus_decreases() {
        let mut pm = PorosityModel::new(200.0e9, 450.0e6);
        pm.update_from_energy_density(10.0e6, 50.0e6, 150.0e6, 0.001);
        let e_eff = pm.effective_young_modulus();
        assert!(e_eff < 200.0e9);
    }

    #[test]
    fn test_porosity_model_zero_in_optimal_window() {
        let mut pm = PorosityModel::new(200.0e9, 450.0e6);
        pm.update_from_energy_density(100.0e6, 50.0e6, 150.0e6, 0.0);
        assert!((pm.lof_fraction + pm.keyhole_fraction).abs() < 1.0e-15);
    }

    // --- DistortionCompensation tests ---

    #[test]
    fn test_distortion_rms_zero_initially() {
        let positions = vec![[0.0; 3], [0.01, 0.0, 0.0], [0.0, 0.01, 0.0]];
        let dc = DistortionCompensation::new(positions, 1.0);
        assert!(dc.rms_distortion() < 1.0e-15);
    }

    #[test]
    fn test_distortion_max_distortion() {
        let nominal = vec![[0.0; 3], [1.0, 0.0, 0.0]];
        let mut dc = DistortionCompensation::new(nominal, 1.0);
        dc.set_distorted_positions(vec![[0.0, 0.0, 0.001], [1.0, 0.0, 0.001]]);
        let max_d = dc.max_distortion();
        assert!((max_d - 0.001).abs() < 1.0e-12);
    }

    #[test]
    fn test_distortion_compensation_direction() {
        let nominal = vec![[0.0; 3]];
        let mut dc = DistortionCompensation::new(nominal, 1.0);
        dc.set_distorted_positions(vec![[0.0, 0.0, 0.005]]);
        // With full compensation, compensated_pos should be opposite
        let comp = dc.compensated_positions[0];
        assert!(comp[2] < 0.0);
    }

    // --- BuildProcessOptimization tests ---

    #[test]
    fn test_bpo_optimise_returns_valid_index() {
        let weights = OptimizationWeights::uniform();
        let candidates = vec![
            ProcessParameters::new(150.0, 0.8, 1.2e-4, 5.0e-5),
            ProcessParameters::new(200.0, 1.0, 1.0e-4, 5.0e-5),
            ProcessParameters::new(250.0, 1.5, 8.0e-5, 5.0e-5),
        ];
        let mut opt = BuildProcessOptimization::new(weights, candidates);
        let best = opt.optimise(5.0e9);
        assert!(best < 3);
    }

    #[test]
    fn test_bpo_best_parameters_not_none() {
        let weights = OptimizationWeights::uniform();
        let candidates = vec![ProcessParameters::new(200.0, 1.0, 1.0e-4, 5.0e-5)];
        let mut opt = BuildProcessOptimization::new(weights, candidates);
        opt.optimise(5.0e9);
        assert!(opt.best_parameters().is_some());
    }

    #[test]
    fn test_process_parameters_energy_density() {
        let p = ProcessParameters::new(200.0, 1.0, 1.0e-4, 5.0e-5);
        let ev = p.energy_density();
        let expected = 200.0 / (1.0 * 1.0e-4 * 5.0e-5);
        assert!((ev - expected).abs() < 1.0);
    }

    #[test]
    fn test_optimization_weights_normalise() {
        let w = OptimizationWeights {
            stress_weight: 2.0,
            distortion_weight: 2.0,
            time_weight: 2.0,
            porosity_weight: 2.0,
        };
        let wn = w.normalise();
        let total = wn.stress_weight + wn.distortion_weight + wn.time_weight + wn.porosity_weight;
        assert!((total - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_temperature_dependent_conductivity() {
        let positions = vec![[0.0; 3]];
        let src = GoldakSource::new(200.0, 0.35, 1.5e-4, 1.0e-4, 1.0e-4, 2.0e-4);
        let thf = ThermalHistoryFem::new(positions, src, 1923.0, 1723.0, 12.0, 620.0, 8000.0);
        // At reference temperature, conductivity = base conductivity
        let k = thf.temperature_dependent_conductivity(298.15, 1.0e-4);
        assert!((k - 12.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_stress_pressure() {
        let mut s = StressState::zero();
        s.stress[0] = 100.0e6;
        s.stress[1] = 200.0e6;
        s.stress[2] = 300.0e6;
        let p = s.pressure();
        assert!((p - 200.0e6).abs() < 1.0);
    }
}
