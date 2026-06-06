// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Biological soft body simulation.
//!
//! Implements:
//! - [`CorticalTension`] — cell cortical tension and internal pressure
//! - [`CellDivision`] — mitotic cell division simulation
//! - `CellMigration` — lamellipodia and filopodia protrusion mechanics
//! - [`CellAdhesion`] — cadherin-based cell–cell adhesion
//! - [`ExtracellularMatrix`] — ECM mechanics (collagen, fibronectin)
//! - [`WoundHealing`] — collective cell migration in wound closure
//! - [`TumorGrowth`] — tumor growth mechanics and solid stress
//! - [`RedBloodCell`] — biconcave RBC membrane, spectrin network
//! - [`PlateletAggregation`] — platelet activation and clot formation
//! - [`BoneRemodeling`] — Wolff's law osteoblast/osteoclast coupling
//!
//! # Conventions
//!
//! * SI units throughout (m, kg, s, N, Pa).
//! * Right-hand coordinate system.
//! * All angles in radians.
//!
//! # References
//!
//! * Alberts et al. (2014) – Molecular Biology of the Cell.
//! * Discher, Boal & Boey (1998) – Simulations of the erythrocyte cytoskeleton.
//! * Wolff (1892) – Das Gesetz der Transformation der Knochen.
//! * Chauvière et al. (2010) – Modeling cell migration mechanics.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

/// 3-D vector type alias.
type Vec3 = [f64; 3];

/// Dot product of two 3-D vectors.
#[inline]
fn dot3(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean norm of a 3-D vector.
#[inline]
fn norm3(v: Vec3) -> f64 {
    dot3(v, v).sqrt()
}

/// Normalise a 3-D vector (returns zero vector if degenerate).
#[inline]
fn normalize3(v: Vec3) -> Vec3 {
    let n = norm3(v);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}

/// Element-wise addition of two 3-D vectors.
#[inline]
fn add3(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Element-wise subtraction of two 3-D vectors.
#[inline]
fn sub3(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-D vector by scalar `s`.
#[inline]
fn scale3(v: Vec3, s: f64) -> Vec3 {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Distance between two 3-D points.
#[inline]
fn dist3(a: Vec3, b: Vec3) -> f64 {
    norm3(sub3(a, b))
}

/// Clamp a value to `[lo, hi]`.
#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    if x < lo {
        lo
    } else if x > hi {
        hi
    } else {
        x
    }
}

// ---------------------------------------------------------------------------
// 1. Cortical Tension and Internal Pressure
// ---------------------------------------------------------------------------

/// Parameters governing the cortical tension and turgor pressure of a cell.
///
/// The Young–Laplace relation `ΔP = 2 γ / R` balances cortical tension `γ`
/// against the internal gauge pressure `ΔP` for a spherical cell of radius `R`.
#[derive(Debug, Clone)]
pub struct CorticalTension {
    /// Cortical surface tension (N m⁻¹).
    pub tension: f64,
    /// Internal hydrostatic pressure relative to external (Pa).
    pub internal_pressure: f64,
    /// Cell radius (m).
    pub radius: f64,
    /// Membrane bending modulus κ (J).
    pub bending_modulus: f64,
    /// Spontaneous curvature C₀ (m⁻¹).
    pub spontaneous_curvature: f64,
}

impl CorticalTension {
    /// Create a new `CorticalTension` model.
    ///
    /// * `tension` — cortical tension γ (N m⁻¹).
    /// * `radius` — cell radius (m).
    /// * `bending_modulus` — Helfrich bending modulus κ (J).
    /// * `spontaneous_curvature` — spontaneous curvature C₀ (m⁻¹).
    pub fn new(
        tension: f64,
        radius: f64,
        bending_modulus: f64,
        spontaneous_curvature: f64,
    ) -> Self {
        let internal_pressure = if radius > 1e-15 {
            2.0 * tension / radius
        } else {
            0.0
        };
        Self {
            tension,
            internal_pressure,
            radius,
            bending_modulus,
            spontaneous_curvature,
        }
    }

    /// Compute the Young–Laplace pressure difference `ΔP = 2 γ / R`.
    pub fn young_laplace_pressure(&self) -> f64 {
        if self.radius > 1e-15 {
            2.0 * self.tension / self.radius
        } else {
            0.0
        }
    }

    /// Helfrich bending energy for a sphere: `E = 8π κ (1 - C₀ R)²`.
    ///
    /// For a pure sphere the mean curvature H = 1/R everywhere.
    pub fn helfrich_energy(&self) -> f64 {
        let h = if self.radius > 1e-15 {
            1.0 / self.radius
        } else {
            0.0
        };
        let dh = h - self.spontaneous_curvature;
        8.0 * PI * self.bending_modulus * dh * dh * self.radius * self.radius
    }

    /// Cortical tension contribution to the net outward force on a hemisphere
    /// cross-section (N): `F = π R γ` (line tension at equator).
    pub fn equatorial_line_force(&self) -> f64 {
        PI * self.radius * self.tension
    }

    /// Update the equilibrium radius given a target pressure `dp` (Pa).
    ///
    /// Solves `dp = 2 γ / R` → `R = 2 γ / dp`.
    pub fn equilibrium_radius(&self, dp: f64) -> f64 {
        if dp > 1e-15 {
            2.0 * self.tension / dp
        } else {
            f64::INFINITY
        }
    }

    /// Effective spring constant for radial perturbations (N m⁻¹).
    ///
    /// Linearising `ΔP(R) = 2γ/R` around `R₀` gives `k = 2γ / R₀²`.
    pub fn radial_spring_constant(&self) -> f64 {
        if self.radius > 1e-15 {
            2.0 * self.tension / (self.radius * self.radius)
        } else {
            0.0
        }
    }
}

// ---------------------------------------------------------------------------
// 2. Cell Division
// ---------------------------------------------------------------------------

/// Phase of the cell cycle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CellCyclePhase {
    /// G1 gap phase — cell growth.
    G1,
    /// S phase — DNA synthesis.
    S,
    /// G2 gap phase — preparation for division.
    G2,
    /// M phase — mitosis.
    M,
}

/// State of a single simulated cell undergoing division.
#[derive(Debug, Clone)]
pub struct CellDivision {
    /// Cell centre position (m).
    pub position: Vec3,
    /// Cell radius (m).
    pub radius: f64,
    /// Current phase of the cell cycle.
    pub phase: CellCyclePhase,
    /// Elapsed time in the current phase (s).
    pub phase_time: f64,
    /// Duration of each phase \[G1, S, G2, M\] (s).
    pub phase_durations: [f64; 4],
    /// Contractile ring tension driving cytokinesis (N m⁻¹).
    pub contractile_ring_tension: f64,
    /// Cleavage furrow ingression fraction in \[0, 1\].
    pub furrow_fraction: f64,
    /// Whether the cell has completed division.
    pub divided: bool,
}

impl CellDivision {
    /// Construct a new `CellDivision` cell starting in G1.
    ///
    /// * `position` — initial centre (m).
    /// * `radius` — initial cell radius (m).
    /// * `phase_durations` — duration \[G1, S, G2, M\] in seconds.
    /// * `contractile_ring_tension` — actomyosin ring tension (N m⁻¹).
    pub fn new(
        position: Vec3,
        radius: f64,
        phase_durations: [f64; 4],
        contractile_ring_tension: f64,
    ) -> Self {
        Self {
            position,
            radius,
            phase: CellCyclePhase::G1,
            phase_time: 0.0,
            phase_durations,
            contractile_ring_tension,
            furrow_fraction: 0.0,
            divided: false,
        }
    }

    /// Advance the cell cycle by `dt` seconds.
    ///
    /// Returns `Some((daughter_a, daughter_b))` positions when division
    /// completes, or `None` otherwise.
    pub fn step(&mut self, dt: f64) -> Option<(Vec3, Vec3)> {
        if self.divided {
            return None;
        }
        self.phase_time += dt;
        let duration = match self.phase {
            CellCyclePhase::G1 => self.phase_durations[0],
            CellCyclePhase::S => self.phase_durations[1],
            CellCyclePhase::G2 => self.phase_durations[2],
            CellCyclePhase::M => self.phase_durations[3],
        };
        if self.phase_time >= duration {
            self.phase_time = 0.0;
            self.phase = match self.phase {
                CellCyclePhase::G1 => CellCyclePhase::S,
                CellCyclePhase::S => CellCyclePhase::G2,
                CellCyclePhase::G2 => CellCyclePhase::M,
                CellCyclePhase::M => {
                    // Division event
                    self.divided = true;
                    let offset = self.radius * 0.6;
                    let da = [
                        self.position[0] + offset,
                        self.position[1],
                        self.position[2],
                    ];
                    let db = [
                        self.position[0] - offset,
                        self.position[1],
                        self.position[2],
                    ];
                    return Some((da, db));
                }
            };
        }
        // During M phase, advance furrow ingression
        if self.phase == CellCyclePhase::M {
            let frac = self.phase_time / duration;
            self.furrow_fraction = clamp(frac, 0.0, 1.0);
        }
        None
    }

    /// Effective equatorial radius during cytokinesis (contractile ring).
    ///
    /// `R_eq(f) = R * (1 - f)` where `f` is the furrow ingression fraction.
    pub fn equatorial_radius(&self) -> f64 {
        self.radius * (1.0 - self.furrow_fraction)
    }

    /// Force generated by the contractile ring at current ingression: `F = 2π R_eq γ`.
    pub fn contractile_force(&self) -> f64 {
        2.0 * PI * self.equatorial_radius() * self.contractile_ring_tension
    }
}

// ---------------------------------------------------------------------------
// 3. Cell Migration (Lamellipodia / Filopodia)
// ---------------------------------------------------------------------------

/// A single filopodium protrusion.
#[derive(Debug, Clone)]
pub struct Filopodium {
    /// Tip position (m).
    pub tip: Vec3,
    /// Base (attachment) position (m).
    pub base: Vec3,
    /// Protrusion speed (m s⁻¹).
    pub protrusion_speed: f64,
    /// Retraction speed (m s⁻¹).
    pub retraction_speed: f64,
    /// Tension force along the filopodium (N).
    pub tension: f64,
    /// Whether the tip is currently adhered.
    pub adhered: bool,
}

impl Filopodium {
    /// Create a new filopodium.
    ///
    /// * `base` — attachment point at cell edge (m).
    /// * `direction` — unit vector pointing outward.
    /// * `initial_length` — initial extension (m).
    /// * `protrusion_speed` — extension speed (m s⁻¹).
    /// * `retraction_speed` — retraction speed (m s⁻¹).
    /// * `tension` — internal actin bundle tension (N).
    pub fn new(
        base: Vec3,
        direction: Vec3,
        initial_length: f64,
        protrusion_speed: f64,
        retraction_speed: f64,
        tension: f64,
    ) -> Self {
        let d = normalize3(direction);
        let tip = add3(base, scale3(d, initial_length));
        Self {
            tip,
            base,
            protrusion_speed,
            retraction_speed,
            tension,
            adhered: false,
        }
    }

    /// Current length of the filopodium (m).
    pub fn length(&self) -> f64 {
        dist3(self.tip, self.base)
    }

    /// Unit vector from base to tip.
    pub fn direction(&self) -> Vec3 {
        normalize3(sub3(self.tip, self.base))
    }

    /// Extend the filopodium by `dt` seconds.  Returns new length.
    pub fn protrude(&mut self, dt: f64) -> f64 {
        let d = self.direction();
        self.tip = add3(self.tip, scale3(d, self.protrusion_speed * dt));
        self.length()
    }

    /// Retract the filopodium by `dt` seconds.  Returns new length.
    pub fn retract(&mut self, dt: f64) -> f64 {
        let len = self.length();
        let retraction = self.retraction_speed * dt;
        if retraction >= len {
            self.tip = self.base;
            return 0.0;
        }
        let d = self.direction();
        self.tip = sub3(self.tip, scale3(d, retraction));
        self.length()
    }

    /// Traction force exerted at the base (N), directed tip→base when adhered.
    pub fn traction_force(&self) -> Vec3 {
        if !self.adhered {
            return [0.0; 3];
        }
        let d = normalize3(sub3(self.base, self.tip));
        scale3(d, self.tension)
    }
}

/// Lamellipodium model — a thin actin-rich sheet at the cell leading edge.
#[derive(Debug, Clone)]
pub struct Lamellipodium {
    /// Half-width of the lamellipodium (m).
    pub half_width: f64,
    /// Protrusion velocity (m s⁻¹).
    pub protrusion_velocity: f64,
    /// Retrograde actin flow speed (m s⁻¹).
    pub retrograde_flow: f64,
    /// Traction stress on the substrate (Pa).
    pub traction_stress: f64,
    /// Leading-edge position (m, 1-D).
    pub leading_edge: f64,
}

impl Lamellipodium {
    /// Construct a new lamellipodium.
    ///
    /// * `initial_edge` — initial leading-edge position (m).
    /// * `half_width` — lateral half-width (m).
    /// * `protrusion_velocity` — forward protrusion speed (m s⁻¹).
    /// * `retrograde_flow` — retrograde actin flow (m s⁻¹).
    /// * `traction_stress` — substrate traction stress (Pa).
    pub fn new(
        initial_edge: f64,
        half_width: f64,
        protrusion_velocity: f64,
        retrograde_flow: f64,
        traction_stress: f64,
    ) -> Self {
        Self {
            half_width,
            protrusion_velocity,
            retrograde_flow,
            traction_stress,
            leading_edge: initial_edge,
        }
    }

    /// Net edge velocity: `v_edge = v_protrusion - v_retrograde`.
    pub fn net_velocity(&self) -> f64 {
        self.protrusion_velocity - self.retrograde_flow
    }

    /// Advance the leading edge by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        self.leading_edge += self.net_velocity() * dt;
    }

    /// Total traction force over the lamellipodium area `F = σ × 2w × Δx`.
    pub fn total_traction(&self, width_depth: f64) -> f64 {
        self.traction_stress * 2.0 * self.half_width * width_depth
    }
}

// ---------------------------------------------------------------------------
// 4. Cell–Cell Adhesion (Cadherin)
// ---------------------------------------------------------------------------

/// Cadherin-mediated cell–cell adhesion bond model.
///
/// Each bond is modelled as a spring with stochastic catch-bond kinetics.
#[derive(Debug, Clone)]
pub struct CadherinBond {
    /// Bond rest length (m).
    pub rest_length: f64,
    /// Bond spring stiffness (N m⁻¹).
    pub stiffness: f64,
    /// Off-rate at zero force (s⁻¹).
    pub k_off_zero: f64,
    /// Characteristic force for catch-bond slip (`Bell` parameter, N).
    pub bell_force: f64,
    /// Whether the bond is currently engaged.
    pub engaged: bool,
}

impl CadherinBond {
    /// Create a new cadherin bond.
    ///
    /// * `rest_length` — equilibrium bond length (m).
    /// * `stiffness` — bond spring constant (N m⁻¹).
    /// * `k_off_zero` — zero-force off-rate (s⁻¹).
    /// * `bell_force` — Bell-model characteristic force (N).
    pub fn new(rest_length: f64, stiffness: f64, k_off_zero: f64, bell_force: f64) -> Self {
        Self {
            rest_length,
            stiffness,
            k_off_zero,
            bell_force,
            engaged: true,
        }
    }

    /// Spring force magnitude for separation distance `d` (N).
    pub fn force(&self, d: f64) -> f64 {
        if !self.engaged {
            return 0.0;
        }
        self.stiffness * (d - self.rest_length).abs()
    }

    /// Bell-model off-rate: `k_off(F) = k_off_0 * exp(F / F_bell)`.
    pub fn off_rate(&self, force: f64) -> f64 {
        self.k_off_zero * (force / self.bell_force).exp()
    }

    /// Probability of bond surviving for time `dt` under force `F`.
    ///
    /// `P_survive = exp(-k_off(F) * dt)`.
    pub fn survival_probability(&self, force: f64, dt: f64) -> f64 {
        (-self.off_rate(force) * dt).exp()
    }

    /// Stochastically update bond state given a random value `r` in \[0, 1).
    ///
    /// The bond breaks if `r > P_survive`.
    pub fn update_stochastic(&mut self, force: f64, dt: f64, r: f64) {
        if self.engaged && r > self.survival_probability(force, dt) {
            self.engaged = false;
        }
    }
}

/// A cluster of cadherin bonds between two cells.
#[derive(Debug, Clone)]
pub struct CellAdhesion {
    /// Position of cell A centre (m).
    pub pos_a: Vec3,
    /// Position of cell B centre (m).
    pub pos_b: Vec3,
    /// Collection of individual cadherin bonds.
    pub bonds: Vec<CadherinBond>,
}

impl CellAdhesion {
    /// Create a new adhesion junction with `n_bonds` identical bonds.
    ///
    /// * `pos_a`, `pos_b` — cell centre positions (m).
    /// * `n_bonds` — initial number of bonds.
    /// * `rest_length` — each bond's rest length (m).
    /// * `stiffness` — each bond's stiffness (N m⁻¹).
    /// * `k_off_zero` — zero-force off-rate (s⁻¹).
    /// * `bell_force` — Bell-model characteristic force (N).
    pub fn new(
        pos_a: Vec3,
        pos_b: Vec3,
        n_bonds: usize,
        rest_length: f64,
        stiffness: f64,
        k_off_zero: f64,
        bell_force: f64,
    ) -> Self {
        let bonds = (0..n_bonds)
            .map(|_| CadherinBond::new(rest_length, stiffness, k_off_zero, bell_force))
            .collect();
        Self {
            pos_a,
            pos_b,
            bonds,
        }
    }

    /// Number of currently engaged bonds.
    pub fn active_bond_count(&self) -> usize {
        self.bonds.iter().filter(|b| b.engaged).count()
    }

    /// Total adhesion force magnitude (N): sum of individual bond forces.
    pub fn total_force(&self) -> f64 {
        let d = dist3(self.pos_a, self.pos_b);
        self.bonds
            .iter()
            .filter(|b| b.engaged)
            .map(|b| b.force(d))
            .sum()
    }

    /// Unit vector from A to B (adhesion direction).
    pub fn adhesion_direction(&self) -> Vec3 {
        normalize3(sub3(self.pos_b, self.pos_a))
    }

    /// Force vector acting on cell A (directed towards B).
    pub fn force_on_a(&self) -> Vec3 {
        scale3(self.adhesion_direction(), self.total_force())
    }
}

// ---------------------------------------------------------------------------
// 5. Extracellular Matrix Mechanics
// ---------------------------------------------------------------------------

/// Extracellular matrix (ECM) mechanical model.
///
/// Implements a fibre-reinforced viscoelastic solid with strain-stiffening.
#[derive(Debug, Clone)]
pub struct ExtracellularMatrix {
    /// Young's modulus at zero strain (Pa).
    pub young_modulus: f64,
    /// Poisson ratio (dimensionless).
    pub poisson_ratio: f64,
    /// Viscous relaxation time (s).
    pub relaxation_time: f64,
    /// Strain-stiffening coefficient β (dimensionless).
    pub strain_stiffening: f64,
    /// Current isotropic strain (dimensionless).
    pub strain: f64,
    /// Collagen fibre volume fraction (0–1).
    pub collagen_fraction: f64,
}

impl ExtracellularMatrix {
    /// Construct a new ECM model.
    ///
    /// * `young_modulus` — linear Young's modulus (Pa).
    /// * `poisson_ratio` — Poisson ratio.
    /// * `relaxation_time` — Maxwell relaxation time (s).
    /// * `strain_stiffening` — exponent for exponential strain stiffening.
    /// * `collagen_fraction` — volume fraction of collagen fibres (0–1).
    pub fn new(
        young_modulus: f64,
        poisson_ratio: f64,
        relaxation_time: f64,
        strain_stiffening: f64,
        collagen_fraction: f64,
    ) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
            relaxation_time,
            strain_stiffening,
            strain: 0.0,
            collagen_fraction,
        }
    }

    /// Effective stiffness with strain stiffening: `E_eff = E₀ * exp(β ε)`.
    pub fn effective_stiffness(&self) -> f64 {
        self.young_modulus * (self.strain_stiffening * self.strain).exp()
    }

    /// Cauchy stress for current strain (Pa): `σ = E_eff * ε`.
    pub fn stress(&self) -> f64 {
        self.effective_stiffness() * self.strain
    }

    /// Shear modulus: `G = E / (2(1 + ν))`.
    pub fn shear_modulus(&self) -> f64 {
        self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Bulk modulus: `K = E / (3(1 - 2ν))`.
    pub fn bulk_modulus(&self) -> f64 {
        self.young_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }

    /// Maxwell viscoelastic stress relaxation over time step `dt`.
    ///
    /// Updates the strain by relaxation: `dε/dt = -ε / τ`.
    pub fn relax(&mut self, dt: f64) {
        let decay = (-dt / self.relaxation_time).exp();
        self.strain *= decay;
    }

    /// Fibre-reinforced stiffness contribution from collagen.
    ///
    /// Uses rule-of-mixtures: `E_composite = (1 - φ) E_matrix + φ E_fibre`.
    /// Collagen fibre modulus is ~1.5 GPa.
    pub fn composite_stiffness(&self) -> f64 {
        let e_fibre = 1.5e9_f64;
        (1.0 - self.collagen_fraction) * self.young_modulus + self.collagen_fraction * e_fibre
    }
}

// ---------------------------------------------------------------------------
// 6. Wound Healing Simulation
// ---------------------------------------------------------------------------

/// A 1-D wound-healing model tracking wound gap closure.
///
/// Models the advancing cell sheet using a Fisher–KPP-inspired phenomenology
/// with contractile cable tension and cell proliferation.
#[derive(Debug, Clone)]
pub struct WoundHealing {
    /// Width of the wound (m).
    pub wound_width: f64,
    /// Cell sheet velocity at the wound edge (m s⁻¹).
    pub edge_velocity: f64,
    /// Actomyosin cable tension at wound edge (N m⁻¹).
    pub cable_tension: f64,
    /// Cell proliferation rate (s⁻¹).
    pub proliferation_rate: f64,
    /// Cell migration speed due to lamellipodia (m s⁻¹).
    pub migration_speed: f64,
    /// Elapsed simulation time (s).
    pub time: f64,
    /// History of wound widths over time (m).
    pub width_history: Vec<f64>,
}

impl WoundHealing {
    /// Create a new wound healing simulation.
    ///
    /// * `initial_width` — initial wound gap width (m).
    /// * `cable_tension` — actomyosin purse-string tension (N m⁻¹).
    /// * `proliferation_rate` — cell division rate at wound edge (s⁻¹).
    /// * `migration_speed` — lamellipodial migration speed (m s⁻¹).
    pub fn new(
        initial_width: f64,
        cable_tension: f64,
        proliferation_rate: f64,
        migration_speed: f64,
    ) -> Self {
        Self {
            wound_width: initial_width,
            edge_velocity: 0.0,
            cable_tension,
            proliferation_rate,
            migration_speed,
            time: 0.0,
            width_history: vec![initial_width],
        }
    }

    /// Advance wound healing by `dt` seconds.
    ///
    /// The wound closes through both active migration and purse-string contraction.
    /// Closure rate: `dW/dt = -2 (v_mig + T_cable / (η W))`.
    pub fn step(&mut self, dt: f64, viscosity: f64) {
        let purse_string_velocity = if self.wound_width > 1e-10 {
            self.cable_tension / (viscosity * self.wound_width)
        } else {
            0.0
        };
        let closure_rate = 2.0 * (self.migration_speed + purse_string_velocity);
        let new_width = (self.wound_width - closure_rate * dt).max(0.0);
        self.edge_velocity = closure_rate;
        self.wound_width = new_width;
        self.time += dt;
        self.width_history.push(new_width);
    }

    /// Whether the wound is fully closed.
    pub fn is_closed(&self) -> bool {
        self.wound_width < 1e-9
    }

    /// Estimate closure time using exponential fit: `τ = W₀ / (2 v_mig)`.
    pub fn estimated_closure_time(&self) -> f64 {
        if self.migration_speed > 1e-15 {
            self.width_history[0] / (2.0 * self.migration_speed)
        } else {
            f64::INFINITY
        }
    }

    /// Fractional closure: `(W₀ - W) / W₀`.
    pub fn closure_fraction(&self) -> f64 {
        let w0 = self.width_history[0];
        if w0 > 1e-15 {
            (w0 - self.wound_width) / w0
        } else {
            1.0
        }
    }
}

// ---------------------------------------------------------------------------
// 7. Tumor Growth Mechanics
// ---------------------------------------------------------------------------

/// Isotropic solid-stress model for a growing spherical tumor.
///
/// Based on the biphasic mixture theory of Ambrosi & Mollica (2002).
#[derive(Debug, Clone)]
pub struct TumorGrowth {
    /// Tumor radius (m).
    pub radius: f64,
    /// Tumor growth rate constant (s⁻¹).
    pub growth_rate: f64,
    /// Stiffness of surrounding tissue (Pa).
    pub tissue_stiffness: f64,
    /// Osmotic pressure inside the tumor (Pa).
    pub osmotic_pressure: f64,
    /// Interstitial fluid pressure (Pa).
    pub interstitial_pressure: f64,
    /// Necrotic core radius (m).
    pub necrotic_radius: f64,
}

impl TumorGrowth {
    /// Create a new tumor growth model.
    ///
    /// * `initial_radius` — initial tumor radius (m).
    /// * `growth_rate` — exponential growth rate (s⁻¹).
    /// * `tissue_stiffness` — Young's modulus of surrounding tissue (Pa).
    /// * `osmotic_pressure` — oncotic pressure (Pa).
    pub fn new(
        initial_radius: f64,
        growth_rate: f64,
        tissue_stiffness: f64,
        osmotic_pressure: f64,
    ) -> Self {
        Self {
            radius: initial_radius,
            growth_rate,
            tissue_stiffness,
            osmotic_pressure,
            interstitial_pressure: 0.0,
            necrotic_radius: 0.0,
        }
    }

    /// Advance growth by `dt` seconds.
    ///
    /// Exponential growth limited by mechanical stress from surrounding tissue:
    /// `dR/dt = k_g R (1 - σ_mech / σ_max)`.
    pub fn step(&mut self, dt: f64, max_stress: f64) {
        let sigma = self.radial_stress_at_surface();
        let inhibition = clamp(1.0 - sigma / max_stress, 0.0, 1.0);
        self.radius += self.growth_rate * self.radius * inhibition * dt;
        // Update interstitial fluid pressure via Darcy flow approximation
        self.interstitial_pressure = self.osmotic_pressure * (self.radius / (self.radius + 1e-6));
        // Necrotic core forms when radius exceeds diffusion limit (~200 µm)
        let diffusion_limit = 2e-4;
        if self.radius > diffusion_limit {
            self.necrotic_radius = self.radius - diffusion_limit;
        }
    }

    /// Radial compressive stress at the tumor surface (Pa).
    ///
    /// Approximation: `σ_r = E_tissue / 3 * (R/R₀)³ - 1` for small deformations.
    pub fn radial_stress_at_surface(&self) -> f64 {
        self.tissue_stiffness * self.radius / 3.0
    }

    /// Tumor volume (m³).
    pub fn volume(&self) -> f64 {
        4.0 / 3.0 * PI * self.radius.powi(3)
    }

    /// Viable rim volume (m³) — excluding necrotic core.
    pub fn viable_volume(&self) -> f64 {
        let v_total = self.volume();
        let v_necrotic = 4.0 / 3.0 * PI * self.necrotic_radius.powi(3);
        (v_total - v_necrotic).max(0.0)
    }

    /// Doubling time based on current growth rate (s).
    pub fn doubling_time(&self) -> f64 {
        if self.growth_rate > 1e-15 {
            2_f64.ln() / self.growth_rate
        } else {
            f64::INFINITY
        }
    }
}

// ---------------------------------------------------------------------------
// 8. Red Blood Cell Membrane
// ---------------------------------------------------------------------------

/// Spectrin network spring connecting two nodes on a RBC membrane.
#[derive(Debug, Clone)]
pub struct SpectrinSpring {
    /// Indices of the two nodes.
    pub nodes: [usize; 2],
    /// Rest length (m).
    pub rest_length: f64,
    /// Worm-like chain persistence length (m).
    pub persistence_length: f64,
    /// Contour length of the spectrin tetramer (m).
    pub contour_length: f64,
}

impl SpectrinSpring {
    /// Create a new spectrin spring.
    ///
    /// * `nodes` — indices of the two connected nodes.
    /// * `rest_length` — equilibrium spring length (m).
    /// * `persistence_length` — WLC persistence length (m).
    /// * `contour_length` — spectrin tetramer contour length (m).
    pub fn new(
        nodes: [usize; 2],
        rest_length: f64,
        persistence_length: f64,
        contour_length: f64,
    ) -> Self {
        Self {
            nodes,
            rest_length,
            persistence_length,
            contour_length,
        }
    }

    /// Worm-like chain force for extension `x` (N).
    ///
    /// `F = kT / (4 Lp) * [1/(1 - x/L)² - 1 + 4x/L]`  (Marko–Siggia).
    /// `kT ≈ 4.1 pN·nm` at 37 °C.
    pub fn wlc_force(&self, extension: f64) -> f64 {
        let kt = 4.28e-21_f64; // kT at 37°C (J)
        let x = clamp(extension, 0.0, self.contour_length * 0.9999);
        let s = x / self.contour_length;
        kt / (4.0 * self.persistence_length) * (1.0 / (1.0 - s).powi(2) - 1.0 + 4.0 * s)
    }

    /// Force based on current separation `d` (N).
    pub fn force_at_length(&self, d: f64) -> f64 {
        self.wlc_force(d)
    }
}

/// Red blood cell membrane model using a coarse-grained spectrin network.
///
/// The biconcave shape is approximated by a flattened sphere with
/// indented poles. The spectrin network provides shear elasticity.
#[derive(Debug, Clone)]
pub struct RedBloodCell {
    /// Node positions on the membrane (m).
    pub nodes: Vec<Vec3>,
    /// Spectrin springs.
    pub springs: Vec<SpectrinSpring>,
    /// Membrane bending modulus κ (J).
    pub bending_modulus: f64,
    /// Area compression modulus KA (N m⁻¹).
    pub area_modulus: f64,
    /// Target surface area (m²).
    pub target_area: f64,
    /// Target cell volume (m³).
    pub target_volume: f64,
}

impl RedBloodCell {
    /// Construct an RBC model.
    ///
    /// * `nodes` — surface node positions.
    /// * `springs` — spectrin spring topology.
    /// * `bending_modulus` — κ ≈ 2×10⁻¹⁹ J for healthy RBC.
    /// * `area_modulus` — KA ≈ 4.8×10⁻⁶ N m⁻¹ for healthy RBC.
    /// * `target_area` — target surface area (m²).
    /// * `target_volume` — target enclosed volume (m³).
    pub fn new(
        nodes: Vec<Vec3>,
        springs: Vec<SpectrinSpring>,
        bending_modulus: f64,
        area_modulus: f64,
        target_area: f64,
        target_volume: f64,
    ) -> Self {
        Self {
            nodes,
            springs,
            bending_modulus,
            area_modulus,
            target_area,
            target_volume,
        }
    }

    /// Compute total spectrin elastic energy (J).
    pub fn elastic_energy(&self) -> f64 {
        self.springs
            .iter()
            .map(|s| {
                let d = dist3(self.nodes[s.nodes[0]], self.nodes[s.nodes[1]]);
                // Harmonic approximation
                let k = 1e-5_f64; // pN/nm scale
                0.5 * k * (d - s.rest_length).powi(2)
            })
            .sum()
    }

    /// Biconcave shape metric: ratio of polar to equatorial diameter.
    ///
    /// For a healthy RBC ≈ 0.4.
    pub fn biconcave_ratio(&self) -> f64 {
        if self.nodes.is_empty() {
            return 0.0;
        }
        let max_r = self
            .nodes
            .iter()
            .map(|p| (p[0] * p[0] + p[2] * p[2]).sqrt())
            .fold(0.0_f64, f64::max);
        let max_z = self
            .nodes
            .iter()
            .map(|p| p[1].abs())
            .fold(0.0_f64, f64::max);
        if max_r > 1e-15 { max_z / max_r } else { 0.0 }
    }

    /// Area constraint energy: `E_A = K_A/2 * (A - A₀)² / A₀`.
    pub fn area_constraint_energy(&self, current_area: f64) -> f64 {
        self.area_modulus * 0.5 * (current_area - self.target_area).powi(2) / self.target_area
    }

    /// Reduced volume: `V* = V / (4π/3 (A/4π)^(3/2))`.
    ///
    /// For a sphere V* = 1; for a healthy RBC V* ≈ 0.64.
    pub fn reduced_volume(&self, current_volume: f64, current_area: f64) -> f64 {
        if current_area < 1e-30 {
            return 0.0;
        }
        let r_eff = (current_area / (4.0 * PI)).sqrt();
        let v_sphere = 4.0 / 3.0 * PI * r_eff.powi(3);
        if v_sphere > 1e-30 {
            current_volume / v_sphere
        } else {
            0.0
        }
    }
}

// ---------------------------------------------------------------------------
// 9. Platelet Aggregation
// ---------------------------------------------------------------------------

/// Activation state of a platelet.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlateletState {
    /// Resting — discoid shape, low adhesion.
    Resting,
    /// Activated — shape change, high adhesion, granule release.
    Activated,
    /// Aggregated — bound into thrombus.
    Aggregated,
}

/// A single platelet agent.
#[derive(Debug, Clone)]
pub struct Platelet {
    /// Position (m).
    pub position: Vec3,
    /// Velocity (m s⁻¹).
    pub velocity: Vec3,
    /// Activation state.
    pub state: PlateletState,
    /// Activation level in \[0, 1\].
    pub activation_level: f64,
    /// GPIbα receptor density (m⁻²).
    pub gpib_density: f64,
    /// αIIbβ3 integrin active fraction (0–1).
    pub integrin_active_fraction: f64,
}

impl Platelet {
    /// Construct a resting platelet.
    ///
    /// * `position` — initial position (m).
    /// * `gpib_density` — vWF receptor density (m⁻²).
    pub fn new(position: Vec3, gpib_density: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            state: PlateletState::Resting,
            activation_level: 0.0,
            gpib_density,
            integrin_active_fraction: 0.05,
        }
    }

    /// Activate the platelet: inside-out signalling increases integrin affinity.
    pub fn activate(&mut self) {
        self.state = PlateletState::Activated;
        self.activation_level = 1.0;
        self.integrin_active_fraction = 0.9;
    }

    /// Integrate platelet position by `dt` under external force `f` and drag
    /// using Stokes law for a disc of radius `r_disc`.
    pub fn step(&mut self, dt: f64, force: Vec3, r_disc: f64, viscosity: f64) {
        // Stokes drag coefficient for a disc: γ ≈ 16 η r
        let gamma = 16.0 * viscosity * r_disc;
        for (i, &fi) in force.iter().enumerate() {
            let a = (fi - gamma * self.velocity[i]) / 1e-14; // platelet mass ~10 pg
            self.velocity[i] += a * dt;
            self.position[i] += self.velocity[i] * dt;
        }
    }

    /// ADP-mediated activation increment for ADP concentration `c` (mol m⁻³).
    ///
    /// Uses a Hill equation: `dA/dt = k_max * c^n / (K_d^n + c^n)`.
    pub fn agonist_activation_increment(&self, c_adp: f64, dt: f64) -> f64 {
        let k_max = 2.0_f64; // s⁻¹
        let k_d = 5e-6_f64; // mol m⁻³ (5 µM)
        let n = 1.5_f64;
        k_max * c_adp.powf(n) / (k_d.powf(n) + c_adp.powf(n)) * dt
    }
}

/// Platelet aggregation model tracking a clot formation patch.
#[derive(Debug, Clone)]
pub struct PlateletAggregation {
    /// Collection of platelets in the simulation.
    pub platelets: Vec<Platelet>,
    /// Clot perimeter radius (m).
    pub clot_radius: f64,
    /// Wall shear rate (s⁻¹).
    pub wall_shear_rate: f64,
    /// ADP agonist concentration (mol m⁻³).
    pub adp_concentration: f64,
}

impl PlateletAggregation {
    /// Construct the aggregation model.
    ///
    /// * `n_platelets` — number of platelet agents.
    /// * `wall_shear_rate` — hydrodynamic shear rate at the wall (s⁻¹).
    /// * `adp_concentration` — local ADP concentration (mol m⁻³).
    pub fn new(n_platelets: usize, wall_shear_rate: f64, adp_concentration: f64) -> Self {
        let platelets = (0..n_platelets)
            .map(|i| {
                let x = (i as f64) * 2e-6;
                Platelet::new([x, 0.0, 0.0], 5e13_f64)
            })
            .collect();
        Self {
            platelets,
            clot_radius: 0.0,
            wall_shear_rate,
            adp_concentration,
        }
    }

    /// Number of activated platelets.
    pub fn activated_count(&self) -> usize {
        self.platelets
            .iter()
            .filter(|p| p.state != PlateletState::Resting)
            .count()
    }

    /// Thrombus coverage area (m²): π r_clot².
    pub fn coverage_area(&self) -> f64 {
        PI * self.clot_radius.powi(2)
    }

    /// Advance the simulation by `dt`.
    ///
    /// Platelets near the activation zone become activated via ADP signalling.
    pub fn step(&mut self, dt: f64) {
        let c = self.adp_concentration;
        for p in &mut self.platelets {
            let inc = p.agonist_activation_increment(c, dt);
            p.activation_level += inc;
            if p.activation_level > 0.5 && p.state == PlateletState::Resting {
                p.activate();
            }
        }
        // Grow clot radius proportional to activated count
        let active = self.activated_count() as f64;
        self.clot_radius = (active * 2e-6 / PI).sqrt();
    }
}

// ---------------------------------------------------------------------------
// 10. Bone Remodeling (Wolff's Law)
// ---------------------------------------------------------------------------

/// Bone remodeling model based on Wolff's law and Frost's mechanostat theory.
///
/// Bone density adapts to mechanical loading through the osteoblast–osteoclast
/// coupling governed by the strain energy density (SED).
#[derive(Debug, Clone)]
pub struct BoneRemodeling {
    /// Local apparent bone density ρ (kg m⁻³).
    pub density: f64,
    /// Reference (target) strain energy density (J m⁻³).
    pub reference_sed: f64,
    /// Remodeling rate constant C (kg m⁻³ s⁻¹ per J m⁻³).
    pub remodeling_rate: f64,
    /// Dead zone half-width about the reference SED (J m⁻³).
    pub dead_zone: f64,
    /// Minimum bone density ρ_min (kg m⁻³).
    pub density_min: f64,
    /// Maximum bone density ρ_max (kg m⁻³).
    pub density_max: f64,
    /// Current von Mises stress (Pa).
    pub von_mises_stress: f64,
    /// Current Young's modulus (Pa) — Currey power law.
    pub young_modulus: f64,
}

impl BoneRemodeling {
    /// Create a new bone remodeling site.
    ///
    /// * `initial_density` — initial apparent density (kg m⁻³). Cortical ≈ 1800.
    /// * `reference_sed` — homeostatic strain energy density (J m⁻³).
    /// * `remodeling_rate` — remodeling coefficient C (kg m⁻³ s⁻¹ / (J m⁻³)).
    /// * `dead_zone` — lazy zone half-width (J m⁻³).
    pub fn new(
        initial_density: f64,
        reference_sed: f64,
        remodeling_rate: f64,
        dead_zone: f64,
    ) -> Self {
        let young_modulus = Self::currey_modulus(initial_density);
        Self {
            density: initial_density,
            reference_sed,
            remodeling_rate,
            dead_zone,
            density_min: 0.01,
            density_max: 2100.0,
            von_mises_stress: 0.0,
            young_modulus,
        }
    }

    /// Currey's power law: `E = c * ρ^n` (Pa).
    /// Uses `c = 7.88×10⁻³` Pa and `n = 3.09` (Currey 1988).
    pub fn currey_modulus(density: f64) -> f64 {
        7.88e-3 * density.powf(3.09)
    }

    /// Strain energy density for applied stress `σ` (J m⁻³): `U = σ²/(2E)`.
    pub fn strain_energy_density(&self, stress: f64) -> f64 {
        if self.young_modulus > 1e-10 {
            stress * stress / (2.0 * self.young_modulus)
        } else {
            0.0
        }
    }

    /// Remodeling stimulus: `e_sed - e_ref` outside the dead zone, else zero.
    pub fn remodeling_stimulus(&self, sed: f64) -> f64 {
        let delta = sed - self.reference_sed;
        if delta.abs() <= self.dead_zone {
            0.0
        } else if delta > 0.0 {
            delta - self.dead_zone
        } else {
            delta + self.dead_zone
        }
    }

    /// Advance bone remodeling by `dt` given an applied `stress` (Pa).
    ///
    /// Updates density, Young's modulus, and von Mises stress.
    pub fn step(&mut self, dt: f64, stress: f64) {
        self.von_mises_stress = stress;
        let sed = self.strain_energy_density(stress);
        let stimulus = self.remodeling_stimulus(sed);
        self.density += self.remodeling_rate * stimulus * dt;
        self.density = clamp(self.density, self.density_min, self.density_max);
        self.young_modulus = Self::currey_modulus(self.density);
    }

    /// Osteoblast activity index (0–1) proportional to anabolic stimulus.
    pub fn osteoblast_activity(&self, sed: f64) -> f64 {
        let s = self.remodeling_stimulus(sed);
        clamp(s / (self.reference_sed + 1.0), 0.0, 1.0)
    }

    /// Osteoclast activity index (0–1) proportional to catabolic stimulus.
    pub fn osteoclast_activity(&self, sed: f64) -> f64 {
        let s = -self.remodeling_stimulus(sed);
        clamp(s / (self.reference_sed + 1.0), 0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // --- CorticalTension ---

    #[test]
    fn test_cortical_tension_young_laplace() {
        let ct = CorticalTension::new(5e-4, 10e-6, 2e-19, 0.0);
        let dp = ct.young_laplace_pressure();
        // ΔP = 2 * 5e-4 / 10e-6 = 100 Pa
        assert!((dp - 100.0).abs() < 1e-6, "dp={dp}");
    }

    #[test]
    fn test_cortical_tension_helfrich_energy_positive() {
        let ct = CorticalTension::new(5e-4, 10e-6, 2e-19, 0.0);
        let e = ct.helfrich_energy();
        assert!(e >= 0.0, "Helfrich energy must be non-negative, got {e}");
    }

    #[test]
    fn test_cortical_tension_equatorial_force() {
        let ct = CorticalTension::new(1.0, 1.0, 1e-19, 0.0);
        let f = ct.equatorial_line_force();
        assert!((f - PI).abs() < EPS, "f={f}");
    }

    #[test]
    fn test_cortical_tension_equilibrium_radius() {
        let ct = CorticalTension::new(2.0, 1.0, 1e-19, 0.0);
        let r = ct.equilibrium_radius(4.0); // R = 2*2/4 = 1
        assert!((r - 1.0).abs() < EPS, "r={r}");
    }

    #[test]
    fn test_cortical_tension_spring_constant() {
        let ct = CorticalTension::new(1.0, 2.0, 1e-19, 0.0);
        let k = ct.radial_spring_constant();
        assert!((k - 0.5).abs() < EPS, "k={k}"); // 2*1/4 = 0.5
    }

    // --- CellDivision ---

    #[test]
    fn test_cell_division_phase_progression() {
        let durations = [100.0, 200.0, 50.0, 50.0];
        let mut cell = CellDivision::new([0.0; 3], 10e-6, durations, 1e-3);
        assert_eq!(cell.phase, CellCyclePhase::G1);
        // Advance past G1
        for _ in 0..110 {
            cell.step(1.0);
        }
        assert_eq!(cell.phase, CellCyclePhase::S);
    }

    #[test]
    fn test_cell_division_produces_daughters() {
        let durations = [0.1, 0.1, 0.1, 0.1];
        let mut cell = CellDivision::new([0.0, 0.0, 0.0], 10e-6, durations, 1e-3);
        let mut daughters: Option<(Vec3, Vec3)> = None;
        for _ in 0..100 {
            if let Some(d) = cell.step(0.05) {
                daughters = Some(d);
                break;
            }
        }
        assert!(daughters.is_some(), "Cell should have divided");
        let (da, db) = daughters.unwrap();
        // Daughters should be separated
        assert!(dist3(da, db) > 1e-10, "Daughters should be separated");
    }

    #[test]
    fn test_cell_division_contractile_force_positive() {
        let mut cell = CellDivision::new([0.0; 3], 10e-6, [0.1, 0.1, 0.1, 1.0], 1e-3);
        // Advance into M phase
        for _ in 0..30 {
            cell.step(0.1);
        }
        let f = cell.contractile_force();
        assert!(f >= 0.0, "Contractile force must be non-negative: {f}");
    }

    // --- Filopodium ---

    #[test]
    fn test_filopodium_protrudes() {
        let mut filo = Filopodium::new([0.0; 3], [1.0, 0.0, 0.0], 1e-6, 1e-7, 5e-8, 1e-11);
        let len_before = filo.length();
        filo.protrude(1.0);
        assert!(filo.length() > len_before, "Filopodium should grow");
    }

    #[test]
    fn test_filopodium_retracts() {
        let mut filo = Filopodium::new([0.0; 3], [1.0, 0.0, 0.0], 1e-6, 1e-7, 5e-8, 1e-11);
        let len_before = filo.length();
        filo.retract(1.0);
        assert!(filo.length() <= len_before, "Filopodium should retract");
    }

    #[test]
    fn test_filopodium_traction_zero_unadhered() {
        let filo = Filopodium::new([0.0; 3], [1.0, 0.0, 0.0], 1e-6, 1e-7, 5e-8, 1e-11);
        let f = filo.traction_force();
        assert!(
            norm3(f) < EPS,
            "Unadhered filopodium traction should be zero"
        );
    }

    // --- Lamellipodium ---

    #[test]
    fn test_lamellipodium_edge_advances() {
        let mut lam = Lamellipodium::new(0.0, 5e-6, 1e-7, 3e-8, 100.0);
        lam.step(10.0);
        assert!(lam.leading_edge > 0.0, "Edge should advance");
    }

    #[test]
    fn test_lamellipodium_net_velocity() {
        let lam = Lamellipodium::new(0.0, 5e-6, 1e-7, 3e-8, 100.0);
        let v = lam.net_velocity();
        assert!((v - (1e-7 - 3e-8)).abs() < EPS, "v={v}");
    }

    // --- CadherinBond ---

    #[test]
    fn test_cadherin_bond_force_at_rest() {
        let bond = CadherinBond::new(10e-9, 1e-3, 0.1, 5e-12);
        let f = bond.force(10e-9);
        assert!(f.abs() < EPS, "At rest length force should be zero: {f}");
    }

    #[test]
    fn test_cadherin_bond_force_stretched() {
        let bond = CadherinBond::new(10e-9, 1e-3, 0.1, 5e-12);
        let f = bond.force(20e-9);
        assert!(f > 0.0, "Stretched bond should have positive force");
    }

    #[test]
    fn test_cadherin_bond_survival_probability_unity_at_zero_force() {
        let bond = CadherinBond::new(10e-9, 1e-3, 1e-3, 5e-12);
        // At zero force, k_off = k_off_zero * e^0 = k_off_zero
        let p = bond.survival_probability(0.0, 0.0); // dt=0 → P=1
        assert!((p - 1.0).abs() < EPS, "P at dt=0 should be 1: {p}");
    }

    // --- CellAdhesion ---

    #[test]
    fn test_cell_adhesion_active_bonds_initial() {
        let adh = CellAdhesion::new([0.0; 3], [1e-5, 0.0, 0.0], 10, 10e-9, 1e-3, 0.1, 5e-12);
        assert_eq!(adh.active_bond_count(), 10);
    }

    #[test]
    fn test_cell_adhesion_force_direction() {
        let adh = CellAdhesion::new([0.0; 3], [1e-5, 0.0, 0.0], 5, 10e-9, 1e-3, 0.1, 5e-12);
        let f = adh.force_on_a();
        // Force on A should be in +x direction
        assert!(
            f[0] > 0.0 || adh.total_force() == 0.0,
            "Force on A should point towards B"
        );
    }

    // --- ExtracellularMatrix ---

    #[test]
    fn test_ecm_effective_stiffness_unity_at_zero_strain() {
        let ecm = ExtracellularMatrix::new(1000.0, 0.3, 10.0, 1.0, 0.1);
        let e = ecm.effective_stiffness();
        assert!((e - 1000.0).abs() < EPS, "At zero strain E_eff = E0: {e}");
    }

    #[test]
    fn test_ecm_strain_stiffening() {
        let mut ecm = ExtracellularMatrix::new(1000.0, 0.3, 10.0, 2.0, 0.0);
        ecm.strain = 0.5;
        let e = ecm.effective_stiffness();
        let expected = 1000.0 * (2.0_f64 * 0.5).exp();
        assert!(
            (e - expected).abs() < 1e-6,
            "Strain stiffening mismatch: {e} vs {expected}"
        );
    }

    #[test]
    fn test_ecm_relaxation_decreases_strain() {
        let mut ecm = ExtracellularMatrix::new(1000.0, 0.3, 10.0, 1.0, 0.0);
        ecm.strain = 0.1;
        let before = ecm.strain;
        ecm.relax(1.0);
        assert!(ecm.strain < before, "Strain should relax");
    }

    #[test]
    fn test_ecm_composite_stiffness_exceeds_matrix() {
        let ecm = ExtracellularMatrix::new(1000.0, 0.3, 10.0, 1.0, 0.3);
        let e_c = ecm.composite_stiffness();
        assert!(
            e_c > ecm.young_modulus,
            "Composite stiffness should exceed matrix: {e_c}"
        );
    }

    // --- WoundHealing ---

    #[test]
    fn test_wound_healing_closes() {
        let mut wh = WoundHealing::new(100e-6, 1e-4, 1e-4, 2e-7);
        for _ in 0..1000 {
            wh.step(1.0, 1e3);
        }
        assert!(
            wh.wound_width < 100e-6,
            "Wound should have partially closed"
        );
    }

    #[test]
    fn test_wound_healing_closure_fraction_increases() {
        let mut wh = WoundHealing::new(100e-6, 1e-4, 1e-4, 2e-7);
        let f0 = wh.closure_fraction();
        wh.step(100.0, 1e3);
        let f1 = wh.closure_fraction();
        assert!(f1 >= f0, "Closure fraction should not decrease");
    }

    #[test]
    fn test_wound_healing_estimated_closure_time_positive() {
        let wh = WoundHealing::new(100e-6, 1e-4, 1e-4, 2e-7);
        let t = wh.estimated_closure_time();
        assert!(
            t > 0.0 && t.is_finite(),
            "Closure time should be positive finite: {t}"
        );
    }

    // --- TumorGrowth ---

    #[test]
    fn test_tumor_growth_volume_positive() {
        let tg = TumorGrowth::new(1e-3, 1e-5, 100.0, 1000.0);
        assert!(tg.volume() > 0.0, "Volume must be positive");
    }

    #[test]
    fn test_tumor_growth_step_increases_radius() {
        let mut tg = TumorGrowth::new(1e-4, 1e-5, 10.0, 100.0);
        let r0 = tg.radius;
        tg.step(3600.0, 1e10); // large max stress → no inhibition
        assert!(tg.radius >= r0, "Tumor should grow");
    }

    #[test]
    fn test_tumor_growth_doubling_time_consistent() {
        let tg = TumorGrowth::new(1e-3, 1.0, 100.0, 0.0);
        let td = tg.doubling_time();
        assert!((td - 2_f64.ln()).abs() < 1e-10, "Doubling time: {td}");
    }

    // --- RedBloodCell ---

    #[test]
    fn test_rbc_elastic_energy_non_negative() {
        let nodes = vec![[0.0_f64; 3], [1e-6, 0.0, 0.0]];
        let springs = vec![SpectrinSpring::new([0, 1], 1e-6, 20e-9, 75e-9)];
        let rbc = RedBloodCell::new(nodes, springs, 2e-19, 4.8e-6, 135e-12, 94e-18);
        let e = rbc.elastic_energy();
        assert!(e >= 0.0, "Elastic energy must be non-negative: {e}");
    }

    #[test]
    fn test_rbc_area_constraint_zero_at_target() {
        let rbc = RedBloodCell::new(vec![], vec![], 2e-19, 4.8e-6, 135e-12, 94e-18);
        let e = rbc.area_constraint_energy(135e-12);
        assert!(
            e.abs() < EPS,
            "Area constraint energy should be zero at target: {e}"
        );
    }

    #[test]
    fn test_rbc_wlc_force_increases_with_extension() {
        let s = SpectrinSpring::new([0, 1], 40e-9, 20e-9, 75e-9);
        let f1 = s.wlc_force(10e-9);
        let f2 = s.wlc_force(60e-9);
        assert!(
            f2 > f1,
            "WLC force should increase with extension: f1={f1}, f2={f2}"
        );
    }

    // --- PlateletAggregation ---

    #[test]
    fn test_platelet_initially_resting() {
        let pa = PlateletAggregation::new(5, 1000.0, 1e-6);
        assert_eq!(pa.activated_count(), 0);
    }

    #[test]
    fn test_platelet_activation_step() {
        let mut pa = PlateletAggregation::new(5, 1000.0, 1e-3);
        pa.step(10.0);
        // At high ADP some platelets should activate
        assert!(
            pa.activated_count() > 0,
            "Some platelets should activate at high ADP"
        );
    }

    #[test]
    fn test_platelet_coverage_area_grows() {
        let mut pa = PlateletAggregation::new(20, 1000.0, 1e-3);
        pa.step(1.0);
        pa.step(5.0);
        let a1 = pa.coverage_area();
        pa.step(20.0);
        let a2 = pa.coverage_area();
        assert!(a2 >= a1, "Coverage area should grow: a1={a1}, a2={a2}");
    }

    // --- BoneRemodeling ---

    #[test]
    fn test_bone_remodeling_modulus_positive() {
        let br = BoneRemodeling::new(1800.0, 5000.0, 1e-3, 500.0);
        assert!(br.young_modulus > 0.0, "Modulus must be positive");
    }

    #[test]
    fn test_bone_remodeling_density_increases_under_load() {
        let mut br = BoneRemodeling::new(1000.0, 1000.0, 1.0, 100.0);
        let rho0 = br.density;
        // Apply high stress to produce SED >> reference_sed
        for _ in 0..100 {
            br.step(1.0, 1e8);
        }
        assert!(
            br.density >= rho0,
            "Density should increase under high load"
        );
    }

    #[test]
    fn test_bone_remodeling_density_bounded() {
        let mut br = BoneRemodeling::new(1800.0, 1000.0, 100.0, 0.0);
        for _ in 0..10000 {
            br.step(1.0, 1e9);
        }
        assert!(
            br.density <= br.density_max,
            "Density must not exceed maximum"
        );
        assert!(
            br.density >= br.density_min,
            "Density must not go below minimum"
        );
    }

    #[test]
    fn test_bone_remodeling_dead_zone_no_change() {
        let mut br = BoneRemodeling::new(1000.0, 5000.0, 1.0, 1e9);
        let rho0 = br.density;
        // Stress giving SED = 5000 exactly, within dead zone
        br.step(1.0, 0.0); // zero stress → SED=0, below dead zone only if dead_zone > reference
        // Since dead_zone = 1e9 >> reference_sed = 5000, zero stress is in dead zone
        // density should be unchanged
        assert!(
            (br.density - rho0).abs() < EPS,
            "In dead zone density should not change"
        );
    }

    #[test]
    fn test_bone_remodeling_stimulus_outside_dead_zone() {
        let br = BoneRemodeling::new(1000.0, 1000.0, 1.0, 100.0);
        let stim = br.remodeling_stimulus(2000.0); // 2000 >> 1000 + 100
        assert!(
            stim > 0.0,
            "Stimulus should be positive above dead zone: {stim}"
        );
    }
}
