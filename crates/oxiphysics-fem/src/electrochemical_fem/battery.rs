// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Pseudo-2D (P2D) Doyle-Fuller-Newman lithium-ion battery model.
//!
//! The cell is discretised with cell-centred finite volumes along the
//! through-thickness coordinate x (negative electrode | separator | positive
//! electrode) and with node-centred spherical finite volumes along the particle
//! radius r at every electrode cell.  One time step of [`BatteryFem::step`]
//! performs a linearly implicit operator split:
//!
//! 1. charge conservation in the solid (φ_s) and electrolyte (φ_e) phases,
//!    coupled through Butler-Volmer kinetics, solved by damped Newton
//!    iteration at the start-of-step concentrations;
//! 2. implicit solid-phase diffusion in every particle, driven by the
//!    interfacial current from step 1;
//! 3. implicit electrolyte diffusion with the (1 − t₊) a_s i_n / F source from
//!    the same interfacial current.
//!
//! Because both concentration updates use the same converged interfacial
//! current and both finite-volume schemes are exactly conservative, the total
//! lithium inventory ([`BatteryFem::total_lithium`]) is conserved up to the
//! Newton residual.

use super::{ButlerVolmerFem, FARADAY, NernstPlanckFem, R_GAS};
use crate::error::{Error, Result};

/// Maximum number of Newton iterations for the potential sub-problem.
const MAX_NEWTON_ITERATIONS: usize = 60;
/// Maximum number of step halvings in the Newton line search.
const MAX_LINE_SEARCH_HALVINGS: usize = 40;

// ============================================================================
// Domain descriptors
// ============================================================================

/// Domain type for the DFN model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DfnDomain {
    /// Negative electrode (anode).
    NegativeElectrode,
    /// Separator.
    Separator,
    /// Positive electrode (cathode).
    PositiveElectrode,
}

/// Parameters for one electrode domain in the P2D model.
#[derive(Debug, Clone)]
pub struct ElectrodeDomain {
    /// Domain type.
    pub domain: DfnDomain,
    /// Thickness (m).
    pub thickness: f64,
    /// Number of macro-scale (x-direction) finite-volume cells.
    pub n_macro: usize,
    /// Number of micro-scale (r-direction) nodes in each particle (≥ 2).
    pub n_micro: usize,
    /// Solid-phase Li diffusivity D_s (m²/s).
    pub diffusivity_solid: f64,
    /// Maximum solid concentration c_s_max (mol/m³).
    pub c_solid_max: f64,
    /// Particle radius r_p (m).
    pub particle_radius: f64,
    /// Electrode (electrolyte-filled) porosity ε_l.
    pub porosity: f64,
    /// Bruggeman exponent p.
    pub bruggeman: f64,
    /// Solid-phase conductivity σ (S/m).
    pub sigma_solid: f64,
    /// Exchange current density j₀_ref (A/m²) at c_s = c_s_max/2 and
    /// c_e = [`BatteryFem::c_electrolyte_ref`].
    pub j0_ref: f64,
}

impl ElectrodeDomain {
    /// Create a default negative electrode (graphite) parameter set.
    pub fn default_negative() -> Self {
        Self {
            domain: DfnDomain::NegativeElectrode,
            thickness: 100e-6,
            n_macro: 10,
            n_micro: 5,
            diffusivity_solid: 3.9e-14,
            c_solid_max: 31370.0,
            particle_radius: 2.0e-6,
            porosity: 0.3,
            bruggeman: 1.5,
            sigma_solid: 100.0,
            j0_ref: 1.0,
        }
    }

    /// Create a default positive electrode (LiCoO₂) parameter set.
    pub fn default_positive() -> Self {
        Self {
            domain: DfnDomain::PositiveElectrode,
            thickness: 80e-6,
            n_macro: 10,
            n_micro: 5,
            diffusivity_solid: 1.0e-14,
            c_solid_max: 51218.0,
            particle_radius: 3.0e-6,
            porosity: 0.3,
            bruggeman: 1.5,
            sigma_solid: 3.8,
            j0_ref: 1.0,
        }
    }

    /// Effective electrolyte diffusivity D_e_eff = D_e * ε^p (m²/s).
    pub fn effective_electrolyte_diffusivity(&self, d_electrolyte: f64) -> f64 {
        d_electrolyte * self.porosity.powf(self.bruggeman)
    }

    /// Effective electrolyte conductivity κ_eff = κ * ε^p (S/m).
    pub fn effective_conductivity(&self, kappa: f64) -> f64 {
        kappa * self.porosity.powf(self.bruggeman)
    }

    /// Active-material (solid) volume fraction ε_s = 1 − ε_l.
    pub fn solid_volume_fraction(&self) -> f64 {
        1.0 - self.porosity
    }

    /// Specific interfacial area a_s = 3 ε_s / r_p (1/m).
    pub fn specific_interfacial_area(&self) -> f64 {
        3.0 * self.solid_volume_fraction() / self.particle_radius
    }

    /// Effective solid-phase conductivity σ_eff = σ ε_s^p (S/m).
    pub fn effective_solid_conductivity(&self) -> f64 {
        self.sigma_solid * self.solid_volume_fraction().powf(self.bruggeman)
    }

    /// Width of one macro-scale finite-volume cell (m).
    pub fn cell_width(&self) -> f64 {
        self.thickness / self.n_macro as f64
    }

    /// Open-circuit potential U(θ) of this electrode at surface stoichiometry
    /// θ = c_s / c_s_max (V vs. Li/Li⁺).
    ///
    /// The empirical correlations of [`ButlerVolmerFem`] are parameterised by
    /// the electrode state of charge: θ for the graphite anode and 1 − θ for
    /// the LiCoO₂ cathode (a delithiated cathode is fully charged).
    pub fn open_circuit_potential(&self, theta: f64) -> f64 {
        match self.domain {
            DfnDomain::PositiveElectrode => ButlerVolmerFem::ocp_licoo2(1.0 - theta),
            DfnDomain::NegativeElectrode | DfnDomain::Separator => {
                ButlerVolmerFem::ocp_graphite(theta)
            }
        }
    }
}

// ============================================================================
// Spherical finite-volume helpers
// ============================================================================

/// Control-volume "volumes" V_i / (4π) of the node-centred spherical grid.
///
/// Nodes sit at r_i = i h (i = 0..n), h = R/(n − 1); node i owns
/// \[max(r_i − h/2, 0), min(r_i + h/2, R)\].  The volumes tile the sphere
/// exactly: Σ V_i = R³/3.
fn spherical_control_volumes(n: usize, radius: f64) -> Vec<f64> {
    let h = radius / (n as f64 - 1.0);
    (0..n)
        .map(|i| {
            let r = i as f64 * h;
            let lo = (r - 0.5 * h).max(0.0);
            let hi = (r + 0.5 * h).min(radius);
            (hi.powi(3) - lo.powi(3)) / 3.0
        })
        .collect()
}

/// Volume-weighted particle-average concentration.
fn particle_average(c: &[f64], volumes: &[f64], radius: f64) -> f64 {
    let total: f64 = c.iter().zip(volumes).map(|(ci, vi)| ci * vi).sum();
    total / (radius.powi(3) / 3.0)
}

/// Dense LU solve with partial pivoting (in place): on return `b` holds x.
fn solve_dense_in_place(a: &mut [f64], b: &mut [f64], n: usize) -> Result<()> {
    if a.len() != n * n || b.len() != n {
        return Err(Error::InvalidInput(format!(
            "dense solve: matrix {} / rhs {} inconsistent with n = {n}",
            a.len(),
            b.len()
        )));
    }
    let scale = a.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    if !(scale.is_finite() && scale > 0.0) {
        return Err(Error::SingularSystem);
    }
    let tiny = scale * f64::EPSILON * n as f64;
    for col in 0..n {
        let mut piv = col;
        let mut piv_abs = a[col * n + col].abs();
        for row in col + 1..n {
            let v = a[row * n + col].abs();
            if v > piv_abs {
                piv = row;
                piv_abs = v;
            }
        }
        if piv_abs <= tiny {
            return Err(Error::SingularSystem);
        }
        if piv != col {
            for j in 0..n {
                a.swap(col * n + j, piv * n + j);
            }
            b.swap(col, piv);
        }
        let pivot = a[col * n + col];
        for row in col + 1..n {
            let factor = a[row * n + col] / pivot;
            if factor == 0.0 {
                continue;
            }
            for j in col..n {
                a[row * n + j] -= factor * a[col * n + j];
            }
            b[row] -= factor * b[col];
        }
    }
    for i in (0..n).rev() {
        let mut sum = b[i];
        for j in i + 1..n {
            sum -= a[i * n + j] * b[j];
        }
        b[i] = sum / a[i * n + i];
    }
    Ok(())
}

/// Conductance of the face between two finite-volume cells (harmonic mean).
fn face_conductance(h_left: f64, k_left: f64, h_right: f64, k_right: f64) -> f64 {
    1.0 / (0.5 * h_left / k_left + 0.5 * h_right / k_right)
}

// ============================================================================
// BatteryFem
// ============================================================================

/// Per-cell electrolyte properties along x (cell-centred).
struct ElectrolyteCells {
    width: Vec<f64>,
    porosity: Vec<f64>,
    kappa_eff: Vec<f64>,
    diff_eff: Vec<f64>,
    /// Specific interfacial area (0 in the separator).
    area: Vec<f64>,
}

/// Pseudo-2D (P2D) Doyle-Fuller-Newman lithium-ion battery model.
///
/// Implements the DFN model with:
/// - macro-scale 1-D finite volumes along the cell thickness (x direction),
///   x = 0 at the negative and x = L at the positive current collector;
/// - micro-scale 1-D spherical finite volumes for solid diffusion inside the
///   particles (r direction), ∂c_s/∂t = D_s r⁻² ∂_r(r² ∂_r c_s) with
///   −D_s ∂_r c_s|_R = i_n / F;
/// - electrolyte diffusion ε ∂c_e/∂t = ∂_x(D_eff ∂_x c_e) + (1 − t₊) a_s i_n / F;
/// - solid charge conservation ∂_x(σ_eff ∂_x φ_s) = a_s i_n with i_s = I at
///   both collectors (φ_s = 0 at the negative collector);
/// - electrolyte charge conservation
///   ∂_x(κ_eff ∂_x φ_e) − ∂_x(κ_D ∂_x ln c_e) = −a_s i_n, κ_D = 2RTκ_eff(1 − t₊)/F;
/// - Butler-Volmer kinetics i_n = i₀ \[exp(α_a f η) − exp(−α_c f η)\] with
///   η = φ_s − φ_e − U(θ_surf) and i₀ = j₀_ref √(c_e/c_e,ref) · 2√(θ(1 − θ)).
///
/// Index 0 of the per-electrode vectors is the negative and index 1 the
/// positive electrode.  The applied current density is positive on discharge.
#[derive(Debug, Clone)]
pub struct BatteryFem {
    /// Negative electrode parameters.
    pub neg: ElectrodeDomain,
    /// Positive electrode parameters.
    pub pos: ElectrodeDomain,
    /// Separator thickness (m).
    pub separator_thickness: f64,
    /// Number of finite-volume cells in the separator.
    pub n_sep: usize,
    /// Separator porosity.
    pub separator_porosity: f64,
    /// Separator Bruggeman exponent.
    pub separator_bruggeman: f64,
    /// Electrolyte diffusivity D_e (m²/s).
    pub d_electrolyte: f64,
    /// Electrolyte conductivity κ (S/m).
    pub kappa_electrolyte: f64,
    /// Li⁺ transference number t₊.
    pub transference: f64,
    /// Reference electrolyte concentration for the exchange current (mol/m³).
    pub c_electrolyte_ref: f64,
    /// Anodic charge-transfer coefficient α_a.
    pub alpha_anodic: f64,
    /// Cathodic charge-transfer coefficient α_c.
    pub alpha_cathodic: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Applied current density (A/m²), positive = discharge.
    pub applied_current: f64,
    /// Solid concentration: `c_solid[electrode][x_cell][r_node]` (mol/m³).
    pub c_solid: Vec<Vec<Vec<f64>>>,
    /// Electrolyte concentration per x cell (mol/m³).
    pub c_electrolyte: Vec<f64>,
    /// Solid-phase potential per electrode x cell (V).
    pub phi_solid: Vec<Vec<f64>>,
    /// Electrolyte potential per x cell (V).
    pub phi_electrolyte: Vec<f64>,
    /// Interfacial (Butler-Volmer) current density per electrode x cell (A/m²),
    /// positive for de-intercalation.
    pub reaction_current: Vec<Vec<f64>>,
    /// Additional area-specific contact resistance of the cell (Ω·m²).
    ///
    /// Terminal voltage = φ_s(L) − φ_s(0) − `internal_resistance` × `applied_current`.
    pub internal_resistance: f64,
    /// Current time (s).
    pub time: f64,
    /// Time step (s).
    pub dt: f64,
}

impl BatteryFem {
    /// Construct a new P2D/DFN battery model with default parameters.
    ///
    /// The cell starts at rest: uniform concentrations, η = 0 everywhere
    /// (φ_s,neg = 0, φ_e = −U_n, φ_s,pos = U_p − U_n) and a uniform reaction
    /// distribution consistent with the applied current.
    pub fn new(applied_current: f64, temperature: f64, dt: f64) -> Self {
        let neg = ElectrodeDomain::default_negative();
        let pos = ElectrodeDomain::default_positive();
        let n_sep = 5_usize;
        let n_total = neg.n_macro + n_sep + pos.n_macro;
        let theta_n = 0.5;
        let theta_p = 0.8;
        let c_solid = vec![
            vec![vec![neg.c_solid_max * theta_n; neg.n_micro]; neg.n_macro],
            vec![vec![pos.c_solid_max * theta_p; pos.n_micro]; pos.n_macro],
        ];
        let u_n = neg.open_circuit_potential(theta_n);
        let u_p = pos.open_circuit_potential(theta_p);
        let i_neg = applied_current / (neg.specific_interfacial_area() * neg.thickness);
        let i_pos = -applied_current / (pos.specific_interfacial_area() * pos.thickness);
        let reaction_current = vec![vec![i_neg; neg.n_macro], vec![i_pos; pos.n_macro]];
        let phi_solid = vec![vec![0.0; neg.n_macro], vec![u_p - u_n; pos.n_macro]];
        Self {
            neg,
            pos,
            separator_thickness: 25e-6,
            n_sep,
            separator_porosity: 0.5,
            separator_bruggeman: 1.5,
            d_electrolyte: 2.7e-10,
            kappa_electrolyte: 1.0,
            transference: 0.363,
            c_electrolyte_ref: 1000.0,
            alpha_anodic: 0.5,
            alpha_cathodic: 0.5,
            temperature,
            applied_current,
            c_solid,
            c_electrolyte: vec![1000.0; n_total],
            phi_solid,
            phi_electrolyte: vec![-u_n; n_total],
            reaction_current,
            internal_resistance: 0.0,
            time: 0.0,
            dt,
        }
    }

    /// Compute total number of macro-scale cells.
    pub fn n_total_nodes(&self) -> usize {
        self.neg.n_macro + self.n_sep + self.pos.n_macro
    }

    fn electrode(&self, electrode_idx: usize) -> Result<&ElectrodeDomain> {
        match electrode_idx {
            0 => Ok(&self.neg),
            1 => Ok(&self.pos),
            _ => Err(Error::InvalidInput(format!(
                "electrode index {electrode_idx} (expected 0 = negative, 1 = positive)"
            ))),
        }
    }

    /// Volume-averaged lithium stoichiometry c̄_s / c_s,max of one electrode.
    fn mean_stoichiometry(&self, electrode_idx: usize) -> f64 {
        let Ok(dom) = self.electrode(electrode_idx) else {
            return 0.0;
        };
        let Some(cells) = self.c_solid.get(electrode_idx) else {
            return 0.0;
        };
        if cells.is_empty() || dom.n_micro < 2 {
            return 0.0;
        }
        let volumes = spherical_control_volumes(dom.n_micro, dom.particle_radius);
        let sum: f64 = cells
            .iter()
            .map(|c| particle_average(c, &volumes, dom.particle_radius))
            .sum();
        sum / cells.len() as f64 / dom.c_solid_max
    }

    /// Volume-averaged lithium stoichiometry of the negative electrode
    /// (c̄_s / c_s,max, equal to its state of charge).
    pub fn soc_negative(&self) -> f64 {
        self.mean_stoichiometry(0)
    }

    /// Volume-averaged lithium stoichiometry c̄_s / c_s,max of the positive
    /// electrode (its state of charge is 1 minus this value).
    pub fn soc_positive(&self) -> f64 {
        self.mean_stoichiometry(1)
    }

    /// Total lithium inventory per unit electrode area (mol/m²): lithium in
    /// the particles of both electrodes plus Li⁺ in the electrolyte.
    pub fn total_lithium(&self) -> f64 {
        let mut total = 0.0;
        for (e, dom) in [&self.neg, &self.pos].into_iter().enumerate() {
            if dom.n_micro < 2 {
                continue;
            }
            let volumes = spherical_control_volumes(dom.n_micro, dom.particle_radius);
            let weight = dom.solid_volume_fraction() * dom.cell_width();
            if let Some(cells) = self.c_solid.get(e) {
                total += cells
                    .iter()
                    .map(|c| weight * particle_average(c, &volumes, dom.particle_radius))
                    .sum::<f64>();
            }
        }
        let props = self.electrolyte_cells();
        total += self
            .c_electrolyte
            .iter()
            .zip(props.porosity.iter().zip(&props.width))
            .map(|(c, (eps, h))| c * eps * h)
            .sum::<f64>();
        total
    }

    /// Butler-Volmer overpotential η = φ_s − φ_e − U(θ_surf) of one electrode
    /// cell (V), or `None` for an out-of-range index.
    pub fn overpotential(&self, electrode_idx: usize, cell: usize) -> Option<f64> {
        let dom = self.electrode(electrode_idx).ok()?;
        let x_cell = self.electrolyte_index(electrode_idx, cell);
        let theta = self.c_solid.get(electrode_idx)?.get(cell)?.last()? / dom.c_solid_max;
        let phi_s = *self.phi_solid.get(electrode_idx)?.get(cell)?;
        let phi_e = *self.phi_electrolyte.get(x_cell)?;
        Some(phi_s - phi_e - dom.open_circuit_potential(theta))
    }

    /// Cell terminal voltage (V).
    ///
    /// V_cell = φ_s(x = L) − φ_s(x = 0) − R_int · I, with φ_s(0) = 0 and the
    /// positive-collector potential extrapolated from the last cell centre.
    /// After a [`Self::step`] / [`Self::solve_potentials`] this includes the
    /// open-circuit potentials at the particle surfaces, both kinetic
    /// overpotentials and the ohmic and concentration polarisation of the
    /// electrolyte and solid phases; before the first solve it reflects the
    /// open-circuit initial state.
    pub fn terminal_voltage(&self) -> f64 {
        let sigma_p = self.pos.effective_solid_conductivity();
        let phi_last = self
            .phi_solid
            .get(1)
            .and_then(|v| v.last())
            .copied()
            .unwrap_or(0.0);
        let phi_collector = phi_last - self.applied_current * 0.5 * self.pos.cell_width() / sigma_p;
        phi_collector - self.applied_current * self.internal_resistance
    }

    /// Index of electrode cell `cell` in the global electrolyte (x) cell list.
    fn electrolyte_index(&self, electrode_idx: usize, cell: usize) -> usize {
        if electrode_idx == 0 {
            cell
        } else {
            self.neg.n_macro + self.n_sep + cell
        }
    }

    fn electrolyte_cells(&self) -> ElectrolyteCells {
        let n_t = self.n_total_nodes();
        let mut props = ElectrolyteCells {
            width: Vec::with_capacity(n_t),
            porosity: Vec::with_capacity(n_t),
            kappa_eff: Vec::with_capacity(n_t),
            diff_eff: Vec::with_capacity(n_t),
            area: Vec::with_capacity(n_t),
        };
        let h_sep = if self.n_sep > 0 {
            self.separator_thickness / self.n_sep as f64
        } else {
            0.0
        };
        let sep_factor = self.separator_porosity.powf(self.separator_bruggeman);
        for _ in 0..self.neg.n_macro {
            props.width.push(self.neg.cell_width());
            props.porosity.push(self.neg.porosity);
            props
                .kappa_eff
                .push(self.neg.effective_conductivity(self.kappa_electrolyte));
            props.diff_eff.push(
                self.neg
                    .effective_electrolyte_diffusivity(self.d_electrolyte),
            );
            props.area.push(self.neg.specific_interfacial_area());
        }
        for _ in 0..self.n_sep {
            props.width.push(h_sep);
            props.porosity.push(self.separator_porosity);
            props.kappa_eff.push(self.kappa_electrolyte * sep_factor);
            props.diff_eff.push(self.d_electrolyte * sep_factor);
            props.area.push(0.0);
        }
        for _ in 0..self.pos.n_macro {
            props.width.push(self.pos.cell_width());
            props.porosity.push(self.pos.porosity);
            props
                .kappa_eff
                .push(self.pos.effective_conductivity(self.kappa_electrolyte));
            props.diff_eff.push(
                self.pos
                    .effective_electrolyte_diffusivity(self.d_electrolyte),
            );
            props.area.push(self.pos.specific_interfacial_area());
        }
        props
    }

    /// Check array shapes and parameter validity before a solve.
    fn validate(&self) -> Result<()> {
        if !(self.dt.is_finite() && self.dt > 0.0) {
            return Err(Error::InvalidInput(format!(
                "time step must be positive, got {}",
                self.dt
            )));
        }
        if !(self.temperature.is_finite() && self.temperature > 0.0) {
            return Err(Error::InvalidInput(format!(
                "temperature must be positive, got {}",
                self.temperature
            )));
        }
        if self.n_sep > 0 && self.separator_thickness <= 0.0 {
            return Err(Error::InvalidInput(
                "separator thickness must be positive".to_owned(),
            ));
        }
        for (e, dom) in [&self.neg, &self.pos].into_iter().enumerate() {
            if dom.n_macro == 0 || dom.n_micro < 2 {
                return Err(Error::InvalidInput(format!(
                    "electrode {e}: need n_macro ≥ 1 and n_micro ≥ 2 (got {} / {})",
                    dom.n_macro, dom.n_micro
                )));
            }
            let positive = [
                dom.thickness,
                dom.particle_radius,
                dom.diffusivity_solid,
                dom.c_solid_max,
                dom.sigma_solid,
                dom.porosity,
            ];
            if positive.iter().any(|v| !(v.is_finite() && *v > 0.0)) || dom.porosity >= 1.0 {
                return Err(Error::InvalidInput(format!(
                    "electrode {e}: non-physical geometry or transport parameters"
                )));
            }
            let cells = self.c_solid.get(e).map(Vec::len).unwrap_or(usize::MAX);
            let phis = self.phi_solid.get(e).map(Vec::len).unwrap_or(usize::MAX);
            let reac = self
                .reaction_current
                .get(e)
                .map(Vec::len)
                .unwrap_or(usize::MAX);
            if cells != dom.n_macro || phis != dom.n_macro || reac != dom.n_macro {
                return Err(Error::InvalidInput(format!(
                    "electrode {e}: state arrays do not match n_macro = {}",
                    dom.n_macro
                )));
            }
            let bad_particle = self.c_solid[e].iter().any(|p| p.len() != dom.n_micro);
            if bad_particle {
                return Err(Error::InvalidInput(format!(
                    "electrode {e}: particle arrays do not match n_micro = {}",
                    dom.n_micro
                )));
            }
        }
        let n_t = self.n_total_nodes();
        if self.c_electrolyte.len() != n_t || self.phi_electrolyte.len() != n_t {
            return Err(Error::InvalidInput(format!(
                "electrolyte arrays must have {n_t} cells"
            )));
        }
        if self
            .c_electrolyte
            .iter()
            .any(|c| !(c.is_finite() && *c > 0.0))
        {
            return Err(Error::General(
                "electrolyte depleted: non-positive Li⁺ concentration".to_owned(),
            ));
        }
        if !(self.kappa_electrolyte > 0.0 && self.d_electrolyte > 0.0) {
            return Err(Error::InvalidInput(
                "electrolyte conductivity and diffusivity must be positive".to_owned(),
            ));
        }
        Ok(())
    }

    /// Exchange current density i₀ (A/m²) and its surface stoichiometry check.
    fn exchange_current(&self, dom: &ElectrodeDomain, c_e: f64, c_surf: f64) -> Result<f64> {
        let theta = c_surf / dom.c_solid_max;
        if !(theta > 0.0 && theta < 1.0) {
            return Err(Error::General(format!(
                "{:?}: particle surface stoichiometry {theta} left (0, 1) — \
                 electrode fully (de)lithiated",
                dom.domain
            )));
        }
        let c_ratio = c_e / self.c_electrolyte_ref;
        Ok(dom.j0_ref * c_ratio.sqrt() * 2.0 * (theta * (1.0 - theta)).sqrt())
    }

    /// Solve the coupled φ_s / φ_e / Butler-Volmer charge balance at the
    /// current concentrations (damped Newton iteration).
    ///
    /// Updates `phi_solid`, `phi_electrolyte` and `reaction_current`.  On
    /// convergence Σ a_s h i_n equals I over the negative and −I over the
    /// positive electrode to within the Newton residual.
    pub fn solve_potentials(&mut self) -> Result<()> {
        self.validate()?;
        let n_n = self.neg.n_macro;
        let n_p = self.pos.n_macro;
        let n_t = self.n_total_nodes();
        let n_dof = n_n + n_t + n_p;
        let props = self.electrolyte_cells();
        let f_rt = FARADAY / (R_GAS * self.temperature);
        let chi = 2.0 * R_GAS * self.temperature * (1.0 - self.transference) / FARADAY;
        let ln_c: Vec<f64> = self.c_electrolyte.iter().map(|c| c.ln()).collect();

        // Per-electrode-cell kinetic data frozen over the step.
        let mut i0 = [vec![0.0; n_n], vec![0.0; n_p]];
        let mut u_ocp = [vec![0.0; n_n], vec![0.0; n_p]];
        for e in 0..2 {
            let dom = if e == 0 { &self.neg } else { &self.pos };
            for (k, particle) in self.c_solid[e].iter().enumerate() {
                let c_surf = particle.last().copied().unwrap_or(0.0);
                let x = self.electrolyte_index(e, k);
                i0[e][k] = self.exchange_current(dom, self.c_electrolyte[x], c_surf)?;
                u_ocp[e][k] = dom.open_circuit_potential(c_surf / dom.c_solid_max);
            }
        }

        let mut u = vec![0.0; n_dof];
        u[..n_n].copy_from_slice(&self.phi_solid[0]);
        u[n_n..n_n + n_t].copy_from_slice(&self.phi_electrolyte);
        u[n_n + n_t..].copy_from_slice(&self.phi_solid[1]);

        let ctx = PotentialSystem {
            n_n,
            n_p,
            n_t,
            pos_offset: n_n + self.n_sep,
            sigma: [
                self.neg.effective_solid_conductivity(),
                self.pos.effective_solid_conductivity(),
            ],
            h_el: [self.neg.cell_width(), self.pos.cell_width()],
            area: [
                self.neg.specific_interfacial_area(),
                self.pos.specific_interfacial_area(),
            ],
            props: &props,
            ln_c: &ln_c,
            chi,
            f_rt,
            alpha_a: self.alpha_anodic,
            alpha_c: self.alpha_cathodic,
            i0: &i0,
            u_ocp: &u_ocp,
            applied: self.applied_current,
        };

        let mut res = vec![0.0; n_dof];
        let mut jac = vec![0.0; n_dof * n_dof];
        ctx.evaluate(&u, &mut res, None);
        let mut res_norm = l2_norm(&res);
        // Residual floor (A/m²) accepted when the line search stalls at round-off.
        let res_floor = 1e-8 * (self.applied_current.abs() + 1.0);
        let mut converged = false;
        let mut iterations = 0;
        let mut trial = vec![0.0; n_dof];
        let mut trial_res = vec![0.0; n_dof];
        while iterations < MAX_NEWTON_ITERATIONS {
            iterations += 1;
            if !res_norm.is_finite() {
                break;
            }
            if res_norm == 0.0 {
                converged = true;
                break;
            }
            ctx.evaluate(&u, &mut res, Some(&mut jac));
            let mut delta: Vec<f64> = res.iter().map(|r| -r).collect();
            solve_dense_in_place(&mut jac, &mut delta, n_dof)?;
            // Backtracking line search on the residual 2-norm.
            let mut lambda = 1.0;
            let mut full_step = true;
            let mut accepted = false;
            for _ in 0..MAX_LINE_SEARCH_HALVINGS {
                for ((t, ui), di) in trial.iter_mut().zip(&u).zip(&delta) {
                    *t = ui + lambda * di;
                }
                ctx.evaluate(&trial, &mut trial_res, None);
                let trial_norm = l2_norm(&trial_res);
                if trial_norm.is_finite() && trial_norm <= (1.0 - 1e-4 * lambda) * res_norm {
                    accepted = true;
                    break;
                }
                lambda *= 0.5;
                full_step = false;
            }
            let step_norm = delta.iter().fold(0.0_f64, |m, d| m.max(d.abs()));
            let u_norm = u.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            let tiny_step = step_norm <= 1e-13 * (1.0 + u_norm);
            if !accepted {
                // At the round-off floor no further descent is possible.
                converged = tiny_step || res_norm <= res_floor;
                break;
            }
            u.copy_from_slice(&trial);
            res.copy_from_slice(&trial_res);
            res_norm = l2_norm(&res);
            if tiny_step || (full_step && step_norm <= 1e-11 * (1.0 + u_norm)) {
                // A full Newton step this small lands on the round-off floor.
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(Error::NotConverged {
                solver: "BatteryFem potential Newton".to_owned(),
                iterations,
                residual: res_norm / (self.applied_current.abs() + 1.0),
            });
        }
        self.phi_solid[0].copy_from_slice(&u[..n_n]);
        self.phi_electrolyte.copy_from_slice(&u[n_n..n_n + n_t]);
        self.phi_solid[1].copy_from_slice(&u[n_n + n_t..]);
        for e in 0..2 {
            let count = if e == 0 { n_n } else { n_p };
            for k in 0..count {
                let (i_n, _) = ctx.kinetics(&u, e, k);
                self.reaction_current[e][k] = i_n;
            }
        }
        Ok(())
    }

    /// Solid-diffusion update of one electrode (returns the new particle
    /// profiles without committing them).
    fn solid_diffusion_update(&self, electrode_idx: usize) -> Result<Vec<Vec<f64>>> {
        let dom = self.electrode(electrode_idx)?;
        let n_r = dom.n_micro;
        if n_r < 2 {
            return Err(Error::InvalidInput(format!(
                "n_micro must be at least 2, got {n_r}"
            )));
        }
        let radius = dom.particle_radius;
        let ds = dom.diffusivity_solid;
        let h_r = radius / (n_r as f64 - 1.0);
        let volumes = spherical_control_volumes(n_r, radius);
        // Face conductances D A_{i+1/2} / h with A = r² (per 4π).
        let face: Vec<f64> = (0..n_r - 1)
            .map(|i| ds * ((i as f64 + 0.5) * h_r).powi(2) / h_r)
            .collect();
        let inv_dt = 1.0 / self.dt;
        let particles = self.c_solid.get(electrode_idx).ok_or_else(|| {
            Error::InvalidInput(format!("missing solid state for electrode {electrode_idx}"))
        })?;
        let currents = self.reaction_current.get(electrode_idx).ok_or_else(|| {
            Error::InvalidInput(format!(
                "missing reaction current for electrode {electrode_idx}"
            ))
        })?;
        if currents.len() != particles.len() {
            return Err(Error::InvalidInput(
                "reaction_current and c_solid lengths differ".to_owned(),
            ));
        }
        let mut updated = Vec::with_capacity(particles.len());
        for (c, &i_n) in particles.iter().zip(currents) {
            if c.len() != n_r {
                return Err(Error::InvalidInput(format!(
                    "particle profile has {} nodes, expected {n_r}",
                    c.len()
                )));
            }
            let molar_flux = i_n / FARADAY;
            let mut lower = vec![0.0_f64; n_r];
            let mut diag = vec![0.0_f64; n_r];
            let mut upper = vec![0.0_f64; n_r];
            let mut rhs = vec![0.0_f64; n_r];
            for i in 0..n_r {
                diag[i] = volumes[i] * inv_dt;
                rhs[i] = volumes[i] * inv_dt * c[i];
                if i > 0 {
                    diag[i] += face[i - 1];
                    lower[i] = -face[i - 1];
                }
                if i + 1 < n_r {
                    diag[i] += face[i];
                    upper[i] = -face[i];
                }
            }
            // Outward molar flux i_n/F through the particle surface (area R²).
            rhs[n_r - 1] -= radius * radius * molar_flux;
            updated.push(NernstPlanckFem::solve_tridiagonal(
                &lower, &diag, &upper, &rhs,
            ));
        }
        Ok(updated)
    }

    /// Advance solid-phase diffusion in the particles of one electrode by one
    /// time step (implicit Euler, exactly conservative spherical finite volumes).
    ///
    /// The surface flux is −D_s ∂c_s/∂r|_R = i_n / F with i_n taken from
    /// `reaction_current` (initialised to the uniform distribution
    /// +I/(a_s L) / −I/(a_s L) and updated by [`Self::solve_potentials`]).
    /// The particle-average concentration therefore changes by exactly
    /// −3 i_n dt / (F r_p) per step.
    pub fn step_solid_diffusion(&mut self, electrode_idx: usize) -> Result<()> {
        if !(self.dt.is_finite() && self.dt > 0.0) {
            return Err(Error::InvalidInput(format!(
                "time step must be positive, got {}",
                self.dt
            )));
        }
        let updated = self.solid_diffusion_update(electrode_idx)?;
        self.c_solid[electrode_idx] = updated;
        Ok(())
    }

    /// Electrolyte diffusion update (returns the new profile without committing).
    fn electrolyte_update(&self) -> Vec<f64> {
        let props = self.electrolyte_cells();
        let n_t = self.n_total_nodes();
        let inv_dt = 1.0 / self.dt;
        let mut lower = vec![0.0_f64; n_t];
        let mut diag = vec![0.0_f64; n_t];
        let mut upper = vec![0.0_f64; n_t];
        let mut rhs = vec![0.0_f64; n_t];
        for k in 0..n_t {
            let storage = props.porosity[k] * props.width[k] * inv_dt;
            diag[k] = storage;
            rhs[k] = storage * self.c_electrolyte[k];
        }
        for k in 0..n_t.saturating_sub(1) {
            let g = face_conductance(
                props.width[k],
                props.diff_eff[k],
                props.width[k + 1],
                props.diff_eff[k + 1],
            );
            diag[k] += g;
            upper[k] -= g;
            lower[k + 1] -= g;
            diag[k + 1] += g;
        }
        let src = (1.0 - self.transference) / FARADAY;
        for e in 0..2 {
            let dom = if e == 0 { &self.neg } else { &self.pos };
            let area = dom.specific_interfacial_area();
            let h = dom.cell_width();
            for (k, &i_n) in self.reaction_current[e].iter().enumerate() {
                rhs[self.electrolyte_index(e, k)] += src * area * h * i_n;
            }
        }
        NernstPlanckFem::solve_tridiagonal(&lower, &diag, &upper, &rhs)
    }

    /// Advance the full DFN model by one time step.
    ///
    /// Solves the potential/kinetics problem, then advances solid and
    /// electrolyte diffusion with the resulting interfacial current.  The
    /// state is only committed when every sub-step succeeds and the new
    /// concentrations are physical; otherwise an error is returned and the
    /// concentrations are left at the start-of-step values.
    pub fn step(&mut self) -> Result<()> {
        self.solve_potentials()?;
        let neg_new = self.solid_diffusion_update(0)?;
        let pos_new = self.solid_diffusion_update(1)?;
        let c_e_new = self.electrolyte_update();
        for (dom, particles) in [(&self.neg, &neg_new), (&self.pos, &pos_new)] {
            let out_of_range = particles
                .iter()
                .flatten()
                .any(|&c| !(c.is_finite() && c > 0.0 && c < dom.c_solid_max));
            if out_of_range {
                return Err(Error::General(format!(
                    "{:?}: solid concentration left (0, c_max) — cut-off reached",
                    dom.domain
                )));
            }
        }
        if c_e_new.iter().any(|c| !(c.is_finite() && *c > 0.0)) {
            return Err(Error::General(
                "electrolyte depleted: non-positive Li⁺ concentration".to_owned(),
            ));
        }
        self.c_solid[0] = neg_new;
        self.c_solid[1] = pos_new;
        self.c_electrolyte = c_e_new;
        self.time += self.dt;
        Ok(())
    }

    /// Compute C-rate for a given capacity (Ah).
    pub fn c_rate(&self, capacity_ah: f64) -> f64 {
        self.applied_current / (capacity_ah * 3600.0)
    }

    /// Compute specific energy (Wh/kg) given electrode mass density (kg/m³).
    pub fn specific_energy(&self, mass_density: f64, volume: f64) -> f64 {
        let energy = self.terminal_voltage() * self.applied_current * self.time;
        let mass = mass_density * volume;
        energy / mass / 3600.0
    }
}

fn l2_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// Frozen data of the potential sub-problem for one Newton solve.
struct PotentialSystem<'a> {
    n_n: usize,
    n_p: usize,
    n_t: usize,
    /// Index of the first positive-electrode cell in the electrolyte list.
    pos_offset: usize,
    sigma: [f64; 2],
    h_el: [f64; 2],
    area: [f64; 2],
    props: &'a ElectrolyteCells,
    ln_c: &'a [f64],
    chi: f64,
    f_rt: f64,
    alpha_a: f64,
    alpha_c: f64,
    i0: &'a [Vec<f64>; 2],
    u_ocp: &'a [Vec<f64>; 2],
    applied: f64,
}

impl PotentialSystem<'_> {
    fn solid_dof(&self, e: usize, k: usize) -> usize {
        if e == 0 { k } else { self.n_n + self.n_t + k }
    }

    fn electrolyte_dof(&self, x: usize) -> usize {
        self.n_n + x
    }

    fn electrolyte_cell(&self, e: usize, k: usize) -> usize {
        if e == 0 { k } else { self.pos_offset + k }
    }

    /// Interfacial current i_n and ∂i_n/∂η for electrode cell (e, k).
    fn kinetics(&self, u: &[f64], e: usize, k: usize) -> (f64, f64) {
        let phi_s = u[self.solid_dof(e, k)];
        let phi_e = u[self.electrolyte_dof(self.electrolyte_cell(e, k))];
        let eta = phi_s - phi_e - self.u_ocp[e][k];
        let ea = (self.alpha_a * self.f_rt * eta).exp();
        let ec = (-self.alpha_c * self.f_rt * eta).exp();
        let i0 = self.i0[e][k];
        (
            i0 * (ea - ec),
            i0 * self.f_rt * (self.alpha_a * ea + self.alpha_c * ec),
        )
    }

    /// Evaluate the residual (A/m²) and optionally the dense Jacobian.
    fn evaluate(&self, u: &[f64], res: &mut [f64], mut jac: Option<&mut Vec<f64>>) {
        let n_dof = res.len();
        res.iter_mut().for_each(|r| *r = 0.0);
        if let Some(j) = jac.as_deref_mut() {
            j.iter_mut().for_each(|v| *v = 0.0);
        }
        // Linear face flux i = −G (u_q − u_p) + extra, leaving cell p, entering q.
        let add_face = |res: &mut [f64],
                        jac: &mut Option<&mut Vec<f64>>,
                        p: usize,
                        q: usize,
                        g: f64,
                        extra: f64| {
            let flux = -g * (u[q] - u[p]) + extra;
            res[p] += flux;
            res[q] -= flux;
            if let Some(j) = jac.as_deref_mut() {
                j[p * n_dof + p] += g;
                j[p * n_dof + q] -= g;
                j[q * n_dof + p] -= g;
                j[q * n_dof + q] += g;
            }
        };
        // Solid phase, both electrodes.
        for e in 0..2 {
            let count = if e == 0 { self.n_n } else { self.n_p };
            let g = self.sigma[e] / self.h_el[e];
            for k in 0..count.saturating_sub(1) {
                let p = self.solid_dof(e, k);
                let q = self.solid_dof(e, k + 1);
                add_face(res, &mut jac, p, q, g, 0.0);
            }
        }
        // Negative collector: Dirichlet φ_s = 0 on the x = 0 face (half cell).
        {
            let p = self.solid_dof(0, 0);
            let g = 2.0 * self.sigma[0] / self.h_el[0];
            // Current entering the cell through the face: −G (φ_p − 0).
            res[p] -= -g * u[p];
            if let Some(j) = jac.as_deref_mut() {
                j[p * n_dof + p] += g;
            }
        }
        // Positive collector: solid current I leaves through the x = L face.
        {
            let p = self.solid_dof(1, self.n_p - 1);
            res[p] += self.applied;
        }
        // Electrolyte phase (zero current at both ends).
        let props = self.props;
        for x in 0..self.n_t.saturating_sub(1) {
            let g = face_conductance(
                props.width[x],
                props.kappa_eff[x],
                props.width[x + 1],
                props.kappa_eff[x + 1],
            );
            let diffusion = self.chi * g * (self.ln_c[x + 1] - self.ln_c[x]);
            let p = self.electrolyte_dof(x);
            let q = self.electrolyte_dof(x + 1);
            add_face(res, &mut jac, p, q, g, diffusion);
        }
        // Butler-Volmer exchange between the phases.
        for e in 0..2 {
            let count = if e == 0 { self.n_n } else { self.n_p };
            let w = self.area[e] * self.h_el[e];
            for k in 0..count {
                let (i_n, di_n) = self.kinetics(u, e, k);
                let s = self.solid_dof(e, k);
                let l = self.electrolyte_dof(self.electrolyte_cell(e, k));
                res[s] += w * i_n;
                res[l] -= w * i_n;
                if let Some(j) = jac.as_deref_mut() {
                    j[s * n_dof + s] += w * di_n;
                    j[s * n_dof + l] -= w * di_n;
                    j[l * n_dof + s] -= w * di_n;
                    j[l * n_dof + l] += w * di_n;
                }
            }
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_battery_fem_creation() {
        let bat = BatteryFem::new(1.0, 298.15, 1.0);
        assert_eq!(bat.neg.domain, DfnDomain::NegativeElectrode);
        assert_eq!(bat.pos.domain, DfnDomain::PositiveElectrode);
        assert_eq!(bat.phi_solid[0].len(), bat.neg.n_macro);
        assert_eq!(bat.phi_solid[1].len(), bat.pos.n_macro);
    }

    #[test]
    fn test_battery_n_total_nodes() {
        let bat = BatteryFem::new(1.0, 298.15, 1.0);
        let expected = bat.neg.n_macro + bat.n_sep + bat.pos.n_macro;
        assert_eq!(bat.n_total_nodes(), expected);
    }

    #[test]
    fn test_battery_soc_initial() {
        let bat = BatteryFem::new(1.0, 298.15, 1.0);
        let soc_n = bat.soc_negative();
        let soc_p = bat.soc_positive();
        assert!((soc_n - 0.5).abs() < 1e-10);
        assert!((soc_p - 0.8).abs() < 1e-10);
    }

    #[test]
    fn test_battery_initial_state_is_at_rest() {
        let bat = BatteryFem::new(1.0, 298.15, 1.0);
        for e in 0..2 {
            for k in 0..bat.phi_solid[e].len() {
                let eta = bat.overpotential(e, k).unwrap_or(f64::NAN);
                assert!(eta.abs() < 1e-12, "η[{e}][{k}] = {eta}");
            }
        }
        assert!(bat.overpotential(2, 0).is_none());
    }

    #[test]
    fn test_battery_terminal_voltage_range() {
        let bat = BatteryFem::new(1.0, 298.15, 1.0);
        let v = bat.terminal_voltage();
        assert!(v > 2.0 && v < 5.0, "V_cell = {v}");
    }

    #[test]
    fn test_battery_ir_drop_terminal_voltage() {
        // R_int = 0.05 Ω·m², I = 2 A/m² → IR drop = 0.1 V exactly.
        let mut bat = BatteryFem::new(2.0, 298.15, 1.0);
        bat.internal_resistance = 0.05;
        let v_no_r = {
            let mut b2 = bat.clone();
            b2.internal_resistance = 0.0;
            b2.terminal_voltage()
        };
        let v_with_r = bat.terminal_voltage();
        assert!(
            v_with_r < v_no_r,
            "IR drop not applied: v_with_r={v_with_r}, v_no_r={v_no_r}"
        );
        assert!(
            (v_no_r - v_with_r - 0.1).abs() < 1e-10,
            "expected IR drop of 0.1 V, got {:.6}",
            v_no_r - v_with_r
        );
    }

    #[test]
    fn test_battery_step() {
        let mut bat = BatteryFem::new(1.0, 298.15, 10.0);
        bat.step().expect("DFN step");
        assert!((bat.time - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_electrode_effective_conductivity() {
        let neg = ElectrodeDomain::default_negative();
        let kappa_eff = neg.effective_conductivity(1.0);
        assert!(kappa_eff < 1.0 && kappa_eff > 0.0);
    }

    #[test]
    fn test_electrode_effective_diffusivity() {
        let pos = ElectrodeDomain::default_positive();
        let d_eff = pos.effective_electrolyte_diffusivity(2.7e-10);
        assert!(d_eff < 2.7e-10 && d_eff > 0.0);
    }

    #[test]
    fn test_spherical_control_volumes_tile_sphere() {
        for n in [2, 3, 5, 17] {
            let r = 2.5e-6;
            let v = spherical_control_volumes(n, r);
            let sum: f64 = v.iter().sum();
            assert!(((sum - r.powi(3) / 3.0) / (r.powi(3) / 3.0)).abs() < 1e-14);
        }
    }

    #[test]
    fn test_solid_diffusion_obeys_faradays_law() {
        // Regression (F1160): the surface row used to be a Dirichlet update
        // with a spurious dt factor and a slab-type centre row.
        let mut bat = BatteryFem::new(5.0, 298.15, 2.0);
        let r_n = bat.neg.particle_radius;
        let vol = spherical_control_volumes(bat.neg.n_micro, r_n);
        let before: Vec<f64> = bat.c_solid[0]
            .iter()
            .map(|c| particle_average(c, &vol, r_n))
            .collect();
        let i_n = bat.reaction_current[0][0];
        bat.step_solid_diffusion(0).expect("solid diffusion");
        let expected = -3.0 * i_n * bat.dt / (FARADAY * r_n);
        for (k, c) in bat.c_solid[0].iter().enumerate() {
            let change = particle_average(c, &vol, r_n) - before[k];
            assert!(
                ((change - expected) / expected).abs() < 1e-10,
                "cell {k}: Δc̄ = {change}, expected {expected}"
            );
        }
        assert!(bat.step_solid_diffusion(2).is_err());
    }

    #[test]
    fn test_solid_diffusion_electrode_signs() {
        // Discharge: Li leaves the negative particles, enters the positive ones.
        let mut bat = BatteryFem::new(5.0, 298.15, 5.0);
        let n0 = bat.soc_negative();
        let p0 = bat.soc_positive();
        bat.step_solid_diffusion(0).expect("neg");
        bat.step_solid_diffusion(1).expect("pos");
        assert!(bat.soc_negative() < n0);
        assert!(bat.soc_positive() > p0);
        // Surface depleted (neg) / enriched (pos) relative to the centre.
        let neg = &bat.c_solid[0][0];
        let pos = &bat.c_solid[1][0];
        assert!(neg[neg.len() - 1] < neg[0]);
        assert!(pos[pos.len() - 1] > pos[0]);
    }

    #[test]
    fn test_solid_diffusion_quasi_steady_profile() {
        // Constant outward flux j: after the transient c_surf − c̄ = −j R / (5 D).
        let mut bat = BatteryFem::new(10.0, 298.15, 1.0);
        bat.neg.n_micro = 41;
        bat.neg.n_macro = 1;
        bat.c_solid[0] = vec![vec![15000.0; 41]];
        bat.reaction_current[0] = vec![0.1];
        let r = bat.neg.particle_radius;
        let d = bat.neg.diffusivity_solid;
        let vol = spherical_control_volumes(41, r);
        for _ in 0..1500 {
            bat.step_solid_diffusion(0).expect("solid diffusion");
        }
        let c = &bat.c_solid[0][0];
        let j = 0.1 / FARADAY;
        let got = c[40] - particle_average(c, &vol, r);
        let expected = -j * r / (5.0 * d);
        assert!(
            ((got - expected) / expected).abs() < 0.02,
            "c_surf − c̄ = {got}, expected {expected}"
        );
    }

    #[test]
    fn test_battery_total_lithium_conserved() {
        let mut bat = BatteryFem::new(15.0, 298.15, 10.0);
        let li0 = bat.total_lithium();
        for _ in 0..40 {
            bat.step().expect("DFN step");
        }
        let li = bat.total_lithium();
        assert!(((li - li0) / li0).abs() < 1e-10, "Li {li0} -> {li}");
    }

    #[test]
    fn test_battery_reaction_current_integrates_to_applied_current() {
        let mut bat = BatteryFem::new(20.0, 298.15, 5.0);
        bat.step().expect("DFN step");
        let sum_neg: f64 = bat.reaction_current[0]
            .iter()
            .map(|i| i * bat.neg.specific_interfacial_area() * bat.neg.cell_width())
            .sum();
        let sum_pos: f64 = bat.reaction_current[1]
            .iter()
            .map(|i| i * bat.pos.specific_interfacial_area() * bat.pos.cell_width())
            .sum();
        assert!(((sum_neg - 20.0) / 20.0).abs() < 1e-9, "Σ_neg = {sum_neg}");
        assert!(((sum_pos + 20.0) / 20.0).abs() < 1e-9, "Σ_pos = {sum_pos}");
        // Anodic in the negative, cathodic in the positive electrode.
        assert!(bat.reaction_current[0].iter().all(|&i| i > 0.0));
        assert!(bat.reaction_current[1].iter().all(|&i| i < 0.0));
        for k in 0..bat.neg.n_macro {
            assert!(bat.overpotential(0, k).unwrap_or(-1.0) > 0.0);
        }
        for k in 0..bat.pos.n_macro {
            assert!(bat.overpotential(1, k).unwrap_or(1.0) < 0.0);
        }
    }

    #[test]
    fn test_battery_lithium_moves_from_negative_to_positive() {
        let current = 10.0;
        let mut bat = BatteryFem::new(current, 298.15, 10.0);
        let lith = |b: &BatteryFem, e: usize| -> f64 {
            let dom = if e == 0 { &b.neg } else { &b.pos };
            let vol = spherical_control_volumes(dom.n_micro, dom.particle_radius);
            b.c_solid[e]
                .iter()
                .map(|c| {
                    dom.solid_volume_fraction()
                        * dom.cell_width()
                        * particle_average(c, &vol, dom.particle_radius)
                })
                .sum()
        };
        let n0 = lith(&bat, 0);
        let p0 = lith(&bat, 1);
        let n_steps = 20;
        for _ in 0..n_steps {
            bat.step().expect("DFN step");
        }
        let transferred = current * bat.dt * n_steps as f64 / FARADAY;
        let dn = lith(&bat, 0) - n0;
        let dp = lith(&bat, 1) - p0;
        assert!(
            ((dn + transferred) / transferred).abs() < 1e-8,
            "Δneg = {dn}"
        );
        assert!(
            ((dp - transferred) / transferred).abs() < 1e-8,
            "Δpos = {dp}"
        );
    }

    #[test]
    fn test_battery_electrolyte_polarisation_during_discharge() {
        let mut bat = BatteryFem::new(20.0, 298.15, 5.0);
        let props = bat.electrolyte_cells();
        let salt = |b: &BatteryFem| -> f64 {
            b.c_electrolyte
                .iter()
                .zip(props.porosity.iter().zip(&props.width))
                .map(|(c, (e, h))| c * e * h)
                .sum()
        };
        let s0 = salt(&bat);
        for _ in 0..30 {
            bat.step().expect("DFN step");
        }
        let n_t = bat.n_total_nodes();
        // Li⁺ is released in the negative and consumed in the positive electrode.
        assert!(bat.c_electrolyte[0] > 1000.0);
        assert!(bat.c_electrolyte[n_t - 1] < 1000.0);
        assert!(((salt(&bat) - s0) / s0).abs() < 1e-10);
        // Ionic current flows towards the positive electrode: φ_e decreases.
        assert!(bat.phi_electrolyte[0] > bat.phi_electrolyte[n_t - 1]);
    }

    #[test]
    fn test_battery_voltage_polarisation_sign() {
        let ocv = BatteryFem::new(0.0, 298.15, 1.0).terminal_voltage();
        let mut dis = BatteryFem::new(20.0, 298.15, 1.0);
        dis.solve_potentials().expect("discharge potentials");
        let mut chg = BatteryFem::new(-20.0, 298.15, 1.0);
        chg.solve_potentials().expect("charge potentials");
        assert!(dis.terminal_voltage() < ocv - 1e-4);
        assert!(chg.terminal_voltage() > ocv + 1e-4);
        // Symmetric kinetics and linear transport: polarisation is odd in I
        // to leading order.
        let dv_dis = ocv - dis.terminal_voltage();
        let dv_chg = chg.terminal_voltage() - ocv;
        assert!((dv_dis - dv_chg).abs() / dv_dis < 0.05);
    }

    #[test]
    fn test_battery_rest_state_is_stationary() {
        let mut bat = BatteryFem::new(0.0, 298.15, 50.0);
        let v0 = bat.terminal_voltage();
        let c_s0 = bat.c_solid.clone();
        for _ in 0..10 {
            bat.step().expect("rest step");
        }
        assert!((bat.terminal_voltage() - v0).abs() < 1e-10);
        for (a, b) in bat
            .c_solid
            .iter()
            .flatten()
            .flatten()
            .zip(c_s0.iter().flatten().flatten())
        {
            assert!((a - b).abs() < 1e-8 * b);
        }
        for &c in &bat.c_electrolyte {
            assert!((c - 1000.0).abs() < 1e-9);
        }
    }

    #[test]
    fn test_battery_voltage_decreases_during_discharge() {
        let mut bat = BatteryFem::new(20.0, 298.15, 20.0);
        bat.step().expect("first step");
        let v_start = bat.terminal_voltage();
        for _ in 0..30 {
            bat.step().expect("DFN step");
        }
        assert!(bat.terminal_voltage() < v_start);
    }

    #[test]
    fn test_battery_rejects_inconsistent_state() {
        let mut bat = BatteryFem::new(1.0, 298.15, 1.0);
        bat.c_electrolyte.pop();
        assert!(bat.step().is_err());
        let mut bat = BatteryFem::new(1.0, 298.15, 0.0);
        assert!(bat.step().is_err());
    }

    #[test]
    fn test_battery_cutoff_reports_error_and_keeps_state() {
        // Drive the negative electrode to full delithiation.
        let mut bat = BatteryFem::new(200.0, 298.15, 200.0);
        let mut failed = false;
        for _ in 0..200 {
            let snapshot = bat.c_solid.clone();
            if bat.step().is_err() {
                assert_eq!(bat.c_solid, snapshot);
                failed = true;
                break;
            }
        }
        assert!(failed, "cut-off never reported");
    }

    #[test]
    fn test_dense_solver_matches_known_solution() {
        let mut a = vec![4.0, -2.0, 1.0, -2.0, 4.0, -2.0, 1.0, -2.0, 4.0];
        let mut b = vec![11.0, -16.0, 17.0];
        solve_dense_in_place(&mut a, &mut b, 3).expect("solve");
        let x = [1.0, -2.0, 3.0];
        for (got, want) in b.iter().zip(x) {
            assert!((got - want).abs() < 1e-12);
        }
        let mut singular = vec![1.0, 2.0, 2.0, 4.0];
        let mut rhs = vec![1.0, 2.0];
        assert!(solve_dense_in_place(&mut singular, &mut rhs, 2).is_err());
    }
}
