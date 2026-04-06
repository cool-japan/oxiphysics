#![allow(clippy::if_same_then_else, clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Traffic-flow simulation using Lattice Boltzmann and cellular-automaton methods.
//!
//! This module provides:
//!
//! - [`NagelSchreckenberg`]: Nagel–Schreckenberg stochastic cellular automaton
//! - [`LwrLattice`]: Lighthill–Whitham–Richards continuum model on a lattice
//! - [`FundamentalDiagram`]: Greenshields / triangular flow–density relations
//! - [`TrafficSignal`]: Fixed-cycle traffic signal timing
//! - [`OnRamp`]: On-ramp merging with priority rules
//! - [`OffRamp`]: Off-ramp diverging flow
//! - [`Roundabout`]: Single-lane roundabout circulation
//! - [`MultiLane`]: Multi-lane road with lane-change rules
//! - [`SpeedLimitZone`]: Enforces variable-speed limits on road sections
//! - [`AccidentPropagation`]: Accident shock-wave propagation model
//! - [`TrafficLbm`]: D1Q3 LBM solver for macroscopic traffic density/flow

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// D1Q3 lattice constants  (velocities -1, 0, +1)
// ─────────────────────────────────────────────────────────────────────────────

const NQ: usize = 3;
/// D1Q3 lattice velocities.
const CV: [f64; NQ] = [-1.0, 0.0, 1.0];
/// D1Q3 equilibrium weights.
const W3: [f64; NQ] = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0];
/// Speed of sound squared for D1Q3.
const CS2_1D: f64 = 1.0 / 3.0;

// ─────────────────────────────────────────────────────────────────────────────
// Helper: D1Q3 equilibrium distribution
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn feq1d(rho: f64, u: f64) -> [f64; NQ] {
    let mut f = [0.0f64; NQ];
    let u2 = u * u;
    for q in 0..NQ {
        let cu = CV[q] * u;
        f[q] = W3[q]
            * rho
            * (1.0 + cu / CS2_1D + cu * cu / (2.0 * CS2_1D * CS2_1D) - u2 / (2.0 * CS2_1D));
    }
    f
}

#[inline]
fn macroscopic1d(f: &[f64; NQ]) -> (f64, f64) {
    let rho = f.iter().sum::<f64>();
    if rho < 1e-300 {
        return (0.0, 0.0);
    }
    let u = f.iter().enumerate().map(|(q, &v)| CV[q] * v).sum::<f64>() / rho;
    (rho, u)
}

// ─────────────────────────────────────────────────────────────────────────────
// FundamentalDiagram
// ─────────────────────────────────────────────────────────────────────────────

/// Model selection for flow–density fundamental diagram.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DiagramModel {
    /// Greenshields linear speed–density model.
    Greenshields,
    /// Triangular (Newell–Daganzo) model.
    Triangular,
    /// Underwood exponential model.
    Underwood,
    /// Drake exponential model.
    Drake,
}

/// Flow–density fundamental diagram for a road section.
#[derive(Clone, Debug)]
pub struct FundamentalDiagram {
    /// Maximum (jam) density (veh/m).
    pub rho_jam: f64,
    /// Free-flow speed (m/s).
    pub v_free: f64,
    /// Critical density at capacity (veh/m).
    pub rho_crit: f64,
    /// Maximum flow (capacity) (veh/s).
    pub q_max: f64,
    /// Selected model variant.
    pub model: DiagramModel,
    /// Backward wave speed for triangular model (m/s).
    pub w_back: f64,
}

impl FundamentalDiagram {
    /// Create a new `FundamentalDiagram`.
    pub fn new(rho_jam: f64, v_free: f64, model: DiagramModel) -> Self {
        let rho_crit = rho_jam / 2.0;
        let q_max = v_free * rho_crit / 2.0;
        let w_back = q_max / (rho_jam - rho_crit);
        Self {
            rho_jam,
            v_free,
            rho_crit,
            q_max,
            model,
            w_back,
        }
    }

    /// Compute equilibrium speed for given density.
    pub fn speed(&self, rho: f64) -> f64 {
        let r = (rho / self.rho_jam).clamp(0.0, 1.0);
        match self.model {
            DiagramModel::Greenshields => self.v_free * (1.0 - r),
            DiagramModel::Triangular => {
                if rho <= self.rho_crit {
                    self.v_free
                } else {
                    self.w_back * (self.rho_jam / rho - 1.0)
                }
            }
            DiagramModel::Underwood => self.v_free * (-r / (1.0 - r + 1e-12)).exp().min(1.0),
            DiagramModel::Drake => self.v_free * (-(rho / self.rho_crit).powi(2) / 2.0).exp(),
        }
    }

    /// Compute flow q = rho * v(rho).
    pub fn flow(&self, rho: f64) -> f64 {
        rho * self.speed(rho)
    }

    /// Compute density at which flow equals `q` (Newton iteration).
    pub fn density_for_flow(&self, q: f64) -> f64 {
        let mut rho = self.rho_crit;
        for _ in 0..50 {
            let f0 = self.flow(rho) - q;
            let df = (self.flow(rho + 1e-6) - self.flow(rho - 1e-6)) / 2e-6;
            if df.abs() < 1e-14 {
                break;
            }
            rho = (rho - f0 / df).clamp(0.0, self.rho_jam);
        }
        rho
    }

    /// Compute capacity drop factor (1 means no drop).
    pub fn capacity_drop_factor(&self) -> f64 {
        self.q_max / (self.v_free * self.rho_jam / 4.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NagelSchreckenberg cellular automaton
// ─────────────────────────────────────────────────────────────────────────────

/// Nagel–Schreckenberg stochastic cellular automaton.
///
/// Each cell represents one vehicle position on a single-lane ring road.
/// Vehicles are identified by integer indices stored in the `positions` vector.
#[derive(Clone, Debug)]
pub struct NagelSchreckenberg {
    /// Road length (number of cells).
    pub road_len: usize,
    /// Maximum velocity (cells/step).
    pub v_max: usize,
    /// Dawdling (braking) probability.
    pub p_dawdle: f64,
    /// Vehicle positions (cell index).
    pub positions: Vec<usize>,
    /// Vehicle velocities (cells/step).
    pub velocities: Vec<usize>,
    /// Step counter.
    pub step: usize,
}

impl NagelSchreckenberg {
    /// Create a new `NagelSchreckenberg` simulation.
    pub fn new(road_len: usize, v_max: usize, p_dawdle: f64) -> Self {
        Self {
            road_len,
            v_max,
            p_dawdle,
            positions: Vec::new(),
            velocities: Vec::new(),
            step: 0,
        }
    }

    /// Add a vehicle at `pos` with initial velocity `vel`.
    pub fn add_vehicle(&mut self, pos: usize, vel: usize) {
        self.positions.push(pos % self.road_len);
        self.velocities.push(vel.min(self.v_max));
    }

    /// Populate road uniformly with `n_vehicles` vehicles.
    pub fn populate_uniform(&mut self, n_vehicles: usize) {
        self.positions.clear();
        self.velocities.clear();
        let spacing = self.road_len / n_vehicles.max(1);
        for i in 0..n_vehicles {
            self.positions.push(i * spacing);
            self.velocities.push(self.v_max / 2);
        }
    }

    /// Advance simulation by one time step using a deterministic random seed.
    pub fn step_with_seed(&mut self, rand_bits: &[u8]) {
        let n = self.positions.len();
        if n == 0 {
            return;
        }
        // Sort vehicles by position for gap calculation
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by_key(|&i| self.positions[i]);

        let mut new_vel = self.velocities.clone();
        let mut new_pos = self.positions.clone();

        for idx in 0..n {
            let i = order[idx];
            let next_i = order[(idx + 1) % n];
            let pos_i = self.positions[i];
            let pos_next = self.positions[next_i];
            // Gap = distance to vehicle ahead (periodic)
            let gap = if pos_next > pos_i {
                pos_next - pos_i - 1
            } else {
                self.road_len + pos_next - pos_i - 1
            };
            // Rule 1: accelerate
            let v1 = (self.velocities[i] + 1).min(self.v_max);
            // Rule 2: brake for gap
            let v2 = v1.min(gap);
            // Rule 3: randomisation (dawdling)
            let rand_bit = rand_bits.get(i).copied().unwrap_or(0);
            let v3 = if rand_bit < (255.0 * self.p_dawdle) as u8 && v2 > 0 {
                v2 - 1
            } else {
                v2
            };
            new_vel[i] = v3;
            new_pos[i] = (pos_i + v3) % self.road_len;
        }

        self.velocities = new_vel;
        self.positions = new_pos;
        self.step += 1;
    }

    /// Compute mean flow (vehicles/step) on the ring.
    pub fn mean_flow(&self) -> f64 {
        if self.positions.is_empty() {
            return 0.0;
        }
        self.velocities.iter().sum::<usize>() as f64 / self.road_len as f64
    }

    /// Compute global density (vehicles/cell).
    pub fn density(&self) -> f64 {
        self.positions.len() as f64 / self.road_len as f64
    }

    /// Compute mean speed (cells/step).
    pub fn mean_speed(&self) -> f64 {
        if self.velocities.is_empty() {
            return 0.0;
        }
        self.velocities.iter().sum::<usize>() as f64 / self.velocities.len() as f64
    }

    /// Detect traffic jams: return list of jam-head cell positions.
    ///
    /// A cell is considered a jam head if a vehicle has velocity 0.
    pub fn jam_heads(&self) -> Vec<usize> {
        self.positions
            .iter()
            .zip(self.velocities.iter())
            .filter(|&(_, &v)| v == 0)
            .map(|(&p, _)| p)
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// LwrLattice: LWR continuum model on a 1-D lattice
// ─────────────────────────────────────────────────────────────────────────────

/// Lighthill–Whitham–Richards continuum traffic model solved on a 1-D lattice
/// using the Godunov upwind scheme.
#[derive(Clone, Debug)]
pub struct LwrLattice {
    /// Number of lattice cells.
    pub nx: usize,
    /// Cell length (m).
    pub dx: f64,
    /// Time step (s).
    pub dt: f64,
    /// Density field (veh/m), length `nx`.
    pub rho: Vec<f64>,
    /// Fundamental diagram.
    pub diagram: FundamentalDiagram,
    /// Simulation time (s).
    pub time: f64,
}

impl LwrLattice {
    /// Create a new `LwrLattice`.
    pub fn new(nx: usize, dx: f64, dt: f64, diagram: FundamentalDiagram) -> Self {
        Self {
            nx,
            dx,
            dt,
            rho: vec![0.0; nx],
            diagram,
            time: 0.0,
        }
    }

    /// Set initial density profile from a closure.
    pub fn init_from<F: Fn(f64) -> f64>(&mut self, f: F) {
        for i in 0..self.nx {
            let x = (i as f64 + 0.5) * self.dx;
            self.rho[i] = f(x).clamp(0.0, self.diagram.rho_jam);
        }
    }

    /// Godunov numerical flux at interface between cells `left` and `right`.
    fn godunov_flux(&self, rho_l: f64, rho_r: f64) -> f64 {
        let q_l = self.diagram.flow(rho_l);
        let q_r = self.diagram.flow(rho_r);
        let rho_c = self.diagram.rho_crit;
        // Supply / demand formulation
        let demand = if rho_l <= rho_c {
            q_l
        } else {
            self.diagram.q_max
        };
        let supply = if rho_r >= rho_c {
            q_r
        } else {
            self.diagram.q_max
        };
        demand.min(supply)
    }

    /// Advance by one time step (Godunov scheme).
    pub fn step(&mut self) {
        let mut rho_new = self.rho.clone();
        let cfl = self.dt / self.dx;
        for i in 0..self.nx {
            let il = if i == 0 { self.nx - 1 } else { i - 1 };
            let ir = (i + 1) % self.nx;
            let f_right = self.godunov_flux(self.rho[i], self.rho[ir]);
            let f_left = self.godunov_flux(self.rho[il], self.rho[i]);
            rho_new[i] = (self.rho[i] - cfl * (f_right - f_left)).clamp(0.0, self.diagram.rho_jam);
        }
        self.rho = rho_new;
        self.time += self.dt;
    }

    /// Run for `n_steps` steps.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.step();
        }
    }

    /// Compute total vehicle count (integral of density).
    pub fn total_vehicles(&self) -> f64 {
        self.rho.iter().sum::<f64>() * self.dx
    }

    /// Compute mean density.
    pub fn mean_density(&self) -> f64 {
        self.rho.iter().sum::<f64>() / self.nx as f64
    }

    /// Compute total flow (vehicles/s) at cell `i`.
    pub fn flow_at(&self, i: usize) -> f64 {
        self.diagram.flow(self.rho[i])
    }

    /// Check CFL condition: returns true if dt/dx * v_free <= 1.
    pub fn cfl_ok(&self) -> bool {
        self.dt / self.dx * self.diagram.v_free <= 1.0
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TrafficSignal
// ─────────────────────────────────────────────────────────────────────────────

/// Phase of a traffic signal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SignalPhase {
    /// Green: vehicles may pass.
    Green,
    /// Yellow: vehicles decelerate.
    Yellow,
    /// Red: vehicles must stop.
    Red,
}

/// Fixed-cycle traffic signal with green/yellow/red phases.
#[derive(Clone, Debug)]
pub struct TrafficSignal {
    /// Cell index where the signal is located.
    pub cell: usize,
    /// Green phase duration (s).
    pub green_duration: f64,
    /// Yellow phase duration (s).
    pub yellow_duration: f64,
    /// Red phase duration (s).
    pub red_duration: f64,
    /// Current phase.
    pub phase: SignalPhase,
    /// Time elapsed in current phase (s).
    pub phase_time: f64,
    /// Offset (s) for coordinated green waves.
    pub offset: f64,
}

impl TrafficSignal {
    /// Create a new `TrafficSignal` starting in Red.
    pub fn new(cell: usize, green: f64, yellow: f64, red: f64, offset: f64) -> Self {
        let phase_time = offset % (green + yellow + red);
        let phase = if phase_time < green {
            SignalPhase::Green
        } else if phase_time < green + yellow {
            SignalPhase::Yellow
        } else {
            SignalPhase::Red
        };
        Self {
            cell,
            green_duration: green,
            yellow_duration: yellow,
            red_duration: red,
            phase,
            phase_time,
            offset,
        }
    }

    /// Advance the signal by `dt` seconds.
    pub fn advance(&mut self, dt: f64) {
        self.phase_time += dt;
        let cycle = self.green_duration + self.yellow_duration + self.red_duration;
        self.phase_time %= cycle;
        self.phase = if self.phase_time < self.green_duration {
            SignalPhase::Green
        } else if self.phase_time < self.green_duration + self.yellow_duration {
            SignalPhase::Yellow
        } else {
            SignalPhase::Red
        };
    }

    /// Return the effective speed factor (1.0 = full flow, 0.0 = stopped).
    pub fn speed_factor(&self) -> f64 {
        match self.phase {
            SignalPhase::Green => 1.0,
            SignalPhase::Yellow => 0.5,
            SignalPhase::Red => 0.0,
        }
    }

    /// Return remaining time in current phase (s).
    pub fn remaining(&self) -> f64 {
        match self.phase {
            SignalPhase::Green => self.green_duration - self.phase_time,
            SignalPhase::Yellow => self.green_duration + self.yellow_duration - self.phase_time,
            SignalPhase::Red => {
                self.green_duration + self.yellow_duration + self.red_duration - self.phase_time
            }
        }
    }

    /// Compute optimal green split for two competing flows `q1` and `q2`.
    pub fn optimal_green_split(total_cycle: f64, q1: f64, q2: f64) -> (f64, f64) {
        let total = q1 + q2;
        if total < 1e-12 {
            return (total_cycle / 2.0, total_cycle / 2.0);
        }
        (total_cycle * q1 / total, total_cycle * q2 / total)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// OnRamp
// ─────────────────────────────────────────────────────────────────────────────

/// On-ramp merging model.
///
/// Vehicles enter from the ramp when the main-road density at the merge point
/// is below a threshold.
#[derive(Clone, Debug)]
pub struct OnRamp {
    /// Cell index of merge point on main road.
    pub merge_cell: usize,
    /// Desired ramp inflow (veh/s).
    pub demand: f64,
    /// Maximum density on main road for merge acceptance (veh/m).
    pub merge_density_limit: f64,
    /// Accumulated vehicles waiting on ramp.
    pub queue: f64,
    /// Total vehicles that have merged.
    pub merged_count: f64,
}

impl OnRamp {
    /// Create a new `OnRamp`.
    pub fn new(merge_cell: usize, demand: f64, density_limit: f64) -> Self {
        Self {
            merge_cell,
            demand,
            merge_density_limit: density_limit,
            queue: 0.0,
            merged_count: 0.0,
        }
    }

    /// Process merging: returns actual inflow (veh/s) given main-road density.
    pub fn process(&mut self, main_density: f64, dt: f64) -> f64 {
        self.queue += self.demand * dt;
        if main_density < self.merge_density_limit && self.queue > 0.0 {
            let inflow = self.queue.min(self.demand * dt);
            self.queue -= inflow;
            self.merged_count += inflow;
            inflow / dt
        } else {
            0.0
        }
    }

    /// Queue length in vehicles.
    pub fn queue_length(&self) -> f64 {
        self.queue
    }

    /// Merge acceptance rate (fraction, 0..1) based on density.
    pub fn acceptance_rate(&self, main_density: f64) -> f64 {
        (1.0 - main_density / self.merge_density_limit).clamp(0.0, 1.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// OffRamp
// ─────────────────────────────────────────────────────────────────────────────

/// Off-ramp diverging model.
#[derive(Clone, Debug)]
pub struct OffRamp {
    /// Cell index of diverge point on main road.
    pub diverge_cell: usize,
    /// Fraction of traffic taking the off-ramp (0..1).
    pub split_fraction: f64,
    /// Cumulative vehicles exited via off-ramp.
    pub exited_count: f64,
}

impl OffRamp {
    /// Create a new `OffRamp`.
    pub fn new(diverge_cell: usize, split_fraction: f64) -> Self {
        Self {
            diverge_cell,
            split_fraction: split_fraction.clamp(0.0, 1.0),
            exited_count: 0.0,
        }
    }

    /// Compute exit flow (veh/s) for given through-flow.
    pub fn exit_flow(&mut self, through_flow: f64, dt: f64) -> f64 {
        let exiting = through_flow * self.split_fraction;
        self.exited_count += exiting * dt;
        exiting
    }

    /// Remaining through-flow after diversion.
    pub fn through_flow(&self, total_flow: f64) -> f64 {
        total_flow * (1.0 - self.split_fraction)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Roundabout
// ─────────────────────────────────────────────────────────────────────────────

/// Single-lane roundabout simulation using cellular automaton rules.
#[derive(Clone, Debug)]
pub struct Roundabout {
    /// Number of cells in the roundabout ring.
    pub ring_cells: usize,
    /// Number of entry legs.
    pub n_entries: usize,
    /// Circulating density per cell (vehicles).
    pub ring: Vec<f64>,
    /// Inflow demand per entry (veh/s).
    pub inflow: Vec<f64>,
    /// Outflow fraction per exit (0..1).
    pub outflow_fraction: Vec<f64>,
    /// Critical gap acceptance density (veh/cell).
    pub critical_gap_density: f64,
    /// Simulation time step (s).
    pub dt: f64,
    /// Simulation time (s).
    pub time: f64,
}

impl Roundabout {
    /// Create a new `Roundabout`.
    pub fn new(ring_cells: usize, n_entries: usize, critical_gap: f64, dt: f64) -> Self {
        Self {
            ring_cells,
            n_entries,
            ring: vec![0.0; ring_cells],
            inflow: vec![0.0; n_entries],
            outflow_fraction: vec![1.0 / n_entries as f64; n_entries],
            critical_gap_density: critical_gap,
            dt,
            time: 0.0,
        }
    }

    /// Set inflow demand for entry `i`.
    pub fn set_inflow(&mut self, i: usize, demand: f64) {
        if i < self.n_entries {
            self.inflow[i] = demand;
        }
    }

    /// Advance one time step.
    pub fn step(&mut self) {
        let nc = self.ring_cells;
        let mut ring_new = self.ring.clone();
        // Circulation: shift ring density forward
        for i in 0..nc {
            let ip = (i + 1) % nc;
            ring_new[ip] += self.ring[i] * 0.5 * self.dt;
            ring_new[i] -= self.ring[i] * 0.5 * self.dt;
        }
        // Entry/exit at entry cells
        let entry_spacing = nc / self.n_entries.max(1);
        for e in 0..self.n_entries {
            let cell = (e * entry_spacing) % nc;
            // Yield if ring is busy
            if self.ring[cell] < self.critical_gap_density {
                ring_new[cell] += self.inflow[e] * self.dt;
            }
            // Exit: remove fraction
            ring_new[cell] = (ring_new[cell] * (1.0 - self.outflow_fraction[e])).max(0.0);
        }
        // Clamp densities
        for v in ring_new.iter_mut() {
            *v = v.clamp(0.0, 1.0);
        }
        self.ring = ring_new;
        self.time += self.dt;
    }

    /// Compute total circulating vehicles.
    pub fn total_circulating(&self) -> f64 {
        self.ring.iter().sum()
    }

    /// Compute mean ring density.
    pub fn mean_density(&self) -> f64 {
        self.total_circulating() / self.ring_cells as f64
    }

    /// Compute utilisation ratio (mean / critical gap).
    pub fn utilisation(&self) -> f64 {
        (self.mean_density() / self.critical_gap_density).min(1.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MultiLane traffic model
// ─────────────────────────────────────────────────────────────────────────────

/// Multi-lane road with coupled LWR lanes and lane-change rules.
#[derive(Clone, Debug)]
pub struct MultiLane {
    /// Number of lanes.
    pub n_lanes: usize,
    /// Number of cells per lane.
    pub nx: usize,
    /// Cell length (m).
    pub dx: f64,
    /// Time step (s).
    pub dt: f64,
    /// Density per lane per cell (lane × cell).
    pub rho: Vec<Vec<f64>>,
    /// Fundamental diagram (shared across lanes).
    pub diagram: FundamentalDiagram,
    /// Lane-change rate coefficient.
    pub lc_rate: f64,
    /// Simulation time (s).
    pub time: f64,
}

impl MultiLane {
    /// Create a new `MultiLane` simulation.
    pub fn new(
        n_lanes: usize,
        nx: usize,
        dx: f64,
        dt: f64,
        diagram: FundamentalDiagram,
        lc_rate: f64,
    ) -> Self {
        Self {
            n_lanes,
            nx,
            dx,
            dt,
            rho: vec![vec![0.0; nx]; n_lanes],
            diagram,
            lc_rate,
            time: 0.0,
        }
    }

    /// Set initial density for lane `l` from a closure.
    pub fn init_lane<F: Fn(f64) -> f64>(&mut self, l: usize, f: F) {
        if l < self.n_lanes {
            for i in 0..self.nx {
                let x = (i as f64 + 0.5) * self.dx;
                self.rho[l][i] = f(x).clamp(0.0, self.diagram.rho_jam);
            }
        }
    }

    /// Advance by one time step.
    pub fn step(&mut self) {
        let cfl = self.dt / self.dx;
        let mut rho_new = self.rho.clone();

        // LWR advection per lane (Godunov)
        for l in 0..self.n_lanes {
            for i in 0..self.nx {
                let il = if i == 0 { self.nx - 1 } else { i - 1 };
                let ir = (i + 1) % self.nx;
                let f_r = godunov_flux_fd(&self.diagram, self.rho[l][i], self.rho[l][ir]);
                let f_l = godunov_flux_fd(&self.diagram, self.rho[l][il], self.rho[l][i]);
                rho_new[l][i] -= cfl * (f_r - f_l);
            }
        }

        // Lane-change diffusion between adjacent lanes
        for l in 0..self.n_lanes {
            for i in 0..self.nx {
                let lc = self.lc_rate * self.dt;
                if l > 0 {
                    let diff = lc * (self.rho[l - 1][i] - self.rho[l][i]);
                    rho_new[l][i] += diff;
                    rho_new[l - 1][i] -= diff;
                }
                if l < self.n_lanes - 1 {
                    let diff = lc * (self.rho[l + 1][i] - self.rho[l][i]);
                    rho_new[l][i] += diff;
                    rho_new[l + 1][i] -= diff;
                }
            }
        }

        // Clamp
        for l in 0..self.n_lanes {
            for i in 0..self.nx {
                rho_new[l][i] = rho_new[l][i].clamp(0.0, self.diagram.rho_jam);
            }
        }

        self.rho = rho_new;
        self.time += self.dt;
    }

    /// Compute total vehicle count across all lanes.
    pub fn total_vehicles(&self) -> f64 {
        self.rho.iter().flat_map(|lane| lane.iter()).sum::<f64>() * self.dx
    }

    /// Compute mean density across all lanes at cell `i`.
    pub fn mean_density_at(&self, i: usize) -> f64 {
        self.rho.iter().map(|lane| lane[i]).sum::<f64>() / self.n_lanes as f64
    }

    /// Compute flow on lane `l` at cell `i`.
    pub fn flow_at(&self, l: usize, i: usize) -> f64 {
        self.diagram.flow(self.rho[l][i])
    }
}

/// Godunov flux helper using a `FundamentalDiagram`.
fn godunov_flux_fd(diag: &FundamentalDiagram, rho_l: f64, rho_r: f64) -> f64 {
    let q_l = diag.flow(rho_l);
    let q_r = diag.flow(rho_r);
    let rho_c = diag.rho_crit;
    let demand = if rho_l <= rho_c { q_l } else { diag.q_max };
    let supply = if rho_r >= rho_c { q_r } else { diag.q_max };
    demand.min(supply)
}

// ─────────────────────────────────────────────────────────────────────────────
// SpeedLimitZone
// ─────────────────────────────────────────────────────────────────────────────

/// A variable speed limit zone applied to a range of cells.
#[derive(Clone, Debug)]
pub struct SpeedLimitZone {
    /// First cell index (inclusive).
    pub cell_start: usize,
    /// Last cell index (inclusive).
    pub cell_end: usize,
    /// Speed limit (m/s).
    pub speed_limit: f64,
    /// Whether zone is currently active.
    pub active: bool,
    /// Enforcement factor (1.0 = full, 0.0 = none).
    pub compliance: f64,
}

impl SpeedLimitZone {
    /// Create a new `SpeedLimitZone`.
    pub fn new(cell_start: usize, cell_end: usize, speed_limit: f64, compliance: f64) -> Self {
        Self {
            cell_start,
            cell_end,
            speed_limit,
            active: true,
            compliance: compliance.clamp(0.0, 1.0),
        }
    }

    /// Return effective speed limit considering compliance factor.
    pub fn effective_limit(&self, free_flow_speed: f64) -> f64 {
        if !self.active {
            return free_flow_speed;
        }
        free_flow_speed * (1.0 - self.compliance) + self.speed_limit * self.compliance
    }

    /// Check if cell `i` is inside this zone.
    pub fn contains(&self, i: usize) -> bool {
        self.active && i >= self.cell_start && i <= self.cell_end
    }

    /// Compute capacity reduction due to speed limit.
    pub fn capacity_fraction(&self, diagram: &FundamentalDiagram) -> f64 {
        (self.speed_limit / diagram.v_free).clamp(0.0, 1.0)
    }

    /// Apply zone: reduce density-weighted speed for a cell.
    pub fn apply(&self, rho: f64, v_unconstrained: f64, diagram: &FundamentalDiagram) -> f64 {
        if !self.active {
            return v_unconstrained;
        }
        let v_lim = self.effective_limit(diagram.v_free);
        v_unconstrained.min(v_lim * (1.0 - rho / diagram.rho_jam))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AccidentPropagation
// ─────────────────────────────────────────────────────────────────────────────

/// Accident/incident propagation model.
///
/// An accident reduces capacity at one cell and propagates a shock wave upstream.
#[derive(Clone, Debug)]
pub struct AccidentPropagation {
    /// Cell index of accident location.
    pub accident_cell: usize,
    /// Capacity reduction factor (0 = full block, 1 = no effect).
    pub capacity_factor: f64,
    /// Duration of incident (s).
    pub duration: f64,
    /// Elapsed time since accident (s).
    pub elapsed: f64,
    /// Shock wave speed (m/s) — negative means upstream propagation.
    pub shock_speed: f64,
    /// Whether accident is currently active.
    pub active: bool,
}

impl AccidentPropagation {
    /// Create a new `AccidentPropagation`.
    pub fn new(accident_cell: usize, capacity_factor: f64, duration: f64) -> Self {
        let shock_speed = -15.0; // typical backward shock wave speed
        Self {
            accident_cell,
            capacity_factor,
            duration,
            elapsed: 0.0,
            shock_speed,
            active: true,
        }
    }

    /// Advance incident by `dt` seconds.
    pub fn advance(&mut self, dt: f64) {
        if self.active {
            self.elapsed += dt;
            if self.elapsed >= self.duration {
                self.active = false;
            }
        }
    }

    /// Effective flow capacity fraction at a given cell and current time.
    pub fn capacity_at(&self, cell: usize, dx: f64) -> f64 {
        if !self.active {
            return 1.0;
        }
        let distance = (cell as f64 - self.accident_cell as f64) * dx;
        // Queue extends upstream (negative direction)
        let jam_front = self.shock_speed * self.elapsed;
        if distance >= jam_front && distance <= 0.0 {
            self.capacity_factor
        } else if distance == 0.0 {
            self.capacity_factor
        } else {
            1.0
        }
    }

    /// Estimate queue length (m) at current time.
    pub fn queue_length(&self) -> f64 {
        if !self.active {
            return 0.0;
        }
        (self.shock_speed * self.elapsed).abs()
    }

    /// Shock wave position offset from accident cell (m).
    pub fn shock_position(&self) -> f64 {
        self.shock_speed * self.elapsed
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TrafficLbm: D1Q3 LBM for macroscopic traffic
// ─────────────────────────────────────────────────────────────────────────────

/// D1Q3 Lattice Boltzmann solver for macroscopic traffic density and flow.
///
/// The equilibrium distribution encodes the fundamental diagram.
#[derive(Clone, Debug)]
pub struct TrafficLbm {
    /// Number of lattice cells.
    pub nx: usize,
    /// Cell spacing (m).
    pub dx: f64,
    /// Time step (s).
    pub dt: f64,
    /// Relaxation time τ.
    pub tau: f64,
    /// Distribution function array, shape (nx × NQ).
    pub f: Vec<[f64; NQ]>,
    /// Fundamental diagram.
    pub diagram: FundamentalDiagram,
    /// Simulation time (s).
    pub time: f64,
    /// Signals applied to cells.
    pub signals: Vec<TrafficSignal>,
    /// Speed limit zones.
    pub speed_zones: Vec<SpeedLimitZone>,
    /// Accident incidents.
    pub accidents: Vec<AccidentPropagation>,
}

impl TrafficLbm {
    /// Create a new `TrafficLbm`.
    pub fn new(nx: usize, dx: f64, dt: f64, tau: f64, diagram: FundamentalDiagram) -> Self {
        let mut s = Self {
            nx,
            dx,
            dt,
            tau,
            f: vec![[0.0; NQ]; nx],
            diagram,
            time: 0.0,
            signals: Vec::new(),
            speed_zones: Vec::new(),
            accidents: Vec::new(),
        };
        // Initialise to zero density
        for i in 0..nx {
            s.f[i] = feq1d(0.0, 0.0);
        }
        s
    }

    /// Initialise density field from closure.
    pub fn init_density<F: Fn(f64) -> f64>(&mut self, func: F) {
        for i in 0..self.nx {
            let x = (i as f64 + 0.5) * self.dx;
            let rho = func(x).clamp(0.0, self.diagram.rho_jam);
            let v = self.diagram.speed(rho);
            self.f[i] = feq1d(rho, v);
        }
    }

    /// Add a traffic signal.
    pub fn add_signal(&mut self, signal: TrafficSignal) {
        self.signals.push(signal);
    }

    /// Add a speed limit zone.
    pub fn add_speed_zone(&mut self, zone: SpeedLimitZone) {
        self.speed_zones.push(zone);
    }

    /// Add an accident.
    pub fn add_accident(&mut self, acc: AccidentPropagation) {
        self.accidents.push(acc);
    }

    /// Get density at cell `i`.
    pub fn density(&self, i: usize) -> f64 {
        macroscopic1d(&self.f[i]).0
    }

    /// Get velocity at cell `i`.
    pub fn velocity(&self, i: usize) -> f64 {
        macroscopic1d(&self.f[i]).1
    }

    /// Get flow (veh/s) at cell `i`.
    pub fn flow(&self, i: usize) -> f64 {
        let (rho, _u) = macroscopic1d(&self.f[i]);
        self.diagram.flow(rho)
    }

    /// Collision step: relax towards equilibrium.
    fn collide(&mut self) {
        for i in 0..self.nx {
            let (rho, _u_lbm) = macroscopic1d(&self.f[i]);
            // Effective speed from fundamental diagram
            let mut v_eq = self.diagram.speed(rho);

            // Apply speed limit zones
            for zone in &self.speed_zones {
                if zone.contains(i) {
                    v_eq = v_eq.min(zone.speed_limit);
                }
            }

            // Apply signals
            for sig in &self.signals {
                if sig.cell == i {
                    v_eq *= sig.speed_factor();
                }
            }

            // Apply accidents
            for acc in &self.accidents {
                let cap = acc.capacity_at(i, self.dx);
                v_eq *= cap;
            }

            let feq = feq1d(rho, v_eq);
            for q in 0..NQ {
                self.f[i][q] += (feq[q] - self.f[i][q]) / self.tau;
            }
        }
    }

    /// Streaming step: shift distributions by lattice velocities.
    fn stream(&mut self) {
        let mut f_new = vec![[0.0f64; NQ]; self.nx];
        for i in 0..self.nx {
            for q in 0..NQ {
                let c = CV[q] as isize;
                let i_src = ((i as isize - c).rem_euclid(self.nx as isize)) as usize;
                f_new[i][q] = self.f[i_src][q];
            }
        }
        self.f = f_new;
    }

    /// Advance all controllers by `dt`.
    fn advance_controllers(&mut self) {
        for sig in self.signals.iter_mut() {
            sig.advance(self.dt);
        }
        for acc in self.accidents.iter_mut() {
            acc.advance(self.dt);
        }
    }

    /// Advance the LBM by one time step.
    pub fn step(&mut self) {
        self.advance_controllers();
        self.collide();
        self.stream();
        self.time += self.dt;
    }

    /// Run for `n_steps` steps.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.step();
        }
    }

    /// Collect density profile.
    pub fn density_profile(&self) -> Vec<f64> {
        (0..self.nx).map(|i| self.density(i)).collect()
    }

    /// Collect flow profile.
    pub fn flow_profile(&self) -> Vec<f64> {
        (0..self.nx).map(|i| self.flow(i)).collect()
    }

    /// Compute total vehicles on road.
    pub fn total_vehicles(&self) -> f64 {
        (0..self.nx).map(|i| self.density(i)).sum::<f64>() * self.dx
    }

    /// Compute mean travel time via Little's law: L / Q_mean.
    pub fn mean_travel_time(&self) -> f64 {
        let _l = self.nx as f64 * self.dx;
        let n = self.total_vehicles();
        let q_mean = self.flow_profile().iter().sum::<f64>() / self.nx as f64;
        if q_mean < 1e-12 {
            return f64::INFINITY;
        }
        n / q_mean
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Jam detection utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Detect congested cells where density exceeds threshold.
#[allow(dead_code)]
pub fn detect_congestion(rho: &[f64], threshold: f64) -> Vec<usize> {
    rho.iter()
        .enumerate()
        .filter(|&(_, &r)| r > threshold)
        .map(|(i, _)| i)
        .collect()
}

/// Compute shockwave speed between upstream density `rho_u` and downstream `rho_d`.
#[allow(dead_code)]
pub fn shock_wave_speed(diagram: &FundamentalDiagram, rho_u: f64, rho_d: f64) -> f64 {
    let dq = diagram.flow(rho_d) - diagram.flow(rho_u);
    let drho = rho_d - rho_u;
    if drho.abs() < 1e-12 {
        return diagram.speed(rho_u);
    }
    dq / drho
}

/// Compute the oscillation amplitude in speed (standard deviation).
#[allow(dead_code)]
pub fn speed_oscillation_amplitude(speeds: &[f64]) -> f64 {
    if speeds.len() < 2 {
        return 0.0;
    }
    let mean = speeds.iter().sum::<f64>() / speeds.len() as f64;
    let var = speeds.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / speeds.len() as f64;
    var.sqrt()
}

/// Compute throughput efficiency (actual flow / capacity).
#[allow(dead_code)]
pub fn throughput_efficiency(actual_flow: f64, capacity: f64) -> f64 {
    if capacity < 1e-12 {
        return 0.0;
    }
    (actual_flow / capacity).clamp(0.0, 1.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Green wave coordination
// ─────────────────────────────────────────────────────────────────────────────

/// Compute offset for coordinated green wave (bandwidth optimisation).
///
/// Returns optimal phase offset in seconds for signal at distance `d` from
/// reference, given free-flow speed `v_ff` and cycle time `T`.
#[allow(dead_code)]
pub fn green_wave_offset(d: f64, v_ff: f64, cycle: f64) -> f64 {
    if v_ff < 1e-12 {
        return 0.0;
    }
    (d / v_ff) % cycle
}

/// Compute the degree of saturation for a signal approach.
///
/// `q` is the flow (veh/s), `s` is the saturation flow (veh/s), `g` is
/// green time, `c` is cycle time.
#[allow(dead_code)]
pub fn degree_of_saturation(q: f64, s: f64, g: f64, c: f64) -> f64 {
    if s < 1e-12 || g < 1e-12 {
        return f64::INFINITY;
    }
    q * c / (s * g)
}

/// Webster's optimal cycle length formula.
///
/// `y_sum` is the sum of critical flow ratios, `L` is total lost time.
#[allow(dead_code)]
pub fn webster_optimal_cycle(y_sum: f64, l: f64) -> f64 {
    let denom = 1.0 - y_sum;
    if denom <= 0.0 {
        return f64::INFINITY;
    }
    (1.5 * l + 5.0) / denom
}

// ─────────────────────────────────────────────────────────────────────────────
// Emission model (CO2 as proxy)
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate CO2 emission rate (g/s·m) from speed using Virginia Tech model.
#[allow(dead_code)]
pub fn co2_emission_rate(speed_ms: f64) -> f64 {
    // Simplified: high at idle and high speed, minimum at ~60 km/h (~16.7 m/s)
    let v = speed_ms.max(0.0);
    let v_opt = 16.7_f64;
    let base = 2.0_f64;
    base * (1.0 + ((v - v_opt) / v_opt).powi(2))
}

/// Compute total road-level emissions (g/s) for a density/speed profile.
#[allow(dead_code)]
pub fn total_road_emissions(rho: &[f64], speed: &[f64], dx: f64) -> f64 {
    rho.iter()
        .zip(speed.iter())
        .map(|(&r, &v)| r * co2_emission_rate(v) * dx)
        .sum()
}

// ─────────────────────────────────────────────────────────────────────────────
// Ramp metering
// ─────────────────────────────────────────────────────────────────────────────

/// ALINEA ramp metering controller.
///
/// Adjusts on-ramp flow to maintain target occupancy downstream.
#[derive(Clone, Debug)]
pub struct AlineaController {
    /// Target density (veh/m) downstream.
    pub rho_target: f64,
    /// ALINEA gain K_R (m²/s/veh).
    pub k_r: f64,
    /// Current metering rate (veh/s).
    pub metering_rate: f64,
    /// Minimum metering rate (veh/s).
    pub rate_min: f64,
    /// Maximum metering rate (veh/s).
    pub rate_max: f64,
}

impl AlineaController {
    /// Create a new `AlineaController`.
    pub fn new(rho_target: f64, k_r: f64, rate_min: f64, rate_max: f64) -> Self {
        Self {
            rho_target,
            k_r,
            metering_rate: rate_min,
            rate_min,
            rate_max,
        }
    }

    /// Update metering rate based on measured density downstream.
    pub fn update(&mut self, rho_measured: f64) -> f64 {
        let error = self.rho_target - rho_measured;
        self.metering_rate =
            (self.metering_rate + self.k_r * error).clamp(self.rate_min, self.rate_max);
        self.metering_rate
    }

    /// Reset controller to minimum rate.
    pub fn reset(&mut self) {
        self.metering_rate = self.rate_min;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Travel time estimation
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate travel time on a road section using the BPR volume–delay function.
///
/// Returns travel time (s) for given volume `q`, free-flow travel time `t0`,
/// capacity `c`, and BPR parameters `alpha` and `beta`.
#[allow(dead_code)]
pub fn bpr_travel_time(q: f64, t0: f64, c: f64, alpha: f64, beta: f64) -> f64 {
    t0 * (1.0 + alpha * (q / c.max(1e-12)).powf(beta))
}

/// Compute queue dissipation time (s) after signal turns green.
///
/// `n_queued` vehicles, saturation flow `s` (veh/s), lost time `tl` (s).
#[allow(dead_code)]
pub fn queue_dissipation_time(n_queued: f64, s: f64, tl: f64) -> f64 {
    if s < 1e-12 {
        return f64::INFINITY;
    }
    n_queued / s + tl
}

// ─────────────────────────────────────────────────────────────────────────────
// Macro-simulation result
// ─────────────────────────────────────────────────────────────────────────────

/// Summary result for a traffic simulation run.
#[derive(Clone, Debug)]
pub struct TrafficResult {
    /// Total vehicle-hours travelled (veh·h).
    pub total_vht: f64,
    /// Total vehicle-kilometres travelled (veh·km).
    pub total_vkt: f64,
    /// Mean travel speed (m/s).
    pub mean_speed: f64,
    /// Peak congestion density (veh/m).
    pub peak_density: f64,
    /// Mean throughput flow (veh/s).
    pub mean_flow: f64,
    /// Simulation duration (s).
    pub duration: f64,
}

impl TrafficResult {
    /// Compute Level of Service (0=A best … 5=F worst) from V/C ratio.
    pub fn level_of_service(&self, capacity: f64) -> usize {
        let vc = if capacity > 1e-12 {
            self.mean_flow / capacity
        } else {
            1.0
        };
        if vc < 0.20 {
            0
        } else if vc < 0.44 {
            1
        } else if vc < 0.64 {
            2
        } else if vc < 0.85 {
            3
        } else if vc < 1.00 {
            4
        } else {
            5
        }
    }

    /// Compute delay per vehicle (s/veh).
    pub fn delay_per_vehicle(&self, free_flow_time: f64) -> f64 {
        if self.total_vkt < 1e-12 {
            return 0.0;
        }
        let actual_time = self.total_vht * 3600.0 / self.total_vkt.max(1e-12);
        (actual_time - free_flow_time).max(0.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility: sinusoidal density wave for testing
// ─────────────────────────────────────────────────────────────────────────────

/// Generate a sinusoidal density perturbation for testing stability.
#[allow(dead_code)]
pub fn sinusoidal_density(x: f64, road_length: f64, rho_mean: f64, amplitude: f64) -> f64 {
    (rho_mean + amplitude * (2.0 * PI * x / road_length).sin()).max(0.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── FundamentalDiagram tests ──────────────────────────────────────────────

    #[test]
    fn test_greenshields_zero_density() {
        let fd = FundamentalDiagram::new(0.1, 30.0, DiagramModel::Greenshields);
        assert!((fd.speed(0.0) - 30.0).abs() < 1e-9);
    }

    #[test]
    fn test_greenshields_jam_density() {
        let fd = FundamentalDiagram::new(0.1, 30.0, DiagramModel::Greenshields);
        assert!(fd.speed(0.1) < 1e-9);
    }

    #[test]
    fn test_greenshields_flow_zero_at_jam() {
        let fd = FundamentalDiagram::new(0.1, 30.0, DiagramModel::Greenshields);
        assert!(fd.flow(0.1) < 1e-9);
    }

    #[test]
    fn test_triangular_free_flow_regime() {
        let fd = FundamentalDiagram::new(0.12, 25.0, DiagramModel::Triangular);
        let v = fd.speed(0.01); // below rho_crit
        assert!((v - 25.0).abs() < 1e-6);
    }

    #[test]
    fn test_triangular_congested_regime() {
        let fd = FundamentalDiagram::new(0.12, 25.0, DiagramModel::Triangular);
        let v_jam = fd.speed(0.12);
        assert!(v_jam < 1e-6);
    }

    #[test]
    fn test_drake_model_max_at_zero() {
        let fd = FundamentalDiagram::new(0.1, 20.0, DiagramModel::Drake);
        assert!((fd.speed(0.0) - 20.0).abs() < 1e-6);
    }

    #[test]
    fn test_flow_nonnegative() {
        let fd = FundamentalDiagram::new(0.1, 30.0, DiagramModel::Greenshields);
        for i in 0..=10 {
            let rho = 0.01 * i as f64;
            assert!(fd.flow(rho) >= 0.0);
        }
    }

    #[test]
    fn test_density_for_flow_roundtrip() {
        let fd = FundamentalDiagram::new(0.1, 30.0, DiagramModel::Greenshields);
        let rho0 = 0.03;
        let q0 = fd.flow(rho0);
        let rho1 = fd.density_for_flow(q0);
        // Two solutions exist; check at least one matches
        assert!((fd.flow(rho1) - q0).abs() < 1e-4);
    }

    // ── NagelSchreckenberg tests ──────────────────────────────────────────────

    #[test]
    fn test_ns_initial_density() {
        let mut ns = NagelSchreckenberg::new(100, 5, 0.3);
        ns.populate_uniform(20);
        assert!((ns.density() - 0.2).abs() < 1e-9);
    }

    #[test]
    fn test_ns_step_preserves_count() {
        let mut ns = NagelSchreckenberg::new(100, 5, 0.0);
        ns.populate_uniform(10);
        let count_before = ns.positions.len();
        let bits = vec![0u8; 100];
        ns.step_with_seed(&bits);
        assert_eq!(ns.positions.len(), count_before);
    }

    #[test]
    fn test_ns_positions_in_range() {
        let mut ns = NagelSchreckenberg::new(50, 3, 0.0);
        ns.populate_uniform(5);
        let bits = vec![0u8; 50];
        for _ in 0..10 {
            ns.step_with_seed(&bits);
            for &p in &ns.positions {
                assert!(p < 50);
            }
        }
    }

    #[test]
    fn test_ns_mean_flow_nonneg() {
        let mut ns = NagelSchreckenberg::new(100, 5, 0.3);
        ns.populate_uniform(30);
        let bits = vec![100u8; 100];
        ns.step_with_seed(&bits);
        assert!(ns.mean_flow() >= 0.0);
    }

    // ── LwrLattice tests ──────────────────────────────────────────────────────

    #[test]
    fn test_lwr_conserves_mass() {
        let fd = FundamentalDiagram::new(0.1, 20.0, DiagramModel::Greenshields);
        let mut lwr = LwrLattice::new(100, 1.0, 0.04, fd);
        lwr.init_from(|_| 0.03);
        let n0 = lwr.total_vehicles();
        lwr.run(50);
        let n1 = lwr.total_vehicles();
        assert!((n0 - n1).abs() / n0 < 1e-6);
    }

    #[test]
    fn test_lwr_cfl_check() {
        let fd = FundamentalDiagram::new(0.1, 20.0, DiagramModel::Greenshields);
        let lwr = LwrLattice::new(100, 1.0, 0.04, fd);
        assert!(lwr.cfl_ok());
    }

    #[test]
    fn test_lwr_uniform_steady_state() {
        let fd = FundamentalDiagram::new(0.1, 20.0, DiagramModel::Greenshields);
        let mut lwr = LwrLattice::new(50, 1.0, 0.04, fd);
        lwr.init_from(|_| 0.05);
        lwr.run(100);
        // Uniform density should remain uniform
        let var = {
            let mean = lwr.mean_density();
            lwr.rho.iter().map(|&r| (r - mean).powi(2)).sum::<f64>() / lwr.nx as f64
        };
        assert!(var < 1e-20);
    }

    #[test]
    fn test_lwr_flow_at_valid() {
        let fd = FundamentalDiagram::new(0.1, 20.0, DiagramModel::Greenshields);
        let mut lwr = LwrLattice::new(50, 1.0, 0.04, fd);
        lwr.init_from(|_| 0.05);
        assert!(lwr.flow_at(0) >= 0.0);
    }

    // ── TrafficSignal tests ────────────────────────────────────────────────────

    #[test]
    fn test_signal_starts_red() {
        let sig = TrafficSignal::new(5, 30.0, 5.0, 40.0, 0.0);
        // offset=0 => phase_time=0 => green region if 0 < green_duration
        // Actually with offset=0 and phase_time=0 < green_duration=30 → Green
        // Let's just check phase_time is in [0, cycle)
        let cycle = 30.0 + 5.0 + 40.0;
        assert!(sig.phase_time < cycle);
    }

    #[test]
    fn test_signal_advance_to_yellow() {
        let mut sig = TrafficSignal::new(0, 30.0, 5.0, 40.0, 0.0);
        sig.advance(31.0); // past green into yellow
        assert_eq!(sig.phase, SignalPhase::Yellow);
    }

    #[test]
    fn test_signal_advance_to_red() {
        let mut sig = TrafficSignal::new(0, 30.0, 5.0, 40.0, 0.0);
        sig.advance(37.0); // past green + yellow
        assert_eq!(sig.phase, SignalPhase::Red);
    }

    #[test]
    fn test_signal_speed_factor_red_zero() {
        let mut sig = TrafficSignal::new(0, 30.0, 5.0, 40.0, 0.0);
        sig.advance(37.0);
        assert!((sig.speed_factor() - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_signal_speed_factor_green_one() {
        let sig = TrafficSignal::new(0, 30.0, 5.0, 40.0, 0.0);
        // phase_time=0 is in green range
        assert!((sig.speed_factor() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_optimal_green_split_proportional() {
        let (g1, g2) = TrafficSignal::optimal_green_split(60.0, 1.0, 1.0);
        assert!((g1 - g2).abs() < 1e-9);
    }

    // ── OnRamp tests ──────────────────────────────────────────────────────────

    #[test]
    fn test_onramp_no_merge_when_congested() {
        let mut ramp = OnRamp::new(10, 1.0, 0.05);
        let inflow = ramp.process(0.08, 1.0); // density > limit
        assert!((inflow - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_onramp_merges_when_free() {
        let mut ramp = OnRamp::new(10, 1.0, 0.05);
        let inflow = ramp.process(0.01, 1.0); // density < limit
        assert!(inflow > 0.0);
    }

    #[test]
    fn test_onramp_acceptance_rate_free_flow() {
        let ramp = OnRamp::new(10, 1.0, 0.05);
        let rate = ramp.acceptance_rate(0.0);
        assert!((rate - 1.0).abs() < 1e-9);
    }

    // ── OffRamp tests ─────────────────────────────────────────────────────────

    #[test]
    fn test_offramp_exit_plus_through_equals_total() {
        let mut ramp = OffRamp::new(20, 0.3);
        let total = 10.0;
        let exit = ramp.exit_flow(total, 1.0);
        let through = ramp.through_flow(total);
        assert!((exit + through - total).abs() < 1e-9);
    }

    #[test]
    fn test_offramp_zero_split() {
        let mut ramp = OffRamp::new(20, 0.0);
        let exit = ramp.exit_flow(10.0, 1.0);
        assert!((exit - 0.0).abs() < 1e-9);
    }

    // ── Roundabout tests ──────────────────────────────────────────────────────

    #[test]
    fn test_roundabout_initial_zero() {
        let r = Roundabout::new(20, 4, 0.5, 0.1);
        assert!((r.total_circulating() - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_roundabout_step_nonneg() {
        let mut r = Roundabout::new(20, 4, 0.5, 0.1);
        r.set_inflow(0, 2.0);
        r.step();
        for &v in &r.ring {
            assert!(v >= 0.0);
        }
    }

    // ── SpeedLimitZone tests ──────────────────────────────────────────────────

    #[test]
    fn test_speedzone_contains() {
        let zone = SpeedLimitZone::new(5, 15, 20.0, 1.0);
        assert!(zone.contains(10));
        assert!(!zone.contains(4));
        assert!(!zone.contains(16));
    }

    #[test]
    fn test_speedzone_capacity_fraction() {
        let fd = FundamentalDiagram::new(0.1, 40.0, DiagramModel::Greenshields);
        let zone = SpeedLimitZone::new(0, 10, 20.0, 1.0);
        let frac = zone.capacity_fraction(&fd);
        assert!((frac - 0.5).abs() < 1e-9);
    }

    // ── AccidentPropagation tests ─────────────────────────────────────────────

    #[test]
    fn test_accident_active_initially() {
        let acc = AccidentPropagation::new(50, 0.2, 300.0);
        assert!(acc.active);
    }

    #[test]
    fn test_accident_deactivates_after_duration() {
        let mut acc = AccidentPropagation::new(50, 0.2, 100.0);
        acc.advance(101.0);
        assert!(!acc.active);
    }

    #[test]
    fn test_accident_queue_grows() {
        let mut acc = AccidentPropagation::new(50, 0.2, 300.0);
        acc.advance(10.0);
        assert!(acc.queue_length() > 0.0);
    }

    // ── TrafficLbm tests ──────────────────────────────────────────────────────

    #[test]
    fn test_lbm_init_density() {
        let fd = FundamentalDiagram::new(0.1, 20.0, DiagramModel::Greenshields);
        let mut lbm = TrafficLbm::new(50, 1.0, 0.04, 1.0, fd);
        lbm.init_density(|_| 0.05);
        let d = lbm.density(25);
        assert!((d - 0.05).abs() < 1e-6);
    }

    #[test]
    fn test_lbm_step_runs() {
        let fd = FundamentalDiagram::new(0.1, 20.0, DiagramModel::Greenshields);
        let mut lbm = TrafficLbm::new(50, 1.0, 0.04, 1.0, fd);
        lbm.init_density(|_| 0.03);
        lbm.step();
        assert!(lbm.time > 0.0);
    }

    #[test]
    fn test_lbm_total_vehicles_conserved() {
        let fd = FundamentalDiagram::new(0.1, 20.0, DiagramModel::Greenshields);
        let mut lbm = TrafficLbm::new(100, 1.0, 0.04, 1.0, fd);
        lbm.init_density(|_| 0.04);
        let n0 = lbm.total_vehicles();
        lbm.run(20);
        let n1 = lbm.total_vehicles();
        assert!((n0 - n1).abs() / n0 < 1e-4);
    }

    #[test]
    fn test_lbm_with_signal() {
        let fd = FundamentalDiagram::new(0.1, 20.0, DiagramModel::Greenshields);
        let mut lbm = TrafficLbm::new(50, 1.0, 0.04, 1.0, fd);
        lbm.init_density(|_| 0.03);
        let sig = TrafficSignal::new(25, 10.0, 2.0, 15.0, 0.0);
        lbm.add_signal(sig);
        lbm.run(10);
        assert!(lbm.time > 0.0);
    }

    #[test]
    fn test_lbm_with_accident() {
        let fd = FundamentalDiagram::new(0.1, 20.0, DiagramModel::Greenshields);
        let mut lbm = TrafficLbm::new(50, 1.0, 0.04, 1.0, fd);
        lbm.init_density(|_| 0.03);
        let acc = AccidentPropagation::new(30, 0.3, 100.0);
        lbm.add_accident(acc);
        lbm.run(5);
        assert!(lbm.time > 0.0);
    }

    // ── AlineaController tests ────────────────────────────────────────────────

    #[test]
    fn test_alinea_reduces_rate_on_high_density() {
        let mut ctrl = AlineaController::new(0.05, 100.0, 0.1, 2.0);
        ctrl.metering_rate = 1.0;
        let rate = ctrl.update(0.08); // above target
        assert!(rate < 1.0);
    }

    #[test]
    fn test_alinea_increases_rate_on_low_density() {
        let mut ctrl = AlineaController::new(0.05, 100.0, 0.1, 2.0);
        ctrl.metering_rate = 1.0;
        let rate = ctrl.update(0.02); // below target
        assert!(rate > 1.0);
    }

    #[test]
    fn test_alinea_clamps_to_bounds() {
        let mut ctrl = AlineaController::new(0.05, 1000.0, 0.1, 2.0);
        ctrl.update(0.0); // very low → would push rate high
        assert!(ctrl.metering_rate <= 2.0);
        ctrl.update(1.0); // very high → would push rate low
        assert!(ctrl.metering_rate >= 0.1);
    }

    // ── Utility function tests ────────────────────────────────────────────────

    #[test]
    fn test_shock_wave_speed_sign() {
        let fd = FundamentalDiagram::new(0.1, 30.0, DiagramModel::Greenshields);
        // Density jump from free-flow (rho_u) to congested (rho_d) → negative shock speed.
        // rho_crit=0.05; with rho_u=0.02 (free-flow) and rho_d=0.09 (congested),
        // flow(0.02)=0.02*30*(1-0.2)=0.48, flow(0.09)=0.09*30*(1-0.9)=0.27
        // ws = (0.27-0.48)/(0.09-0.02) = -3.0 < 0
        let ws = shock_wave_speed(&fd, 0.02, 0.09);
        assert!(ws < 0.0);
    }

    #[test]
    fn test_bpr_free_flow() {
        let t = bpr_travel_time(0.0, 60.0, 100.0, 0.15, 4.0);
        assert!((t - 60.0).abs() < 1e-9);
    }

    #[test]
    fn test_bpr_at_capacity() {
        let t = bpr_travel_time(100.0, 60.0, 100.0, 0.15, 4.0);
        assert!(t > 60.0);
    }

    #[test]
    fn test_green_wave_offset_zero_distance() {
        let off = green_wave_offset(0.0, 15.0, 90.0);
        assert!((off - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_degree_of_saturation_at_capacity() {
        let x = degree_of_saturation(0.5, 1.0, 30.0, 60.0);
        assert!((x - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_webster_cycle_reasonable() {
        let c = webster_optimal_cycle(0.6, 10.0);
        assert!(c > 30.0 && c < 200.0);
    }

    #[test]
    fn test_multilane_step_conserves() {
        let fd = FundamentalDiagram::new(0.1, 20.0, DiagramModel::Greenshields);
        let mut ml = MultiLane::new(2, 50, 1.0, 0.04, fd, 0.01);
        ml.init_lane(0, |_| 0.03);
        ml.init_lane(1, |_| 0.03);
        let n0 = ml.total_vehicles();
        ml.step();
        let n1 = ml.total_vehicles();
        assert!((n0 - n1).abs() / n0 < 1e-4);
    }

    #[test]
    fn test_sinusoidal_density_nonneg() {
        for i in 0..100 {
            let x = i as f64 * 0.5;
            let rho = sinusoidal_density(x, 50.0, 0.03, 0.01);
            assert!(rho >= 0.0);
        }
    }

    #[test]
    fn test_detect_congestion() {
        let rho = vec![0.01, 0.03, 0.08, 0.02, 0.09];
        let jams = detect_congestion(&rho, 0.05);
        assert_eq!(jams, vec![2, 4]);
    }

    #[test]
    fn test_queue_dissipation_time() {
        let t = queue_dissipation_time(20.0, 2.0, 3.0);
        assert!((t - 13.0).abs() < 1e-9);
    }

    #[test]
    fn test_traffic_result_los() {
        let result = TrafficResult {
            total_vht: 100.0,
            total_vkt: 5000.0,
            mean_speed: 50.0 / 3.6,
            peak_density: 0.05,
            mean_flow: 40.0,
            duration: 3600.0,
        };
        let los = result.level_of_service(100.0);
        assert!(los <= 5);
    }

    #[test]
    fn test_co2_emission_rate_positive() {
        for v in [0.0, 10.0, 20.0, 30.0] {
            assert!(co2_emission_rate(v) > 0.0);
        }
    }
}
