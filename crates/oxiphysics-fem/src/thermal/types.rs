//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Convection boundary condition: q = h (T - T_inf).
///
/// Adds a convective heat transfer term at a boundary node.
#[derive(Debug, Clone)]
pub struct ConvectionBC {
    /// Node index where convection is applied.
    pub node: usize,
    /// Convection coefficient h \[W/(m^2 K)\].
    pub h: f64,
    /// Area associated with this BC \[m^2\].
    pub area: f64,
    /// Ambient (fluid) temperature T_inf \[K\].
    pub t_inf: f64,
}
impl ConvectionBC {
    /// Create a new convection BC.
    pub fn new(node: usize, h: f64, area: f64, t_inf: f64) -> Self {
        Self {
            node,
            h,
            area,
            t_inf,
        }
    }
    /// Compute the convective heat flux at the current temperature.
    ///
    /// q = h * A * (T - T_inf)
    pub fn heat_flux(&self, temperature: f64) -> f64 {
        self.h * self.area * (temperature - self.t_inf)
    }
    /// Contribution to the conductance matrix (diagonal entry).
    pub fn stiffness_contribution(&self) -> f64 {
        self.h * self.area
    }
    /// Contribution to the load vector.
    pub fn load_contribution(&self) -> f64 {
        self.h * self.area * self.t_inf
    }
}
/// Radiation boundary condition: q = epsilon * sigma * A * (T^4 - T_surr^4).
///
/// Uses the Stefan-Boltzmann law for thermal radiation.
#[derive(Debug, Clone)]
pub struct RadiationBC {
    /// Node index where radiation is applied.
    pub node: usize,
    /// Surface emissivity epsilon (0 to 1).
    pub emissivity: f64,
    /// Stefan-Boltzmann constant \[W/(m^2 K^4)\].
    pub stefan_boltzmann: f64,
    /// Area associated with this BC \[m^2\].
    pub area: f64,
    /// Surrounding temperature T_surr \[K\].
    pub t_surr: f64,
}
impl RadiationBC {
    /// Create a new radiation BC with the standard Stefan-Boltzmann constant.
    pub fn new(node: usize, emissivity: f64, area: f64, t_surr: f64) -> Self {
        Self {
            node,
            emissivity,
            stefan_boltzmann: 5.670374419e-8,
            area,
            t_surr,
        }
    }
    /// Compute the radiative heat flux at the current temperature.
    ///
    /// q = epsilon * sigma * A * (T^4 - T_surr^4)
    pub fn heat_flux(&self, temperature: f64) -> f64 {
        self.emissivity
            * self.stefan_boltzmann
            * self.area
            * (temperature.powi(4) - self.t_surr.powi(4))
    }
    /// Linearized radiation coefficient (for Newton-Raphson).
    ///
    /// h_rad = 4 * epsilon * sigma * T^3
    pub fn linearized_coefficient(&self, temperature: f64) -> f64 {
        4.0 * self.emissivity * self.stefan_boltzmann * temperature.powi(3)
    }
    /// Linearized stiffness contribution at temperature T.
    pub fn linearized_stiffness(&self, temperature: f64) -> f64 {
        self.linearized_coefficient(temperature) * self.area
    }
    /// Linearized load contribution.
    pub fn linearized_load(&self, temperature: f64) -> f64 {
        let h_rad = self.linearized_coefficient(temperature);
        h_rad * self.area * self.t_surr
    }
}
/// Material/geometry properties for a single thermal element.
pub struct ThermalElement {
    /// Thermal conductivity k \[W/(m K)\].
    pub conductivity: f64,
    /// Mass density rho \[kg/m^3\].
    pub density: f64,
    /// Specific heat capacity cp \[J/(kg K)\].
    pub specific_heat: f64,
    /// Element volume V \[m^3\].
    pub volume: f64,
}
impl ThermalElement {
    /// Create a new `ThermalElement`.
    pub fn new(k: f64, rho: f64, cp: f64, vol: f64) -> Self {
        Self {
            conductivity: k,
            density: rho,
            specific_heat: cp,
            volume: vol,
        }
    }
    /// Thermal diffusivity alpha = k / (rho cp) \[m^2/s\].
    pub fn diffusivity(&self) -> f64 {
        self.conductivity / (self.density * self.specific_heat)
    }
    /// Lumped thermal capacitance C = rho cp V \[J/K\].
    pub fn capacitance(&self) -> f64 {
        self.density * self.specific_heat * self.volume
    }
    /// 1-D rod conductance G = k A / L \[W/K\].
    pub fn conductance(&self, area: f64, length: f64) -> f64 {
        self.conductivity * area / length
    }
}
/// Thermal contact resistance between two surfaces.
///
/// The heat flux across the contact is:
/// q = (T1 - T2) / R_tc
///
/// where R_tc = 1 / (h_c * A) is the thermal contact resistance.
#[derive(Debug, Clone)]
pub struct ThermalContactResistance {
    /// Node index on surface 1.
    pub node_1: usize,
    /// Node index on surface 2.
    pub node_2: usize,
    /// Contact conductance h_c \[W/(m^2 K)\].
    pub contact_conductance: f64,
    /// Contact area \[m^2\].
    pub area: f64,
}
impl ThermalContactResistance {
    /// Create a new thermal contact resistance.
    pub fn new(node_1: usize, node_2: usize, h_c: f64, area: f64) -> Self {
        Self {
            node_1,
            node_2,
            contact_conductance: h_c,
            area,
        }
    }
    /// Compute the heat flux from node_1 to node_2.
    pub fn heat_flux(&self, t1: f64, t2: f64) -> f64 {
        self.contact_conductance * self.area * (t1 - t2)
    }
    /// Contact conductance G_c = h_c * A \[W/K\].
    pub fn conductance(&self) -> f64 {
        self.contact_conductance * self.area
    }
    /// Thermal resistance R_tc = 1 / G_c \[K/W\].
    pub fn resistance(&self) -> f64 {
        1.0 / self.conductance()
    }
}
/// Anisotropic thermal conductivity tensor (3x3 symmetric, W/(m K)).
///
/// For orthotropic materials: k_xy = k_xz = k_yz = 0.
#[derive(Debug, Clone)]
pub struct ThermalConductivityTensor {
    /// Conductivity matrix k\[3\]\[3\] (W/(m K)).
    pub k: [[f64; 3]; 3],
}
impl ThermalConductivityTensor {
    /// Create an isotropic conductivity tensor: k * I.
    pub fn isotropic(k: f64) -> Self {
        let mut mat = [[0.0_f64; 3]; 3];
        for (i, row) in mat.iter_mut().enumerate() {
            row[i] = k;
        }
        Self { k: mat }
    }
    /// Create an orthotropic conductivity tensor.
    ///
    /// # Arguments
    /// * `kx` - conductivity along x (W/(m K))
    /// * `ky` - conductivity along y (W/(m K))
    /// * `kz` - conductivity along z (W/(m K))
    pub fn orthotropic(kx: f64, ky: f64, kz: f64) -> Self {
        let k = [[kx, 0.0, 0.0], [0.0, ky, 0.0], [0.0, 0.0, kz]];
        Self { k }
    }
    /// Compute heat flux vector q = -k * grad_T.
    ///
    /// # Arguments
    /// * `grad_t` - temperature gradient \[dT/dx, dT/dy, dT/dz\]
    pub fn heat_flux(&self, grad_t: &[f64; 3]) -> [f64; 3] {
        let mut q = [0.0_f64; 3];
        for (q_i, k_row) in q.iter_mut().zip(self.k.iter()) {
            *q_i -= k_row
                .iter()
                .zip(grad_t.iter())
                .map(|(&k, &g)| k * g)
                .sum::<f64>();
        }
        q
    }
    /// Maximum (dominant) eigenvalue of the conductivity tensor (approximate, for diagonal).
    pub fn max_conductivity(&self) -> f64 {
        self.k[0][0].max(self.k[1][1]).max(self.k[2][2])
    }
    /// Effective isotropic conductivity: geometric mean of diagonal terms.
    pub fn effective_isotropic(&self) -> f64 {
        let product = self.k[0][0] * self.k[1][1] * self.k[2][2];
        if product <= 0.0 {
            return 0.0;
        }
        product.powf(1.0 / 3.0)
    }
}
/// 1-D FEM thermal mesh.
pub struct ThermalMesh1D {
    /// Number of nodes.
    pub n_nodes: usize,
    /// Current nodal temperatures \[K\].
    pub temperatures: Vec<f64>,
    /// One element per inter-node gap (n_nodes - 1 elements).
    pub elements: Vec<ThermalElement>,
    /// Physical position of each node along the rod \[m\].
    pub node_positions: Vec<f64>,
    /// Dirichlet boundary conditions: (node_index, temperature).
    pub(super) bc: Vec<(usize, f64)>,
    /// Convection BCs.
    pub(super) convection_bcs: Vec<ConvectionBC>,
    /// Radiation BCs.
    pub(super) radiation_bcs: Vec<RadiationBC>,
    /// Thermal contact resistances.
    pub(super) contact_resistances: Vec<ThermalContactResistance>,
}
impl ThermalMesh1D {
    /// Construct a uniform 1-D mesh.
    pub fn new_uniform(n_nodes: usize, length: f64, conductivity: f64, rho: f64, cp: f64) -> Self {
        assert!(n_nodes >= 2, "need at least 2 nodes");
        let n_elem = n_nodes - 1;
        let dx = length / n_elem as f64;
        let vol = dx;
        let node_positions: Vec<f64> = (0..n_nodes).map(|i| i as f64 * dx).collect();
        let temperatures = vec![0.0_f64; n_nodes];
        let elements: Vec<ThermalElement> = (0..n_elem)
            .map(|_| ThermalElement::new(conductivity, rho, cp, vol))
            .collect();
        Self {
            n_nodes,
            temperatures,
            elements,
            node_positions,
            bc: Vec::new(),
            convection_bcs: Vec::new(),
            radiation_bcs: Vec::new(),
            contact_resistances: Vec::new(),
        }
    }
    /// Assemble the global conductance (stiffness) matrix K.
    pub fn assemble_conductance_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.n_nodes;
        let mut k_global = vec![vec![0.0_f64; n]; n];
        for (e, elem) in self.elements.iter().enumerate() {
            let i = e;
            let j = e + 1;
            let dx = self.node_positions[j] - self.node_positions[i];
            let g = elem.conductance(1.0, dx);
            k_global[i][i] += g;
            k_global[i][j] -= g;
            k_global[j][i] -= g;
            k_global[j][j] += g;
        }
        for conv in &self.convection_bcs {
            k_global[conv.node][conv.node] += conv.stiffness_contribution();
        }
        for tcr in &self.contact_resistances {
            let g = tcr.conductance();
            k_global[tcr.node_1][tcr.node_1] += g;
            k_global[tcr.node_1][tcr.node_2] -= g;
            k_global[tcr.node_2][tcr.node_1] -= g;
            k_global[tcr.node_2][tcr.node_2] += g;
        }
        k_global
    }
    /// Assemble the load vector including convection and radiation contributions.
    pub fn assemble_load_vector(&self, heat_sources: &[f64]) -> Vec<f64> {
        let mut rhs = heat_sources.to_vec();
        for conv in &self.convection_bcs {
            rhs[conv.node] += conv.load_contribution();
        }
        for rad in &self.radiation_bcs {
            rhs[rad.node] += rad.linearized_load(self.temperatures[rad.node]);
        }
        rhs
    }
    /// Assemble the lumped capacitance vector C.
    pub fn assemble_capacitance_vector(&self) -> Vec<f64> {
        let mut c_vec = vec![0.0_f64; self.n_nodes];
        for (e, elem) in self.elements.iter().enumerate() {
            let cap = elem.capacitance();
            c_vec[e] += cap * 0.5;
            c_vec[e + 1] += cap * 0.5;
        }
        c_vec
    }
    /// Impose a Dirichlet (fixed-temperature) boundary condition at `node`.
    pub fn set_temperature_bc(&mut self, node: usize, temp: f64) {
        assert!(node < self.n_nodes);
        self.bc.retain(|&(n, _)| n != node);
        self.bc.push((node, temp));
        self.temperatures[node] = temp;
    }
    /// Add a convection boundary condition.
    pub fn add_convection_bc(&mut self, node: usize, h: f64, area: f64, t_inf: f64) {
        assert!(node < self.n_nodes);
        self.convection_bcs
            .push(ConvectionBC::new(node, h, area, t_inf));
    }
    /// Add a radiation boundary condition.
    pub fn add_radiation_bc(&mut self, node: usize, emissivity: f64, area: f64, t_surr: f64) {
        assert!(node < self.n_nodes);
        self.radiation_bcs
            .push(RadiationBC::new(node, emissivity, area, t_surr));
    }
    /// Add a thermal contact resistance between two nodes.
    pub fn add_contact_resistance(&mut self, node_1: usize, node_2: usize, h_c: f64, area: f64) {
        assert!(node_1 < self.n_nodes);
        assert!(node_2 < self.n_nodes);
        self.contact_resistances
            .push(ThermalContactResistance::new(node_1, node_2, h_c, area));
    }
    /// Explicit Euler step: C*(T_new - T)/dt = -K*T + Q.
    pub fn step_explicit(&mut self, dt: f64, heat_sources: &[f64]) {
        assert_eq!(heat_sources.len(), self.n_nodes);
        let k = self.assemble_conductance_matrix();
        let c = self.assemble_capacitance_vector();
        let mut t_new = self.temperatures.clone();
        for i in 0..self.n_nodes {
            if self.bc.iter().any(|&(n, _)| n == i) {
                continue;
            }
            let kt_i: f64 = (0..self.n_nodes)
                .map(|j| k[i][j] * self.temperatures[j])
                .sum();
            let mut q_extra = 0.0;
            for conv in &self.convection_bcs {
                if conv.node == i {
                    q_extra += conv.load_contribution()
                        - conv.stiffness_contribution() * self.temperatures[i];
                }
            }
            for rad in &self.radiation_bcs {
                if rad.node == i {
                    q_extra -= rad.heat_flux(self.temperatures[i]);
                }
            }
            t_new[i] = self.temperatures[i] + dt / c[i] * (-kt_i + heat_sources[i] + q_extra);
        }
        self.temperatures = t_new;
    }
    /// Backward Euler (implicit) step: (C/dt + K)*T_new = (C/dt)*T + Q.
    ///
    /// Uses the Thomas algorithm on the tridiagonal system.
    pub fn step_implicit(&mut self, dt: f64, heat_sources: &[f64]) {
        assert_eq!(heat_sources.len(), self.n_nodes);
        let n = self.n_nodes;
        let k = self.assemble_conductance_matrix();
        let c = self.assemble_capacitance_vector();
        let load = self.assemble_load_vector(heat_sources);
        let mut a = vec![0.0_f64; n];
        let mut b = vec![0.0_f64; n];
        let mut cv = vec![0.0_f64; n];
        let mut rhs = vec![0.0_f64; n];
        for i in 0..n {
            b[i] = c[i] / dt + k[i][i];
            if i > 0 {
                a[i] = k[i][i - 1];
            }
            if i < n - 1 {
                cv[i] = k[i][i + 1];
            }
            rhs[i] = c[i] / dt * self.temperatures[i] + load[i];
        }
        for rad in &self.radiation_bcs {
            let i = rad.node;
            b[i] += rad.linearized_stiffness(self.temperatures[i]);
        }
        for &(node, temp) in &self.bc {
            a[node] = 0.0;
            b[node] = 1.0;
            cv[node] = 0.0;
            rhs[node] = temp;
        }
        self.temperatures = thomas_algorithm(&a, &b, &cv, &rhs);
    }
    /// Steady-state solve: K*T = Q.
    pub fn steady_state(&mut self, heat_sources: &[f64]) {
        assert_eq!(heat_sources.len(), self.n_nodes);
        let n = self.n_nodes;
        let k = self.assemble_conductance_matrix();
        let load = self.assemble_load_vector(heat_sources);
        let mut a = vec![0.0_f64; n];
        let mut b = vec![0.0_f64; n];
        let mut cv = vec![0.0_f64; n];
        let mut rhs = load;
        for i in 0..n {
            b[i] = k[i][i];
            if i > 0 {
                a[i] = k[i][i - 1];
            }
            if i < n - 1 {
                cv[i] = k[i][i + 1];
            }
        }
        for &(node, temp) in &self.bc {
            a[node] = 0.0;
            b[node] = 1.0;
            cv[node] = 0.0;
            rhs[node] = temp;
        }
        self.temperatures = thomas_algorithm(&a, &b, &cv, &rhs);
    }
    /// Heat flux at element `elem_idx`: q = -k * dT/dx \[W/m^2\].
    pub fn heat_flux_at(&self, elem_idx: usize) -> f64 {
        assert!(elem_idx < self.elements.len());
        let i = elem_idx;
        let j = elem_idx + 1;
        let dx = self.node_positions[j] - self.node_positions[i];
        let dt_dx = (self.temperatures[j] - self.temperatures[i]) / dx;
        -self.elements[elem_idx].conductivity * dt_dx
    }
    /// Total heat stored in the mesh: Q_total = sum(rho * cp * V * T).
    pub fn total_heat_content(&self) -> f64 {
        let c = self.assemble_capacitance_vector();
        c.iter()
            .zip(self.temperatures.iter())
            .map(|(ci, ti)| ci * ti)
            .sum()
    }
    /// Maximum temperature in the mesh.
    pub fn max_temperature(&self) -> f64 {
        self.temperatures
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Minimum temperature in the mesh.
    pub fn min_temperature(&self) -> f64 {
        self.temperatures
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min)
    }
    /// Compute the CFL stability limit for explicit time integration.
    ///
    /// dt_crit = dx^2 / (2 * alpha) for 1D.
    pub fn cfl_limit(&self) -> f64 {
        let mut dt_min = f64::INFINITY;
        for (e, elem) in self.elements.iter().enumerate() {
            let dx = self.node_positions[e + 1] - self.node_positions[e];
            let alpha = elem.diffusivity();
            let dt = dx * dx / (2.0 * alpha);
            if dt < dt_min {
                dt_min = dt;
            }
        }
        dt_min
    }
}
/// Crank-Nicolson time integration result.
pub struct CrankNicolsonResult {
    /// Temperature history: each inner Vec is temperatures at one time step.
    pub history: Vec<Vec<f64>>,
    /// Final temperatures.
    pub final_temperatures: Vec<f64>,
}
