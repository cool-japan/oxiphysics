//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::Material;

/// Extended material with full elastic constants.
#[derive(Debug, Clone)]
pub struct ExtendedMaterial {
    /// Base material (density, friction, restitution).
    pub base: Material,
    /// Young's modulus (Pa).
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Yield strength (Pa), 0 if not applicable.
    pub yield_strength: f64,
    /// Ultimate tensile strength (Pa), 0 if not applicable.
    pub uts: f64,
    /// Thermal conductivity (W/(m*K)).
    pub thermal_conductivity: f64,
    /// Specific heat capacity (J/(kg*K)).
    pub specific_heat: f64,
    /// Coefficient of thermal expansion (1/K).
    pub thermal_expansion: f64,
}
impl ExtendedMaterial {
    /// Create an extended material.
    pub fn new(
        base: Material,
        young_modulus: f64,
        poisson_ratio: f64,
        yield_strength: f64,
        uts: f64,
        thermal_conductivity: f64,
        specific_heat: f64,
        thermal_expansion: f64,
    ) -> Self {
        Self {
            base,
            young_modulus,
            poisson_ratio,
            yield_strength,
            uts,
            thermal_conductivity,
            specific_heat,
            thermal_expansion,
        }
    }
    /// Bulk modulus K = E / (3*(1-2*nu)).
    pub fn bulk_modulus(&self) -> f64 {
        let nu = self.poisson_ratio;
        self.young_modulus / (3.0 * (1.0 - 2.0 * nu))
    }
    /// Shear modulus G = E / (2*(1+nu)).
    pub fn shear_modulus(&self) -> f64 {
        let nu = self.poisson_ratio;
        self.young_modulus / (2.0 * (1.0 + nu))
    }
    /// Specific stiffness E / rho (m^2/s^2).
    pub fn specific_stiffness(&self) -> f64 {
        self.young_modulus / self.base.density
    }
    /// Specific strength sigma_y / rho (m^2/s^2).
    pub fn specific_strength(&self) -> f64 {
        self.yield_strength / self.base.density
    }
    /// Thermal diffusivity k / (rho * cp) (m^2/s).
    pub fn thermal_diffusivity(&self) -> f64 {
        self.thermal_conductivity / (self.base.density * self.specific_heat)
    }
}
/// Material category for filtering.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MaterialCategory {
    /// Metals and alloys.
    Metal,
    /// Polymers and plastics.
    Polymer,
    /// Ceramics.
    Ceramic,
    /// Composite materials.
    Composite,
    /// Biomaterials.
    Biomaterial,
    /// Building and construction materials.
    Building,
    /// Fluids (liquids and gases).
    Fluid,
    /// Sports equipment materials.
    Sports,
}
/// A preset material with its category.
pub struct CategorisedMaterial {
    /// The material.
    pub material: Material,
    /// The category.
    pub category: MaterialCategory,
}
