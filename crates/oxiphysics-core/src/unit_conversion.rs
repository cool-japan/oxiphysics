// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit conversion and fundamental physical constants.
//!
//! Provides enumerations for common measurement units and conversion functions
//! between them, as well as a struct containing widely-used physical constants
//! in SI units.

// ──────────────────────────────────────────────────────────────────────────────
// PhysicalConstants
// ──────────────────────────────────────────────────────────────────────────────

/// Fundamental physical constants in SI units (CODATA 2018 values).
#[derive(Debug, Clone, PartialEq)]
pub struct PhysicalConstants {
    /// Speed of light in vacuum (m/s).
    pub c: f64,
    /// Planck constant (J·s).
    pub h: f64,
    /// Boltzmann constant (J/K).
    pub k_b: f64,
    /// Avogadro constant (mol⁻¹).
    pub n_a: f64,
    /// Molar gas constant (J/(mol·K)).
    pub r: f64,
    /// Newtonian gravitational constant (m³/(kg·s²)).
    pub g: f64,
    /// Electric constant / permittivity of free space (F/m).
    pub epsilon_0: f64,
    /// Magnetic constant / permeability of free space (H/m).
    pub mu_0: f64,
}

impl Default for PhysicalConstants {
    fn default() -> Self {
        Self {
            c: 2.997_924_58e8,
            h: 6.626_070_15e-34,
            k_b: 1.380_649e-23,
            n_a: 6.022_140_76e23,
            r: 8.314_462_618,
            g: 6.674_30e-11,
            epsilon_0: 8.854_187_812_8e-12,
            mu_0: 1.256_637_062_12e-6,
        }
    }
}

impl PhysicalConstants {
    /// Construct constants with the default (CODATA 2018) values.
    pub fn new() -> Self {
        Self::default()
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Length
// ──────────────────────────────────────────────────────────────────────────────

/// Units of length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LengthUnit {
    /// Meter (SI base unit).
    Meter,
    /// Centimeter (1e-2 m).
    Centimeter,
    /// Millimeter (1e-3 m).
    Millimeter,
    /// Micrometre / micron (1e-6 m).
    Micrometer,
    /// Nanometer (1e-9 m).
    Nanometer,
    /// Ångström (1e-10 m).
    Angstrom,
    /// Foot (0.3048 m).
    Foot,
    /// Inch (0.0254 m).
    Inch,
    /// Light-year (9.460_730_472_580_8e15 m).
    LightYear,
    /// Astronomical Unit (1.495_978_707e11 m).
    AstronomicalUnit,
    /// Parsec (3.085_677_581_49e16 m).
    Parsec,
}

impl LengthUnit {
    /// Return the conversion factor to metres.
    fn to_meters(self) -> f64 {
        match self {
            LengthUnit::Meter => 1.0,
            LengthUnit::Centimeter => 1e-2,
            LengthUnit::Millimeter => 1e-3,
            LengthUnit::Micrometer => 1e-6,
            LengthUnit::Nanometer => 1e-9,
            LengthUnit::Angstrom => 1e-10,
            LengthUnit::Foot => 0.3048,
            LengthUnit::Inch => 0.0254,
            LengthUnit::LightYear => 9.460_730_472_580_8e15,
            LengthUnit::AstronomicalUnit => 1.495_978_707e11,
            LengthUnit::Parsec => 3.085_677_581_49e16,
        }
    }
}

/// Convert a length value from `from` units to `to` units.
///
/// # Examples
/// ```no_run
/// use oxiphysics_core::unit_conversion::{LengthUnit, convert_length};
/// let inches = convert_length(1.0, LengthUnit::Foot, LengthUnit::Inch);
/// assert!((inches - 12.0).abs() < 1e-10);
/// ```
pub fn convert_length(value: f64, from: LengthUnit, to: LengthUnit) -> f64 {
    value * from.to_meters() / to.to_meters()
}

// ──────────────────────────────────────────────────────────────────────────────
// Mass
// ──────────────────────────────────────────────────────────────────────────────

/// Units of mass.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MassUnit {
    /// Kilogram (SI base unit).
    Kilogram,
    /// Gram (1e-3 kg).
    Gram,
    /// Milligram (1e-6 kg).
    Milligram,
    /// Microgram (1e-9 kg).
    Microgram,
    /// Pound-mass (0.453_592_37 kg).
    Pound,
    /// Ounce (0.028_349_523_125 kg).
    Ounce,
    /// Atomic mass unit / dalton (1.660_539_066_60e-27 kg).
    AtomicMassUnit,
    /// Dalton — identical to `AtomicMassUnit`.
    Dalton,
}

impl MassUnit {
    /// Return the conversion factor to kilograms.
    fn to_kg(self) -> f64 {
        match self {
            MassUnit::Kilogram => 1.0,
            MassUnit::Gram => 1e-3,
            MassUnit::Milligram => 1e-6,
            MassUnit::Microgram => 1e-9,
            MassUnit::Pound => 0.453_592_37,
            MassUnit::Ounce => 0.028_349_523_125,
            MassUnit::AtomicMassUnit | MassUnit::Dalton => 1.660_539_066_60e-27,
        }
    }
}

/// Convert a mass value from `from` units to `to` units.
///
/// # Examples
/// ```no_run
/// use oxiphysics_core::unit_conversion::{MassUnit, convert_mass};
/// let grams = convert_mass(1.0, MassUnit::Kilogram, MassUnit::Gram);
/// assert!((grams - 1000.0).abs() < 1e-10);
/// ```
pub fn convert_mass(value: f64, from: MassUnit, to: MassUnit) -> f64 {
    value * from.to_kg() / to.to_kg()
}

// ──────────────────────────────────────────────────────────────────────────────
// Energy
// ──────────────────────────────────────────────────────────────────────────────

/// Units of energy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyUnit {
    /// Joule (SI).
    Joule,
    /// Kilojoule (1e3 J).
    Kilojoule,
    /// Electronvolt (1.602_176_634e-19 J).
    Electronvolt,
    /// Kilocalorie (4184 J, thermochemical).
    Kilocalorie,
    /// Calorie (4.184 J, thermochemical).
    Calorie,
    /// Erg (1e-7 J, CGS).
    Erg,
    /// British thermal unit (1055.05585262 J).
    Btu,
    /// Hartree energy (4.359_744_650_98e-18 J).
    Hartree,
}

impl EnergyUnit {
    /// Return the conversion factor to joules.
    fn to_joules(self) -> f64 {
        match self {
            EnergyUnit::Joule => 1.0,
            EnergyUnit::Kilojoule => 1e3,
            EnergyUnit::Electronvolt => 1.602_176_634e-19,
            EnergyUnit::Kilocalorie => 4184.0,
            EnergyUnit::Calorie => 4.184,
            EnergyUnit::Erg => 1e-7,
            EnergyUnit::Btu => 1055.05585262,
            EnergyUnit::Hartree => 4.359_744_650_98e-18,
        }
    }
}

/// Convert an energy value from `from` units to `to` units.
///
/// # Examples
/// ```no_run
/// use oxiphysics_core::unit_conversion::{EnergyUnit, convert_energy};
/// let kj = convert_energy(1000.0, EnergyUnit::Joule, EnergyUnit::Kilojoule);
/// assert!((kj - 1.0).abs() < 1e-10);
/// ```
pub fn convert_energy(value: f64, from: EnergyUnit, to: EnergyUnit) -> f64 {
    value * from.to_joules() / to.to_joules()
}

// ──────────────────────────────────────────────────────────────────────────────
// Temperature
// ──────────────────────────────────────────────────────────────────────────────

/// Units of temperature.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TemperatureUnit {
    /// Kelvin (absolute thermodynamic temperature).
    Kelvin,
    /// Degrees Celsius.
    Celsius,
    /// Degrees Fahrenheit.
    Fahrenheit,
    /// Rankine (absolute Fahrenheit scale).
    Rankine,
}

/// Convert a temperature value from `from` units to `to` units.
///
/// Temperature conversion is not a simple multiplicative factor, so this
/// function applies the appropriate affine transformation.
///
/// # Examples
/// ```no_run
/// use oxiphysics_core::unit_conversion::{TemperatureUnit, convert_temperature};
/// let f = convert_temperature(100.0, TemperatureUnit::Celsius, TemperatureUnit::Fahrenheit);
/// assert!((f - 212.0).abs() < 1e-8);
/// ```
pub fn convert_temperature(value: f64, from: TemperatureUnit, to: TemperatureUnit) -> f64 {
    // First convert to Kelvin.
    let kelvin = match from {
        TemperatureUnit::Kelvin => value,
        TemperatureUnit::Celsius => value + 273.15,
        TemperatureUnit::Fahrenheit => (value - 32.0) * 5.0 / 9.0 + 273.15,
        TemperatureUnit::Rankine => value * 5.0 / 9.0,
    };
    // Then convert from Kelvin to target.
    match to {
        TemperatureUnit::Kelvin => kelvin,
        TemperatureUnit::Celsius => kelvin - 273.15,
        TemperatureUnit::Fahrenheit => (kelvin - 273.15) * 9.0 / 5.0 + 32.0,
        TemperatureUnit::Rankine => kelvin * 9.0 / 5.0,
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Pressure
// ──────────────────────────────────────────────────────────────────────────────

/// Units of pressure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PressureUnit {
    /// Pascal (SI, N/m²).
    Pascal,
    /// Kilopascal (1e3 Pa).
    Kilopascal,
    /// Megapascal (1e6 Pa).
    Megapascal,
    /// Gigapascal (1e9 Pa).
    Gigapascal,
    /// Bar (1e5 Pa).
    Bar,
    /// Standard atmosphere (101_325 Pa).
    Atmosphere,
    /// Pound-force per square inch (6894.757_293_168_36 Pa).
    Psi,
    /// Torr / mmHg (101_325/760 Pa).
    Torr,
}

impl PressureUnit {
    /// Return the conversion factor to pascals.
    fn to_pascals(self) -> f64 {
        match self {
            PressureUnit::Pascal => 1.0,
            PressureUnit::Kilopascal => 1e3,
            PressureUnit::Megapascal => 1e6,
            PressureUnit::Gigapascal => 1e9,
            PressureUnit::Bar => 1e5,
            PressureUnit::Atmosphere => 101_325.0,
            PressureUnit::Psi => 6_894.757_293_168_361,
            PressureUnit::Torr => 101_325.0 / 760.0,
        }
    }
}

/// Convert a pressure value from `from` units to `to` units.
///
/// # Examples
/// ```no_run
/// use oxiphysics_core::unit_conversion::{PressureUnit, convert_pressure};
/// let pa = convert_pressure(1.0, PressureUnit::Atmosphere, PressureUnit::Pascal);
/// assert!((pa - 101_325.0).abs() < 1e-4);
/// ```
pub fn convert_pressure(value: f64, from: PressureUnit, to: PressureUnit) -> f64 {
    value * from.to_pascals() / to.to_pascals()
}

// ──────────────────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── PhysicalConstants ────────────────────────────────────────────────────

    #[test]
    fn test_speed_of_light_positive() {
        let c = PhysicalConstants::new();
        assert!(c.c > 0.0);
    }

    #[test]
    fn test_planck_constant_order_of_magnitude() {
        let pc = PhysicalConstants::new();
        assert!(pc.h > 1e-35 && pc.h < 1e-33);
    }

    #[test]
    fn test_boltzmann_constant_order() {
        let pc = PhysicalConstants::new();
        assert!(pc.k_b > 1e-24 && pc.k_b < 1e-22);
    }

    #[test]
    fn test_avogadro_order() {
        let pc = PhysicalConstants::new();
        assert!(pc.n_a > 6.0e23 && pc.n_a < 6.1e23);
    }

    #[test]
    fn test_gas_constant_relation() {
        // R = N_A * k_B
        let pc = PhysicalConstants::new();
        let computed = pc.n_a * pc.k_b;
        assert!((computed - pc.r).abs() / pc.r < 1e-6);
    }

    #[test]
    fn test_gravitational_constant_order() {
        let pc = PhysicalConstants::new();
        assert!(pc.g > 6.6e-11 && pc.g < 6.7e-11);
    }

    // ── Length ───────────────────────────────────────────────────────────────

    #[test]
    fn test_length_identity() {
        assert!((convert_length(5.0, LengthUnit::Meter, LengthUnit::Meter) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_foot_to_inches() {
        let inches = convert_length(1.0, LengthUnit::Foot, LengthUnit::Inch);
        assert!((inches - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_meter_to_centimeter() {
        let cm = convert_length(1.0, LengthUnit::Meter, LengthUnit::Centimeter);
        assert!((cm - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_meter_to_millimeter() {
        let mm = convert_length(1.0, LengthUnit::Meter, LengthUnit::Millimeter);
        assert!((mm - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_nanometer_to_angstrom() {
        let ang = convert_length(1.0, LengthUnit::Nanometer, LengthUnit::Angstrom);
        assert!((ang - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_length_roundtrip() {
        let val = 42.0;
        let cm = convert_length(val, LengthUnit::Meter, LengthUnit::Centimeter);
        let back = convert_length(cm, LengthUnit::Centimeter, LengthUnit::Meter);
        assert!((back - val).abs() < 1e-10);
    }

    #[test]
    fn test_au_larger_than_ly() {
        // 1 AU < 1 ly, so converting 1 ly → AU yields a large number.
        let aus = convert_length(1.0, LengthUnit::LightYear, LengthUnit::AstronomicalUnit);
        assert!(aus > 60_000.0);
    }

    #[test]
    fn test_inch_to_meter() {
        let m = convert_length(1.0, LengthUnit::Inch, LengthUnit::Meter);
        assert!((m - 0.0254).abs() < 1e-10);
    }

    // ── Mass ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_mass_identity() {
        assert!((convert_mass(3.0, MassUnit::Kilogram, MassUnit::Kilogram) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_kg_to_gram() {
        let g = convert_mass(1.0, MassUnit::Kilogram, MassUnit::Gram);
        assert!((g - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_pound_to_kg() {
        let kg = convert_mass(1.0, MassUnit::Pound, MassUnit::Kilogram);
        assert!((kg - 0.453_592_37).abs() < 1e-8);
    }

    #[test]
    fn test_ounce_to_gram() {
        // 1 lb = 16 oz, 1 lb = 453.59237 g → 1 oz ≈ 28.3495 g
        let g = convert_mass(1.0, MassUnit::Ounce, MassUnit::Gram);
        assert!((g - 28.349_523_125).abs() < 1e-6);
    }

    #[test]
    fn test_amu_dalton_equal() {
        let a = convert_mass(1.0, MassUnit::AtomicMassUnit, MassUnit::Kilogram);
        let b = convert_mass(1.0, MassUnit::Dalton, MassUnit::Kilogram);
        assert!((a - b).abs() < 1e-40);
    }

    #[test]
    fn test_mass_roundtrip() {
        let val = 7.5;
        let lb = convert_mass(val, MassUnit::Kilogram, MassUnit::Pound);
        let back = convert_mass(lb, MassUnit::Pound, MassUnit::Kilogram);
        assert!((back - val).abs() < 1e-10);
    }

    // ── Energy ───────────────────────────────────────────────────────────────

    #[test]
    fn test_energy_identity() {
        assert!((convert_energy(2.0, EnergyUnit::Joule, EnergyUnit::Joule) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_kj_to_j() {
        let j = convert_energy(1.0, EnergyUnit::Kilojoule, EnergyUnit::Joule);
        assert!((j - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_calorie_to_joule() {
        let j = convert_energy(1.0, EnergyUnit::Calorie, EnergyUnit::Joule);
        assert!((j - 4.184).abs() < 1e-10);
    }

    #[test]
    fn test_kcal_to_cal() {
        let cal = convert_energy(1.0, EnergyUnit::Kilocalorie, EnergyUnit::Calorie);
        assert!((cal - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_erg_to_joule() {
        let j = convert_energy(1.0, EnergyUnit::Erg, EnergyUnit::Joule);
        assert!((j - 1e-7).abs() < 1e-18);
    }

    #[test]
    fn test_energy_roundtrip() {
        let val = 100.0;
        let ev = convert_energy(val, EnergyUnit::Joule, EnergyUnit::Electronvolt);
        let back = convert_energy(ev, EnergyUnit::Electronvolt, EnergyUnit::Joule);
        assert!((back - val).abs() / val < 1e-10);
    }

    // ── Temperature ──────────────────────────────────────────────────────────

    #[test]
    fn test_celsius_to_kelvin_freezing() {
        let k = convert_temperature(0.0, TemperatureUnit::Celsius, TemperatureUnit::Kelvin);
        assert!((k - 273.15).abs() < 1e-8);
    }

    #[test]
    fn test_celsius_to_fahrenheit_boiling() {
        let f = convert_temperature(100.0, TemperatureUnit::Celsius, TemperatureUnit::Fahrenheit);
        assert!((f - 212.0).abs() < 1e-8);
    }

    #[test]
    fn test_fahrenheit_freezing_to_celsius() {
        let c = convert_temperature(32.0, TemperatureUnit::Fahrenheit, TemperatureUnit::Celsius);
        assert!(c.abs() < 1e-8);
    }

    #[test]
    fn test_rankine_to_kelvin() {
        // 0 R = 0 K
        let k = convert_temperature(0.0, TemperatureUnit::Rankine, TemperatureUnit::Kelvin);
        assert!(k.abs() < 1e-10);
    }

    #[test]
    fn test_rankine_459_67_is_zero_fahrenheit() {
        // 459.67 R = 0 °F
        let f = convert_temperature(
            459.67,
            TemperatureUnit::Rankine,
            TemperatureUnit::Fahrenheit,
        );
        assert!(f.abs() < 1e-6);
    }

    #[test]
    fn test_temperature_roundtrip_k_c() {
        let val = 500.0_f64;
        let c = convert_temperature(val, TemperatureUnit::Kelvin, TemperatureUnit::Celsius);
        let back = convert_temperature(c, TemperatureUnit::Celsius, TemperatureUnit::Kelvin);
        assert!((back - val).abs() < 1e-8);
    }

    #[test]
    fn test_temperature_identity_kelvin() {
        let val = 300.0_f64;
        let k = convert_temperature(val, TemperatureUnit::Kelvin, TemperatureUnit::Kelvin);
        assert!((k - val).abs() < 1e-10);
    }

    // ── Pressure ─────────────────────────────────────────────────────────────

    #[test]
    fn test_pressure_identity() {
        assert!(
            (convert_pressure(5.0, PressureUnit::Pascal, PressureUnit::Pascal) - 5.0).abs() < 1e-12
        );
    }

    #[test]
    fn test_atm_to_pascal() {
        let pa = convert_pressure(1.0, PressureUnit::Atmosphere, PressureUnit::Pascal);
        assert!((pa - 101_325.0).abs() < 1e-4);
    }

    #[test]
    fn test_bar_to_kpa() {
        let kpa = convert_pressure(1.0, PressureUnit::Bar, PressureUnit::Kilopascal);
        assert!((kpa - 100.0).abs() < 1e-8);
    }

    #[test]
    fn test_torr_to_pa() {
        let pa = convert_pressure(760.0, PressureUnit::Torr, PressureUnit::Pascal);
        assert!((pa - 101_325.0).abs() < 1e-4);
    }

    #[test]
    fn test_psi_to_pa() {
        let pa = convert_pressure(1.0, PressureUnit::Psi, PressureUnit::Pascal);
        assert!((pa - 6_894.757_293_168_361).abs() < 1e-4);
    }

    #[test]
    fn test_pressure_roundtrip() {
        let val = 3.5;
        let atm = convert_pressure(val, PressureUnit::Megapascal, PressureUnit::Atmosphere);
        let back = convert_pressure(atm, PressureUnit::Atmosphere, PressureUnit::Megapascal);
        assert!((back - val).abs() / val < 1e-10);
    }

    #[test]
    fn test_gpa_to_mpa() {
        let mpa = convert_pressure(1.0, PressureUnit::Gigapascal, PressureUnit::Megapascal);
        assert!((mpa - 1000.0).abs() < 1e-8);
    }
}
