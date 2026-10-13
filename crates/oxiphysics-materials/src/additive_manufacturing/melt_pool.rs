// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Isotherm geometry of the Rosenthal moving point-source solution.
//!
//! The steady temperature field of a point source moving at speed `v` over the
//! surface of a semi-infinite solid (see [`RosenthalSolution::temperature`]) is
//!
//! ```text
//! T(ξ, y, z) = T₀ + q / (2π k R) · exp(-v (R + ξ) / (2 α)),    R = √(ξ² + y² + z²)
//! ```
//!
//! with `ξ = x − v t` the coordinate along the scan direction relative to the
//! source (negative behind the source, positive ahead of it), `q` the absorbed
//! power, `k` the thermal conductivity and `α` the thermal diffusivity.
//!
//! `T` depends on `y` and `z` only through `y² + z²`, so every isotherm is a
//! surface of revolution about the scan axis: its depth below the free surface
//! equals its half-width on the surface.
//!
//! For an isotherm `T = T_iso` put `L = q / (2π k (T_iso − T₀))`.  The isotherm
//! is the set of points with `ln(L / R) = v (R + ξ) / (2 α)`.  Behind the source
//! on the scan axis (`ξ < 0`, `y = z = 0`, `R = −ξ`) the exponential equals one,
//! so the isotherm crosses the axis at exactly `ξ = −L`.  Ahead of the source
//! the crossing `ξ = +F` is the root of `ln(L / ξ) = v ξ / α` and satisfies
//! `F <= L`.  Between the two crossings the half-width `y(ξ)` of the isotherm
//! rises to a single maximum and falls back to zero, which is located with a
//! golden-section search on `ξ` around a bisection for `y`.

use std::f64::consts::PI;

use super::types::{MeltPoolGeometry, RosenthalSolution};

/// Bisection steps of every root search.  The bracket shrinks by `2^-64`,
/// which is below the resolution of `f64` for the brackets used here.
const BISECTION_STEPS: usize = 64;

/// Golden-section steps of the width maximisation.  The bracket shrinks by a
/// factor of about `1e-13`; the maximum itself converges quadratically faster.
const GOLDEN_SECTION_STEPS: usize = 64;

/// `1/φ = (√5 − 1)/2`, the golden-section ratio.
const INVERSE_GOLDEN_RATIO: f64 = 0.618_033_988_749_894_9;

/// Validated, dimension-reduced parameters of one isotherm.
#[derive(Debug, Clone, Copy)]
struct IsothermParameters {
    /// Trailing crossing distance `L = q / (2π k (T_iso − T₀))` (m).
    rear: f64,
    /// Decay rate `v / (2 α)` of the leading flank (1/m); zero for a
    /// stationary source.
    decay: f64,
}

impl IsothermParameters {
    /// `ln(L / R) − decay·(R + ξ)`: positive inside the isotherm, zero on it
    /// and negative outside.  `r` is the distance from the source.
    fn excess(&self, xi: f64, r: f64) -> f64 {
        (self.rear / r).ln() - self.decay * (r + xi)
    }

    /// Distance `F` ahead of the source at which the isotherm crosses the scan
    /// axis: the root of `ln(L / ξ) = 2·decay·ξ` in `(0, L]`.
    fn front_extent(&self) -> f64 {
        let (mut lo, mut hi) = (0.0_f64, self.rear);
        for _ in 0..BISECTION_STEPS {
            let mid = 0.5 * (lo + hi);
            if self.excess(mid, mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Half-width `y` of the isotherm on the free surface at the position `xi`
    /// along the scan axis, or exactly zero when `xi` lies outside the pool.
    fn half_width_at(&self, xi: f64) -> f64 {
        let axis_distance = xi.abs();
        if axis_distance > 0.0 && self.excess(xi, axis_distance) <= 0.0 {
            return 0.0;
        }
        let (mut lo, mut hi) = (0.0_f64, self.rear);
        for _ in 0..BISECTION_STEPS {
            let mid = 0.5 * (lo + hi);
            if self.excess(xi, xi.hypot(mid)) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Maximum of [`Self::half_width_at`] over the pool extent
    /// `[−L, front]` by golden-section search (the width profile is unimodal).
    fn max_half_width(&self, front: f64) -> f64 {
        let (mut a, mut b) = (-self.rear, front);
        let mut c = b - INVERSE_GOLDEN_RATIO * (b - a);
        let mut d = a + INVERSE_GOLDEN_RATIO * (b - a);
        let mut f_c = self.half_width_at(c);
        let mut f_d = self.half_width_at(d);
        for _ in 0..GOLDEN_SECTION_STEPS {
            if f_c > f_d {
                b = d;
                d = c;
                f_d = f_c;
                c = b - INVERSE_GOLDEN_RATIO * (b - a);
                f_c = self.half_width_at(c);
            } else {
                a = c;
                c = d;
                f_c = f_d;
                d = a + INVERSE_GOLDEN_RATIO * (b - a);
                f_d = self.half_width_at(d);
            }
        }
        f_c.max(f_d)
    }
}

impl RosenthalSolution {
    /// Validate the solution parameters for the isotherm `t_iso` (K).
    ///
    /// Returns `None` unless the absorbed power, the conductivity and the
    /// temperature rise `t_iso − T₀` are positive and finite, the scan speed
    /// is finite and non-negative, and — for a moving source — the diffusivity
    /// is positive and finite.
    fn isotherm_parameters(&self, t_iso: f64) -> Option<IsothermParameters> {
        let absorbed_power = self.absorptivity * self.power;
        let rise = t_iso - self.ambient_temp;
        let conductivity = self.thermal_conductivity;
        let speed = self.scan_speed;
        let positive_finite = |x: f64| x.is_finite() && x > 0.0;
        if !(positive_finite(absorbed_power)
            && positive_finite(conductivity)
            && positive_finite(rise)
            && speed.is_finite()
            && speed >= 0.0)
        {
            return None;
        }
        let decay = if speed == 0.0 {
            0.0
        } else if positive_finite(self.diffusivity) {
            speed / (2.0 * self.diffusivity)
        } else {
            return None;
        };
        let rear = absorbed_power / (2.0 * PI * conductivity * rise);
        (positive_finite(rear) && decay.is_finite()).then_some(IsothermParameters { rear, decay })
    }

    /// Distance (m) behind the source at which the isotherm `t_iso` (K)
    /// crosses the scan axis on the trailing side,
    /// `L = q / (2π k (t_iso − T₀))`.
    ///
    /// On that half-axis the exponential factor of the Rosenthal field is one,
    /// so the result is exact: `temperature(-L, 0, 0) == t_iso`.  Returns
    /// `None` when `t_iso` does not exceed the ambient temperature or the
    /// solution parameters are not physical (non-positive power or
    /// conductivity, negative scan speed, non-positive diffusivity for a
    /// moving source).
    pub fn isotherm_rear_extent(&self, t_iso: f64) -> Option<f64> {
        self.isotherm_parameters(t_iso).map(|p| p.rear)
    }

    /// Distance (m) ahead of the source at which the isotherm `t_iso` (K)
    /// crosses the scan axis on the leading side: the root of
    /// `temperature(ξ, 0, 0) = t_iso` for `ξ > 0`, found by bisection.
    ///
    /// It never exceeds [`Self::isotherm_rear_extent`] and equals it for a
    /// stationary source.  Returns `None` under the conditions listed there.
    pub fn isotherm_front_extent(&self, t_iso: f64) -> Option<f64> {
        self.isotherm_parameters(t_iso).map(|p| p.front_extent())
    }

    /// Half-width `y` (m) of the isotherm `t_iso` (K) on the free surface at
    /// the position `xi` (m) along the scan axis (negative behind the source),
    /// i.e. the `y > 0` with `temperature(xi, y, 0) = t_iso`, found by
    /// bisection.
    ///
    /// The result is exactly zero for positions outside the pool.  Returns
    /// `None` under the conditions listed at [`Self::isotherm_rear_extent`].
    pub fn isotherm_half_width_at(&self, xi: f64, t_iso: f64) -> Option<f64> {
        if !xi.is_finite() {
            return None;
        }
        self.isotherm_parameters(t_iso).map(|p| p.half_width_at(xi))
    }

    /// Maximum half-width (m) of the isotherm `t_iso` (K) on the free
    /// surface: the maximum over the scan-axis position `ξ` of
    /// [`Self::isotherm_half_width_at`], located by golden-section search.
    ///
    /// For a fast source (`v·L/α ≫ 1`) it approaches the classical thick-plate
    /// value `√(2 q α / (π e k v (t_iso − T₀)))`.  Returns `None` under the
    /// conditions listed at [`Self::isotherm_rear_extent`].
    pub fn isotherm_max_half_width(&self, t_iso: f64) -> Option<f64> {
        let parameters = self.isotherm_parameters(t_iso)?;
        Some(parameters.max_half_width(parameters.front_extent()))
    }

    /// Melt-pool geometry for the solidus/melting temperature `t_melt` (K):
    /// the `t_melt` isotherm of the Rosenthal field.
    ///
    /// * `half_length` – half of the total pool length along the scan axis,
    ///   `(L_rear + L_front) / 2`, with the exact trailing extent
    ///   `L_rear = q / (2π k (t_melt − T₀))` and the leading extent
    ///   `L_front` from [`Self::isotherm_front_extent`];
    /// * `half_width` – the maximum surface half-width of the isotherm
    ///   ([`Self::isotherm_max_half_width`]);
    /// * `depth` – equal to `half_width`, because the field is axisymmetric
    ///   about the scan axis.
    ///
    /// The point-source solution is singular at the source, so the geometry is
    /// only meaningful for pools that are large compared with the beam.
    /// Returns `None` when `t_melt` does not exceed the ambient temperature or
    /// the parameters are not physical.
    pub fn melt_pool_geometry(&self, t_melt: f64) -> Option<MeltPoolGeometry> {
        let parameters = self.isotherm_parameters(t_melt)?;
        let front = parameters.front_extent();
        let half_width = parameters.max_half_width(front);
        Some(MeltPoolGeometry {
            half_length: 0.5 * (parameters.rear + front),
            half_width,
            depth: half_width,
        })
    }

    /// Temperature gradient magnitude (K/m) along the trailing centre line at
    /// the rear tip of the `t_melt` isotherm,
    /// `G = 2π k (t_melt − T₀)² / q`.
    ///
    /// On that half-axis `T = T₀ + q / (2π k |ξ|)`, whose derivative at the
    /// crossing `|ξ| = q / (2π k (t_melt − T₀))` is `G`.  The rear tip
    /// advances at the scan speed, so `v·G` is the cooling rate of the
    /// solidifying material.  Returns `None` under the conditions listed at
    /// [`Self::isotherm_rear_extent`].
    pub fn trailing_edge_gradient(&self, t_melt: f64) -> Option<f64> {
        let parameters = self.isotherm_parameters(t_melt)?;
        // With L = q / (2π k rise) the slope q / (2π k L²) equals
        // 2π k rise² / q.
        let gradient = self.absorptivity * self.power
            / (2.0 * PI * self.thermal_conductivity * parameters.rear * parameters.rear);
        gradient.is_finite().then_some(gradient)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::additive_manufacturing::{LayerModel, ScanStrategy, ThermalHistory};

    const POWER: f64 = 200.0;
    const ABSORPTIVITY: f64 = 0.35;
    const CONDUCTIVITY: f64 = 6.7;
    const DENSITY: f64 = 4430.0;
    const SPECIFIC_HEAT: f64 = 560.0;
    const MELTING_POINT: f64 = 1933.0;
    const AMBIENT: f64 = 300.0;

    fn diffusivity() -> f64 {
        CONDUCTIVITY / (DENSITY * SPECIFIC_HEAT)
    }

    fn solution(speed: f64) -> RosenthalSolution {
        RosenthalSolution::new(
            POWER,
            ABSORPTIVITY,
            speed,
            CONDUCTIVITY,
            diffusivity(),
            AMBIENT,
        )
    }

    fn layer(speed: f64, power: f64, spot_radius: f64) -> LayerModel {
        LayerModel::new(
            30.0e-6,
            120.0e-6,
            spot_radius,
            speed,
            power,
            ScanStrategy::Rotating67,
        )
    }

    fn history(speed: f64) -> ThermalHistory {
        ThermalHistory::from_process(
            &layer(speed, POWER, 55.0e-6),
            ABSORPTIVITY,
            CONDUCTIVITY,
            DENSITY,
            SPECIFIC_HEAT,
            MELTING_POINT,
            AMBIENT,
        )
    }

    fn relative_error(value: f64, reference: f64) -> f64 {
        ((value - reference) / reference).abs()
    }

    /// Independent evaluation of the maximum half-width.  The isotherm is
    /// parametrised by the distance `R` from the source: with `u = ln(L/R)` and
    /// `c = v/(2α)` it satisfies `R + ξ = u/c`, so the squared half-width is
    /// `R² − ξ² = (u/c)(2R − u/c)`.  That function of `R` is strictly concave on
    /// `[front, L]`; its stationary point `u − 1 + u/(c R) = 0` is found by
    /// bisection.
    fn parametric_max_half_width(rear: f64, front: f64, decay: f64) -> f64 {
        let derivative = |r: f64| {
            let u = (rear / r).ln();
            u - 1.0 + u / (decay * r)
        };
        let (mut lo, mut hi) = (front, rear);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if derivative(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let r = 0.5 * (lo + hi);
        let u = (rear / r).ln();
        ((u / decay) * (2.0 * r - u / decay)).sqrt()
    }

    #[test]
    fn rear_extent_is_the_exact_trailing_isotherm() {
        let rs = solution(1.0);
        let rear = rs.isotherm_rear_extent(MELTING_POINT).unwrap();
        let expected = ABSORPTIVITY * POWER / (2.0 * PI * CONDUCTIVITY * (MELTING_POINT - AMBIENT));
        assert!(relative_error(rear, expected) < 1.0e-14, "rear = {rear}");
        let t = rs.temperature(-rear, 0.0, 0.0);
        assert!(
            (t - MELTING_POINT).abs() < 1.0e-9,
            "T(-L, 0, 0) = {t}, expected {MELTING_POINT}"
        );
    }

    #[test]
    fn front_extent_is_the_leading_isotherm_crossing() {
        for speed in [0.05, 0.25, 1.0, 4.0] {
            let rs = solution(speed);
            let rear = rs.isotherm_rear_extent(MELTING_POINT).unwrap();
            let front = rs.isotherm_front_extent(MELTING_POINT).unwrap();
            assert!(front > 0.0 && front < rear, "v = {speed}: front = {front}");
            let t = rs.temperature(front, 0.0, 0.0);
            assert!(
                (t - MELTING_POINT).abs() < 1.0e-8,
                "v = {speed}: T(front, 0, 0) = {t}"
            );
        }
    }

    #[test]
    fn isotherm_points_lie_on_the_melting_temperature() {
        let rs = solution(1.0);
        let rear = rs.isotherm_rear_extent(MELTING_POINT).unwrap();
        let front = rs.isotherm_front_extent(MELTING_POINT).unwrap();
        for fraction in [-0.95, -0.6, -0.37, -0.1, 0.0, 0.5, 0.9] {
            let xi = if fraction < 0.0 {
                fraction * rear
            } else {
                fraction * front
            };
            let y = rs.isotherm_half_width_at(xi, MELTING_POINT).unwrap();
            assert!(y > 0.0, "xi = {xi}: y = {y}");
            let t = rs.temperature(xi, y, 0.0);
            assert!(
                (t - MELTING_POINT).abs() < 1.0e-8 * (MELTING_POINT - AMBIENT),
                "xi = {xi}: T = {t}"
            );
        }
    }

    #[test]
    fn half_width_is_zero_outside_the_pool() {
        let rs = solution(1.0);
        let rear = rs.isotherm_rear_extent(MELTING_POINT).unwrap();
        let front = rs.isotherm_front_extent(MELTING_POINT).unwrap();
        assert_eq!(
            rs.isotherm_half_width_at(-1.01 * rear, MELTING_POINT),
            Some(0.0)
        );
        assert_eq!(
            rs.isotherm_half_width_at(1.01 * front, MELTING_POINT),
            Some(0.0)
        );
        assert_eq!(rs.isotherm_half_width_at(f64::NAN, MELTING_POINT), None);
    }

    #[test]
    fn max_half_width_matches_a_brute_force_scan_of_the_temperature_field() {
        let rs = solution(1.0);
        let rear = rs.isotherm_rear_extent(MELTING_POINT).unwrap();
        let front = rs.isotherm_front_extent(MELTING_POINT).unwrap();
        // Bisection on y with the public temperature field only.
        let width_at = |xi: f64| {
            let (mut lo, mut hi) = (0.0_f64, rear);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if rs.temperature(xi, mid, 0.0) > MELTING_POINT {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            0.5 * (lo + hi)
        };
        let samples = 4000;
        let mut best = 0.0_f64;
        let mut best_xi = 0.0_f64;
        for i in 1..samples {
            let xi = -rear + (rear + front) * f64::from(i) / f64::from(samples);
            let w = width_at(xi);
            if w > best {
                best = w;
                best_xi = xi;
            }
        }
        let found = rs.isotherm_max_half_width(MELTING_POINT).unwrap();
        assert!(
            found >= best * (1.0 - 1.0e-12),
            "search {found} below scan {best}"
        );
        assert!(
            relative_error(found, best) < 1.0e-5,
            "search {found} vs scan {best} (at xi = {best_xi})"
        );
        // The widest section lies behind the source.
        assert!(best_xi < 0.0 && best_xi > -rear, "xi* = {best_xi}");
    }

    #[test]
    fn max_half_width_matches_the_parametric_solution() {
        for speed in [0.01, 0.1, 0.5, 1.0, 2.0, 8.0] {
            let rs = solution(speed);
            let rear = rs.isotherm_rear_extent(MELTING_POINT).unwrap();
            let front = rs.isotherm_front_extent(MELTING_POINT).unwrap();
            let decay = speed / (2.0 * diffusivity());
            let reference = parametric_max_half_width(rear, front, decay);
            let found = rs.isotherm_max_half_width(MELTING_POINT).unwrap();
            assert!(
                relative_error(found, reference) < 1.0e-9,
                "v = {speed}: {found} vs {reference}"
            );
        }
    }

    #[test]
    fn fast_source_approaches_the_classical_thick_plate_width() {
        let speed = 20.0;
        let rs = solution(speed);
        let q = ABSORPTIVITY * POWER;
        let classical = (2.0 * q * diffusivity()
            / (PI * std::f64::consts::E * CONDUCTIVITY * speed * (MELTING_POINT - AMBIENT)))
            .sqrt();
        let found = rs.isotherm_max_half_width(MELTING_POINT).unwrap();
        assert!(
            relative_error(found, classical) < 1.0e-3,
            "{found} vs classical {classical}"
        );
    }

    #[test]
    fn stationary_source_gives_a_hemispherical_pool() {
        let rs = solution(0.0);
        let rear = rs.isotherm_rear_extent(MELTING_POINT).unwrap();
        let pool = rs.melt_pool_geometry(MELTING_POINT).unwrap();
        assert!(relative_error(pool.half_length, rear) < 1.0e-12);
        assert!(relative_error(pool.half_width, rear) < 1.0e-9);
        assert!(relative_error(pool.depth, rear) < 1.0e-9);
        // The diffusivity is irrelevant without motion.
        let no_diffusivity =
            RosenthalSolution::new(POWER, ABSORPTIVITY, 0.0, CONDUCTIVITY, 0.0, AMBIENT);
        assert!(no_diffusivity.melt_pool_geometry(MELTING_POINT).is_some());
    }

    #[test]
    fn pool_shrinks_with_scan_speed() {
        let speeds = [0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0];
        let pools: Vec<MeltPoolGeometry> = speeds
            .iter()
            .map(|&v| solution(v).melt_pool_geometry(MELTING_POINT).unwrap())
            .collect();
        for pair in pools.windows(2) {
            assert!(pair[1].half_width < pair[0].half_width, "{pair:?}");
            assert!(pair[1].depth < pair[0].depth, "{pair:?}");
            assert!(pair[1].half_length < pair[0].half_length, "{pair:?}");
        }
        // Reference values of an independent evaluation (Ti-6Al-4V, 70 W
        // absorbed, 1 m/s).
        let pool = solution(1.0).melt_pool_geometry(MELTING_POINT).unwrap();
        assert!(relative_error(pool.half_length, 5.151257250847936e-4) < 1.0e-8);
        assert!(relative_error(pool.half_width, 6.338753027789476e-5) < 1.0e-8);
        assert_eq!(pool.depth, pool.half_width);
    }

    #[test]
    fn pool_depends_on_diffusivity_and_melting_point() {
        let base = solution(1.0).melt_pool_geometry(MELTING_POINT).unwrap();
        let diffusive = RosenthalSolution::new(
            POWER,
            ABSORPTIVITY,
            1.0,
            CONDUCTIVITY,
            4.0 * diffusivity(),
            AMBIENT,
        )
        .melt_pool_geometry(MELTING_POINT)
        .unwrap();
        // Faster heat spreading widens the pool and lengthens its leading flank.
        assert!(diffusive.half_width > 1.5 * base.half_width);
        assert!(diffusive.half_length > base.half_length);
        let hotter = solution(1.0)
            .melt_pool_geometry(MELTING_POINT + 400.0)
            .unwrap();
        assert!(hotter.half_width < base.half_width);
        assert!(hotter.half_length < base.half_length);
    }

    #[test]
    fn pool_grows_with_absorbed_power() {
        let low = solution(1.0).melt_pool_geometry(MELTING_POINT).unwrap();
        let high = RosenthalSolution::new(
            2.0 * POWER,
            ABSORPTIVITY,
            1.0,
            CONDUCTIVITY,
            diffusivity(),
            AMBIENT,
        )
        .melt_pool_geometry(MELTING_POINT)
        .unwrap();
        // Fast-source scaling: width ~ sqrt(q), trailing extent ~ q.
        let width_ratio = high.half_width / low.half_width;
        assert!(
            (width_ratio - 2.0_f64.sqrt()).abs() < 5.0e-3,
            "width ratio = {width_ratio}"
        );
        assert!(high.half_length > 1.9 * low.half_length);
    }

    #[test]
    fn depth_equals_half_width_by_axisymmetry() {
        let rs = solution(1.0);
        let pool = rs.melt_pool_geometry(MELTING_POINT).unwrap();
        assert_eq!(pool.depth, pool.half_width);
        let xi = -0.37 * rs.isotherm_rear_extent(MELTING_POINT).unwrap();
        let w = rs.isotherm_half_width_at(xi, MELTING_POINT).unwrap();
        let lateral = rs.temperature(xi, w, 0.0);
        let downward = rs.temperature(xi, 0.0, w);
        assert!((lateral - downward).abs() < 1.0e-9);
    }

    #[test]
    fn trailing_edge_gradient_is_the_temperature_slope_at_the_rear_tip() {
        let rs = solution(1.0);
        let gradient = rs.trailing_edge_gradient(MELTING_POINT).unwrap();
        let q = ABSORPTIVITY * POWER;
        let expected = 2.0 * PI * CONDUCTIVITY * (MELTING_POINT - AMBIENT).powi(2) / q;
        assert!(relative_error(gradient, expected) < 1.0e-14);
        let rear = rs.isotherm_rear_extent(MELTING_POINT).unwrap();
        let h = 1.0e-3 * rear;
        let slope =
            (rs.temperature(-rear + h, 0.0, 0.0) - rs.temperature(-rear - h, 0.0, 0.0)) / (2.0 * h);
        assert!(
            relative_error(slope, gradient) < 1.0e-5,
            "finite-difference slope {slope} vs {gradient}"
        );
    }

    #[test]
    fn unphysical_parameters_are_rejected() {
        let t = MELTING_POINT;
        let with = |power: f64, absorptivity: f64, speed: f64, k: f64, alpha: f64| {
            RosenthalSolution::new(power, absorptivity, speed, k, alpha, AMBIENT)
        };
        let alpha = diffusivity();
        assert!(
            with(0.0, ABSORPTIVITY, 1.0, CONDUCTIVITY, alpha)
                .melt_pool_geometry(t)
                .is_none()
        );
        assert!(
            with(POWER, 0.0, 1.0, CONDUCTIVITY, alpha)
                .melt_pool_geometry(t)
                .is_none()
        );
        assert!(
            with(POWER, ABSORPTIVITY, 1.0, 0.0, alpha)
                .melt_pool_geometry(t)
                .is_none()
        );
        assert!(
            with(POWER, ABSORPTIVITY, -1.0, CONDUCTIVITY, alpha)
                .melt_pool_geometry(t)
                .is_none()
        );
        assert!(
            with(POWER, ABSORPTIVITY, 1.0, CONDUCTIVITY, 0.0)
                .melt_pool_geometry(t)
                .is_none()
        );
        assert!(
            with(POWER, ABSORPTIVITY, f64::NAN, CONDUCTIVITY, alpha)
                .melt_pool_geometry(t)
                .is_none()
        );
        let rs = solution(1.0);
        assert!(rs.melt_pool_geometry(AMBIENT).is_none());
        assert!(rs.melt_pool_geometry(AMBIENT - 10.0).is_none());
        assert!(rs.melt_pool_geometry(f64::NAN).is_none());
        assert!(rs.trailing_edge_gradient(AMBIENT).is_none());
        assert!(rs.isotherm_rear_extent(AMBIENT).is_none());
        assert!(rs.isotherm_front_extent(f64::INFINITY).is_none());
        assert!(rs.isotherm_max_half_width(AMBIENT).is_none());
    }

    // ---- ThermalHistory::from_process ------------------------------------

    #[test]
    fn from_process_pool_is_not_a_fixed_fraction_of_the_spot_radius() {
        // Regression: the pool used to be (r, 0.6 r, 0.4 r) whatever the inputs.
        let th = history(1.0);
        let pool = th.melt_pool;
        let spot = 55.0e-6;
        assert!((pool.half_length - spot).abs() > 100.0e-6, "{pool:?}");
        assert!((pool.half_width / pool.half_length - 0.6).abs() > 0.1);
        assert!((pool.depth / pool.half_length - 0.4).abs() > 0.1);
        assert!(relative_error(pool.half_length, 5.151257250847936e-4) < 1.0e-8);
        assert!(relative_error(pool.half_width, 6.338753027789476e-5) < 1.0e-8);
        assert_eq!(pool.depth, pool.half_width);
    }

    #[test]
    fn from_process_higher_scan_speed_gives_a_smaller_pool() {
        let slow = history(0.5).melt_pool;
        let fast = history(2.0).melt_pool;
        assert!(fast.half_width < slow.half_width, "{fast:?} vs {slow:?}");
        assert!(fast.depth < slow.depth);
        assert!(fast.half_length < slow.half_length);
        assert!(relative_error(fast.half_width, 4.49013364461862e-5) < 1.0e-8);
        assert!(relative_error(slow.half_width, 8.933083306184714e-5) < 1.0e-8);
    }

    #[test]
    fn from_process_depends_on_diffusivity_and_melting_point() {
        let scan_layer = layer(1.0, POWER, 55.0e-6);
        let build = |density: f64, melting_point: f64| {
            ThermalHistory::from_process(
                &scan_layer,
                ABSORPTIVITY,
                CONDUCTIVITY,
                density,
                SPECIFIC_HEAT,
                melting_point,
                AMBIENT,
            )
        };
        let base = build(DENSITY, MELTING_POINT);
        // A lighter material has a larger diffusivity and a wider pool.
        assert!(
            build(0.5 * DENSITY, MELTING_POINT).melt_pool.half_width > base.melt_pool.half_width
        );
        // A higher melting point gives a smaller pool.
        assert!(build(DENSITY, 2300.0).melt_pool.half_width < base.melt_pool.half_width);
    }

    #[test]
    fn from_process_gradient_and_cooling_rate_follow_the_rear_tip() {
        let speed = 1.2;
        let th = history(speed);
        let q = ABSORPTIVITY * POWER;
        let gradient = 2.0 * PI * CONDUCTIVITY * (MELTING_POINT - AMBIENT).powi(2) / q;
        assert!(relative_error(th.thermal_gradient, gradient) < 1.0e-14);
        assert!(relative_error(th.cooling_rate, speed * gradient) < 1.0e-14);
        // Faster scanning cools faster at an unchanged gradient.
        assert!(history(2.4).cooling_rate > 1.9 * th.cooling_rate);
        // The solidification rate is the scan speed.
        assert!(relative_error(th.solidification_rate(), speed) < 1.0e-12);
        assert!(th.primary_dendrite_arm_spacing() > 0.0);
    }

    #[test]
    fn from_process_peak_temperature_is_the_trailing_axis_value_at_one_beam_radius() {
        let th = history(1.0);
        let rs = solution(1.0);
        let at_beam_edge = rs.temperature(-55.0e-6, 0.0, 0.0);
        assert!(relative_error(th.peak_temperature, at_beam_edge) < 1.0e-12);
        assert!(th.recrystallised);
    }

    #[test]
    fn from_process_haz_is_the_half_rise_isotherm_outside_the_pool() {
        let th = history(1.0);
        assert!(
            relative_error(th.haz_width, 2.6415142614477638e-5) < 1.0e-7,
            "haz = {}",
            th.haz_width
        );
        let rs = solution(1.0);
        let haz_isotherm = 0.5 * (MELTING_POINT + AMBIENT);
        let outer = rs.isotherm_max_half_width(haz_isotherm).unwrap();
        assert!(relative_error(th.melt_pool.half_width + th.haz_width, outer) < 1.0e-12);
    }

    #[test]
    fn from_process_without_melting_has_no_pool() {
        // 1 W absorbed stays far below the melting point at one beam radius.
        let weak = ThermalHistory::from_process(
            &layer(1.0, 1.0, 55.0e-6),
            ABSORPTIVITY,
            CONDUCTIVITY,
            DENSITY,
            SPECIFIC_HEAT,
            MELTING_POINT,
            AMBIENT,
        );
        assert!(weak.peak_temperature < MELTING_POINT);
        assert_eq!(weak.melt_pool.half_length, 0.0);
        assert_eq!(weak.melt_pool.half_width, 0.0);
        assert_eq!(weak.melt_pool.depth, 0.0);
        assert_eq!(weak.thermal_gradient, 0.0);
        assert_eq!(weak.cooling_rate, 0.0);
        assert_eq!(weak.haz_width, 0.0);
        assert_eq!(weak.solidification_rate(), 0.0);
    }

    #[test]
    fn from_process_pool_threshold_is_the_beam_radius() {
        // The pool exists exactly when its trailing extent exceeds the beam
        // radius, i.e. when the peak temperature exceeds the melting point.
        let spot = 55.0e-6;
        let threshold = 2.0 * PI * CONDUCTIVITY * spot * (MELTING_POINT - AMBIENT) / ABSORPTIVITY;
        let build = |power: f64| {
            ThermalHistory::from_process(
                &layer(1.0, power, spot),
                ABSORPTIVITY,
                CONDUCTIVITY,
                DENSITY,
                SPECIFIC_HEAT,
                MELTING_POINT,
                AMBIENT,
            )
        };
        assert_eq!(build(0.99 * threshold).melt_pool.half_width, 0.0);
        let just_above = build(1.01 * threshold).melt_pool;
        assert!(just_above.half_width > 0.0);
        assert!(just_above.half_length > 0.0);
    }

    #[test]
    fn from_process_stationary_beam_has_zero_cooling_rate() {
        let th = history(0.0);
        assert!(th.melt_pool.half_width > 0.0);
        assert!(th.thermal_gradient > 0.0);
        assert_eq!(th.cooling_rate, 0.0);
    }
}
