// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Static elastic catenary of a mooring line resting on a frictionless seabed.
//!
//! The line has unstretched length `L`, submerged weight per unit length `w`
//! (N/m) and axial stiffness `EA` (N).  Its anchor sits on the seabed at the
//! origin and its fairlead at horizontal distance `X` and height `Z` above the
//! anchor.  The unknowns are the horizontal tension `H` (constant along the
//! line because the seabed is frictionless) and the vertical tension `V` at the
//! fairlead.  With `L_B = max(0, L - V/w)` the unstretched length lying on the
//! seabed and `V_A = max(0, V - w·L)` the vertical tension at the anchor
//! (non-zero only when the whole line is lifted off the seabed):
//!
//! ```text
//! X = L_B + (H/w)·[asinh(V/H) - asinh(V_A/H)] + H·L/EA
//! Z = (sqrt(H² + V²) - sqrt(H² + V_A²))/w + (V² - V_A²)/(2·w·EA)
//! ```
//!
//! These are the classical elastic catenary equations (Faltinsen 1990;
//! Jonkman 2007, NREL/TP-500-41958); for `EA → ∞` they reduce to `L_s = sqrt(Z·(Z + 2a))`,
//! `X = L - L_s + a·asinh(L_s/a)` and `V = w·L_s` with `a = H/w` while part of
//! the line rests on the seabed.
//!
//! The solver brackets the root in `H` and in `V` and refines both with
//! Brent's method, so it converges for every geometry that admits an
//! equilibrium.

use thiserror::Error;

/// Failure modes of [`solve_static_catenary`].
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum CatenaryError {
    /// A line parameter or span is not finite, or violates its sign constraint.
    #[error(
        "catenary parameter `{name}` must be finite and non-negative (positive for length and stiffness), got {value}"
    )]
    InvalidParameter {
        /// Name of the offending parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },
    /// An inextensible line is shorter than (or as long as) the distance it has
    /// to span, so the tension would be infinite.
    #[error(
        "inextensible line of length {length} m cannot span the {chord} m between anchor and fairlead"
    )]
    LineTooShort {
        /// Unstretched line length (m).
        length: f64,
        /// Straight-line distance between anchor and fairlead (m).
        chord: f64,
    },
    /// The root search failed to bracket or converge (extreme parameters).
    #[error("catenary equilibrium could not be bracketed or converged")]
    NoConvergence,
}

/// Static equilibrium of a catenary mooring line, see [`solve_static_catenary`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CatenarySolution {
    /// Horizontal tension `H` (N); uniform along a line on a frictionless seabed.
    pub horizontal_tension: f64,
    /// Vertical tension component `V` at the fairlead (N), equal to the weight
    /// of the suspended part when the line touches the seabed.
    pub vertical_tension: f64,
    /// Vertical tension component at the anchor (N); zero while part of the
    /// line rests on the seabed.
    pub anchor_vertical_tension: f64,
    /// Unstretched length of line resting on the seabed (m).
    pub grounded_length: f64,
    /// Unstretched line length `L` (m).
    pub length: f64,
    /// Submerged weight per unit length `w` (N/m).
    pub weight_per_length: f64,
    /// Axial stiffness `EA` (N); `f64::INFINITY` for an inextensible line.
    pub axial_stiffness: f64,
    /// Horizontal distance between anchor and fairlead (m).
    pub horizontal_span: f64,
    /// Fairlead height above the anchor (m).
    pub vertical_span: f64,
}

impl CatenarySolution {
    /// Magnitude of the tension at the fairlead, `sqrt(H² + V²)` (N).
    pub fn fairlead_tension(&self) -> f64 {
        self.horizontal_tension.hypot(self.vertical_tension)
    }

    /// Magnitude of the tension at the anchor, `sqrt(H² + V_A²)` (N).
    pub fn anchor_tension(&self) -> f64 {
        self.horizontal_tension.hypot(self.anchor_vertical_tension)
    }

    /// Unstretched length of the part of the line that is lifted off the seabed (m).
    pub fn suspended_length(&self) -> f64 {
        self.length - self.grounded_length
    }

    /// Catenary parameter `a = H/w` (m); `0` for a weightless line.
    pub fn catenary_parameter(&self) -> f64 {
        if self.weight_per_length > 0.0 {
            self.horizontal_tension / self.weight_per_length
        } else {
            0.0
        }
    }

    /// `1/EA`, zero for an inextensible line.
    fn inverse_stiffness(&self) -> f64 {
        if self.axial_stiffness.is_finite() && self.axial_stiffness > 0.0 {
            1.0 / self.axial_stiffness
        } else {
            0.0
        }
    }

    /// Position `[x, z]` (m) of the line point with unstretched arc length `s`
    /// measured from the anchor (`s` is clamped to `[0, L]`).
    ///
    /// The anchor is at `[0, 0]` and `s = L` reproduces the fairlead
    /// `[horizontal_span, vertical_span]`.  For a slack line (`H = 0`) the
    /// excess line piled up on the seabed is compressed into the horizontal
    /// span, and the rest rises vertically at the fairlead.
    pub fn position_at(&self, s: f64) -> [f64; 2] {
        let s = s.clamp(0.0, self.length);
        let h = self.horizontal_tension;
        let w = self.weight_per_length;
        let inv_ea = self.inverse_stiffness();

        if w <= 0.0 {
            // Weightless line: straight between anchor and fairlead.
            let t = if self.length > 0.0 {
                s / self.length
            } else {
                0.0
            };
            return [t * self.horizontal_span, t * self.vertical_span];
        }
        if h <= 0.0 {
            // Slack line: lying part compressed into the span, then a vertical riser.
            if s <= self.grounded_length {
                let t = if self.grounded_length > 0.0 {
                    s / self.grounded_length
                } else {
                    0.0
                };
                return [t * self.horizontal_span, 0.0];
            }
            let rise = s - self.grounded_length;
            let v = self.anchor_vertical_tension + w * rise;
            let v_a = self.anchor_vertical_tension;
            let z =
                (v * v - v_a * v_a) * (1.0 / (v + v_a).max(f64::MIN_POSITIVE) + 0.5 * inv_ea) / w;
            return [self.horizontal_span, z];
        }
        if s <= self.grounded_length {
            return [s * (1.0 + h * inv_ea), 0.0];
        }
        let suspended = s - self.grounded_length;
        let v_a = self.anchor_vertical_tension;
        let v = v_a + w * suspended;
        let t_a = h.hypot(v_a);
        let t_s = h.hypot(v);
        let x = self.grounded_length * (1.0 + h * inv_ea)
            + h / w * ((v / h).asinh() - (v_a / h).asinh())
            + h * suspended * inv_ea;
        let z = (v * v - v_a * v_a) * (1.0 / (t_s + t_a) + 0.5 * inv_ea) / w;
        [x, z]
    }
}

/// Line properties shared by the geometry functions of the solver.
#[derive(Clone, Copy)]
struct Line {
    length: f64,
    weight: f64,
    inv_ea: f64,
}

impl Line {
    /// Vertical tension at the anchor for fairlead vertical tension `v`.
    fn anchor_vertical(&self, v: f64) -> f64 {
        (v - self.weight * self.length).max(0.0)
    }

    /// Unstretched length resting on the seabed for fairlead vertical tension `v`.
    fn grounded_length(&self, v: f64) -> f64 {
        (self.length - v / self.weight).max(0.0)
    }

    /// Fairlead height above the anchor for the tensions `(h, v)`.
    fn height(&self, h: f64, v: f64) -> f64 {
        let v_a = self.anchor_vertical(v);
        let delta_sq = v * v - v_a * v_a;
        let t_sum = h.hypot(v) + h.hypot(v_a);
        let geometric = if t_sum > 0.0 { delta_sq / t_sum } else { 0.0 };
        (geometric + 0.5 * self.inv_ea * delta_sq) / self.weight
    }

    /// Horizontal anchor-to-fairlead distance for the tensions `(h, v)`.
    fn span(&self, h: f64, v: f64) -> f64 {
        let v_a = self.anchor_vertical(v);
        let catenary = if h > 0.0 {
            h / self.weight * ((v / h).asinh() - (v_a / h).asinh())
        } else {
            0.0
        };
        self.grounded_length(v) + catenary + h * self.length * self.inv_ea
    }

    /// Fairlead vertical tension that lifts the fairlead to height `z` for a
    /// given horizontal tension `h` (`height` is increasing in `v`).
    fn vertical_tension_for_height(&self, h: f64, z: f64) -> f64 {
        if z <= 0.0 {
            return 0.0;
        }
        let residual = |v: f64| self.height(h, v) - z;
        let mut lo = 0.0;
        let mut hi = self.weight * z;
        let mut f_lo = residual(lo);
        let mut f_hi = residual(hi);
        for _ in 0..MAX_BRACKET_STEPS {
            if f_hi >= 0.0 {
                break;
            }
            lo = hi;
            f_lo = f_hi;
            hi *= 2.0;
            f_hi = residual(hi);
        }
        if f_hi < 0.0 {
            // Unreachable height (inextensible line): best effort.
            return hi;
        }
        brent(residual, lo, hi, f_lo, f_hi)
    }
}

/// Maximum number of interval doublings when bracketing a root.
const MAX_BRACKET_STEPS: usize = 200;
/// Maximum number of Brent iterations.
const MAX_BRENT_ITERATIONS: usize = 200;

/// Brent's method on the bracket `[a, b]` with `f(a) = fa`, `f(b) = fb` of
/// opposite sign (or one of them zero).  Converges to a relative accuracy of
/// a few machine epsilons.
fn brent<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, fa: f64, fb: f64) -> f64 {
    let (mut a, mut b, mut fa, mut fb) = (a, b, fa, fb);
    if fa == 0.0 {
        return a;
    }
    if fb == 0.0 {
        return b;
    }
    let mut c = b;
    let mut fc = fb;
    let mut d = b - a;
    let mut e = d;
    for _ in 0..MAX_BRENT_ITERATIONS {
        if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol = 2.0 * f64::EPSILON * b.abs();
        let xm = 0.5 * (c - b);
        if xm.abs() <= tol || fb == 0.0 {
            return b;
        }
        if e.abs() >= tol && fa.abs() > fb.abs() {
            let s = fb / fa;
            let (mut p, mut q);
            if a == c {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                let q0 = fa / fc;
                let r = fb / fc;
                p = s * (2.0 * xm * q0 * (q0 - r) - (b - a) * (r - 1.0));
                q = (q0 - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if p > 0.0 {
                q = -q;
            }
            p = p.abs();
            let bound_1 = 3.0 * xm * q - (tol * q).abs();
            let bound_2 = (e * q).abs();
            if 2.0 * p < bound_1.min(bound_2) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        if d.abs() > tol {
            b += d;
        } else {
            b += tol.copysign(xm);
        }
        fb = f(b);
    }
    b
}

fn require_finite(name: &'static str, value: f64) -> Result<(), CatenaryError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(CatenaryError::InvalidParameter { name, value })
    }
}

fn require_non_negative(name: &'static str, value: f64) -> Result<(), CatenaryError> {
    require_finite(name, value)?;
    if value >= 0.0 {
        Ok(())
    } else {
        Err(CatenaryError::InvalidParameter { name, value })
    }
}

/// Solve the static equilibrium of an elastic catenary mooring line.
///
/// * `length` — unstretched line length `L > 0` (m).
/// * `weight_per_length` — submerged weight per unit length `w >= 0` (N/m).
///   (Multiply a mass density in kg/m by `g`.)  `w = 0` is the weightless
///   taut-string limit.
/// * `axial_stiffness` — `EA > 0` (N); pass `f64::INFINITY` for an
///   inextensible line.
/// * `horizontal_span` — horizontal anchor-to-fairlead distance `X >= 0` (m).
/// * `vertical_span` — fairlead height above the anchor `Z >= 0` (m); the
///   anchor rests on the (horizontal, frictionless) seabed.
///
/// When `X` is smaller than the span at which the line first needs any
/// horizontal tension, the line is slack: the result has `H = 0` and the
/// vertical tension is the weight of the line hanging from the fairlead.
///
/// # Errors
///
/// * [`CatenaryError::InvalidParameter`] for non-finite or negative inputs;
/// * [`CatenaryError::LineTooShort`] for an inextensible line with
///   `L <= sqrt(X² + Z²)`;
/// * [`CatenaryError::NoConvergence`] for parameters so extreme that the root
///   cannot be bracketed.
pub fn solve_static_catenary(
    length: f64,
    weight_per_length: f64,
    axial_stiffness: f64,
    horizontal_span: f64,
    vertical_span: f64,
) -> Result<CatenarySolution, CatenaryError> {
    require_finite("length", length)?;
    if length <= 0.0 {
        return Err(CatenaryError::InvalidParameter {
            name: "length",
            value: length,
        });
    }
    require_non_negative("weight_per_length", weight_per_length)?;
    if axial_stiffness.is_nan() || axial_stiffness <= 0.0 {
        return Err(CatenaryError::InvalidParameter {
            name: "axial_stiffness",
            value: axial_stiffness,
        });
    }
    require_non_negative("horizontal_span", horizontal_span)?;
    require_non_negative("vertical_span", vertical_span)?;

    let inv_ea = if axial_stiffness.is_infinite() {
        0.0
    } else {
        1.0 / axial_stiffness
    };
    let chord = horizontal_span.hypot(vertical_span);
    if inv_ea == 0.0 && chord >= length {
        return Err(CatenaryError::LineTooShort { length, chord });
    }

    let assemble = |h: f64, v: f64, line: &Line| CatenarySolution {
        horizontal_tension: h,
        vertical_tension: v,
        anchor_vertical_tension: line.anchor_vertical(v),
        grounded_length: line.grounded_length(v),
        length,
        weight_per_length,
        axial_stiffness,
        horizontal_span,
        vertical_span,
    };

    // Weightless limit: a straight elastic string along the chord.
    if weight_per_length == 0.0 {
        let tension = if inv_ea > 0.0 {
            (chord / length - 1.0).max(0.0) / inv_ea
        } else {
            0.0
        };
        let (h, v) = if chord > 0.0 {
            (
                tension * horizontal_span / chord,
                tension * vertical_span / chord,
            )
        } else {
            (0.0, 0.0)
        };
        return Ok(CatenarySolution {
            horizontal_tension: h,
            vertical_tension: v,
            anchor_vertical_tension: v,
            grounded_length: if tension > 0.0 { 0.0 } else { length },
            length,
            weight_per_length,
            axial_stiffness,
            horizontal_span,
            vertical_span,
        });
    }

    let line = Line {
        length,
        weight: weight_per_length,
        inv_ea,
    };

    // Slack line: even zero horizontal tension reaches the requested span.
    let v_slack = line.vertical_tension_for_height(0.0, vertical_span);
    let span_slack = line.span(0.0, v_slack);
    if horizontal_span <= span_slack {
        return Ok(assemble(0.0, v_slack, &line));
    }

    let residual = |h: f64| {
        let v = line.vertical_tension_for_height(h, vertical_span);
        line.span(h, v) - horizontal_span
    };
    let scale = weight_per_length * length;
    let mut h_lo = 1e-12 * scale;
    let mut f_lo = residual(h_lo);
    if f_lo >= 0.0 {
        // The root lies below 1e-12 of the line weight: effectively slack.
        return Ok(assemble(0.0, v_slack, &line));
    }
    let mut h_hi = scale;
    let mut f_hi = residual(h_hi);
    let mut steps = 0;
    while f_hi < 0.0 {
        h_lo = h_hi;
        f_lo = f_hi;
        h_hi *= 2.0;
        f_hi = residual(h_hi);
        steps += 1;
        if steps > MAX_BRACKET_STEPS || !f_hi.is_finite() {
            return Err(CatenaryError::NoConvergence);
        }
    }
    let h = brent(residual, h_lo, h_hi, f_lo, f_hi);
    let v = line.vertical_tension_for_height(h, vertical_span);
    let span_error = (line.span(h, v) - horizontal_span).abs();
    if !h.is_finite() || !v.is_finite() || span_error > 1e-6 * horizontal_span.max(1.0) {
        return Err(CatenaryError::NoConvergence);
    }
    Ok(assemble(h, v, &line))
}
