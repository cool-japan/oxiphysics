// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Scientific visualization: tensor fields, phase diagrams, uncertainty bands.
//!
//! Provides data structures and algorithms for visualizing:
//! - 3D tensor fields via ellipsoid/superquadric glyphs
//! - Phase diagrams (T-P or T-x) with coexistence lines
//! - Error bands (mean ± std, bootstrap confidence intervals)
//! - Bifurcation diagrams (parameter vs steady state)
//! - Fractal sets (Mandelbrot/Julia) with multiple coloring schemes
//! - Force-directed network layouts (Fruchterman-Reingold)
//! - Phase portraits (vector field + trajectories)
//! - Hodographs and Lissajous figures

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// TensorGlyph — 3D tensor visualization
// ---------------------------------------------------------------------------

/// Shape type for tensor glyphs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GlyphShape {
    /// Ellipsoid glyph — axes scaled by eigenvalues.
    Ellipsoid,
    /// Superquadric glyph with shape parameter gamma.
    Superquadric,
    /// Cuboid glyph aligned with eigenvectors.
    Cuboid,
    /// Cylinder glyph (aligned with major eigenvector).
    Cylinder,
}

/// A rendered tensor glyph at a point in space.
#[derive(Debug, Clone)]
pub struct TensorGlyph {
    /// Center position `[x, y, z]`.
    pub position: [f64; 3],
    /// Eigenvalues `[lambda1, lambda2, lambda3]` (principal stresses/strains).
    pub eigenvalues: [f64; 3],
    /// Eigenvectors stored row-major: `[v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z]`.
    pub eigenvectors: [f64; 9],
    /// Selected glyph shape.
    pub shape: GlyphShape,
    /// Superquadric shape parameter (gamma). Only used for [`GlyphShape::Superquadric`].
    pub superquadric_gamma: f64,
    /// Scale factor applied uniformly to eigenvalue axes.
    pub scale: f64,
    /// RGBA color `[r, g, b, a]` in `[0, 1]`.
    pub color: [f32; 4],
}

impl TensorGlyph {
    /// Construct a new tensor glyph at `position` from eigendecomposition data.
    pub fn new(
        position: [f64; 3],
        eigenvalues: [f64; 3],
        eigenvectors: [f64; 9],
        shape: GlyphShape,
        scale: f64,
    ) -> Self {
        Self {
            position,
            eigenvalues,
            eigenvectors,
            shape,
            superquadric_gamma: 2.0,
            scale,
            color: [0.7, 0.3, 0.1, 1.0],
        }
    }

    /// Return the three semi-axis lengths of the ellipsoid (absolute eigenvalues × scale).
    pub fn ellipsoid_axes(&self) -> [f64; 3] {
        [
            self.eigenvalues[0].abs() * self.scale,
            self.eigenvalues[1].abs() * self.scale,
            self.eigenvalues[2].abs() * self.scale,
        ]
    }

    /// Compute the fractional anisotropy (FA) from the eigenvalues.
    ///
    /// FA = 0 for isotropic, FA → 1 for highly anisotropic.
    pub fn fractional_anisotropy(&self) -> f64 {
        let [l1, l2, l3] = self.eigenvalues;
        let mean = (l1 + l2 + l3) / 3.0;
        let num = ((l1 - mean).powi(2) + (l2 - mean).powi(2) + (l3 - mean).powi(2)).sqrt();
        let den = (l1 * l1 + l2 * l2 + l3 * l3).sqrt();
        if den < 1e-14 {
            0.0
        } else {
            (3.0 / 2.0_f64).sqrt() * num / den
        }
    }

    /// Compute the linear anisotropy metric cl = (λ1 − λ2) / λ1.
    pub fn linear_anisotropy(&self) -> f64 {
        let l1 = self.eigenvalues[0].abs().max(1e-14);
        let l2 = self.eigenvalues[1].abs();
        (l1 - l2) / l1
    }

    /// Compute the planar anisotropy metric cp = (λ2 − λ3) / λ1.
    pub fn planar_anisotropy(&self) -> f64 {
        let l1 = self.eigenvalues[0].abs().max(1e-14);
        let l2 = self.eigenvalues[1].abs();
        let l3 = self.eigenvalues[2].abs();
        (l2 - l3) / l1
    }

    /// Generate a simple point cloud approximating the ellipsoid surface (lon × lat grid).
    ///
    /// Returns `(n_lon * n_lat)` points as `[x, y, z]` arrays.
    pub fn ellipsoid_surface_points(&self, n_lon: usize, n_lat: usize) -> Vec<[f64; 3]> {
        let [a, b, c] = self.ellipsoid_axes();
        let mut pts = Vec::with_capacity(n_lon * n_lat);
        for i in 0..n_lat {
            let phi = PI * (i as f64 + 0.5) / n_lat as f64; // [0, π]
            for j in 0..n_lon {
                let theta = 2.0 * PI * j as f64 / n_lon as f64;
                let x = self.position[0] + a * phi.sin() * theta.cos();
                let y = self.position[1] + b * phi.sin() * theta.sin();
                let z = self.position[2] + c * phi.cos();
                pts.push([x, y, z]);
            }
        }
        pts
    }
}

/// Build tensor glyphs for a field of symmetric 3×3 matrices (Voigt notation).
///
/// `tensors`: each entry is `[Txx, Tyy, Tzz, Txy, Txz, Tyz]`.
/// `positions`: matching list of `[x, y, z]` positions.
/// Returns one `TensorGlyph` per entry after 3×3 eigen-decomposition (power iteration).
pub fn build_tensor_glyphs(
    tensors: &[[f64; 6]],
    positions: &[[f64; 3]],
    shape: GlyphShape,
    scale: f64,
) -> Vec<TensorGlyph> {
    assert_eq!(tensors.len(), positions.len());
    tensors
        .iter()
        .zip(positions.iter())
        .map(|(t, p)| {
            let (evals, evecs) = eigen3_sym(*t);
            TensorGlyph::new(*p, evals, evecs, shape, scale)
        })
        .collect()
}

/// Approximate eigendecomposition for a 3×3 symmetric matrix in Voigt form.
///
/// Uses Jacobi iteration (10 sweeps) for correctness.
fn eigen3_sym(voigt: [f64; 6]) -> ([f64; 3], [f64; 9]) {
    let [txx, tyy, tzz, txy, txz, tyz] = voigt;
    // Build 3×3 matrix
    let mut a = [[txx, txy, txz], [txy, tyy, tyz], [txz, tyz, tzz]];
    // Identity eigenvector matrix
    let mut v = [[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    // Jacobi iterations
    for _ in 0..50 {
        let mut max_off = 0.0_f64;
        let mut p = 0;
        let mut q = 1;
        for (i, a_row) in a.iter().enumerate() {
            for (j, a_ij) in a_row.iter().enumerate().skip(i + 1) {
                if a_ij.abs() > max_off {
                    max_off = a_ij.abs();
                    p = i;
                    q = j;
                }
            }
        }
        if max_off < 1e-14 {
            break;
        }
        let theta = 0.5 * (a[q][q] - a[p][p]) / a[p][q];
        let t_val = theta.signum() / (theta.abs() + (1.0 + theta * theta).sqrt());
        let c = 1.0 / (1.0 + t_val * t_val).sqrt();
        let s = t_val * c;

        // Update matrix
        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];
        a[p][p] = c * c * app - 2.0 * c * s * apq + s * s * aqq;
        a[q][q] = s * s * app + 2.0 * c * s * apq + c * c * aqq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;
        let sym_updates: Vec<(usize, f64, f64)> = (0..3)
            .filter(|&r| r != p && r != q)
            .map(|r| (r, c * a[p][r] - s * a[q][r], s * a[p][r] + c * a[q][r]))
            .collect();
        for (r, new_pr, new_qr) in sym_updates {
            a[p][r] = new_pr;
            a[r][p] = new_pr;
            a[q][r] = new_qr;
            a[r][q] = new_qr;
        }
        // Update eigenvectors
        for v_row in v.iter_mut() {
            let vp = v_row[p];
            let vq = v_row[q];
            v_row[p] = c * vp - s * vq;
            v_row[q] = s * vp + c * vq;
        }
    }

    let evals = [a[0][0], a[1][1], a[2][2]];
    let evecs = [
        v[0][0], v[0][1], v[0][2], v[1][0], v[1][1], v[1][2], v[2][0], v[2][1], v[2][2],
    ];
    (evals, evecs)
}

// ---------------------------------------------------------------------------
// PhaseDiagramViz — 2D phase diagram (T-P or T-x)
// ---------------------------------------------------------------------------

/// A point on a phase boundary (temperature, pressure).
#[derive(Debug, Clone, Copy)]
pub struct PhasePoint {
    /// Temperature in Kelvin.
    pub temperature: f64,
    /// Pressure in Pascals.
    pub pressure: f64,
}

impl PhasePoint {
    /// Construct a new phase point.
    pub fn new(temperature: f64, pressure: f64) -> Self {
        Self {
            temperature,
            pressure,
        }
    }
}

/// Named phase region in a phase diagram.
#[derive(Debug, Clone)]
pub struct PhaseRegion {
    /// Human-readable phase name (e.g., "Liquid", "Vapor", "Solid").
    pub name: String,
    /// RGBA color for the region.
    pub color: [f32; 4],
    /// Polygon vertices `[T, P]` outlining the region.
    pub boundary: Vec<[f64; 2]>,
}

/// A tie-line connecting two coexisting phases at fixed temperature.
#[derive(Debug, Clone, Copy)]
pub struct TieLine {
    /// Temperature of the tie-line (K).
    pub temperature: f64,
    /// Composition of phase 1 (mole fraction).
    pub x1: f64,
    /// Composition of phase 2 (mole fraction).
    pub x2: f64,
}

/// 2D phase diagram visualization data.
#[derive(Debug, Clone, Default)]
pub struct PhaseDiagramViz {
    /// Phase regions (filled polygons).
    pub regions: Vec<PhaseRegion>,
    /// Coexistence lines (list of connected `[T, P]` points).
    pub coexistence_lines: Vec<Vec<[f64; 2]>>,
    /// Triple point (if applicable).
    pub triple_point: Option<PhasePoint>,
    /// Critical point (if applicable).
    pub critical_point: Option<PhasePoint>,
    /// Tie-lines for two-phase regions.
    pub tie_lines: Vec<TieLine>,
}

impl PhaseDiagramViz {
    /// Create an empty phase diagram.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a coexistence line (e.g., Clausius-Clapeyron vapor pressure curve).
    pub fn add_coexistence_line(&mut self, points: Vec<[f64; 2]>) {
        self.coexistence_lines.push(points);
    }

    /// Add a phase region.
    pub fn add_region(&mut self, region: PhaseRegion) {
        self.regions.push(region);
    }

    /// Set the triple point.
    pub fn set_triple_point(&mut self, t: f64, p: f64) {
        self.triple_point = Some(PhasePoint::new(t, p));
    }

    /// Set the critical point.
    pub fn set_critical_point(&mut self, t: f64, p: f64) {
        self.critical_point = Some(PhasePoint::new(t, p));
    }

    /// Construct tie-lines at the given temperatures for a binary T-x diagram.
    ///
    /// `liquidus(T)` returns the liquid-phase composition; `solidus(T)` the solid.
    pub fn compute_tie_lines(
        &mut self,
        temperatures: &[f64],
        liquidus: impl Fn(f64) -> f64,
        solidus: impl Fn(f64) -> f64,
    ) {
        self.tie_lines = temperatures
            .iter()
            .map(|&t| TieLine {
                temperature: t,
                x1: solidus(t),
                x2: liquidus(t),
            })
            .collect();
    }

    /// Generate Clausius-Clapeyron vapor pressure curve points.
    ///
    /// `ln(P/P0) = -ΔH_vap/R * (1/T - 1/T0)`
    pub fn clausius_clapeyron_line(
        t_range: [f64; 2],
        n_pts: usize,
        delta_h_vap: f64,
        t0: f64,
        p0: f64,
    ) -> Vec<[f64; 2]> {
        let r = 8.314_f64;
        (0..n_pts)
            .map(|i| {
                let t = t_range[0] + (t_range[1] - t_range[0]) * i as f64 / (n_pts - 1) as f64;
                let ln_p = -delta_h_vap / r * (1.0 / t - 1.0 / t0);
                let p = p0 * ln_p.exp();
                [t, p]
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// ErrorBandPlot — mean ± std or confidence intervals
// ---------------------------------------------------------------------------

/// A single error band sample series.
#[derive(Debug, Clone)]
pub struct ErrorBandPlot {
    /// X-axis values.
    pub x: Vec<f64>,
    /// Mean values.
    pub mean: Vec<f64>,
    /// Lower bound of the band (mean − k·std or lower CI).
    pub lower: Vec<f64>,
    /// Upper bound of the band (mean + k·std or upper CI).
    pub upper: Vec<f64>,
    /// Label for this series.
    pub label: String,
    /// RGBA fill color for the uncertainty band.
    pub band_color: [f32; 4],
    /// RGBA line color for the mean curve.
    pub line_color: [f32; 4],
}

impl ErrorBandPlot {
    /// Construct from raw data samples (`rows × columns` where columns = x-positions).
    ///
    /// Each row is one bootstrap/Monte-Carlo sample. Computes mean and
    /// `k`-sigma bands across samples.
    pub fn from_samples(x: Vec<f64>, samples: &[Vec<f64>], k_sigma: f64) -> Self {
        let n = x.len();
        let m = samples.len();
        let mut mean = vec![0.0; n];
        let mut std = vec![0.0; n];

        if m == 0 {
            return Self {
                x,
                mean,
                lower: std.clone(),
                upper: std,
                label: String::new(),
                band_color: [0.3, 0.5, 0.9, 0.3],
                line_color: [0.1, 0.3, 0.8, 1.0],
            };
        }

        for j in 0..n {
            let sum: f64 = samples.iter().map(|s| s[j]).sum();
            mean[j] = sum / m as f64;
        }
        for j in 0..n {
            let var: f64 = samples
                .iter()
                .map(|s| (s[j] - mean[j]).powi(2))
                .sum::<f64>()
                / m as f64;
            std[j] = var.sqrt();
        }

        let lower = mean
            .iter()
            .zip(std.iter())
            .map(|(m, s)| m - k_sigma * s)
            .collect();
        let upper = mean
            .iter()
            .zip(std.iter())
            .map(|(m, s)| m + k_sigma * s)
            .collect();

        Self {
            x,
            mean,
            lower,
            upper,
            label: String::new(),
            band_color: [0.3, 0.5, 0.9, 0.3],
            line_color: [0.1, 0.3, 0.8, 1.0],
        }
    }

    /// Compute bootstrap confidence interval at `confidence_level` (e.g., 0.95).
    ///
    /// Performs `n_bootstrap` resamples and returns the percentile-based CI.
    pub fn bootstrap_ci(
        x: Vec<f64>,
        data: &[Vec<f64>],
        confidence_level: f64,
        n_bootstrap: usize,
    ) -> Self {
        let n = x.len();
        let m = data.len();
        if m == 0 || n == 0 {
            return Self {
                x,
                mean: vec![0.0; n],
                lower: vec![0.0; n],
                upper: vec![0.0; n],
                label: String::new(),
                band_color: [0.3, 0.5, 0.9, 0.3],
                line_color: [0.1, 0.3, 0.8, 1.0],
            };
        }

        let alpha = 1.0 - confidence_level;
        let lo_idx = ((alpha / 2.0) * n_bootstrap as f64) as usize;
        let hi_idx =
            ((1.0 - alpha / 2.0) * n_bootstrap as f64).min((n_bootstrap - 1) as f64) as usize;

        let mean: Vec<f64> = (0..n)
            .map(|j| data.iter().map(|s| s[j]).sum::<f64>() / m as f64)
            .collect();

        // Simple deterministic bootstrap via circular resampling
        let mut lower = vec![0.0; n];
        let mut upper = vec![0.0; n];
        for j in 0..n {
            let mut boot_means: Vec<f64> = (0..n_bootstrap)
                .map(|b| {
                    let sum: f64 = (0..m).map(|i| data[(i + b) % m][j]).sum();
                    sum / m as f64
                })
                .collect();
            boot_means.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            lower[j] = boot_means[lo_idx.min(boot_means.len() - 1)];
            upper[j] = boot_means[hi_idx.min(boot_means.len() - 1)];
        }

        Self {
            x,
            mean,
            lower,
            upper,
            label: String::new(),
            band_color: [0.3, 0.5, 0.9, 0.3],
            line_color: [0.1, 0.3, 0.8, 1.0],
        }
    }

    /// Check that the mean value is within the band at every x position.
    pub fn mean_within_band(&self) -> bool {
        self.mean
            .iter()
            .zip(self.lower.iter().zip(self.upper.iter()))
            .all(|(m, (lo, hi))| m >= lo && m <= hi)
    }
}

// ---------------------------------------------------------------------------
// BifurcationDiagram — parameter vs steady-state plot
// ---------------------------------------------------------------------------

/// Stability class for a bifurcation branch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BranchStability {
    /// Stable fixed point or limit cycle.
    Stable,
    /// Unstable fixed point.
    Unstable,
    /// Chaotic attractor region.
    Chaotic,
}

/// A point on a bifurcation diagram.
#[derive(Debug, Clone, Copy)]
pub struct BifurcationPoint {
    /// Control parameter value.
    pub parameter: f64,
    /// Steady-state value (e.g., x after transients).
    pub x: f64,
    /// Stability classification.
    pub stability: BranchStability,
}

/// Bifurcation diagram data.
#[derive(Debug, Clone, Default)]
pub struct BifurcationDiagram {
    /// All diagram points.
    pub points: Vec<BifurcationPoint>,
    /// Label for the control parameter axis.
    pub parameter_label: String,
    /// Label for the state variable axis.
    pub state_label: String,
}

impl BifurcationDiagram {
    /// Create an empty bifurcation diagram.
    pub fn new() -> Self {
        Self::default()
    }

    /// Compute a bifurcation diagram for the logistic map `x_{n+1} = r * x(1-x)`.
    ///
    /// - `r_range`: range of the control parameter r.
    /// - `n_r`: number of r values to sample.
    /// - `n_transient`: number of iterations to discard (transients).
    /// - `n_plot`: number of subsequent iterations to record.
    pub fn logistic_map(r_range: [f64; 2], n_r: usize, n_transient: usize, n_plot: usize) -> Self {
        let mut diag = Self::new();
        diag.parameter_label = "r".to_string();
        diag.state_label = "x".to_string();

        for i in 0..n_r {
            let r = r_range[0] + (r_range[1] - r_range[0]) * i as f64 / (n_r - 1) as f64;
            let mut x = 0.5;
            // Discard transients
            for _ in 0..n_transient {
                x = r * x * (1.0 - x);
            }
            // Collect attractor
            let mut attractor: Vec<f64> = Vec::with_capacity(n_plot);
            for _ in 0..n_plot {
                x = r * x * (1.0 - x);
                attractor.push(x);
            }
            // Estimate stability via number of distinct values
            attractor.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            attractor.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
            let stability = if attractor.len() == 1 {
                BranchStability::Stable
            } else if attractor.len() > n_plot / 4 {
                BranchStability::Chaotic
            } else {
                BranchStability::Stable
            };

            for xv in attractor {
                diag.points.push(BifurcationPoint {
                    parameter: r,
                    x: xv,
                    stability,
                });
            }
        }

        diag
    }

    /// Return all stable points.
    pub fn stable_points(&self) -> Vec<&BifurcationPoint> {
        self.points
            .iter()
            .filter(|p| p.stability == BranchStability::Stable)
            .collect()
    }

    /// Return all chaotic points.
    pub fn chaotic_points(&self) -> Vec<&BifurcationPoint> {
        self.points
            .iter()
            .filter(|p| p.stability == BranchStability::Chaotic)
            .collect()
    }
}

// ---------------------------------------------------------------------------
// FractalViz — Mandelbrot/Julia set renderer
// ---------------------------------------------------------------------------

/// Coloring algorithm for fractal escape-time rendering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FractalColoring {
    /// Raw iteration count mapped linearly to color.
    EscapeTime,
    /// Smooth (continuous) coloring using fractional escape count.
    Smooth,
    /// Orbit-trap coloring based on minimum distance to a trap shape.
    OrbitTrap,
}

/// Mandelbrot/Julia set renderer.
#[derive(Debug, Clone)]
pub struct FractalViz {
    /// Width in pixels.
    pub width: usize,
    /// Height in pixels.
    pub height: usize,
    /// Maximum number of iterations.
    pub max_iter: u32,
    /// Escape radius (typically 2.0).
    pub escape_radius: f64,
    /// Coloring scheme.
    pub coloring: FractalColoring,
    /// Julia set constant `c = (c_re, c_im)`. If `None`, renders Mandelbrot.
    pub julia_c: Option<(f64, f64)>,
}

impl FractalViz {
    /// Create a default Mandelbrot renderer.
    pub fn mandelbrot(width: usize, height: usize, max_iter: u32) -> Self {
        Self {
            width,
            height,
            max_iter,
            escape_radius: 2.0,
            coloring: FractalColoring::Smooth,
            julia_c: None,
        }
    }

    /// Create a Julia set renderer with constant `c`.
    pub fn julia(width: usize, height: usize, max_iter: u32, c_re: f64, c_im: f64) -> Self {
        Self {
            width,
            height,
            max_iter,
            escape_radius: 2.0,
            coloring: FractalColoring::Smooth,
            julia_c: Some((c_re, c_im)),
        }
    }

    /// Compute the escape time (or smooth value) for a single complex point.
    ///
    /// Returns a value in `[0, max_iter]` where `max_iter` means "did not escape"
    /// (interior point). For smooth coloring returns a fractional value.
    pub fn escape_value(&self, c_re: f64, c_im: f64) -> f64 {
        let er2 = self.escape_radius * self.escape_radius;
        let (mut zr, mut zi, cr, ci) = if let Some((jcr, jci)) = self.julia_c {
            (c_re, c_im, jcr, jci)
        } else {
            (0.0, 0.0, c_re, c_im)
        };

        for i in 0..self.max_iter {
            let zr2 = zr * zr;
            let zi2 = zi * zi;
            if zr2 + zi2 > er2 {
                if self.coloring == FractalColoring::Smooth {
                    // Smooth coloring: nu = n - log2(log2(|z|))
                    let log_zn = (zr2 + zi2).ln() / 2.0;
                    let nu = (log_zn / 2_f64.ln()).ln() / 2_f64.ln();
                    return (i as f64 + 1.0 - nu).max(0.0);
                }
                return i as f64;
            }
            let new_zr = zr2 - zi2 + cr;
            zi = 2.0 * zr * zi + ci;
            zr = new_zr;
        }
        self.max_iter as f64
    }

    /// Compute the orbit-trap value (minimum |z| over the orbit).
    pub fn orbit_trap_value(&self, c_re: f64, c_im: f64) -> f64 {
        let er2 = self.escape_radius * self.escape_radius;
        let (mut zr, mut zi, cr, ci) = if let Some((jcr, jci)) = self.julia_c {
            (c_re, c_im, jcr, jci)
        } else {
            (0.0, 0.0, c_re, c_im)
        };
        let mut min_dist = f64::INFINITY;
        for _ in 0..self.max_iter {
            let d = (zr * zr + zi * zi).sqrt();
            if d < min_dist {
                min_dist = d;
            }
            if zr * zr + zi * zi > er2 {
                break;
            }
            let new_zr = zr * zr - zi * zi + cr;
            zi = 2.0 * zr * zi + ci;
            zr = new_zr;
        }
        min_dist
    }

    /// Render the fractal to a flat array of escape values (row-major).
    ///
    /// Maps pixel `(i, j)` to complex plane via the viewport `[re_min, re_max] × [im_min, im_max]`.
    pub fn render(&self, re_min: f64, re_max: f64, im_min: f64, im_max: f64) -> Vec<f64> {
        let mut buf = vec![0.0_f64; self.width * self.height];
        for row in 0..self.height {
            let im = im_max - (im_max - im_min) * row as f64 / (self.height - 1).max(1) as f64;
            for col in 0..self.width {
                let re = re_min + (re_max - re_min) * col as f64 / (self.width - 1).max(1) as f64;
                let v = if self.coloring == FractalColoring::OrbitTrap {
                    self.orbit_trap_value(re, im)
                } else {
                    self.escape_value(re, im)
                };
                buf[row * self.width + col] = v;
            }
        }
        buf
    }
}

// ---------------------------------------------------------------------------
// NetworkViz — force-directed graph layout (Fruchterman-Reingold)
// ---------------------------------------------------------------------------

/// A node in the network graph.
#[derive(Debug, Clone)]
pub struct NetworkNode {
    /// Node identifier.
    pub id: usize,
    /// 2D position `[x, y]`.
    pub position: [f64; 2],
    /// Node degree (number of edges).
    pub degree: usize,
    /// Optional label.
    pub label: Option<String>,
    /// RGBA color.
    pub color: [f32; 4],
    /// Radius for rendering.
    pub radius: f64,
}

/// An edge in the network graph.
#[derive(Debug, Clone, Copy)]
pub struct NetworkEdge {
    /// Source node index.
    pub source: usize,
    /// Target node index.
    pub target: usize,
    /// Edge weight (thickness/strength).
    pub weight: f64,
}

/// Force-directed graph layout using Fruchterman-Reingold algorithm.
#[derive(Debug, Clone)]
pub struct NetworkViz {
    /// All nodes.
    pub nodes: Vec<NetworkNode>,
    /// All edges.
    pub edges: Vec<NetworkEdge>,
    /// Canvas width.
    pub width: f64,
    /// Canvas height.
    pub height: f64,
}

impl NetworkViz {
    /// Create a new network with the given canvas dimensions.
    pub fn new(width: f64, height: f64) -> Self {
        Self {
            nodes: Vec::new(),
            edges: Vec::new(),
            width,
            height,
        }
    }

    /// Add a node at position `[x, y]`.
    pub fn add_node(&mut self, id: usize, x: f64, y: f64) {
        self.nodes.push(NetworkNode {
            id,
            position: [x, y],
            degree: 0,
            label: None,
            color: [0.3, 0.6, 0.9, 1.0],
            radius: 5.0,
        });
    }

    /// Add an edge between node indices `src` and `tgt`.
    pub fn add_edge(&mut self, src: usize, tgt: usize, weight: f64) {
        self.edges.push(NetworkEdge {
            source: src,
            target: tgt,
            weight,
        });
        if src < self.nodes.len() {
            self.nodes[src].degree += 1;
        }
        if tgt < self.nodes.len() {
            self.nodes[tgt].degree += 1;
        }
    }

    /// Run Fruchterman-Reingold layout for `n_iter` iterations.
    ///
    /// `k` is the ideal edge length (typically `sqrt(area / N)`).
    pub fn layout_fruchterman_reingold(&mut self, n_iter: usize, k: f64) {
        let n = self.nodes.len();
        if n == 0 {
            return;
        }
        let mut disp: Vec<[f64; 2]> = vec![[0.0, 0.0]; n];
        let mut temperature = self.width / 10.0;
        let cooling = temperature / n_iter.max(1) as f64;

        for _ in 0..n_iter {
            // Reset displacements
            for d in disp.iter_mut() {
                *d = [0.0, 0.0];
            }

            // Repulsive forces (all pairs)
            for (i, disp_i) in disp.iter_mut().enumerate() {
                for j in 0..n {
                    if i == j {
                        continue;
                    }
                    let dx = self.nodes[i].position[0] - self.nodes[j].position[0];
                    let dy = self.nodes[i].position[1] - self.nodes[j].position[1];
                    let dist = (dx * dx + dy * dy).sqrt().max(1e-6);
                    let force = k * k / dist;
                    disp_i[0] += (dx / dist) * force;
                    disp_i[1] += (dy / dist) * force;
                }
            }

            // Attractive forces (edges)
            for edge in &self.edges {
                let s = edge.source.min(n - 1);
                let t = edge.target.min(n - 1);
                if s == t {
                    continue;
                }
                let dx = self.nodes[s].position[0] - self.nodes[t].position[0];
                let dy = self.nodes[s].position[1] - self.nodes[t].position[1];
                let dist = (dx * dx + dy * dy).sqrt().max(1e-6);
                let force = dist * dist / k;
                let fx = (dx / dist) * force;
                let fy = (dy / dist) * force;
                disp[s][0] -= fx;
                disp[s][1] -= fy;
                disp[t][0] += fx;
                disp[t][1] += fy;
            }

            // Apply displacement capped at temperature
            for (i, disp_i) in disp.iter().enumerate() {
                let disp_len = (disp_i[0] * disp_i[0] + disp_i[1] * disp_i[1])
                    .sqrt()
                    .max(1e-10);
                let capped = disp_len.min(temperature);
                self.nodes[i].position[0] += (disp_i[0] / disp_len) * capped;
                self.nodes[i].position[1] += (disp_i[1] / disp_len) * capped;
                // Clamp to canvas
                self.nodes[i].position[0] = self.nodes[i].position[0].clamp(0.0, self.width);
                self.nodes[i].position[1] = self.nodes[i].position[1].clamp(0.0, self.height);
            }

            temperature -= cooling;
        }
    }

    /// Color nodes by their degree (low = blue, high = red).
    pub fn color_by_degree(&mut self) {
        let max_deg = self
            .nodes
            .iter()
            .map(|n| n.degree)
            .max()
            .unwrap_or(1)
            .max(1);
        for node in self.nodes.iter_mut() {
            let t = node.degree as f32 / max_deg as f32;
            node.color = [t, 0.2, 1.0 - t, 1.0];
        }
    }
}

// ---------------------------------------------------------------------------
// PhasePortrait — vector field + trajectories
// ---------------------------------------------------------------------------

/// A trajectory in the phase plane (or 3D phase space).
#[derive(Debug, Clone)]
pub struct Trajectory {
    /// Sequence of state vectors `[x, y]` (or `[x, y, z]` for 3D).
    pub states: Vec<[f64; 2]>,
    /// Integration step size used.
    pub dt: f64,
}

/// A fixed point with its Jacobian-derived stability type.
#[derive(Debug, Clone, Copy)]
pub struct FixedPoint {
    /// Location `[x, y]`.
    pub location: [f64; 2],
    /// Stability type.
    pub stability: FixedPointType,
}

/// Classification of a fixed point based on eigenvalues of the Jacobian.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FixedPointType {
    /// Both eigenvalues have negative real part (stable node or focus).
    StableNode,
    /// Both eigenvalues have positive real part (unstable node or focus).
    UnstableNode,
    /// Eigenvalues are purely imaginary (center).
    Center,
    /// One positive, one negative eigenvalue (saddle point).
    Saddle,
    /// Stable spiral focus.
    StableFocus,
    /// Unstable spiral focus.
    UnstableFocus,
}

/// Phase portrait: vector field and trajectory data for a 2D ODE.
#[derive(Debug, Clone)]
pub struct PhasePortrait {
    /// Background vector field on a grid.
    pub vector_field: Vec<([f64; 2], [f64; 2])>, // (position, velocity)
    /// Integrated trajectories.
    pub trajectories: Vec<Trajectory>,
    /// Detected fixed points.
    pub fixed_points: Vec<FixedPoint>,
}

impl PhasePortrait {
    /// Compute a phase portrait for `f(x, y) -> (dx/dt, dy/dt)` on a grid.
    ///
    /// Uses 4th-order Runge-Kutta to integrate trajectories from `seeds`.
    pub fn compute<F>(
        f: &F,
        x_range: [f64; 2],
        y_range: [f64; 2],
        grid_n: usize,
        seeds: &[[f64; 2]],
        dt: f64,
        n_steps: usize,
    ) -> Self
    where
        F: Fn(f64, f64) -> (f64, f64),
    {
        // Build vector field on grid
        let mut vector_field = Vec::with_capacity(grid_n * grid_n);
        for ix in 0..grid_n {
            let x = x_range[0] + (x_range[1] - x_range[0]) * ix as f64 / grid_n.max(1) as f64;
            for iy in 0..grid_n {
                let y = y_range[0] + (y_range[1] - y_range[0]) * iy as f64 / grid_n.max(1) as f64;
                let (vx, vy) = f(x, y);
                vector_field.push(([x, y], [vx, vy]));
            }
        }

        // Integrate trajectories using RK4
        let trajectories = seeds
            .iter()
            .map(|&seed| {
                let mut states = Vec::with_capacity(n_steps + 1);
                let mut s = seed;
                states.push(s);
                for _ in 0..n_steps {
                    let [x, y] = s;
                    let (k1x, k1y) = f(x, y);
                    let (k2x, k2y) = f(x + 0.5 * dt * k1x, y + 0.5 * dt * k1y);
                    let (k3x, k3y) = f(x + 0.5 * dt * k2x, y + 0.5 * dt * k2y);
                    let (k4x, k4y) = f(x + dt * k3x, y + dt * k3y);
                    s[0] = x + dt / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
                    s[1] = y + dt / 6.0 * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
                    states.push(s);
                }
                Trajectory { states, dt }
            })
            .collect();

        Self {
            vector_field,
            trajectories,
            fixed_points: Vec::new(),
        }
    }

    /// Detect and classify fixed points near `candidates` using finite-difference Jacobian.
    pub fn classify_fixed_points<F>(&mut self, candidates: &[[f64; 2]], f: &F, eps: f64)
    where
        F: Fn(f64, f64) -> (f64, f64),
    {
        for &[x0, y0] in candidates {
            // Numerical Jacobian
            let (fp, _) = f(x0, y0);
            let (fm, _) = f(x0 - eps, y0);
            let (fpe, _) = f(x0, y0);
            let (_, fq) = f(x0, y0);
            let (_, fqm) = f(x0, y0 - eps);
            let a11 = (f(x0 + eps, y0).0 - fm) / (2.0 * eps);
            let a12 = (f(x0, y0 + eps).0 - fp) / (2.0 * eps);
            let a21 = (f(x0 + eps, y0).1 - f(x0 - eps, y0).1) / (2.0 * eps);
            let a22 = (fq - fqm) / eps;
            let _ = fpe; // silence unused

            let trace = a11 + a22;
            let det = a11 * a22 - a12 * a21;
            let disc = trace * trace - 4.0 * det;

            let stability = if det < 0.0 {
                FixedPointType::Saddle
            } else if disc >= 0.0 {
                if trace < 0.0 {
                    FixedPointType::StableNode
                } else {
                    FixedPointType::UnstableNode
                }
            } else if trace.abs() < 1e-8 {
                FixedPointType::Center
            } else if trace < 0.0 {
                FixedPointType::StableFocus
            } else {
                FixedPointType::UnstableFocus
            };

            self.fixed_points.push(FixedPoint {
                location: [x0, y0],
                stability,
            });
        }
    }
}

// ---------------------------------------------------------------------------
// Hodograph — velocity hodograph (v_x vs v_y)
// ---------------------------------------------------------------------------

/// A velocity hodograph: plots the velocity vector tip's trajectory.
#[derive(Debug, Clone, Default)]
pub struct Hodograph {
    /// Sequence of velocity vectors `[vx, vy]`.
    pub velocities: Vec<[f64; 2]>,
    /// Corresponding times.
    pub times: Vec<f64>,
    /// Label for the series.
    pub label: String,
}

impl Hodograph {
    /// Create a new empty hodograph.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a velocity sample.
    pub fn push(&mut self, time: f64, vx: f64, vy: f64) {
        self.times.push(time);
        self.velocities.push([vx, vy]);
    }

    /// Compute the arc length of the hodograph curve.
    pub fn arc_length(&self) -> f64 {
        let mut total = 0.0;
        for i in 1..self.velocities.len() {
            let dvx = self.velocities[i][0] - self.velocities[i - 1][0];
            let dvy = self.velocities[i][1] - self.velocities[i - 1][1];
            total += (dvx * dvx + dvy * dvy).sqrt();
        }
        total
    }

    /// Compute mean speed.
    pub fn mean_speed(&self) -> f64 {
        if self.velocities.is_empty() {
            return 0.0;
        }
        let sum: f64 = self
            .velocities
            .iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1]).sqrt())
            .sum();
        sum / self.velocities.len() as f64
    }

    /// Generate a hodograph from projectile motion (no drag).
    ///
    /// `v0x`, `v0y`: initial velocity components; `g`: gravitational acceleration.
    pub fn projectile(v0x: f64, v0y: f64, g: f64, dt: f64, n_steps: usize) -> Self {
        let mut h = Hodograph::new();
        for i in 0..=n_steps {
            let t = i as f64 * dt;
            h.push(t, v0x, v0y - g * t);
        }
        h
    }
}

// ---------------------------------------------------------------------------
// LissajousFigures — parametric curves for signal analysis
// ---------------------------------------------------------------------------

/// A single Lissajous figure.
#[derive(Debug, Clone)]
pub struct LissajousFigure {
    /// X-axis frequency.
    pub freq_x: f64,
    /// Y-axis frequency.
    pub freq_y: f64,
    /// Phase offset (radians).
    pub phase: f64,
    /// Amplitude for x.
    pub amp_x: f64,
    /// Amplitude for y.
    pub amp_y: f64,
    /// Sampled points `[x, y]`.
    pub points: Vec<[f64; 2]>,
}

impl LissajousFigure {
    /// Generate a Lissajous figure with `n_pts` samples over `t_range`.
    pub fn new(
        freq_x: f64,
        freq_y: f64,
        phase: f64,
        amp_x: f64,
        amp_y: f64,
        t_range: [f64; 2],
        n_pts: usize,
    ) -> Self {
        let points = (0..n_pts)
            .map(|i| {
                let t =
                    t_range[0] + (t_range[1] - t_range[0]) * i as f64 / (n_pts - 1).max(1) as f64;
                [
                    amp_x * (freq_x * t).sin(),
                    amp_y * (freq_y * t + phase).sin(),
                ]
            })
            .collect();
        Self {
            freq_x,
            freq_y,
            phase,
            amp_x,
            amp_y,
            points,
        }
    }

    /// Return true if the figure is a closed Lissajous curve (rational ratio).
    pub fn is_closed(&self, tol: f64) -> bool {
        if self.freq_y.abs() < 1e-14 {
            return false;
        }
        let ratio = self.freq_x / self.freq_y;
        // Check if ratio is close to a rational number p/q with small denominator
        for q in 1..=20_u32 {
            let p = (ratio * q as f64).round() as u32;
            if (ratio - p as f64 / q as f64).abs() < tol {
                return true;
            }
        }
        false
    }

    /// Compute the bounding box `([x_min, x_max], [y_min, y_max])`.
    pub fn bounding_box(&self) -> ([f64; 2], [f64; 2]) {
        let mut xmin = f64::INFINITY;
        let mut xmax = f64::NEG_INFINITY;
        let mut ymin = f64::INFINITY;
        let mut ymax = f64::NEG_INFINITY;
        for &[x, y] in &self.points {
            if x < xmin {
                xmin = x;
            }
            if x > xmax {
                xmax = x;
            }
            if y < ymin {
                ymin = y;
            }
            if y > ymax {
                ymax = y;
            }
        }
        ([xmin, xmax], [ymin, ymax])
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- TensorGlyph ----

    #[test]
    fn tensor_glyph_ellipsoid_axes_match_eigenvalues() {
        let evals = [3.0, 2.0, 1.0];
        let evecs = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let g = TensorGlyph::new([0.0, 0.0, 0.0], evals, evecs, GlyphShape::Ellipsoid, 1.0);
        let axes = g.ellipsoid_axes();
        assert!((axes[0] - 3.0).abs() < 1e-10);
        assert!((axes[1] - 2.0).abs() < 1e-10);
        assert!((axes[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn tensor_glyph_scale_applied_to_axes() {
        let evals = [2.0, 1.0, 0.5];
        let evecs = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let g = TensorGlyph::new([0.0, 0.0, 0.0], evals, evecs, GlyphShape::Ellipsoid, 2.0);
        let axes = g.ellipsoid_axes();
        assert!((axes[0] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn tensor_glyph_isotropic_fa_is_zero() {
        let evals = [1.0, 1.0, 1.0];
        let evecs = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let g = TensorGlyph::new([0.0, 0.0, 0.0], evals, evecs, GlyphShape::Ellipsoid, 1.0);
        assert!(g.fractional_anisotropy() < 1e-10);
    }

    #[test]
    fn tensor_glyph_anisotropic_fa_near_one() {
        let evals = [100.0, 0.01, 0.01];
        let evecs = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let g = TensorGlyph::new([0.0, 0.0, 0.0], evals, evecs, GlyphShape::Ellipsoid, 1.0);
        let fa = g.fractional_anisotropy();
        assert!(
            fa > 0.9,
            "FA={fa:.6} should be close to 1 for very anisotropic tensor"
        );
    }

    #[test]
    fn tensor_glyph_surface_points_count() {
        let evals = [2.0, 1.5, 1.0];
        let evecs = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let g = TensorGlyph::new([0.0, 0.0, 0.0], evals, evecs, GlyphShape::Ellipsoid, 1.0);
        let pts = g.ellipsoid_surface_points(16, 8);
        assert_eq!(pts.len(), 16 * 8);
    }

    #[test]
    fn build_tensor_glyphs_from_diagonal_tensor() {
        let tensors = vec![[4.0, 1.0, 0.25, 0.0, 0.0, 0.0]];
        let positions = vec![[0.0, 0.0, 0.0]];
        let glyphs = build_tensor_glyphs(&tensors, &positions, GlyphShape::Ellipsoid, 1.0);
        assert_eq!(glyphs.len(), 1);
        // All eigenvalues should be real and match the diagonal entries approximately
        let evals = glyphs[0].eigenvalues;
        let sorted: Vec<f64> = {
            let mut v = evals.to_vec();
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            v
        };
        assert!(
            (sorted[2] - 4.0).abs() < 0.1,
            "largest eigenval≈4: {:.6}",
            sorted[2]
        );
    }

    #[test]
    fn tensor_glyph_linear_anisotropy_extreme() {
        let evals = [10.0, 10.0, 0.001];
        let evecs = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let g = TensorGlyph::new([0.0, 0.0, 0.0], evals, evecs, GlyphShape::Ellipsoid, 1.0);
        assert!(g.planar_anisotropy() > 0.5);
    }

    // ---- PhaseDiagramViz ----

    #[test]
    fn phase_diagram_clausius_clapeyron_points_count() {
        let pts =
            PhaseDiagramViz::clausius_clapeyron_line([300.0, 400.0], 50, 40_000.0, 373.0, 101325.0);
        assert_eq!(pts.len(), 50);
    }

    #[test]
    fn phase_diagram_clausius_clapeyron_pressure_increases() {
        let pts =
            PhaseDiagramViz::clausius_clapeyron_line([300.0, 500.0], 10, 40_000.0, 373.0, 101325.0);
        // Vapor pressure should increase with temperature
        assert!(
            pts[9][1] > pts[0][1],
            "vapor pressure should increase with T"
        );
    }

    #[test]
    fn phase_diagram_tie_line_at_correct_temperature() {
        let mut pd = PhaseDiagramViz::new();
        let temps = vec![400.0, 450.0, 500.0];
        pd.compute_tie_lines(&temps, |t| 0.3 + 0.001 * t, |t| 0.1 + 0.0005 * t);
        assert_eq!(pd.tie_lines.len(), 3);
        assert!((pd.tie_lines[1].temperature - 450.0).abs() < 1e-10);
    }

    #[test]
    fn phase_diagram_triple_and_critical_point() {
        let mut pd = PhaseDiagramViz::new();
        pd.set_triple_point(273.16, 611.73);
        pd.set_critical_point(647.1, 22.064e6);
        assert!(pd.triple_point.is_some());
        assert!(pd.critical_point.is_some());
        let tp = pd.triple_point.unwrap();
        assert!((tp.temperature - 273.16).abs() < 1e-6);
    }

    #[test]
    fn phase_diagram_add_region_and_line() {
        let mut pd = PhaseDiagramViz::new();
        pd.add_region(PhaseRegion {
            name: "Liquid".to_string(),
            color: [0.2, 0.5, 1.0, 0.4],
            boundary: vec![[200.0, 1e3], [300.0, 1e5]],
        });
        pd.add_coexistence_line(vec![[273.0, 1e5], [373.0, 2e5]]);
        assert_eq!(pd.regions.len(), 1);
        assert_eq!(pd.coexistence_lines.len(), 1);
    }

    // ---- ErrorBandPlot ----

    #[test]
    fn error_band_mean_within_band() {
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let samples: Vec<Vec<f64>> = (0..100)
            .map(|i| {
                x.iter()
                    .map(|&xi| xi * xi + (i as f64 - 50.0) * 0.01)
                    .collect()
            })
            .collect();
        let band = ErrorBandPlot::from_samples(x, &samples, 2.0);
        assert!(band.mean_within_band(), "mean should lie within ±2σ band");
    }

    #[test]
    fn error_band_upper_greater_than_lower() {
        let x = vec![0.0, 1.0, 2.0];
        let samples: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64; 3]).collect();
        let band = ErrorBandPlot::from_samples(x, &samples, 1.0);
        for (lo, hi) in band.lower.iter().zip(band.upper.iter()) {
            assert!(hi >= lo, "upper >= lower must hold");
        }
    }

    #[test]
    fn error_band_empty_samples() {
        let x = vec![0.0, 1.0];
        let band = ErrorBandPlot::from_samples(x, &[], 1.0);
        assert_eq!(band.mean.len(), 2);
        assert!(band.mean_within_band());
    }

    #[test]
    fn error_band_bootstrap_ci_structure() {
        let x = vec![0.0, 1.0, 2.0];
        let data: Vec<Vec<f64>> = (0..20).map(|i| vec![i as f64 * 0.1; 3]).collect();
        let band = ErrorBandPlot::bootstrap_ci(x, &data, 0.95, 50);
        assert_eq!(band.x.len(), 3);
        assert!(band.mean_within_band());
    }

    // ---- BifurcationDiagram ----

    #[test]
    fn bifurcation_logistic_period1_for_small_r() {
        // r < 3: logistic map has a single stable fixed point
        let diag = BifurcationDiagram::logistic_map([2.5, 2.9], 5, 500, 100);
        let stable = diag.stable_points();
        assert!(
            !stable.is_empty(),
            "should have stable points for r ∈ [2.5, 2.9]"
        );
    }

    #[test]
    fn bifurcation_logistic_chaotic_for_large_r() {
        // r ≈ 4.0: logistic map is fully chaotic
        let diag = BifurcationDiagram::logistic_map([3.9, 4.0], 5, 500, 200);
        let chaotic = diag.chaotic_points();
        assert!(!chaotic.is_empty(), "should have chaotic points for r ≈ 4");
    }

    #[test]
    fn bifurcation_diagram_points_in_unit_interval() {
        let diag = BifurcationDiagram::logistic_map([2.5, 4.0], 10, 100, 50);
        for p in &diag.points {
            assert!(p.x >= 0.0 && p.x <= 1.0, "x={:.6} not in [0,1]", p.x);
        }
    }

    // ---- FractalViz ----

    #[test]
    fn mandelbrot_origin_does_not_escape() {
        let f = FractalViz::mandelbrot(10, 10, 256);
        let v = f.escape_value(0.0, 0.0);
        // Origin z=0 maps to c=0, stays bounded → max_iter
        assert!(
            (v - 256.0).abs() < 10.0,
            "origin escape value={v:.6} should be near max_iter"
        );
    }

    #[test]
    fn mandelbrot_far_point_escapes_fast() {
        let f = FractalViz::mandelbrot(10, 10, 256);
        let v = f.escape_value(10.0, 10.0);
        assert!(v < 10.0, "far-away point should escape quickly, got {v:.6}");
    }

    #[test]
    fn julia_set_render_dimensions() {
        let f = FractalViz::julia(8, 6, 64, -0.7, 0.27);
        let buf = f.render(-1.5, 1.5, -1.0, 1.0);
        assert_eq!(buf.len(), 48);
    }

    #[test]
    fn mandelbrot_render_dimensions() {
        let f = FractalViz::mandelbrot(4, 4, 32);
        let buf = f.render(-2.0, 1.0, -1.5, 1.5);
        assert_eq!(buf.len(), 16);
    }

    #[test]
    fn fractal_orbit_trap_positive() {
        let f = FractalViz {
            width: 4,
            height: 4,
            max_iter: 64,
            escape_radius: 2.0,
            coloring: FractalColoring::OrbitTrap,
            julia_c: None,
        };
        let v = f.orbit_trap_value(0.0, 0.0);
        assert!(v >= 0.0, "orbit trap should be non-negative: {v:.6}");
    }

    #[test]
    fn fractal_escape_time_coloring() {
        let mut f = FractalViz::mandelbrot(4, 4, 100);
        f.coloring = FractalColoring::EscapeTime;
        let v = f.escape_value(2.5, 0.0);
        assert!(v < 100.0);
    }

    // ---- NetworkViz ----

    #[test]
    fn network_layout_nodes_on_canvas() {
        let mut net = NetworkViz::new(100.0, 100.0);
        for i in 0..5 {
            net.add_node(i, 50.0, 50.0);
        }
        net.add_edge(0, 1, 1.0);
        net.add_edge(1, 2, 1.0);
        let k = (100.0 * 100.0 / 5.0_f64).sqrt();
        net.layout_fruchterman_reingold(50, k);
        for node in &net.nodes {
            assert!(node.position[0] >= 0.0 && node.position[0] <= 100.0);
            assert!(node.position[1] >= 0.0 && node.position[1] <= 100.0);
        }
    }

    #[test]
    fn network_degree_updated_on_add_edge() {
        let mut net = NetworkViz::new(100.0, 100.0);
        net.add_node(0, 10.0, 10.0);
        net.add_node(1, 90.0, 90.0);
        net.add_edge(0, 1, 1.0);
        assert_eq!(net.nodes[0].degree, 1);
        assert_eq!(net.nodes[1].degree, 1);
    }

    #[test]
    fn network_color_by_degree_changes_colors() {
        let mut net = NetworkViz::new(100.0, 100.0);
        for i in 0..3 {
            net.add_node(i, 50.0, 50.0);
        }
        net.add_edge(0, 1, 1.0);
        net.add_edge(0, 2, 1.0);
        let orig_color = net.nodes[0].color;
        net.color_by_degree();
        // Node 0 has degree 2 vs nodes 1,2 with degree 1 → different color
        assert!(
            (net.nodes[0].color[0] - orig_color[0]).abs() > 0.01
                || (net.nodes[1].color[0] - orig_color[0]).abs() > 0.01
        );
    }

    // ---- PhasePortrait ----

    #[test]
    fn phase_portrait_fixed_point_at_zero_velocity() {
        // f(x,y) = (-x, -y): origin is a stable node
        // Integrate long enough with a small step to reach dist < 0.1
        let f = |x: f64, y: f64| (-x, -y);
        let seeds = [[-0.5_f64, 0.5], [0.5, -0.5]];
        let pp = PhasePortrait::compute(&f, [-1.0, 1.0], [-1.0, 1.0], 5, &seeds, 0.05, 200);
        // Trajectories should approach origin (e^{-10} ≈ 4.5e-5 << 0.1)
        for traj in &pp.trajectories {
            let last = traj.states.last().unwrap();
            let dist = (last[0] * last[0] + last[1] * last[1]).sqrt();
            assert!(
                dist < 0.1,
                "trajectory should converge to origin: dist={dist:.6}"
            );
        }
    }

    #[test]
    fn phase_portrait_vector_field_count() {
        let f = |x: f64, y: f64| (y, -x);
        let pp = PhasePortrait::compute(&f, [-1.0, 1.0], [-1.0, 1.0], 8, &[], 0.01, 10);
        assert_eq!(pp.vector_field.len(), 64);
    }

    #[test]
    fn phase_portrait_trajectory_length() {
        let f = |_x: f64, _y: f64| (0.0, 0.0);
        let seeds = [[0.5_f64, 0.5]];
        let pp = PhasePortrait::compute(&f, [-1.0, 1.0], [-1.0, 1.0], 4, &seeds, 0.1, 20);
        assert_eq!(pp.trajectories[0].states.len(), 21);
    }

    #[test]
    fn phase_portrait_classify_stable_node() {
        let f = |x: f64, y: f64| (-2.0 * x, -3.0 * y);
        let mut pp = PhasePortrait::compute(&f, [-1.0, 1.0], [-1.0, 1.0], 4, &[], 0.01, 10);
        pp.classify_fixed_points(&[[0.0, 0.0]], &f, 1e-5);
        assert_eq!(pp.fixed_points.len(), 1);
        let fp = pp.fixed_points[0];
        assert!(
            fp.stability == FixedPointType::StableNode
                || fp.stability == FixedPointType::StableFocus,
            "origin of (-2x,-3y) should be stable: {:?}",
            fp.stability
        );
    }

    // ---- Hodograph ----

    #[test]
    fn hodograph_projectile_linear_vy() {
        let h = Hodograph::projectile(10.0, 20.0, 9.81, 0.1, 10);
        // vx stays constant
        for v in &h.velocities {
            assert!((v[0] - 10.0).abs() < 1e-10, "vx should be constant");
        }
        // vy decreases linearly
        assert!(h.velocities.last().unwrap()[1] < h.velocities[0][1]);
    }

    #[test]
    fn hodograph_arc_length_positive() {
        let h = Hodograph::projectile(5.0, 10.0, 9.81, 0.1, 20);
        assert!(h.arc_length() > 0.0);
    }

    #[test]
    fn hodograph_mean_speed_correct() {
        let mut h = Hodograph::new();
        h.push(0.0, 3.0, 4.0); // speed = 5
        h.push(1.0, 0.0, 5.0); // speed = 5
        let ms = h.mean_speed();
        assert!((ms - 5.0).abs() < 1e-10, "mean speed={ms:.6}");
    }

    // ---- LissajousFigures ----

    #[test]
    fn lissajous_point_count() {
        let fig = LissajousFigure::new(1.0, 2.0, 0.0, 1.0, 1.0, [0.0, 2.0 * PI], 200);
        assert_eq!(fig.points.len(), 200);
    }

    #[test]
    fn lissajous_amplitude_bounds() {
        let amp_x = 3.0;
        let amp_y = 2.0;
        let fig = LissajousFigure::new(1.0, 2.0, PI / 4.0, amp_x, amp_y, [0.0, 4.0 * PI], 500);
        for &[x, y] in &fig.points {
            assert!(
                x.abs() <= amp_x + 1e-10,
                "x={x:.6} exceeds amplitude {amp_x}"
            );
            assert!(
                y.abs() <= amp_y + 1e-10,
                "y={y:.6} exceeds amplitude {amp_y}"
            );
        }
    }

    #[test]
    fn lissajous_integer_ratio_is_closed() {
        // 1:2 ratio → closed figure
        let fig = LissajousFigure::new(1.0, 2.0, 0.0, 1.0, 1.0, [0.0, 2.0 * PI], 100);
        assert!(
            fig.is_closed(1e-6),
            "1:2 frequency ratio should be a closed curve"
        );
    }

    #[test]
    fn lissajous_bounding_box_within_amplitude() {
        let fig = LissajousFigure::new(3.0, 2.0, 0.5, 2.0, 1.5, [0.0, 10.0 * PI], 1000);
        let (xb, yb) = fig.bounding_box();
        assert!(xb[0] >= -2.0 - 1e-10 && xb[1] <= 2.0 + 1e-10);
        assert!(yb[0] >= -1.5 - 1e-10 && yb[1] <= 1.5 + 1e-10);
    }

    #[test]
    fn lissajous_circle_at_equal_freqs_and_phase() {
        // freq_x = freq_y, phase = π/2 → circle
        let fig = LissajousFigure::new(1.0, 1.0, PI / 2.0, 1.0, 1.0, [0.0, 2.0 * PI], 360);
        // All points should lie on unit circle
        for &[x, y] in &fig.points {
            let r = (x * x + y * y).sqrt();
            assert!((r - 1.0).abs() < 0.01, "r={r:.6} should be 1 for a circle");
        }
    }
}
