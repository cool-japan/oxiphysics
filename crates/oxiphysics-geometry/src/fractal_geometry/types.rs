//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use rand::RngExt;
use std::collections::HashMap;

/// Lacunarity computation for fractal point sets.
///
/// Lacunarity measures the "gappiness" or texture of a fractal. A higher
/// lacunarity indicates more heterogeneity (more gaps).
pub struct Lacunarity;
impl Lacunarity {
    /// Compute the gliding-box lacunarity of a 2D binary image (grid).
    ///
    /// `grid` is a 2D boolean array where `true` indicates the fractal is present.
    /// `box_size` is the side length of the gliding box.
    ///
    /// Returns the lacunarity value Lambda = Var(mass) / Mean(mass)^2 + 1.
    pub fn gliding_box(grid: &[Vec<bool>], box_size: usize) -> f64 {
        let rows = grid.len();
        if rows == 0 || box_size == 0 {
            return 0.0;
        }
        let cols = grid[0].len();
        if box_size > rows || box_size > cols {
            return 0.0;
        }
        let mut masses = Vec::new();
        for r in 0..=(rows - box_size) {
            for c in 0..=(cols - box_size) {
                let mut mass = 0u64;
                for dr in 0..box_size {
                    for dc in 0..box_size {
                        if grid[r + dr][c + dc] {
                            mass += 1;
                        }
                    }
                }
                masses.push(mass as f64);
            }
        }
        if masses.is_empty() {
            return 0.0;
        }
        let n = masses.len() as f64;
        let mean = masses.iter().sum::<f64>() / n;
        if mean.abs() < 1e-15 {
            return 0.0;
        }
        let variance = masses.iter().map(|&m| (m - mean).powi(2)).sum::<f64>() / n;
        variance / (mean * mean) + 1.0
    }
    /// Compute lacunarity for a point set by first rasterizing onto a grid.
    ///
    /// `resolution` is the grid side length. Points are mapped to grid cells.
    pub fn from_points(points: &[Point2], resolution: usize, box_size: usize) -> f64 {
        if points.is_empty() || resolution == 0 {
            return 0.0;
        }
        let (min_x, max_x, min_y, max_y) = FractalDimension::bounding_box(points);
        let dx = (max_x - min_x).max(1e-15);
        let dy = (max_y - min_y).max(1e-15);
        let mut grid = vec![vec![false; resolution]; resolution];
        for p in points {
            let ix = (((p.x - min_x) / dx) * (resolution - 1) as f64)
                .round()
                .clamp(0.0, (resolution - 1) as f64) as usize;
            let iy = (((p.y - min_y) / dy) * (resolution - 1) as f64)
                .round()
                .clamp(0.0, (resolution - 1) as f64) as usize;
            grid[iy][ix] = true;
        }
        Self::gliding_box(&grid, box_size)
    }
    /// Multi-scale lacunarity: compute lacunarity for a range of box sizes.
    ///
    /// Returns `(box_size, lacunarity)` pairs.
    pub fn multiscale(grid: &[Vec<bool>], box_sizes: &[usize]) -> Vec<(usize, f64)> {
        box_sizes
            .iter()
            .map(|&bs| (bs, Self::gliding_box(grid, bs)))
            .collect()
    }
}
/// Sierpinski arrowhead curve L-system preset.
pub struct SierpinskiArrowhead;
impl SierpinskiArrowhead {
    /// Create the L-system for a Sierpinski arrowhead curve.
    pub fn lsystem() -> LSystem {
        LSystem::new(
            vec!['F', 'G', '+', '-'],
            vec![
                ProductionRule::new('F', "G-F-G"),
                ProductionRule::new('G', "F+G+F"),
            ],
            "F",
        )
    }
    /// Theoretical dimension: log(3)/log(2) ≈ 1.585.
    pub fn theoretical_dimension() -> f64 {
        (3.0_f64).ln() / (2.0_f64).ln()
    }
    /// Generate Sierpinski arrowhead segments.
    pub fn generate(generations: usize, step_size: f64) -> Vec<(Point2, Point2)> {
        let ls = Self::lsystem();
        ls.turtle_interpret(generations, step_size, 60.0)
    }
}
/// An Iterated Function System consisting of affine transforms with associated probabilities.
#[derive(Debug, Clone)]
pub struct IteratedFunctionSystem {
    /// The affine transforms.
    pub transforms: Vec<AffineTransform2D>,
    /// Probability weights for each transform (must sum to ~1.0).
    pub probabilities: Vec<f64>,
}
impl IteratedFunctionSystem {
    /// Create a new IFS with given transforms and probabilities.
    ///
    /// Handles empty transforms gracefully (returns empty results).
    pub fn new(transforms: Vec<AffineTransform2D>, probabilities: Vec<f64>) -> Self {
        if transforms.is_empty() {
            return Self {
                transforms: Vec::new(),
                probabilities: Vec::new(),
            };
        }
        assert_eq!(
            transforms.len(),
            probabilities.len(),
            "transforms and probabilities must have equal length"
        );
        Self {
            transforms,
            probabilities,
        }
    }
    /// Create an IFS with equal probabilities for all transforms.
    pub fn with_equal_probabilities(transforms: Vec<AffineTransform2D>) -> Self {
        let n = transforms.len();
        assert!(n > 0, "IFS must have at least one transform");
        let prob = 1.0 / n as f64;
        let probabilities = vec![prob; n];
        Self {
            transforms,
            probabilities,
        }
    }
    /// Run the chaos game algorithm to generate attractor points.
    ///
    /// Starting from `start`, iterates `iterations` times selecting transforms
    /// according to their probabilities. Returns the generated points after
    /// discarding the first `skip` transient points.
    pub fn chaos_game(&self, start: Point2, iterations: usize, skip: usize) -> Vec<Point2> {
        use rand::RngExt;
        let mut rng = rand::rng();
        let mut point = start;
        let mut points = Vec::with_capacity(iterations.saturating_sub(skip));
        let cumulative = self.cumulative_distribution();
        for i in 0..iterations {
            let r: f64 = rng.random_range(0.0..1.0);
            let idx = cumulative
                .iter()
                .position(|&c| r < c)
                .unwrap_or(self.transforms.len() - 1);
            point = self.transforms[idx].apply(&point);
            if i >= skip {
                points.push(point);
            }
        }
        points
    }
    /// Deterministic iteration: apply all transforms to every point in the set.
    ///
    /// Starting from `initial_points`, applies all transforms to produce
    /// the next generation. Repeats for `generations` steps.
    pub fn deterministic_iterate(
        &self,
        initial_points: &[Point2],
        generations: usize,
    ) -> Vec<Point2> {
        let mut current = initial_points.to_vec();
        for _ in 0..generations {
            let mut next = Vec::with_capacity(current.len() * self.transforms.len());
            for t in &self.transforms {
                for p in &current {
                    next.push(t.apply(p));
                }
            }
            current = next;
        }
        current
    }
    /// Number of transforms in this IFS.
    pub fn num_transforms(&self) -> usize {
        self.transforms.len()
    }
    /// Create the Barnsley Fern IFS.
    pub fn barnsley_fern() -> Self {
        let transforms = vec![
            AffineTransform2D::new(0.0, 0.0, 0.0, 0.16, 0.0, 0.0),
            AffineTransform2D::new(0.85, 0.04, -0.04, 0.85, 0.0, 1.6),
            AffineTransform2D::new(0.20, -0.26, 0.23, 0.22, 0.0, 1.6),
            AffineTransform2D::new(-0.15, 0.28, 0.26, 0.24, 0.0, 0.44),
        ];
        let probabilities = vec![0.01, 0.85, 0.07, 0.07];
        Self::new(transforms, probabilities)
    }
    /// Generate `n` points using a deterministic seeded chaos game.
    ///
    /// Returns points as `[x, y]` arrays.
    pub fn generate(&self, n: usize, seed: u64) -> Vec<[f64; 2]> {
        if self.transforms.is_empty() {
            return Vec::new();
        }
        let cumulative = self.cumulative_distribution();
        let mut x = 0.0_f64;
        let mut y = 0.0_f64;
        let mut rng_state = seed;
        let mut points = Vec::with_capacity(n);
        // Simple LCG pseudo-random for reproducibility
        for _ in 0..n {
            rng_state = rng_state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let r = ((rng_state >> 33) as f64) / (u32::MAX as f64);
            let idx = cumulative
                .iter()
                .position(|&c| r < c)
                .unwrap_or(self.transforms.len() - 1);
            let t = &self.transforms[idx];
            let nx = t.a * x + t.b * y + t.e;
            let ny = t.c * x + t.d * y + t.f;
            x = nx;
            y = ny;
            points.push([x, y]);
        }
        points
    }
    /// Build cumulative distribution from probabilities.
    fn cumulative_distribution(&self) -> Vec<f64> {
        let mut cum = Vec::with_capacity(self.probabilities.len());
        let mut acc = 0.0;
        for &p in &self.probabilities {
            acc += p;
            cum.push(acc);
        }
        cum
    }
}
/// A production rule for an L-system.
#[derive(Debug, Clone)]
pub struct ProductionRule {
    /// The predecessor character.
    pub predecessor: char,
    /// The successor string.
    pub successor: String,
}
impl ProductionRule {
    /// Create a new production rule.
    pub fn new(predecessor: char, successor: impl Into<String>) -> Self {
        Self {
            predecessor,
            successor: successor.into(),
        }
    }
}
/// The Mandelbrot set: z_{n+1} = z_n^2 + c, starting from z_0 = 0.
///
/// Points `c` in the complex plane are classified by the escape time
/// of the orbit {z_n}.
pub struct MandelbrotSet {
    /// Maximum number of iterations before declaring a point in the set.
    pub max_iter: u32,
    /// Escape radius squared.
    pub escape_radius_sq: f64,
}
impl MandelbrotSet {
    /// Create a new Mandelbrot set with given parameters.
    pub fn new(max_iter: u32, escape_radius: f64) -> Self {
        Self {
            max_iter,
            escape_radius_sq: escape_radius * escape_radius,
        }
    }
    /// Default Mandelbrot set with max_iter=256, escape_radius=2.
    pub fn default_set() -> Self {
        Self::new(256, 2.0)
    }
    /// Compute escape time for a complex number c = (cr, ci).
    ///
    /// Returns `None` if the orbit does not escape within `max_iter` iterations,
    /// or `Some(iterations)` with the escape iteration count.
    pub fn escape_time(&self, cr: f64, ci: f64) -> Option<u32> {
        let mut zr = 0.0;
        let mut zi = 0.0;
        for i in 0..self.max_iter {
            let zr2 = zr * zr;
            let zi2 = zi * zi;
            if zr2 + zi2 > self.escape_radius_sq {
                return Some(i);
            }
            zi = 2.0 * zr * zi + ci;
            zr = zr2 - zi2 + cr;
        }
        None
    }
    /// Smooth escape time using renormalization.
    ///
    /// Returns a continuous (smooth) iteration count for coloring.
    pub fn smooth_escape_time(&self, cr: f64, ci: f64) -> f64 {
        let mut zr = 0.0;
        let mut zi = 0.0;
        for i in 0..self.max_iter {
            let zr2 = zr * zr;
            let zi2 = zi * zi;
            if zr2 + zi2 > self.escape_radius_sq {
                let modulus = (zr2 + zi2).sqrt();
                let mu = i as f64 - modulus.ln().ln() / (2.0_f64).ln();
                return mu;
            }
            zi = 2.0 * zr * zi + ci;
            zr = zr2 - zi2 + cr;
        }
        self.max_iter as f64
    }
    /// Check if a point is in the Mandelbrot set.
    pub fn is_in_set(&self, cr: f64, ci: f64) -> bool {
        self.escape_time(cr, ci).is_none()
    }
    /// Detect if a point is near the boundary of the set.
    ///
    /// A point is near the boundary if it is in the set but has a neighbor
    /// that is not (or vice versa), within the given `delta` distance.
    pub fn is_boundary(&self, cr: f64, ci: f64, delta: f64) -> bool {
        let center = self.is_in_set(cr, ci);
        let offsets = [(delta, 0.0), (-delta, 0.0), (0.0, delta), (0.0, -delta)];
        for (dr, di) in offsets {
            if self.is_in_set(cr + dr, ci + di) != center {
                return true;
            }
        }
        false
    }
    /// Compute the area of the Mandelbrot set using Monte Carlo sampling.
    ///
    /// Samples random points in the rectangle \[x_min, x_max\] x \[y_min, y_max\]
    /// and estimates the area as the fraction of points in the set times the
    /// rectangle area.
    pub fn area_monte_carlo(
        &self,
        x_min: f64,
        x_max: f64,
        y_min: f64,
        y_max: f64,
        num_samples: usize,
    ) -> f64 {
        let mut rng = rand::rng();
        let mut count = 0usize;
        for _ in 0..num_samples {
            let cr = rng.random_range(x_min..x_max);
            let ci = rng.random_range(y_min..y_max);
            if self.is_in_set(cr, ci) {
                count += 1;
            }
        }
        let rect_area = (x_max - x_min) * (y_max - y_min);
        rect_area * count as f64 / num_samples as f64
    }
    /// Generate a grid of escape times for rendering.
    ///
    /// Returns a row-major `width x height` grid of escape times.
    pub fn render_grid(
        &self,
        x_min: f64,
        x_max: f64,
        y_min: f64,
        y_max: f64,
        width: usize,
        height: usize,
    ) -> Vec<Vec<Option<u32>>> {
        let dx = (x_max - x_min) / width as f64;
        let dy = (y_max - y_min) / height as f64;
        (0..height)
            .map(|j| {
                let ci = y_min + (j as f64 + 0.5) * dy;
                (0..width)
                    .map(|i| {
                        let cr = x_min + (i as f64 + 0.5) * dx;
                        self.escape_time(cr, ci)
                    })
                    .collect()
            })
            .collect()
    }
}
/// Hilbert curve coordinate mapping (non-L-system approach).
///
/// Provides index-to-coordinate and coordinate-to-index conversions
/// for Hilbert curves of order `n` (grid size `2^n x 2^n`).
pub struct HilbertCurveMap;
impl HilbertCurveMap {
    /// Convert a Hilbert curve index `d` to (x, y) coordinates for a curve of order `n`.
    pub fn d2xy(n: u32, d: u32) -> (u32, u32) {
        let mut rx: u32;
        let mut ry: u32;
        let mut x = 0_u32;
        let mut y = 0_u32;
        let mut d = d;
        let mut s = 1_u32;
        while s < n {
            rx = if (d & 2) != 0 { 1 } else { 0 };
            ry = if (d & 1) != 0 {
                if rx == 0 { 1 } else { 0 }
            } else {
                0
            };
            Self::rot(s, &mut x, &mut y, rx, ry);
            x += s * rx;
            y += s * ry;
            d /= 4;
            s *= 2;
        }
        (x, y)
    }
    /// Convert (x, y) coordinates to a Hilbert curve index for a curve of order `n`.
    pub fn xy2d(n: u32, x: u32, y: u32) -> u32 {
        let mut rx: u32;
        let mut ry: u32;
        let mut d = 0_u32;
        let mut x = x;
        let mut y = y;
        let mut s = n / 2;
        while s > 0 {
            rx = if (x & s) > 0 { 1 } else { 0 };
            ry = if (y & s) > 0 { 1 } else { 0 };
            d += s * s * ((3 * rx) ^ ry);
            Self::rot(s, &mut x, &mut y, rx, ry);
            s /= 2;
        }
        d
    }
    /// Rotate/flip a quadrant.
    fn rot(n: u32, x: &mut u32, y: &mut u32, rx: u32, ry: u32) {
        if ry == 0 {
            if rx == 1 {
                *x = n.wrapping_sub(1).wrapping_sub(*x);
                *y = n.wrapping_sub(1).wrapping_sub(*y);
            }
            std::mem::swap(x, y);
        }
    }
    /// Generate the full Hilbert curve path as a sequence of (x, y) for order `n`.
    ///
    /// `n` must be a power of 2. Returns `n*n` points.
    pub fn generate_path(n: u32) -> Vec<(u32, u32)> {
        let total = n * n;
        (0..total).map(|d| Self::d2xy(n, d)).collect()
    }
}
/// Fractal terrain generator using the diamond-square algorithm.
///
/// Produces a heightmap on a (2^n + 1) x (2^n + 1) grid.
pub struct FractalTerrain;
impl FractalTerrain {
    /// Generate a terrain heightmap using the diamond-square algorithm.
    ///
    /// `power` determines the grid size: side = 2^power + 1.
    /// `roughness` controls how quickly the random displacement decreases (0..1 typical).
    /// `seed_corners` sets the four corner heights `[top_left, top_right, bot_left, bot_right]`.
    pub fn diamond_square(power: u32, roughness: f64, seed_corners: [f64; 4]) -> Vec<Vec<f64>> {
        let size = (1 << power) + 1;
        let mut grid = vec![vec![0.0_f64; size]; size];
        let last = size - 1;
        grid[0][0] = seed_corners[0];
        grid[0][last] = seed_corners[1];
        grid[last][0] = seed_corners[2];
        grid[last][last] = seed_corners[3];
        let mut rng = rand::rng();
        let mut step = last;
        let mut scale = 1.0_f64;
        while step > 1 {
            let half = step / 2;
            let mut y = 0;
            while y < last {
                let mut x = 0;
                while x < last {
                    let avg = (grid[y][x]
                        + grid[y][x + step]
                        + grid[y + step][x]
                        + grid[y + step][x + step])
                        / 4.0;
                    grid[y + half][x + half] = avg + rng.random_range(-scale..scale);
                    x += step;
                }
                y += step;
            }
            let mut y = 0;
            while y <= last {
                let x_start = if (y / half).is_multiple_of(2) {
                    half
                } else {
                    0
                };
                let mut x = x_start;
                while x <= last {
                    let mut sum = 0.0;
                    let mut count = 0.0;
                    if y >= half {
                        sum += grid[y - half][x];
                        count += 1.0;
                    }
                    if y + half <= last {
                        sum += grid[y + half][x];
                        count += 1.0;
                    }
                    if x >= half {
                        sum += grid[y][x - half];
                        count += 1.0;
                    }
                    if x + half <= last {
                        sum += grid[y][x + half];
                        count += 1.0;
                    }
                    grid[y][x] = sum / count + rng.random_range(-scale..scale);
                    x += step;
                }
                y += half;
            }
            step = half;
            scale *= roughness;
        }
        grid
    }
    /// Generate terrain using midpoint displacement (1D variant).
    ///
    /// Produces a heightmap profile (1D array) of size 2^power + 1.
    /// `roughness` is the Hurst exponent factor (0 < roughness < 1).
    pub fn midpoint_displacement_1d(
        power: u32,
        roughness: f64,
        h_left: f64,
        h_right: f64,
    ) -> Vec<f64> {
        let size = (1 << power) + 1;
        let mut heights = vec![0.0_f64; size];
        heights[0] = h_left;
        heights[size - 1] = h_right;
        let mut rng = rand::rng();
        let mut step = size - 1;
        let mut scale = 1.0_f64;
        while step > 1 {
            let half = step / 2;
            let mut i = half;
            while i < size {
                let left = heights[i - half];
                let right = if i + half < size {
                    heights[i + half]
                } else {
                    heights[size - 1]
                };
                heights[i] = (left + right) / 2.0 + rng.random_range(-scale..scale);
                i += step;
            }
            step = half;
            scale *= roughness;
        }
        heights
    }
    /// Compute the RMS (root mean square) roughness of a heightmap.
    pub fn rms_roughness(grid: &[Vec<f64>]) -> f64 {
        let mut sum = 0.0;
        let mut count = 0usize;
        let mut mean = 0.0;
        for row in grid {
            for &h in row {
                mean += h;
                count += 1;
            }
        }
        if count == 0 {
            return 0.0;
        }
        mean /= count as f64;
        for row in grid {
            for &h in row {
                sum += (h - mean).powi(2);
            }
        }
        (sum / count as f64).sqrt()
    }
    /// Find the min and max height in a grid.
    pub fn height_range(grid: &[Vec<f64>]) -> (f64, f64) {
        let mut lo = f64::MAX;
        let mut hi = f64::MIN;
        for row in grid {
            for &h in row {
                if h < lo {
                    lo = h;
                }
                if h > hi {
                    hi = h;
                }
            }
        }
        (lo, hi)
    }
}
/// Simple fractal mesh: triangle mesh from a fractal terrain heightmap.
///
/// Converts a 2D heightmap grid into a triangle mesh with vertices and triangle
/// index lists suitable for rendering.
pub struct FractalMesh;
impl FractalMesh {
    /// Convert a heightmap to a triangle mesh.
    ///
    /// Returns `(vertices, triangles)` where vertices are `[f64; 3]` positions
    /// and triangles are index triples into the vertex array.
    pub fn from_heightmap(grid: &[Vec<f64>], spacing: f64) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let rows = grid.len();
        if rows == 0 {
            return (Vec::new(), Vec::new());
        }
        let cols = grid[0].len();
        let mut vertices = Vec::with_capacity(rows * cols);
        for (r, row) in grid.iter().enumerate() {
            for (c, &h) in row.iter().enumerate() {
                vertices.push([c as f64 * spacing, h, r as f64 * spacing]);
            }
        }
        let mut triangles = Vec::new();
        for r in 0..(rows - 1) {
            for c in 0..(cols - 1) {
                let i0 = r * cols + c;
                let i1 = i0 + 1;
                let i2 = i0 + cols;
                let i3 = i2 + 1;
                triangles.push([i0, i2, i1]);
                triangles.push([i1, i2, i3]);
            }
        }
        (vertices, triangles)
    }
    /// Compute the total surface area of a triangle mesh.
    pub fn surface_area(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> f64 {
        let mut area = 0.0;
        for tri in triangles {
            let a = vertices[tri[0]];
            let b = vertices[tri[1]];
            let c = vertices[tri[2]];
            let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let cross = [
                ab[1] * ac[2] - ab[2] * ac[1],
                ab[2] * ac[0] - ab[0] * ac[2],
                ab[0] * ac[1] - ab[1] * ac[0],
            ];
            area += (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt() / 2.0;
        }
        area
    }
    /// Number of vertices and triangles expected from a heightmap of given size.
    pub fn expected_counts(rows: usize, cols: usize) -> (usize, usize) {
        let verts = rows * cols;
        let tris = if rows > 1 && cols > 1 {
            2 * (rows - 1) * (cols - 1)
        } else {
            0
        };
        (verts, tris)
    }
}
/// Hilbert curve L-system preset.
pub struct HilbertCurve;
impl HilbertCurve {
    /// Create the L-system for a Hilbert curve.
    pub fn lsystem() -> LSystem {
        LSystem::new(
            vec!['A', 'B', 'F', '+', '-'],
            vec![
                ProductionRule::new('A', "+BF-AFA-FB+"),
                ProductionRule::new('B', "-AF+BFB+FA-"),
            ],
            "A",
        )
    }
    /// Generate Hilbert curve segments.
    pub fn generate(generations: usize, step_size: f64) -> Vec<(Point2, Point2)> {
        let ls = Self::lsystem();
        ls.turtle_interpret(generations, step_size, 90.0)
    }
}
/// Methods for estimating the fractal dimension of a point set.
pub struct FractalDimension;
impl FractalDimension {
    /// Estimate the box-counting (Minkowski) dimension of a 2D point set.
    ///
    /// Counts the number of boxes of side `epsilon` needed to cover the point set,
    /// for geometrically decreasing `epsilon` values. Returns `(dimension, r_squared)`
    /// where `r_squared` is the coefficient of determination of the linear fit.
    pub fn box_counting(points: &[Point2], num_scales: usize) -> (f64, f64) {
        if points.is_empty() || num_scales < 2 {
            return (0.0, 0.0);
        }
        let (min_x, max_x, min_y, max_y) = Self::bounding_box(points);
        let extent = (max_x - min_x).max(max_y - min_y).max(1e-15);
        let mut log_eps = Vec::with_capacity(num_scales);
        let mut log_n = Vec::with_capacity(num_scales);
        for i in 0..num_scales {
            let epsilon = extent / (2.0_f64).powi(i as i32 + 1);
            if epsilon < 1e-15 {
                break;
            }
            let mut occupied = std::collections::HashSet::new();
            for p in points {
                let ix = ((p.x - min_x) / epsilon).floor() as i64;
                let iy = ((p.y - min_y) / epsilon).floor() as i64;
                occupied.insert((ix, iy));
            }
            let count = occupied.len() as f64;
            if count > 0.0 {
                log_eps.push(epsilon.ln());
                log_n.push(count.ln());
            }
        }
        if log_eps.len() < 2 {
            return (0.0, 0.0);
        }
        let (slope, r_sq) = Self::linear_regression(&log_eps, &log_n);
        (-slope, r_sq)
    }
    /// Estimate the correlation dimension of a point set.
    ///
    /// Uses the Grassberger-Procaccia algorithm. For various radii `r`,
    /// counts pairs of points within distance `r`, then fits `log C(r)` vs `log r`.
    pub fn correlation_dimension(points: &[Point2], num_scales: usize) -> (f64, f64) {
        let n = points.len();
        if n < 10 || num_scales < 2 {
            return (0.0, 0.0);
        }
        let mut max_dist = 0.0_f64;
        let mut min_dist = f64::MAX;
        let num_pairs = n * (n - 1) / 2;
        let mut distances = Vec::with_capacity(num_pairs);
        for i in 0..n {
            for j in (i + 1)..n {
                let d = points[i].distance(&points[j]);
                distances.push(d);
                if d > max_dist {
                    max_dist = d;
                }
                if d > 0.0 && d < min_dist {
                    min_dist = d;
                }
            }
        }
        if max_dist < 1e-15 {
            return (0.0, 0.0);
        }
        distances.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut log_r = Vec::with_capacity(num_scales);
        let mut log_c = Vec::with_capacity(num_scales);
        let r_min = min_dist.max(max_dist * 0.001);
        let r_max = max_dist * 0.5;
        let ratio = (r_max / r_min).powf(1.0 / (num_scales as f64 - 1.0));
        for i in 0..num_scales {
            let r = r_min * ratio.powi(i as i32);
            let count = distances.partition_point(|&d| d <= r);
            if count > 0 {
                let c = count as f64 / num_pairs as f64;
                log_r.push(r.ln());
                log_c.push(c.ln());
            }
        }
        if log_r.len() < 2 {
            return (0.0, 0.0);
        }
        let (slope, r_sq) = Self::linear_regression(&log_r, &log_c);
        (slope, r_sq)
    }
    /// Estimate Hausdorff dimension using box-counting with Richardson extrapolation.
    ///
    /// This provides a more accurate estimate by using multiple scale ranges.
    pub fn hausdorff_estimate(points: &[Point2], num_scales: usize) -> f64 {
        if points.is_empty() || num_scales < 4 {
            return 0.0;
        }
        let (dim1, _) = Self::box_counting(points, num_scales);
        let half = num_scales / 2;
        let (dim2, _) = Self::box_counting(points, half);
        0.6 * dim1 + 0.4 * dim2
    }
    /// Compute bounding box of a point set.
    pub(super) fn bounding_box(points: &[Point2]) -> (f64, f64, f64, f64) {
        let mut min_x = f64::MAX;
        let mut max_x = f64::MIN;
        let mut min_y = f64::MAX;
        let mut max_y = f64::MIN;
        for p in points {
            if p.x < min_x {
                min_x = p.x;
            }
            if p.x > max_x {
                max_x = p.x;
            }
            if p.y < min_y {
                min_y = p.y;
            }
            if p.y > max_y {
                max_y = p.y;
            }
        }
        (min_x, max_x, min_y, max_y)
    }
    /// Simple linear regression returning (slope, r_squared).
    fn linear_regression(x: &[f64], y: &[f64]) -> (f64, f64) {
        let n = x.len() as f64;
        let sum_x: f64 = x.iter().sum();
        let sum_y: f64 = y.iter().sum();
        let sum_xy: f64 = x.iter().zip(y.iter()).map(|(xi, yi)| xi * yi).sum();
        let sum_x2: f64 = x.iter().map(|xi| xi * xi).sum();
        let sum_y2: f64 = y.iter().map(|yi| yi * yi).sum();
        let denom = n * sum_x2 - sum_x * sum_x;
        if denom.abs() < 1e-30 {
            return (0.0, 0.0);
        }
        let slope = (n * sum_xy - sum_x * sum_y) / denom;
        let ss_tot = sum_y2 - sum_y * sum_y / n;
        let ss_res = sum_y2 - slope * sum_xy - (sum_y - slope * sum_x) * sum_y / n;
        let r_sq = if ss_tot.abs() > 1e-30 {
            1.0 - ss_res / ss_tot
        } else {
            0.0
        };
        (slope, r_sq.clamp(0.0, 1.0))
    }
}
/// Sierpinski triangle IFS preset.
///
/// The Sierpinski triangle is constructed from 3 affine transforms,
/// each scaling by 1/2 and translating to a corner of the triangle.
pub struct SierpinskiTriangle;
impl SierpinskiTriangle {
    /// Create the IFS for a Sierpinski triangle.
    pub fn ifs() -> IteratedFunctionSystem {
        let transforms = vec![
            AffineTransform2D::new(0.5, 0.0, 0.0, 0.5, 0.0, 0.0),
            AffineTransform2D::new(0.5, 0.0, 0.0, 0.5, 0.5, 0.0),
            AffineTransform2D::new(0.5, 0.0, 0.0, 0.5, 0.25, 0.433),
        ];
        IteratedFunctionSystem::with_equal_probabilities(transforms)
    }
    /// Theoretical Hausdorff dimension: log(3)/log(2) ≈ 1.585.
    pub fn theoretical_dimension() -> f64 {
        (3.0_f64).ln() / (2.0_f64).ln()
    }
    /// Generate Sierpinski triangle points via chaos game.
    pub fn generate(num_points: usize) -> Vec<Point2> {
        let ifs = Self::ifs();
        ifs.chaos_game(Point2::new(0.0, 0.0), num_points + 100, 100)
    }
}
/// A 2D affine transformation: 2x2 matrix \[\[a, b\\], \[c, d\]] + translation (e, f).
///
/// Applies as: x' = a*x + b*y + e, y' = c*x + d*y + f
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AffineTransform2D {
    /// Matrix element (0,0).
    pub a: f64,
    /// Matrix element (0,1).
    pub b: f64,
    /// Matrix element (1,0).
    pub c: f64,
    /// Matrix element (1,1).
    pub d: f64,
    /// Translation x.
    pub e: f64,
    /// Translation y.
    pub f: f64,
}
impl AffineTransform2D {
    /// Create a new affine transform.
    pub fn new(a: f64, b: f64, c: f64, d: f64, e: f64, f: f64) -> Self {
        Self { a, b, c, d, e, f }
    }
    /// Identity transform.
    pub fn identity() -> Self {
        Self::new(1.0, 0.0, 0.0, 1.0, 0.0, 0.0)
    }
    /// Create a rotation transform.
    pub fn rotation(angle_rad: f64) -> Self {
        let cos = angle_rad.cos();
        let sin = angle_rad.sin();
        Self::new(cos, -sin, sin, cos, 0.0, 0.0)
    }
    /// Create a scaling transform.
    pub fn scaling(sx: f64, sy: f64) -> Self {
        Self::new(sx, 0.0, 0.0, sy, 0.0, 0.0)
    }
    /// Create a translation transform.
    pub fn translation(tx: f64, ty: f64) -> Self {
        Self::new(1.0, 0.0, 0.0, 1.0, tx, ty)
    }
    /// Apply this transform to a point.
    pub fn apply(&self, p: &Point2) -> Point2 {
        Point2 {
            x: self.a * p.x + self.b * p.y + self.e,
            y: self.c * p.x + self.d * p.y + self.f,
        }
    }
    /// Compose two transforms: self followed by other.
    ///
    /// The result transform applies self first, then other.
    pub fn compose(&self, other: &AffineTransform2D) -> AffineTransform2D {
        AffineTransform2D {
            a: other.a * self.a + other.b * self.c,
            b: other.a * self.b + other.b * self.d,
            c: other.c * self.a + other.d * self.c,
            d: other.c * self.b + other.d * self.d,
            e: other.a * self.e + other.b * self.f + other.e,
            f: other.c * self.e + other.d * self.f + other.f,
        }
    }
    /// Determinant of the linear part.
    pub fn determinant(&self) -> f64 {
        self.a * self.d - self.b * self.c
    }
    /// Inverse transform, if the determinant is nonzero.
    pub fn inverse(&self) -> Option<AffineTransform2D> {
        let det = self.determinant();
        if det.abs() < 1e-15 {
            return None;
        }
        let inv_det = 1.0 / det;
        let a = self.d * inv_det;
        let b = -self.b * inv_det;
        let c = -self.c * inv_det;
        let d = self.a * inv_det;
        let e = -(a * self.e + b * self.f);
        let f = -(c * self.e + d * self.f);
        Some(AffineTransform2D { a, b, c, d, e, f })
    }
    /// Contractivity factor (spectral norm estimate of the linear part).
    pub fn contractivity(&self) -> f64 {
        (self.a * self.a + self.b * self.b + self.c * self.c + self.d * self.d).sqrt()
    }
}
/// Classification of an orbit in the Julia/Mandelbrot iteration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OrbitType {
    /// Orbit escapes to infinity.
    Escaping,
    /// Orbit converges to a fixed point.
    FixedPoint,
    /// Orbit is periodic with the given period.
    Periodic(usize),
    /// Orbit is bounded but no periodicity detected.
    Bounded,
}
/// Koch snowflake: three Koch curves forming a closed triangle.
///
/// The Koch snowflake starts from an equilateral triangle and applies the Koch
/// curve construction to each side, producing a fractal with dimension log(4)/log(3).
pub struct KochSnowflake;
impl KochSnowflake {
    /// Theoretical Hausdorff dimension: same as Koch curve, log(4)/log(3).
    pub fn theoretical_dimension() -> f64 {
        (4.0_f64).ln() / (3.0_f64).ln()
    }
    /// Generate the Koch snowflake as line segments by recursively subdividing
    /// each side of an equilateral triangle.
    pub fn generate(generations: usize, side_length: f64) -> Vec<(Point2, Point2)> {
        let h = side_length * (3.0_f64).sqrt() / 2.0;
        let vertices = [
            Point2::new(0.0, 0.0),
            Point2::new(side_length, 0.0),
            Point2::new(side_length / 2.0, h),
        ];
        let mut segments = Vec::new();
        for i in 0..3 {
            let a = vertices[i];
            let b = vertices[(i + 1) % 3];
            Self::subdivide(a, b, generations, &mut segments);
        }
        segments
    }
    /// Recursively subdivide one edge of the snowflake.
    fn subdivide(a: Point2, b: Point2, depth: usize, out: &mut Vec<(Point2, Point2)>) {
        if depth == 0 {
            out.push((a, b));
            return;
        }
        let dx = b.x - a.x;
        let dy = b.y - a.y;
        let p1 = Point2::new(a.x + dx / 3.0, a.y + dy / 3.0);
        let p3 = Point2::new(a.x + 2.0 * dx / 3.0, a.y + 2.0 * dy / 3.0);
        let cos60 = 0.5_f64;
        let sin60 = (3.0_f64).sqrt() / 2.0;
        let mx = dx / 3.0;
        let my = dy / 3.0;
        let apex = Point2::new(
            p1.x + mx * cos60 - my * sin60,
            p1.y + mx * sin60 + my * cos60,
        );
        Self::subdivide(a, p1, depth - 1, out);
        Self::subdivide(p1, apex, depth - 1, out);
        Self::subdivide(apex, p3, depth - 1, out);
        Self::subdivide(p3, b, depth - 1, out);
    }
    /// Count the number of segments after `n` generations.
    pub fn segment_count(n: usize) -> usize {
        3 * 4_usize.pow(n as u32)
    }
    /// Perimeter of the snowflake after `n` generations with initial side `s`.
    pub fn perimeter(n: usize, side_length: f64) -> f64 {
        3.0 * side_length * (4.0 / 3.0_f64).powi(n as i32)
    }
    /// Area of the Koch snowflake after `n` generations.
    pub fn area(n: usize, side_length: f64) -> f64 {
        let a0 = (3.0_f64).sqrt() / 4.0 * side_length * side_length;
        let mut area = a0;
        for k in 1..=n {
            let num_new = 3 * 4_usize.pow((k - 1) as u32);
            let triangle_side = side_length / 3.0_f64.powi(k as i32);
            let triangle_area = (3.0_f64).sqrt() / 4.0 * triangle_side * triangle_side;
            area += num_new as f64 * triangle_area;
        }
        area
    }
}
/// Multifractal analysis tools.
///
/// Multifractal analysis examines the distribution of a measure over a fractal
/// set, computing the spectrum of generalized dimensions D(q) and the
/// singularity spectrum f(alpha).
pub struct Multifractal;
impl Multifractal {
    /// Compute the generalized fractal dimension D(q) for a 2D point set
    /// using the box-counting approach.
    ///
    /// `q` is the moment order (q=0 gives box-counting dim, q=1 info dim, q=2 correlation dim).
    /// Returns `(dimension, r_squared)`.
    pub fn generalized_dimension(points: &[Point2], q: f64, num_scales: usize) -> (f64, f64) {
        if points.is_empty() || num_scales < 2 {
            return (0.0, 0.0);
        }
        let (min_x, max_x, min_y, max_y) = FractalDimension::bounding_box(points);
        let extent = (max_x - min_x).max(max_y - min_y).max(1e-15);
        let n_total = points.len() as f64;
        let mut log_eps = Vec::with_capacity(num_scales);
        let mut log_partition = Vec::with_capacity(num_scales);
        for i in 0..num_scales {
            let epsilon = extent / (2.0_f64).powi(i as i32 + 1);
            if epsilon < 1e-15 {
                break;
            }
            let mut box_counts: HashMap<(i64, i64), usize> = HashMap::new();
            for p in points {
                let ix = ((p.x - min_x) / epsilon).floor() as i64;
                let iy = ((p.y - min_y) / epsilon).floor() as i64;
                *box_counts.entry((ix, iy)).or_insert(0) += 1;
            }
            if (q - 1.0).abs() < 1e-10 {
                let mut entropy = 0.0;
                for &count in box_counts.values() {
                    let pi = count as f64 / n_total;
                    if pi > 0.0 {
                        entropy -= pi * pi.ln();
                    }
                }
                if entropy > 0.0 {
                    log_eps.push(epsilon.ln());
                    log_partition.push(entropy);
                }
            } else {
                let partition: f64 = box_counts
                    .values()
                    .map(|&count| (count as f64 / n_total).powf(q))
                    .sum();
                if partition > 0.0 {
                    log_eps.push(epsilon.ln());
                    log_partition.push(partition.ln());
                }
            }
        }
        if log_eps.len() < 2 {
            return (0.0, 0.0);
        }
        let (slope, r_sq) = FractalDimension::linear_regression(&log_eps, &log_partition);
        if (q - 1.0).abs() < 1e-10 {
            (-slope, r_sq)
        } else {
            (slope / (q - 1.0), r_sq)
        }
    }
    /// Compute the multifractal spectrum D(q) for a range of q values.
    ///
    /// Returns a vector of `(q, D_q)` pairs.
    pub fn spectrum(points: &[Point2], q_values: &[f64], num_scales: usize) -> Vec<(f64, f64)> {
        q_values
            .iter()
            .map(|&q| {
                let (dq, _) = Self::generalized_dimension(points, q, num_scales);
                (q, dq)
            })
            .collect()
    }
    /// Compute the singularity spectrum f(alpha) from the D(q) spectrum.
    ///
    /// Uses the Legendre transform: alpha(q) = d(tau(q))/dq, f(alpha) = q*alpha - tau(q),
    /// where tau(q) = (q-1)*D(q).
    ///
    /// Returns a vector of `(alpha, f_alpha)` pairs.
    pub fn singularity_spectrum(dq_spectrum: &[(f64, f64)]) -> Vec<(f64, f64)> {
        if dq_spectrum.len() < 3 {
            return Vec::new();
        }
        let tau: Vec<(f64, f64)> = dq_spectrum
            .iter()
            .map(|&(q, dq)| (q, (q - 1.0) * dq))
            .collect();
        let mut result = Vec::new();
        for i in 1..(tau.len() - 1) {
            let dq = tau[i + 1].0 - tau[i - 1].0;
            if dq.abs() < 1e-15 {
                continue;
            }
            let alpha = (tau[i + 1].1 - tau[i - 1].1) / dq;
            let f_alpha = tau[i].0 * alpha - tau[i].1;
            result.push((alpha, f_alpha));
        }
        result
    }
}
/// Dragon curve L-system preset.
pub struct DragonCurve;
impl DragonCurve {
    /// Create the L-system for a Dragon curve.
    pub fn lsystem() -> LSystem {
        LSystem::new(
            vec!['F', 'G', '+', '-'],
            vec![
                ProductionRule::new('F', "F+G"),
                ProductionRule::new('G', "F-G"),
            ],
            "F",
        )
    }
    /// Generate Dragon curve segments via turtle graphics.
    pub fn generate(generations: usize, step_size: f64) -> Vec<(Point2, Point2)> {
        let ls = Self::lsystem();
        ls.turtle_interpret(generations, step_size, 90.0)
    }
}
/// Barnsley fern IFS preset.
///
/// The classic Barnsley fern uses 4 affine transforms with specific probabilities
/// to generate a realistic fern shape.
pub struct BarnsleyFern;
impl BarnsleyFern {
    /// Create the IFS for a Barnsley fern.
    pub fn ifs() -> IteratedFunctionSystem {
        let transforms = vec![
            AffineTransform2D::new(0.0, 0.0, 0.0, 0.16, 0.0, 0.0),
            AffineTransform2D::new(0.85, 0.04, -0.04, 0.85, 0.0, 1.6),
            AffineTransform2D::new(0.20, -0.26, 0.23, 0.22, 0.0, 1.6),
            AffineTransform2D::new(-0.15, 0.28, 0.26, 0.24, 0.0, 0.44),
        ];
        let probabilities = vec![0.01, 0.85, 0.07, 0.07];
        IteratedFunctionSystem::new(transforms, probabilities)
    }
    /// Generate Barnsley fern points.
    pub fn generate(num_points: usize) -> Vec<Point2> {
        let ifs = Self::ifs();
        ifs.chaos_game(Point2::new(0.0, 0.0), num_points + 200, 200)
    }
}
/// A 3D point for fractal geometry (no nalgebra).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3 {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
    /// Z coordinate.
    pub z: f64,
}
impl Point3 {
    /// Create a new 3D point.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    /// Euclidean distance to another point.
    pub fn distance(&self, other: &Point3) -> f64 {
        ((self.x - other.x).powi(2) + (self.y - other.y).powi(2) + (self.z - other.z).powi(2))
            .sqrt()
    }
}
/// A deterministic context-free L-system (D0L-system).
///
/// An L-system generates strings by parallel rewriting: at each step,
/// every character in the current string is replaced according to the
/// production rules simultaneously.
#[derive(Debug, Clone)]
pub struct LSystem {
    /// The alphabet of valid symbols.
    pub alphabet: Vec<char>,
    /// Production rules mapping symbols to replacement strings.
    pub rules: HashMap<char, String>,
    /// The starting string (axiom).
    pub axiom: String,
    /// Default turning angle in degrees (for turtle interpretation).
    pub angle_deg: f64,
    /// Default step size (for turtle interpretation).
    pub step_size: f64,
}
/// State of an L-system after evolution.
#[derive(Debug, Clone)]
pub struct LSystemEvolution {
    /// The generated string.
    pub string: String,
    /// The generation number.
    pub generation: usize,
}
impl LSystem {
    /// Create a new L-system.
    pub fn new(alphabet: Vec<char>, rules: Vec<ProductionRule>, axiom: impl Into<String>) -> Self {
        let mut rule_map = HashMap::new();
        for rule in rules {
            rule_map.insert(rule.predecessor, rule.successor);
        }
        Self {
            alphabet,
            rules: rule_map,
            axiom: axiom.into(),
            angle_deg: 90.0,
            step_size: 1.0,
        }
    }
    /// Koch curve L-system (angle = 60°).
    pub fn koch_curve() -> Self {
        let mut ls = Self::new(vec!['F'], vec![ProductionRule::new('F', "F+F-F-F+F")], "F");
        ls.angle_deg = 90.0;
        ls
    }
    /// Sierpinski triangle L-system (angle = 120°).
    pub fn sierpinski() -> Self {
        let mut ls = Self::new(
            vec!['A', 'B'],
            vec![
                ProductionRule::new('A', "B-A-B"),
                ProductionRule::new('B', "A+B+A"),
            ],
            "A",
        );
        ls.angle_deg = 60.0;
        ls
    }
    /// Dragon curve L-system (angle = 90°).
    pub fn dragon_curve() -> Self {
        let mut ls = Self::new(
            vec!['X', 'Y'],
            vec![
                ProductionRule::new('X', "X+YF+"),
                ProductionRule::new('Y', "-FX-Y"),
            ],
            "FX",
        );
        ls.angle_deg = 90.0;
        ls
    }
    /// Plant branching L-system (angle = 25°).
    pub fn plant_branching() -> Self {
        let mut ls = Self::new(
            vec!['X', 'F'],
            vec![
                ProductionRule::new('X', "F+[[X]-X]-F[-FX]+X"),
                ProductionRule::new('F', "FF"),
            ],
            "X",
        );
        ls.angle_deg = 25.0;
        ls
    }
    /// Evolve the L-system for `n` generations, returning the state.
    pub fn evolve(&self, n: usize) -> LSystemEvolution {
        LSystemEvolution {
            string: self.iterate(n),
            generation: n,
        }
    }
    /// Interpret a previously evolved state as turtle graphics segments.
    pub fn interpret(&self, state: &LSystemEvolution) -> Vec<(Point2, Point2)> {
        let angle_rad = self.angle_deg * std::f64::consts::PI / 180.0;
        let step = self.step_size;
        let mut x = 0.0_f64;
        let mut y = 0.0_f64;
        let mut heading = 0.0_f64;
        let mut stack: Vec<(f64, f64, f64)> = Vec::new();
        let mut segments = Vec::new();
        for ch in state.string.chars() {
            match ch {
                'F' | 'G' | 'A' | 'B' => {
                    let nx = x + step * heading.cos();
                    let ny = y + step * heading.sin();
                    segments.push((Point2::new(x, y), Point2::new(nx, ny)));
                    x = nx;
                    y = ny;
                }
                '+' => heading += angle_rad,
                '-' => heading -= angle_rad,
                '[' => stack.push((x, y, heading)),
                ']' => {
                    if let Some((sx, sy, sh)) = stack.pop() {
                        x = sx;
                        y = sy;
                        heading = sh;
                    }
                }
                _ => {}
            }
        }
        segments
    }
    /// Generate the string after `n` iterations.
    pub fn iterate(&self, n: usize) -> String {
        let mut current = self.axiom.clone();
        for _ in 0..n {
            let mut next = String::with_capacity(current.len() * 2);
            for ch in current.chars() {
                if let Some(replacement) = self.rules.get(&ch) {
                    next.push_str(replacement);
                } else {
                    next.push(ch);
                }
            }
            current = next;
        }
        current
    }
    /// Count the length of the string after `n` iterations without building it.
    pub fn iterate_length(&self, n: usize) -> usize {
        let mut lengths: HashMap<char, usize> = HashMap::new();
        for &ch in &self.alphabet {
            lengths.insert(ch, 1);
        }
        for _ in 0..n {
            let mut new_lengths: HashMap<char, usize> = HashMap::new();
            for &ch in &self.alphabet {
                if let Some(replacement) = self.rules.get(&ch) {
                    let len: usize = replacement
                        .chars()
                        .map(|c| lengths.get(&c).copied().unwrap_or(1))
                        .sum();
                    new_lengths.insert(ch, len);
                } else {
                    new_lengths.insert(ch, *lengths.get(&ch).unwrap_or(&1));
                }
            }
            lengths = new_lengths;
        }
        self.axiom
            .chars()
            .map(|c| lengths.get(&c).copied().unwrap_or(1))
            .sum()
    }
    /// Interpret the generated string as turtle graphics commands.
    ///
    /// - `F`, `G`: move forward by `step_size`
    /// - `+`: turn left by `angle_deg` degrees
    /// - `-`: turn right by `angle_deg` degrees
    /// - `[`: push state
    /// - `]`: pop state
    ///
    /// Returns a list of line segments as (start, end) pairs.
    pub fn turtle_interpret(
        &self,
        n: usize,
        step_size: f64,
        angle_deg: f64,
    ) -> Vec<(Point2, Point2)> {
        let string = self.iterate(n);
        let angle_rad = angle_deg * std::f64::consts::PI / 180.0;
        let mut x = 0.0_f64;
        let mut y = 0.0_f64;
        let mut heading = 0.0_f64;
        let mut stack: Vec<(f64, f64, f64)> = Vec::new();
        let mut segments = Vec::new();
        for ch in string.chars() {
            match ch {
                'F' | 'G' => {
                    let nx = x + step_size * heading.cos();
                    let ny = y + step_size * heading.sin();
                    segments.push((Point2::new(x, y), Point2::new(nx, ny)));
                    x = nx;
                    y = ny;
                }
                '+' => {
                    heading += angle_rad;
                }
                '-' => {
                    heading -= angle_rad;
                }
                '[' => {
                    stack.push((x, y, heading));
                }
                ']' => {
                    if let Some((sx, sy, sh)) = stack.pop() {
                        x = sx;
                        y = sy;
                        heading = sh;
                    }
                }
                _ => {}
            }
        }
        segments
    }
    /// Add a production rule.
    pub fn add_rule(&mut self, predecessor: char, successor: impl Into<String>) {
        let succ: String = successor.into();
        self.rules.insert(predecessor, succ);
        if !self.alphabet.contains(&predecessor) {
            self.alphabet.push(predecessor);
        }
    }
}
/// A Julia set for the quadratic map z_{n+1} = z_n^2 + c with fixed c.
///
/// Points z_0 in the complex plane are classified by orbit behavior.
pub struct JuliaSet {
    /// Real part of the constant c.
    pub cr: f64,
    /// Imaginary part of the constant c.
    pub ci: f64,
    /// Maximum iterations.
    pub max_iter: u32,
    /// Escape radius squared.
    pub escape_radius_sq: f64,
}
impl JuliaSet {
    /// Create a new Julia set with given c parameter.
    pub fn new(cr: f64, ci: f64, max_iter: u32, escape_radius: f64) -> Self {
        Self {
            cr,
            ci,
            max_iter,
            escape_radius_sq: escape_radius * escape_radius,
        }
    }
    /// Classic Julia set c = -0.7 + 0.27015i.
    pub fn classic() -> Self {
        Self::new(-0.7, 0.27015, 256, 2.0)
    }
    /// Dendrite Julia set c = 0 + i.
    pub fn dendrite() -> Self {
        Self::new(0.0, 1.0, 256, 2.0)
    }
    /// Compute escape time for an initial point z0 = (zr, zi).
    pub fn escape_time(&self, zr: f64, zi: f64) -> Option<u32> {
        let mut zr = zr;
        let mut zi = zi;
        for i in 0..self.max_iter {
            let zr2 = zr * zr;
            let zi2 = zi * zi;
            if zr2 + zi2 > self.escape_radius_sq {
                return Some(i);
            }
            let new_zi = 2.0 * zr * zi + self.ci;
            zr = zr2 - zi2 + self.cr;
            zi = new_zi;
        }
        None
    }
    /// Check if a point is in the Julia set (orbit does not escape).
    pub fn is_in_set(&self, zr: f64, zi: f64) -> bool {
        self.escape_time(zr, zi).is_none()
    }
    /// Classify the orbit of z0: returns `OrbitType`.
    pub fn classify_orbit(&self, zr: f64, zi: f64) -> OrbitType {
        let mut zr = zr;
        let mut zi = zi;
        let mut history_r = Vec::with_capacity(self.max_iter as usize);
        let mut history_i = Vec::with_capacity(self.max_iter as usize);
        for _ in 0..self.max_iter {
            let zr2 = zr * zr;
            let zi2 = zi * zi;
            if zr2 + zi2 > self.escape_radius_sq {
                return OrbitType::Escaping;
            }
            history_r.push(zr);
            history_i.push(zi);
            let new_zi = 2.0 * zr * zi + self.ci;
            zr = zr2 - zi2 + self.cr;
            zi = new_zi;
        }
        let n = history_r.len();
        if n < 2 {
            return OrbitType::Bounded;
        }
        let last_r = history_r[n - 1];
        let last_i = history_i[n - 1];
        for period in 1..=16.min(n - 1) {
            let prev_r = history_r[n - 1 - period];
            let prev_i = history_i[n - 1 - period];
            let dr = (last_r - prev_r).abs();
            let di = (last_i - prev_i).abs();
            if dr < 1e-10 && di < 1e-10 {
                if period == 1 {
                    return OrbitType::FixedPoint;
                }
                return OrbitType::Periodic(period);
            }
        }
        OrbitType::Bounded
    }
    /// Smooth escape time for coloring.
    pub fn smooth_escape_time(&self, zr: f64, zi: f64) -> f64 {
        let mut zr = zr;
        let mut zi = zi;
        for i in 0..self.max_iter {
            let zr2 = zr * zr;
            let zi2 = zi * zi;
            if zr2 + zi2 > self.escape_radius_sq {
                let modulus = (zr2 + zi2).sqrt();
                let mu = i as f64 - modulus.ln().ln() / (2.0_f64).ln();
                return mu;
            }
            let new_zi = 2.0 * zr * zi + self.ci;
            zr = zr2 - zi2 + self.cr;
            zi = new_zi;
        }
        self.max_iter as f64
    }
    /// Render a grid of escape times.
    pub fn render_grid(
        &self,
        x_min: f64,
        x_max: f64,
        y_min: f64,
        y_max: f64,
        width: usize,
        height: usize,
    ) -> Vec<Vec<Option<u32>>> {
        let dx = (x_max - x_min) / width as f64;
        let dy = (y_max - y_min) / height as f64;
        (0..height)
            .map(|j| {
                let zi = y_min + (j as f64 + 0.5) * dy;
                (0..width)
                    .map(|i| {
                        let zr = x_min + (i as f64 + 0.5) * dx;
                        self.escape_time(zr, zi)
                    })
                    .collect()
            })
            .collect()
    }
}
/// An axis-aligned box in 3D, used for Menger sponge cubes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb3 {
    /// Minimum corner.
    pub min: Point3,
    /// Maximum corner.
    pub max: Point3,
}
impl Aabb3 {
    /// Create a new AABB.
    pub fn new(min: Point3, max: Point3) -> Self {
        Self { min, max }
    }
    /// Side length (assumes cubic).
    pub fn side(&self) -> f64 {
        self.max.x - self.min.x
    }
    /// Volume.
    pub fn volume(&self) -> f64 {
        (self.max.x - self.min.x) * (self.max.y - self.min.y) * (self.max.z - self.min.z)
    }
    /// Center point.
    pub fn center(&self) -> Point3 {
        Point3::new(
            (self.min.x + self.max.x) / 2.0,
            (self.min.y + self.max.y) / 2.0,
            (self.min.z + self.max.z) / 2.0,
        )
    }
}
/// Morton (Z-order) curve utilities.
///
/// Maps 2D coordinates to a 1D index using bit-interleaving, producing
/// a space-filling Z-shaped traversal pattern.
pub struct MortonCurve;
impl MortonCurve {
    /// Encode 2D coordinates into a Morton code.
    ///
    /// Interleaves the bits of `x` and `y` (each up to 16 bits).
    pub fn encode(x: u32, y: u32) -> u64 {
        Self::spread_bits(x as u64) | (Self::spread_bits(y as u64) << 1)
    }
    /// Decode a Morton code back to 2D coordinates.
    pub fn decode(code: u64) -> (u32, u32) {
        let x = Self::compact_bits(code) as u32;
        let y = Self::compact_bits(code >> 1) as u32;
        (x, y)
    }
    /// Spread bits: insert a zero between each bit.
    fn spread_bits(mut v: u64) -> u64 {
        v &= 0x0000_0000_FFFF_FFFF;
        v = (v | (v << 16)) & 0x0000_FFFF_0000_FFFF;
        v = (v | (v << 8)) & 0x00FF_00FF_00FF_00FF;
        v = (v | (v << 4)) & 0x0F0F_0F0F_0F0F_0F0F;
        v = (v | (v << 2)) & 0x3333_3333_3333_3333;
        v = (v | (v << 1)) & 0x5555_5555_5555_5555;
        v
    }
    /// Compact bits: remove every other bit.
    fn compact_bits(mut v: u64) -> u64 {
        v &= 0x5555_5555_5555_5555;
        v = (v | (v >> 1)) & 0x3333_3333_3333_3333;
        v = (v | (v >> 2)) & 0x0F0F_0F0F_0F0F_0F0F;
        v = (v | (v >> 4)) & 0x00FF_00FF_00FF_00FF;
        v = (v | (v >> 8)) & 0x0000_FFFF_0000_FFFF;
        v = (v | (v >> 16)) & 0x0000_0000_FFFF_FFFF;
        v
    }
    /// Generate the Z-order traversal as a list of 2D points for an `n x n` grid.
    ///
    /// `n` should be a power of 2.
    pub fn generate_path(n: u32) -> Vec<(u32, u32)> {
        let total = n * n;
        let mut pairs: Vec<(u64, u32, u32)> = Vec::with_capacity(total as usize);
        for y in 0..n {
            for x in 0..n {
                pairs.push((Self::encode(x, y), x, y));
            }
        }
        pairs.sort_by_key(|&(code, _, _)| code);
        pairs.iter().map(|&(_, x, y)| (x, y)).collect()
    }
    /// 3D Morton code encoding.
    pub fn encode_3d(x: u32, y: u32, z: u32) -> u64 {
        Self::spread_bits_3d(x as u64)
            | (Self::spread_bits_3d(y as u64) << 1)
            | (Self::spread_bits_3d(z as u64) << 2)
    }
    /// Spread bits for 3D: insert two zeros between each bit.
    fn spread_bits_3d(mut v: u64) -> u64 {
        v &= 0x001F_FFFF;
        v = (v | (v << 32)) & 0x001F_0000_0000_FFFF;
        v = (v | (v << 16)) & 0x001F_0000_FF00_00FF;
        v = (v | (v << 8)) & 0x100F_00F0_0F00_F00F;
        v = (v | (v << 4)) & 0x10C3_0C30_C30C_30C3;
        v = (v | (v << 2)) & 0x1249_2492_4924_9249;
        v
    }
}
/// A 2D point used in fractal geometry computations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2 {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
}
impl Point2 {
    /// Create a new 2D point.
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
    /// Euclidean distance to another point.
    pub fn distance(&self, other: &Point2) -> f64 {
        ((self.x - other.x).powi(2) + (self.y - other.y).powi(2)).sqrt()
    }
    /// Squared distance to another point.
    pub fn distance_sq(&self, other: &Point2) -> f64 {
        (self.x - other.x).powi(2) + (self.y - other.y).powi(2)
    }
}
/// Menger sponge generator.
///
/// The Menger sponge is a 3D fractal formed by recursively removing the center
/// and face-center sub-cubes from a cube. Its Hausdorff dimension is
/// log(20)/log(3) ≈ 2.727.
pub struct MengerSponge;
impl MengerSponge {
    /// Theoretical Hausdorff dimension: log(20)/log(3).
    pub fn theoretical_dimension() -> f64 {
        (20.0_f64).ln() / (3.0_f64).ln()
    }
    /// Generate the Menger sponge as a list of cubes (AABBs) at the given recursion level.
    ///
    /// `level` 0 returns a single unit cube, level 1 returns 20 sub-cubes, etc.
    pub fn generate(level: usize, origin: Point3, size: f64) -> Vec<Aabb3> {
        if level == 0 {
            return vec![Aabb3::new(
                origin,
                Point3::new(origin.x + size, origin.y + size, origin.z + size),
            )];
        }
        let sub = size / 3.0;
        let mut cubes = Vec::new();
        for ix in 0..3_u8 {
            for iy in 0..3_u8 {
                for iz in 0..3_u8 {
                    let center_count = u8::from(ix == 1) + u8::from(iy == 1) + u8::from(iz == 1);
                    if center_count >= 2 {
                        continue;
                    }
                    let o = Point3::new(
                        origin.x + ix as f64 * sub,
                        origin.y + iy as f64 * sub,
                        origin.z + iz as f64 * sub,
                    );
                    cubes.extend(Self::generate(level - 1, o, sub));
                }
            }
        }
        cubes
    }
    /// Number of cubes at a given level.
    pub fn cube_count(level: usize) -> usize {
        20_usize.pow(level as u32)
    }
    /// Total volume at a given level (unit cube start).
    pub fn total_volume(level: usize) -> f64 {
        (20.0 / 27.0_f64).powi(level as i32)
    }
}
/// Koch curve L-system preset.
///
/// The Koch curve replaces each line segment with 4 segments arranged
/// in a triangular bump, creating a fractal curve of dimension log(4)/log(3).
pub struct KochCurve;
impl KochCurve {
    /// Create the L-system for a Koch curve.
    pub fn lsystem() -> LSystem {
        LSystem::new(
            vec!['F', '+', '-'],
            vec![ProductionRule::new('F', "F+F--F+F")],
            "F",
        )
    }
    /// Theoretical dimension: log(4)/log(3) ≈ 1.2619.
    pub fn theoretical_dimension() -> f64 {
        (4.0_f64).ln() / (3.0_f64).ln()
    }
    /// Generate Koch curve segments via turtle graphics.
    pub fn generate(generations: usize, step_size: f64) -> Vec<(Point2, Point2)> {
        let ls = Self::lsystem();
        ls.turtle_interpret(generations, step_size, 60.0)
    }
    /// Compute the number of segments after `n` generations.
    pub fn segment_count(n: usize) -> usize {
        4_usize.pow(n as u32)
    }
}
