//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
/// Parametric surface types: tensor product, Coons patch.
#[derive(Debug, Clone)]
pub struct SurfaceParametric {
    /// Control net: `rows × cols` array of 3D control points
    pub control_net: Vec<Vec<[f64; 3]>>,
}
impl SurfaceParametric {
    /// Create a new parametric surface from a control net.
    pub fn new(control_net: Vec<Vec<[f64; 3]>>) -> Self {
        Self { control_net }
    }
    /// Evaluate a tensor-product Bézier surface at parameters (u, v).
    pub fn tensor_bezier(&self, u: f64, v: f64) -> [f64; 3] {
        let rows = self.control_net.len();
        if rows == 0 {
            return [0.0, 0.0, 0.0];
        }
        let row_pts: Vec<[f64; 3]> = self
            .control_net
            .iter()
            .map(|row| {
                let curve = CurveParametric::new(row.clone());
                curve.bezier_de_casteljau(u)
            })
            .collect();
        let col_curve = CurveParametric::new(row_pts);
        col_curve.bezier_de_casteljau(v)
    }
    /// Evaluate a bilinear Coons patch at parameters (u, v).
    ///
    /// Requires exactly a 2×2 control net (four corner points).
    pub fn coons_patch(&self, u: f64, v: f64) -> [f64; 3] {
        if self.control_net.len() < 2 || self.control_net[0].len() < 2 {
            return [0.0, 0.0, 0.0];
        }
        let p00 = self.control_net[0][0];
        let p01 = self.control_net[0][1];
        let p10 = self.control_net[1][0];
        let p11 = self.control_net[1][1];
        let w00 = (1.0 - u) * (1.0 - v);
        let w01 = (1.0 - u) * v;
        let w10 = u * (1.0 - v);
        let w11 = u * v;
        [
            w00 * p00[0] + w01 * p01[0] + w10 * p10[0] + w11 * p11[0],
            w00 * p00[1] + w01 * p01[1] + w10 * p10[1] + w11 * p11[1],
            w00 * p00[2] + w01 * p01[2] + w10 * p10[2] + w11 * p11[2],
        ]
    }
    /// Compute the approximate surface normal at (u, v) using finite differences.
    pub fn normal(&self, u: f64, v: f64) -> [f64; 3] {
        let eps = 1e-4;
        let pu = self.tensor_bezier((u + eps).min(1.0), v);
        let pu_neg = self.tensor_bezier((u - eps).max(0.0), v);
        let pv = self.tensor_bezier(u, (v + eps).min(1.0));
        let pv_neg = self.tensor_bezier(u, (v - eps).max(0.0));
        let du = [
            (pu[0] - pu_neg[0]) / (2.0 * eps),
            (pu[1] - pu_neg[1]) / (2.0 * eps),
            (pu[2] - pu_neg[2]) / (2.0 * eps),
        ];
        let dv = [
            (pv[0] - pv_neg[0]) / (2.0 * eps),
            (pv[1] - pv_neg[1]) / (2.0 * eps),
            (pv[2] - pv_neg[2]) / (2.0 * eps),
        ];
        let n = [
            du[1] * dv[2] - du[2] * dv[1],
            du[2] * dv[0] - du[0] * dv[2],
            du[0] * dv[1] - du[1] * dv[0],
        ];
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if len < 1e-12 {
            [0.0, 0.0, 1.0]
        } else {
            [n[0] / len, n[1] / len, n[2] / len]
        }
    }
    /// Sample the surface at `nu × nv` grid points.
    pub fn sample(&self, nu: usize, nv: usize) -> Vec<Vec<[f64; 3]>> {
        (0..nu)
            .map(|i| {
                let u = i as f64 / (nu - 1).max(1) as f64;
                (0..nv)
                    .map(|j| {
                        let v = j as f64 / (nv - 1).max(1) as f64;
                        self.tensor_bezier(u, v)
                    })
                    .collect()
            })
            .collect()
    }
}
/// Turtle graphics state for L-system interpretation.
#[derive(Debug, Clone, Copy)]
pub struct TurtleState {
    /// Current x position
    pub x: f64,
    /// Current y position
    pub y: f64,
    /// Current heading angle in radians
    pub angle: f64,
}
impl TurtleState {
    /// Create a new TurtleState at the origin pointing right.
    pub fn new() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            angle: 0.0,
        }
    }
}
/// 2D Perlin noise generator with permutation table.
#[derive(Debug, Clone)]
pub struct PerlinNoise2d {
    /// Permutation table (512 entries)
    pub(super) perm: Vec<u8>,
}
impl PerlinNoise2d {
    /// Create a new Perlin noise generator with a given seed.
    pub fn new(seed: u64) -> Self {
        let mut perm = (0u8..=255).collect::<Vec<_>>();
        let mut s = seed;
        for i in (1..256).rev() {
            s = s
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let j = (s >> 33) as usize % (i + 1);
            perm.swap(i, j);
        }
        let mut full = perm.clone();
        full.extend_from_slice(&perm);
        Self { perm: full }
    }
    fn fade(t: f64) -> f64 {
        t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
    }
    fn lerp(a: f64, b: f64, t: f64) -> f64 {
        a + t * (b - a)
    }
    fn grad2(hash: u8, x: f64, y: f64) -> f64 {
        match hash & 3 {
            0 => x + y,
            1 => -x + y,
            2 => x - y,
            _ => -x - y,
        }
    }
    /// Evaluate 2D Perlin noise at (x, y). Returns value in approximately \[-1, 1\].
    pub fn noise(&self, x: f64, y: f64) -> f64 {
        let xi = x.floor() as i32;
        let yi = y.floor() as i32;
        let xf = x - xi as f64;
        let yf = y - yi as f64;
        let u = Self::fade(xf);
        let v = Self::fade(yf);
        let aa = self.perm[(self.perm[(xi & 255) as usize] as i32 + (yi & 255)) as usize & 255];
        let ab =
            self.perm[(self.perm[(xi & 255) as usize] as i32 + ((yi + 1) & 255)) as usize & 255];
        let ba =
            self.perm[(self.perm[((xi + 1) & 255) as usize] as i32 + (yi & 255)) as usize & 255];
        let bb = self.perm
            [(self.perm[((xi + 1) & 255) as usize] as i32 + ((yi + 1) & 255)) as usize & 255];
        let x1 = Self::lerp(Self::grad2(aa, xf, yf), Self::grad2(ba, xf - 1.0, yf), u);
        let x2 = Self::lerp(
            Self::grad2(ab, xf, yf - 1.0),
            Self::grad2(bb, xf - 1.0, yf - 1.0),
            u,
        );
        Self::lerp(x1, x2, v)
    }
    /// Fractional Brownian Motion (fBm): sum of octaves of Perlin noise.
    ///
    /// Higher octave count produces more detail (roughness).
    pub fn fbm(&self, x: f64, y: f64, octaves: u32, lacunarity: f64, gain: f64) -> f64 {
        let mut value = 0.0;
        let mut amplitude = 1.0;
        let mut frequency = 1.0;
        let mut max_val = 0.0;
        for _ in 0..octaves {
            value += amplitude * self.noise(x * frequency, y * frequency);
            max_val += amplitude;
            amplitude *= gain;
            frequency *= lacunarity;
        }
        if max_val > 0.0 { value / max_val } else { 0.0 }
    }
    /// Turbulence: sum of absolute values of noise octaves.
    pub fn turbulence(&self, x: f64, y: f64, octaves: u32, lacunarity: f64, gain: f64) -> f64 {
        let mut value = 0.0;
        let mut amplitude = 1.0;
        let mut frequency = 1.0;
        for _ in 0..octaves {
            value += amplitude * self.noise(x * frequency, y * frequency).abs();
            amplitude *= gain;
            frequency *= lacunarity;
        }
        value
    }
}
/// Worley (cellular) noise generator.
#[derive(Debug, Clone)]
pub struct WorleyNoise {
    /// Feature points (2D)
    pub(super) points: Vec<[f64; 2]>,
}
impl WorleyNoise {
    /// Create Worley noise with randomly distributed feature points in \[0, 1\]^2.
    pub fn new(num_points: usize, seed: u64) -> Self {
        let mut s = seed;
        let mut points = Vec::with_capacity(num_points);
        for _ in 0..num_points {
            s = s
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let x = (s >> 11) as f64 / (1u64 << 53) as f64;
            s = s
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let y = (s >> 11) as f64 / (1u64 << 53) as f64;
            points.push([x, y]);
        }
        Self { points }
    }
    /// Evaluate F1 (distance to nearest feature point) at (x, y).
    pub fn f1(&self, x: f64, y: f64) -> f64 {
        self.points
            .iter()
            .map(|&[px, py]| {
                let dx = x - px;
                let dy = y - py;
                (dx * dx + dy * dy).sqrt()
            })
            .fold(f64::MAX, f64::min)
    }
    /// Evaluate F2 - F1 (distance to second nearest minus nearest).
    pub fn f2_minus_f1(&self, x: f64, y: f64) -> f64 {
        let mut d1 = f64::MAX;
        let mut d2 = f64::MAX;
        for &[px, py] in &self.points {
            let d = ((x - px).powi(2) + (y - py).powi(2)).sqrt();
            if d < d1 {
                d2 = d1;
                d1 = d;
            } else if d < d2 {
                d2 = d;
            }
        }
        d2 - d1
    }
}
/// Lindenmayer system (L-system) for generating fractal geometry strings.
#[derive(Debug, Clone)]
pub struct LSystem {
    /// Axiom (initial string)
    pub axiom: String,
    /// Production rules
    pub rules: Vec<LSystemRule>,
    /// Angle increment in radians for turtle interpretation
    pub angle_increment: f64,
    /// Step size for turtle movement
    pub step_size: f64,
}
impl LSystem {
    /// Create a new L-system.
    pub fn new(
        axiom: &str,
        rules: Vec<LSystemRule>,
        angle_increment_deg: f64,
        step_size: f64,
    ) -> Self {
        Self {
            axiom: axiom.to_string(),
            rules,
            angle_increment: angle_increment_deg.to_radians(),
            step_size,
        }
    }
    /// Apply the production rules for `n` generations.
    pub fn evolve(&self, generations: u32) -> LSystemState {
        let mut current = self.axiom.clone();
        for _ in 0..generations {
            let mut next = String::new();
            for ch in current.chars() {
                if let Some(rule) = self.rules.iter().find(|r| r.predecessor == ch) {
                    next.push_str(&rule.successor);
                } else {
                    next.push(ch);
                }
            }
            current = next;
        }
        LSystemState {
            string: current,
            generation: generations,
        }
    }
    /// Interpret the L-system string as turtle graphics commands.
    ///
    /// Returns a vector of line segments as pairs of 2D points.
    /// Symbols: `F` = move forward (draw), `f` = move forward (no draw),
    /// `+` = turn left, `-` = turn right, `[` = push state, `]` = pop state.
    pub fn interpret(&self, state: &LSystemState) -> Vec<([f64; 2], [f64; 2])> {
        let mut segments = Vec::new();
        let mut turtle = TurtleState::new();
        let mut stack: Vec<TurtleState> = Vec::new();
        for ch in state.string.chars() {
            match ch {
                'F' | 'G' => {
                    let nx = turtle.x + self.step_size * turtle.angle.cos();
                    let ny = turtle.y + self.step_size * turtle.angle.sin();
                    segments.push(([turtle.x, turtle.y], [nx, ny]));
                    turtle.x = nx;
                    turtle.y = ny;
                }
                'f' => {
                    turtle.x += self.step_size * turtle.angle.cos();
                    turtle.y += self.step_size * turtle.angle.sin();
                }
                '+' => turtle.angle += self.angle_increment,
                '-' => turtle.angle -= self.angle_increment,
                '[' => stack.push(turtle),
                ']' => {
                    if let Some(s) = stack.pop() {
                        turtle = s;
                    }
                }
                _ => {}
            }
        }
        segments
    }
    /// Create an L-system for the Koch snowflake.
    pub fn koch_curve() -> Self {
        Self::new(
            "F--F--F",
            vec![LSystemRule::new('F', "F+F--F+F")],
            60.0,
            1.0,
        )
    }
    /// Create an L-system for the Sierpinski triangle.
    pub fn sierpinski() -> Self {
        Self::new(
            "F-G-G",
            vec![
                LSystemRule::new('F', "F-G+F+G-F"),
                LSystemRule::new('G', "GG"),
            ],
            120.0,
            1.0,
        )
    }
    /// Create an L-system for the dragon curve.
    pub fn dragon_curve() -> Self {
        Self::new(
            "FX",
            vec![
                LSystemRule::new('X', "X+YF+"),
                LSystemRule::new('Y', "-FX-Y"),
            ],
            90.0,
            1.0,
        )
    }
    /// Create an L-system for plant branching.
    pub fn plant_branching() -> Self {
        Self::new(
            "X",
            vec![
                LSystemRule::new('X', "F+[[X]-X]-F[-FX]+X"),
                LSystemRule::new('F', "FF"),
            ],
            25.0,
            1.0,
        )
    }
}
/// A single affine transformation for an IFS.
#[derive(Debug, Clone, Copy)]
pub struct IfsTransform {
    /// 2×2 linear part: a coefficient
    pub a: f64,
    /// 2×2 linear part: b coefficient
    pub b: f64,
    /// 2×2 linear part: c coefficient
    pub c: f64,
    /// 2×2 linear part: d coefficient
    pub d: f64,
    /// Translation x
    pub e: f64,
    /// Translation y
    pub f: f64,
    /// Probability weight
    pub prob: f64,
}
impl IfsTransform {
    /// Create a new IFS affine transform.
    pub fn new(a: f64, b: f64, c: f64, d: f64, e: f64, f_val: f64, prob: f64) -> Self {
        Self {
            a,
            b,
            c,
            d,
            e,
            f: f_val,
            prob,
        }
    }
    fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        (
            self.a * x + self.b * y + self.e,
            self.c * x + self.d * y + self.f,
        )
    }
}
/// Reaction-diffusion simulation using the Gray-Scott model.
///
/// The Gray-Scott model produces Turing patterns by simulating two chemical
/// species U and V with different diffusion rates and reaction kinetics.
#[derive(Debug, Clone)]
pub struct ReactionDiffusion {
    /// Grid width
    pub width: usize,
    /// Grid height
    pub height: usize,
    /// Concentration of species U
    pub u: Vec<f64>,
    /// Concentration of species V
    pub v: Vec<f64>,
    /// Diffusion rate of U
    pub d_u: f64,
    /// Diffusion rate of V
    pub d_v: f64,
    /// Feed rate
    pub feed: f64,
    /// Kill rate
    pub kill: f64,
}
impl ReactionDiffusion {
    /// Create a new Gray-Scott reaction-diffusion grid.
    ///
    /// # Parameters
    /// - `width`, `height`: grid dimensions
    /// - `d_u`: diffusion rate of U (typically 0.2)
    /// - `d_v`: diffusion rate of V (typically 0.1)
    /// - `feed`: feed rate (controls pattern type)
    /// - `kill`: kill rate (controls pattern type)
    pub fn new(width: usize, height: usize, d_u: f64, d_v: f64, feed: f64, kill: f64) -> Self {
        let size = width * height;
        let u = vec![1.0; size];
        let v = vec![0.0; size];
        Self {
            width,
            height,
            u,
            v,
            d_u,
            d_v,
            feed,
            kill,
        }
    }
    /// Seed with a square perturbation at the center.
    pub fn seed_center(&mut self) {
        let cx = self.width / 2;
        let cy = self.height / 2;
        let r = 5;
        for y in cy.saturating_sub(r)..=(cy + r).min(self.height - 1) {
            for x in cx.saturating_sub(r)..=(cx + r).min(self.width - 1) {
                let idx = y * self.width + x;
                self.u[idx] = 0.5;
                self.v[idx] = 0.25;
            }
        }
    }
    fn idx(&self, x: usize, y: usize) -> usize {
        y * self.width + x
    }
    fn laplacian_u(&self, x: usize, y: usize) -> f64 {
        let i = self.idx(x, y);
        let left = if x > 0 {
            self.u[self.idx(x - 1, y)]
        } else {
            self.u[self.idx(self.width - 1, y)]
        };
        let right = if x + 1 < self.width {
            self.u[self.idx(x + 1, y)]
        } else {
            self.u[self.idx(0, y)]
        };
        let up = if y > 0 {
            self.u[self.idx(x, y - 1)]
        } else {
            self.u[self.idx(x, self.height - 1)]
        };
        let down = if y + 1 < self.height {
            self.u[self.idx(x, y + 1)]
        } else {
            self.u[self.idx(x, 0)]
        };
        left + right + up + down - 4.0 * self.u[i]
    }
    fn laplacian_v(&self, x: usize, y: usize) -> f64 {
        let i = self.idx(x, y);
        let left = if x > 0 {
            self.v[self.idx(x - 1, y)]
        } else {
            self.v[self.idx(self.width - 1, y)]
        };
        let right = if x + 1 < self.width {
            self.v[self.idx(x + 1, y)]
        } else {
            self.v[self.idx(0, y)]
        };
        let up = if y > 0 {
            self.v[self.idx(x, y - 1)]
        } else {
            self.v[self.idx(x, self.height - 1)]
        };
        let down = if y + 1 < self.height {
            self.v[self.idx(x, y + 1)]
        } else {
            self.v[self.idx(x, 0)]
        };
        left + right + up + down - 4.0 * self.v[i]
    }
    /// Advance the simulation by `steps` time steps with given `dt`.
    pub fn step(&mut self, steps: u32, dt: f64) {
        for _ in 0..steps {
            let mut du = vec![0.0f64; self.width * self.height];
            let mut dv = vec![0.0f64; self.width * self.height];
            for y in 0..self.height {
                for x in 0..self.width {
                    let i = self.idx(x, y);
                    let u = self.u[i];
                    let v = self.v[i];
                    let uvv = u * v * v;
                    du[i] = self.d_u * self.laplacian_u(x, y) - uvv + self.feed * (1.0 - u);
                    dv[i] = self.d_v * self.laplacian_v(x, y) + uvv - (self.feed + self.kill) * v;
                }
            }
            for i in 0..self.u.len() {
                self.u[i] = (self.u[i] + dt * du[i]).clamp(0.0, 1.0);
                self.v[i] = (self.v[i] + dt * dv[i]).clamp(0.0, 1.0);
            }
        }
    }
}
/// Iterated Function System for generating fractal attractors (e.g., Barnsley fern).
#[derive(Debug, Clone)]
pub struct IteratedFunctionSystem {
    /// The set of affine transformations
    pub transforms: Vec<IfsTransform>,
}
impl IteratedFunctionSystem {
    /// Create a new IFS from a set of transforms.
    pub fn new(transforms: Vec<IfsTransform>) -> Self {
        Self { transforms }
    }
    /// Create the classic Barnsley fern IFS.
    pub fn barnsley_fern() -> Self {
        Self::new(vec![
            IfsTransform::new(0.0, 0.0, 0.0, 0.16, 0.0, 0.0, 0.01),
            IfsTransform::new(0.85, 0.04, -0.04, 0.85, 0.0, 1.6, 0.85),
            IfsTransform::new(0.2, -0.26, 0.23, 0.22, 0.0, 1.6, 0.07),
            IfsTransform::new(-0.15, 0.28, 0.26, 0.24, 0.0, 0.44, 0.07),
        ])
    }
    /// Generate `n` points of the IFS attractor using the chaos game.
    pub fn generate(&self, n: usize, seed: u64) -> Vec<[f64; 2]> {
        if self.transforms.is_empty() {
            return Vec::new();
        }
        let total: f64 = self.transforms.iter().map(|t| t.prob).sum();
        let mut cum = Vec::with_capacity(self.transforms.len());
        let mut acc = 0.0;
        for t in &self.transforms {
            acc += t.prob / total;
            cum.push(acc);
        }
        let mut x = 0.0f64;
        let mut y = 0.0f64;
        let mut s = seed;
        for _ in 0..20 {
            s = s
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let r = (s >> 11) as f64 / (1u64 << 53) as f64;
            let idx = cum
                .iter()
                .position(|&c| r <= c)
                .unwrap_or(self.transforms.len() - 1);
            let (nx, ny) = self.transforms[idx].apply(x, y);
            x = nx;
            y = ny;
        }
        let mut points = Vec::with_capacity(n);
        for _ in 0..n {
            s = s
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let r = (s >> 11) as f64 / (1u64 << 53) as f64;
            let idx = cum
                .iter()
                .position(|&c| r <= c)
                .unwrap_or(self.transforms.len() - 1);
            let (nx, ny) = self.transforms[idx].apply(x, y);
            x = nx;
            y = ny;
            points.push([x, y]);
        }
        points
    }
}
/// Aggregates noise generation methods.
#[derive(Debug, Clone)]
pub struct NoiseGeometry;
impl NoiseGeometry {
    /// Create a 2D Perlin noise instance.
    pub fn perlin2d(seed: u64) -> PerlinNoise2d {
        PerlinNoise2d::new(seed)
    }
    /// Create a 3D Perlin noise instance.
    pub fn perlin3d(seed: u64) -> PerlinNoise3d {
        PerlinNoise3d::new(seed)
    }
    /// Create a Worley noise instance.
    pub fn worley(num_points: usize, seed: u64) -> WorleyNoise {
        WorleyNoise::new(num_points, seed)
    }
    /// Generate a 2D heightmap using fBm noise.
    ///
    /// Returns a grid of `width × height` values in approximately \[-1, 1\].
    pub fn heightmap_fbm(
        width: usize,
        height: usize,
        scale: f64,
        octaves: u32,
        seed: u64,
    ) -> Vec<Vec<f64>> {
        let noise = PerlinNoise2d::new(seed);
        (0..height)
            .map(|row| {
                (0..width)
                    .map(|col| {
                        let x = col as f64 / width as f64 * scale;
                        let y = row as f64 / height as f64 * scale;
                        noise.fbm(x, y, octaves, 2.0, 0.5)
                    })
                    .collect()
            })
            .collect()
    }
}
/// Represents the current state of an L-system after string rewriting.
#[derive(Debug, Clone)]
pub struct LSystemState {
    /// The current L-system string
    pub string: String,
    /// Current generation number
    pub generation: u32,
}
/// Space colonization algorithm for tree and vascular structure generation.
#[derive(Debug, Clone)]
pub struct SpaceColonization {
    /// Attractor points defining the target space
    pub attractors: Vec<[f64; 3]>,
    /// Current tree nodes
    pub nodes: Vec<ColonizationNode>,
    /// Segment length for each growth step
    pub segment_length: f64,
    /// Influence radius of attractors
    pub influence_radius: f64,
    /// Kill radius (attractor removed when node is within this distance)
    pub kill_radius: f64,
}
impl SpaceColonization {
    /// Create a new space colonization system.
    pub fn new(
        attractors: Vec<[f64; 3]>,
        root: [f64; 3],
        segment_length: f64,
        influence_radius: f64,
        kill_radius: f64,
    ) -> Self {
        let root_node = ColonizationNode {
            position: root,
            parent: None,
            direction: [0.0, 1.0, 0.0],
            num_attractors: 0,
        };
        Self {
            attractors,
            nodes: vec![root_node],
            segment_length,
            influence_radius,
            kill_radius,
        }
    }
    fn dist3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        let dz = a[2] - b[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    fn normalize3(v: [f64; 3]) -> [f64; 3] {
        let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if len < 1e-12 {
            [0.0, 1.0, 0.0]
        } else {
            [v[0] / len, v[1] / len, v[2] / len]
        }
    }
    /// Perform one growth step.
    ///
    /// Returns the number of new nodes added.
    pub fn grow(&mut self) -> usize {
        if self.attractors.is_empty() {
            return 0;
        }
        let n = self.nodes.len();
        let mut directions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0]; n];
        let mut influenced: Vec<bool> = vec![false; n];
        for attr in &self.attractors {
            let mut closest_idx = None;
            let mut closest_dist = self.influence_radius;
            for (i, node) in self.nodes.iter().enumerate() {
                let d = Self::dist3(attr, &node.position);
                if d < closest_dist {
                    closest_dist = d;
                    closest_idx = Some(i);
                }
            }
            if let Some(idx) = closest_idx {
                let node_pos = self.nodes[idx].position;
                let dir = [
                    attr[0] - node_pos[0],
                    attr[1] - node_pos[1],
                    attr[2] - node_pos[2],
                ];
                let norm = Self::normalize3(dir);
                directions[idx][0] += norm[0];
                directions[idx][1] += norm[1];
                directions[idx][2] += norm[2];
                influenced[idx] = true;
            }
        }
        let mut new_nodes = Vec::new();
        for (i, &inf) in influenced.iter().enumerate() {
            if inf {
                let dir = Self::normalize3(directions[i]);
                let pos = self.nodes[i].position;
                let new_pos = [
                    pos[0] + dir[0] * self.segment_length,
                    pos[1] + dir[1] * self.segment_length,
                    pos[2] + dir[2] * self.segment_length,
                ];
                new_nodes.push(ColonizationNode {
                    position: new_pos,
                    parent: Some(i),
                    direction: dir,
                    num_attractors: 0,
                });
            }
        }
        let added = new_nodes.len();
        self.nodes.extend(new_nodes);
        let kill_radius = self.kill_radius;
        let nodes = &self.nodes;
        self.attractors.retain(|attr| {
            nodes
                .iter()
                .all(|node| Self::dist3(attr, &node.position) > kill_radius)
        });
        added
    }
    /// Grow until no more attractors remain or max_steps reached.
    pub fn grow_until_done(&mut self, max_steps: u32) {
        for _ in 0..max_steps {
            if self.attractors.is_empty() {
                break;
            }
            self.grow();
        }
    }
}
/// A 2D Voronoi tessellation.
#[derive(Debug, Clone)]
pub struct VoronoiTessellation2d {
    /// The seed points
    pub seeds: Vec<[f64; 2]>,
    /// Domain bounds \[x_min, x_max, y_min, y_max\]
    pub bounds: [f64; 4],
}
impl VoronoiTessellation2d {
    /// Create a new Voronoi tessellation over the given domain.
    pub fn new(seeds: Vec<[f64; 2]>, bounds: [f64; 4]) -> Self {
        Self { seeds, bounds }
    }
    /// Find the index of the nearest seed to point (x, y).
    pub fn nearest_seed(&self, x: f64, y: f64) -> usize {
        self.seeds
            .iter()
            .enumerate()
            .map(|(i, &[sx, sy])| (i, (x - sx).powi(2) + (y - sy).powi(2)))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0)
    }
    /// Compute the approximate area of each Voronoi cell by sampling the domain.
    ///
    /// Uses a raster scan with `resolution × resolution` sample points.
    pub fn cell_areas(&self, resolution: usize) -> Vec<f64> {
        let [x_min, x_max, y_min, y_max] = self.bounds;
        let total_area = (x_max - x_min) * (y_max - y_min);
        let n = self.seeds.len();
        let mut counts = vec![0usize; n];
        let total_samples = resolution * resolution;
        for row in 0..resolution {
            for col in 0..resolution {
                let x = x_min + (col as f64 + 0.5) / resolution as f64 * (x_max - x_min);
                let y = y_min + (row as f64 + 0.5) / resolution as f64 * (y_max - y_min);
                let idx = self.nearest_seed(x, y);
                counts[idx] += 1;
            }
        }
        counts
            .iter()
            .map(|&c| c as f64 / total_samples as f64 * total_area)
            .collect()
    }
    /// Find neighboring cells (cells that share an approximate boundary).
    pub fn neighbors(&self, resolution: usize) -> Vec<Vec<usize>> {
        let [x_min, x_max, y_min, y_max] = self.bounds;
        let n = self.seeds.len();
        let mut neighbor_sets: Vec<std::collections::HashSet<usize>> =
            (0..n).map(|_| std::collections::HashSet::new()).collect();
        for row in 0..resolution {
            for col in 0..resolution {
                let x = x_min + (col as f64 + 0.5) / resolution as f64 * (x_max - x_min);
                let y = y_min + (row as f64 + 0.5) / resolution as f64 * (y_max - y_min);
                let idx = self.nearest_seed(x, y);
                for (dr, dc) in [(-1i32, 0), (1, 0), (0, -1i32), (0, 1)] {
                    let nr = row as i32 + dr;
                    let nc = col as i32 + dc;
                    if nr >= 0 && nr < resolution as i32 && nc >= 0 && nc < resolution as i32 {
                        let nx = x_min + (nc as f64 + 0.5) / resolution as f64 * (x_max - x_min);
                        let ny = y_min + (nr as f64 + 0.5) / resolution as f64 * (y_max - y_min);
                        let nidx = self.nearest_seed(nx, ny);
                        if nidx != idx {
                            neighbor_sets[idx].insert(nidx);
                            neighbor_sets[nidx].insert(idx);
                        }
                    }
                }
            }
        }
        neighbor_sets
            .into_iter()
            .map(|s| s.into_iter().collect())
            .collect()
    }
    /// Generate random seed points within the domain using LCG.
    pub fn random_seeds(n: usize, bounds: [f64; 4], seed: u64) -> Self {
        let [x_min, x_max, y_min, y_max] = bounds;
        let mut s = seed;
        let mut seeds = Vec::with_capacity(n);
        for _ in 0..n {
            s = s
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let x = x_min + (s >> 11) as f64 / (1u64 << 53) as f64 * (x_max - x_min);
            s = s
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let y = y_min + (s >> 11) as f64 / (1u64 << 53) as f64 * (y_max - y_min);
            seeds.push([x, y]);
        }
        Self::new(seeds, bounds)
    }
}
/// Fractal geometry analysis: box-counting dimension and Hausdorff dimension estimation.
#[derive(Debug, Clone)]
pub struct FractalGeometry {
    /// Set of 2D points forming the fractal set
    pub points: Vec<[f64; 2]>,
}
impl FractalGeometry {
    /// Create a new FractalGeometry from a set of 2D points.
    pub fn new(points: Vec<[f64; 2]>) -> Self {
        Self { points }
    }
    /// Estimate fractal dimension via box-counting method.
    ///
    /// Covers the point set with boxes of decreasing size and measures
    /// the scaling of the count N(ε) with box size ε. The slope of
    /// log(N) vs log(1/ε) gives the box-counting dimension.
    pub fn box_counting_dimension(&self, num_scales: usize) -> f64 {
        if self.points.is_empty() || num_scales < 2 {
            return 0.0;
        }
        let (mut x_min, mut x_max) = (f64::MAX, f64::MIN);
        let (mut y_min, mut y_max) = (f64::MAX, f64::MIN);
        for &[x, y] in &self.points {
            x_min = x_min.min(x);
            x_max = x_max.max(x);
            y_min = y_min.min(y);
            y_max = y_max.max(y);
        }
        let extent = (x_max - x_min).max(y_max - y_min).max(1e-10);
        let mut log_inv_eps = Vec::with_capacity(num_scales);
        let mut log_count = Vec::with_capacity(num_scales);
        for k in 1..=num_scales {
            let eps = extent / (2_f64.powi(k as i32));
            if eps < 1e-12 {
                break;
            }
            let inv_eps = 1.0 / eps;
            let mut occupied: std::collections::HashSet<(i64, i64)> =
                std::collections::HashSet::new();
            for &[x, y] in &self.points {
                let ix = ((x - x_min) / eps).floor() as i64;
                let iy = ((y - y_min) / eps).floor() as i64;
                occupied.insert((ix, iy));
            }
            let n = occupied.len() as f64;
            if n > 0.0 {
                log_inv_eps.push(inv_eps.ln());
                log_count.push(n.ln());
            }
        }
        if log_inv_eps.len() < 2 {
            return 0.0;
        }
        let n = log_inv_eps.len() as f64;
        let sum_x: f64 = log_inv_eps.iter().sum();
        let sum_y: f64 = log_count.iter().sum();
        let sum_xy: f64 = log_inv_eps.iter().zip(&log_count).map(|(x, y)| x * y).sum();
        let sum_x2: f64 = log_inv_eps.iter().map(|x| x * x).sum();
        let denom = n * sum_x2 - sum_x * sum_x;
        if denom.abs() < 1e-12 {
            return 0.0;
        }
        (n * sum_xy - sum_x * sum_y) / denom
    }
    /// Estimate Hausdorff dimension (approximated via box-counting on fine scale).
    ///
    /// For self-similar fractals, the Hausdorff dimension equals the box-counting
    /// dimension. This method uses more scales for a more accurate estimate.
    pub fn hausdorff_dimension(&self) -> f64 {
        self.box_counting_dimension(12)
    }
    /// Generate Koch snowflake points up to a given iteration depth.
    pub fn koch_snowflake(iterations: u32, num_points_per_segment: usize) -> Self {
        let mut segments: Vec<([f64; 2], [f64; 2])> = vec![
            ([0.0, 0.0], [1.0, 0.0]),
            ([1.0, 0.0], [0.5, 0.866_025_4]),
            ([0.5, 0.866_025_4], [0.0, 0.0]),
        ];
        for _ in 0..iterations {
            let mut new_segments = Vec::new();
            for (a, b) in &segments {
                let p1 = *a;
                let p2 = [a[0] + (b[0] - a[0]) / 3.0, a[1] + (b[1] - a[1]) / 3.0];
                let p4 = [
                    a[0] + 2.0 * (b[0] - a[0]) / 3.0,
                    a[1] + 2.0 * (b[1] - a[1]) / 3.0,
                ];
                let p5 = *b;
                let mx = (p2[0] + p4[0]) / 2.0;
                let my = (p2[1] + p4[1]) / 2.0;
                let dx = p4[0] - p2[0];
                let dy = p4[1] - p2[1];
                let p3 = [mx - dy * 3_f64.sqrt() / 2.0, my + dx * 3_f64.sqrt() / 2.0];
                new_segments.push((p1, p2));
                new_segments.push((p2, p3));
                new_segments.push((p3, p4));
                new_segments.push((p4, p5));
            }
            segments = new_segments;
        }
        let n = num_points_per_segment.max(2);
        let mut points = Vec::new();
        for (a, b) in &segments {
            for i in 0..n {
                let t = i as f64 / (n - 1) as f64;
                points.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
            }
        }
        Self { points }
    }
    /// Generate Sierpinski triangle points up to a given iteration depth.
    pub fn sierpinski_triangle(iterations: u32) -> Self {
        let mut triangles: Vec<([f64; 2], [f64; 2], [f64; 2])> =
            vec![([0.0, 0.0], [1.0, 0.0], [0.5, 0.866_025_4])];
        for _ in 0..iterations {
            let mut new_triangles = Vec::new();
            for (a, b, c) in &triangles {
                let ab = [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0];
                let bc = [(b[0] + c[0]) / 2.0, (b[1] + c[1]) / 2.0];
                let ca = [(c[0] + a[0]) / 2.0, (c[1] + a[1]) / 2.0];
                new_triangles.push((*a, ab, ca));
                new_triangles.push((ab, *b, bc));
                new_triangles.push((ca, bc, *c));
            }
            triangles = new_triangles;
        }
        let mut points = Vec::new();
        for (a, b, c) in &triangles {
            points.push(*a);
            points.push(*b);
            points.push(*c);
        }
        Self { points }
    }
}
/// A single production rule in an L-system: maps a predecessor character to a successor string.
#[derive(Debug, Clone)]
pub struct LSystemRule {
    /// The predecessor symbol
    pub predecessor: char,
    /// The successor string
    pub successor: String,
}
impl LSystemRule {
    /// Create a new L-system production rule.
    pub fn new(predecessor: char, successor: &str) -> Self {
        Self {
            predecessor,
            successor: successor.to_string(),
        }
    }
}
/// A triangle in the Delaunay triangulation.
#[derive(Debug, Clone, Copy)]
pub struct DelaunayTri {
    /// Indices into the point array
    pub indices: [usize; 3],
}
impl DelaunayTri {
    /// Create a new Delaunay triangle.
    pub fn new(a: usize, b: usize, c: usize) -> Self {
        Self { indices: [a, b, c] }
    }
    /// Compute the circumcircle center and radius squared of this triangle.
    pub fn circumcircle(&self, points: &[[f64; 2]]) -> ([f64; 2], f64) {
        let [ax, ay] = points[self.indices[0]];
        let [bx, by] = points[self.indices[1]];
        let [cx, cy] = points[self.indices[2]];
        let d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
        if d.abs() < 1e-12 {
            return ([0.0, 0.0], f64::MAX);
        }
        let ux = ((ax * ax + ay * ay) * (by - cy)
            + (bx * bx + by * by) * (cy - ay)
            + (cx * cx + cy * cy) * (ay - by))
            / d;
        let uy = ((ax * ax + ay * ay) * (cx - bx)
            + (bx * bx + by * by) * (ax - cx)
            + (cx * cx + cy * cy) * (bx - ax))
            / d;
        let r2 = (ax - ux).powi(2) + (ay - uy).powi(2);
        ([ux, uy], r2)
    }
    /// Check whether point `p` lies inside the circumcircle of this triangle.
    pub fn in_circumcircle(&self, points: &[[f64; 2]], p: [f64; 2]) -> bool {
        let ([cx, cy], r2) = self.circumcircle(points);
        (p[0] - cx).powi(2) + (p[1] - cy).powi(2) < r2 - 1e-10
    }
}
/// Delaunay triangulation using the Bowyer-Watson algorithm.
#[derive(Debug, Clone)]
pub struct DelaunayTriangulation2d {
    /// Input points
    pub points: Vec<[f64; 2]>,
    /// Output triangles
    pub triangles: Vec<DelaunayTri>,
}
impl DelaunayTriangulation2d {
    /// Compute the Delaunay triangulation of the given points using Bowyer-Watson.
    pub fn new(points: Vec<[f64; 2]>) -> Self {
        let mut tris = Self::bowyer_watson(&points);
        let n = points.len();
        tris.retain(|t| t.indices.iter().all(|&i| i < n));
        Self {
            points,
            triangles: tris,
        }
    }
    fn bowyer_watson(points: &[[f64; 2]]) -> Vec<DelaunayTri> {
        if points.is_empty() {
            return Vec::new();
        }
        let (mut x_min, mut x_max) = (f64::MAX, f64::MIN);
        let (mut y_min, mut y_max) = (f64::MAX, f64::MIN);
        for &[x, y] in points {
            x_min = x_min.min(x);
            x_max = x_max.max(x);
            y_min = y_min.min(y);
            y_max = y_max.max(y);
        }
        let dx = x_max - x_min;
        let dy = y_max - y_min;
        let delta_max = dx.max(dy) * 10.0;
        let mid_x = (x_min + x_max) / 2.0;
        let mid_y = (y_min + y_max) / 2.0;
        let n = points.len();
        let mut all_points = points.to_vec();
        all_points.push([mid_x - 20.0 * delta_max, mid_y - delta_max]);
        all_points.push([mid_x, mid_y + 20.0 * delta_max]);
        all_points.push([mid_x + 20.0 * delta_max, mid_y - delta_max]);
        let mut triangles = vec![DelaunayTri::new(n, n + 1, n + 2)];
        for (pi, &point) in points.iter().enumerate() {
            let bad: Vec<usize> = triangles
                .iter()
                .enumerate()
                .filter(|(_, t)| t.in_circumcircle(&all_points, point))
                .map(|(i, _)| i)
                .collect();
            let mut boundary: Vec<(usize, usize)> = Vec::new();
            for &bi in &bad {
                let t = triangles[bi];
                let edges = [
                    (t.indices[0], t.indices[1]),
                    (t.indices[1], t.indices[2]),
                    (t.indices[2], t.indices[0]),
                ];
                for (e0, e1) in edges {
                    let shared = bad.iter().filter(|&&bj| bj != bi).any(|&bj| {
                        let other = triangles[bj];
                        let other_edges = [
                            (other.indices[0], other.indices[1]),
                            (other.indices[1], other.indices[2]),
                            (other.indices[2], other.indices[0]),
                        ];
                        other_edges
                            .iter()
                            .any(|&(a, b)| (a == e0 && b == e1) || (a == e1 && b == e0))
                    });
                    if !shared {
                        boundary.push((e0, e1));
                    }
                }
            }
            let mut keep = vec![true; triangles.len()];
            for &bi in &bad {
                keep[bi] = false;
            }
            triangles = triangles
                .into_iter()
                .zip(keep)
                .filter(|(_, k)| *k)
                .map(|(t, _)| t)
                .collect();
            for (e0, e1) in boundary {
                triangles.push(DelaunayTri::new(e0, e1, pi));
            }
        }
        triangles
    }
    /// Verify the Delaunay property: no point lies inside the circumcircle of any triangle.
    pub fn verify_delaunay(&self) -> bool {
        for tri in &self.triangles {
            for (i, &p) in self.points.iter().enumerate() {
                if tri.indices.contains(&i) {
                    continue;
                }
                if tri.in_circumcircle(&self.points, p) {
                    return false;
                }
            }
        }
        true
    }
}
/// Parametric curve types: Bézier, B-spline, NURBS.
#[derive(Debug, Clone)]
pub struct CurveParametric {
    /// Control points in 3D
    pub control_points: Vec<[f64; 3]>,
}
impl CurveParametric {
    /// Create a new parametric curve from control points.
    pub fn new(control_points: Vec<[f64; 3]>) -> Self {
        Self { control_points }
    }
    /// Evaluate a Bézier curve at parameter t ∈ \[0, 1\] using de Casteljau algorithm.
    ///
    /// Works for any degree (number of control points - 1).
    pub fn bezier_de_casteljau(&self, t: f64) -> [f64; 3] {
        let n = self.control_points.len();
        if n == 0 {
            return [0.0, 0.0, 0.0];
        }
        let mut pts = self.control_points.clone();
        for r in 1..n {
            for i in 0..(n - r) {
                pts[i] = [
                    (1.0 - t) * pts[i][0] + t * pts[i + 1][0],
                    (1.0 - t) * pts[i][1] + t * pts[i + 1][1],
                    (1.0 - t) * pts[i][2] + t * pts[i + 1][2],
                ];
            }
        }
        pts[0]
    }
    /// Evaluate a cubic Bézier curve using the Bernstein polynomial formula.
    ///
    /// Requires exactly 4 control points. Falls back to de Casteljau for other counts.
    pub fn bezier_bernstein(&self, t: f64) -> [f64; 3] {
        if self.control_points.len() != 4 {
            return self.bezier_de_casteljau(t);
        }
        let [p0, p1, p2, p3] = [
            self.control_points[0],
            self.control_points[1],
            self.control_points[2],
            self.control_points[3],
        ];
        let mt = 1.0 - t;
        let b0 = mt * mt * mt;
        let b1 = 3.0 * mt * mt * t;
        let b2 = 3.0 * mt * t * t;
        let b3 = t * t * t;
        [
            b0 * p0[0] + b1 * p1[0] + b2 * p2[0] + b3 * p3[0],
            b0 * p0[1] + b1 * p1[1] + b2 * p2[1] + b3 * p3[1],
            b0 * p0[2] + b1 * p1[2] + b2 * p2[2] + b3 * p3[2],
        ]
    }
    /// Evaluate a B-spline at parameter t using the de Boor algorithm.
    ///
    /// `knots` must have length = `control_points.len() + degree + 1`.
    pub fn bspline_de_boor(&self, t: f64, degree: usize, knots: &[f64]) -> [f64; 3] {
        let n = self.control_points.len();
        if n == 0 || knots.len() < n + degree + 1 {
            return [0.0, 0.0, 0.0];
        }
        let t = t.clamp(knots[degree], knots[n]);
        let mut k = degree;
        for i in degree..n {
            if t >= knots[i] && t < knots[i + 1] {
                k = i;
                break;
            }
        }
        if t >= knots[n] {
            k = n - 1;
        }
        let mut d: Vec<[f64; 3]> = (0..=degree)
            .map(|j| self.control_points[k - degree + j])
            .collect();
        for r in 1..=degree {
            for j in (r..=degree).rev() {
                let i = k - degree + j;
                let denom = knots[i + degree - r + 1] - knots[i];
                let alpha = if denom.abs() < 1e-12 {
                    0.0
                } else {
                    (t - knots[i]) / denom
                };
                d[j] = [
                    (1.0 - alpha) * d[j - 1][0] + alpha * d[j][0],
                    (1.0 - alpha) * d[j - 1][1] + alpha * d[j][1],
                    (1.0 - alpha) * d[j - 1][2] + alpha * d[j][2],
                ];
            }
        }
        d[degree]
    }
    /// Evaluate a NURBS curve at parameter t.
    ///
    /// `weights` must have the same length as `control_points`.
    pub fn nurbs(&self, t: f64, degree: usize, knots: &[f64], weights: &[f64]) -> [f64; 3] {
        let n = self.control_points.len();
        if n == 0 || weights.len() != n {
            return [0.0, 0.0, 0.0];
        }
        let homo_xyz: Vec<[f64; 3]> = self
            .control_points
            .iter()
            .zip(weights)
            .map(|(&[x, y, z], &w)| [x * w, y * w, z * w])
            .collect();
        let w_pts: Vec<[f64; 3]> = weights.iter().map(|&w| [w, 0.0, 0.0]).collect();
        let homo_curve = CurveParametric::new(homo_xyz);
        let w_curve = CurveParametric::new(w_pts);
        let p = homo_curve.bspline_de_boor(t, degree, knots);
        let w = w_curve.bspline_de_boor(t, degree, knots)[0];
        if w.abs() < 1e-12 {
            return [0.0, 0.0, 0.0];
        }
        [p[0] / w, p[1] / w, p[2] / w]
    }
    /// Sample the Bézier curve at `n` uniformly spaced parameter values.
    pub fn sample_bezier(&self, n: usize) -> Vec<[f64; 3]> {
        (0..n)
            .map(|i| {
                let t = i as f64 / (n - 1).max(1) as f64;
                self.bezier_de_casteljau(t)
            })
            .collect()
    }
    /// Compute the approximate arc length of the Bézier curve by sampling.
    pub fn bezier_arc_length(&self, samples: usize) -> f64 {
        let pts = self.sample_bezier(samples);
        pts.windows(2)
            .map(|w| {
                let [ax, ay, az] = w[0];
                let [bx, by, bz] = w[1];
                ((ax - bx).powi(2) + (ay - by).powi(2) + (az - bz).powi(2)).sqrt()
            })
            .sum()
    }
}
/// 3D Perlin noise generator.
#[derive(Debug, Clone)]
pub struct PerlinNoise3d {
    pub(super) perm: Vec<u8>,
}
impl PerlinNoise3d {
    /// Create a new 3D Perlin noise generator.
    pub fn new(seed: u64) -> Self {
        let mut perm = (0u8..=255).collect::<Vec<_>>();
        let mut s = seed;
        for i in (1..256).rev() {
            s = s
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let j = (s >> 33) as usize % (i + 1);
            perm.swap(i, j);
        }
        let mut full = perm.clone();
        full.extend_from_slice(&perm);
        Self { perm: full }
    }
    fn fade(t: f64) -> f64 {
        t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
    }
    fn lerp(a: f64, b: f64, t: f64) -> f64 {
        a + t * (b - a)
    }
    fn grad3(hash: u8, x: f64, y: f64, z: f64) -> f64 {
        match hash & 15 {
            0 => x + y,
            1 => -x + y,
            2 => x - y,
            3 => -x - y,
            4 => x + z,
            5 => -x + z,
            6 => x - z,
            7 => -x - z,
            8 => y + z,
            9 => -y + z,
            10 => y - z,
            11 => -y - z,
            12 => y + x,
            13 => -y + z,
            14 => y - x,
            _ => -y - z,
        }
    }
    fn p(&self, i: usize) -> u8 {
        self.perm[i & 511]
    }
    /// Evaluate 3D Perlin noise at (x, y, z). Returns value in approximately \[-1, 1\].
    pub fn noise(&self, x: f64, y: f64, z: f64) -> f64 {
        let xi = x.floor() as i32 & 255;
        let yi = y.floor() as i32 & 255;
        let zi = z.floor() as i32 & 255;
        let xf = x - x.floor();
        let yf = y - y.floor();
        let zf = z - z.floor();
        let u = Self::fade(xf);
        let v = Self::fade(yf);
        let w = Self::fade(zf);
        let a = self.p(xi as usize).wrapping_add(yi as u8);
        let aa = self.p(a as usize).wrapping_add(zi as u8);
        let ab = self.p((a.wrapping_add(1)) as usize).wrapping_add(zi as u8);
        let b = self.p((xi + 1) as usize).wrapping_add(yi as u8);
        let ba = self.p(b as usize).wrapping_add(zi as u8);
        let bb = self.p((b.wrapping_add(1)) as usize).wrapping_add(zi as u8);
        Self::lerp(
            Self::lerp(
                Self::lerp(
                    Self::grad3(self.p(aa as usize), xf, yf, zf),
                    Self::grad3(self.p(ba as usize), xf - 1.0, yf, zf),
                    u,
                ),
                Self::lerp(
                    Self::grad3(self.p(ab as usize), xf, yf - 1.0, zf),
                    Self::grad3(self.p(bb as usize), xf - 1.0, yf - 1.0, zf),
                    u,
                ),
                v,
            ),
            Self::lerp(
                Self::lerp(
                    Self::grad3(self.p((aa.wrapping_add(1)) as usize), xf, yf, zf - 1.0),
                    Self::grad3(
                        self.p((ba.wrapping_add(1)) as usize),
                        xf - 1.0,
                        yf,
                        zf - 1.0,
                    ),
                    u,
                ),
                Self::lerp(
                    Self::grad3(
                        self.p((ab.wrapping_add(1)) as usize),
                        xf,
                        yf - 1.0,
                        zf - 1.0,
                    ),
                    Self::grad3(
                        self.p((bb.wrapping_add(1)) as usize),
                        xf - 1.0,
                        yf - 1.0,
                        zf - 1.0,
                    ),
                    u,
                ),
                v,
            ),
            w,
        )
    }
    /// 3D Fractional Brownian Motion.
    pub fn fbm(&self, x: f64, y: f64, z: f64, octaves: u32, lacunarity: f64, gain: f64) -> f64 {
        let mut value = 0.0;
        let mut amplitude = 1.0;
        let mut frequency = 1.0;
        let mut max_val = 0.0;
        for _ in 0..octaves {
            value += amplitude * self.noise(x * frequency, y * frequency, z * frequency);
            max_val += amplitude;
            amplitude *= gain;
            frequency *= lacunarity;
        }
        if max_val > 0.0 { value / max_val } else { 0.0 }
    }
}
/// A node in a space colonization tree.
#[derive(Debug, Clone)]
pub struct ColonizationNode {
    /// 3D position of the node
    pub position: [f64; 3],
    /// Index of the parent node (None for root)
    pub parent: Option<usize>,
    /// Accumulated growth direction
    pub direction: [f64; 3],
    /// Number of attractors influencing this node
    pub num_attractors: usize,
}
