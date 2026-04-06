//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use rand::Rng;
use rand::RngExt;

/// Random sequential adsorption (RSA) packing configuration.
#[derive(Debug, Clone)]
pub struct RsaConfig {
    /// Sphere radius.
    pub radius: f64,
    /// Domain dimensions \[lx, ly, lz\].
    pub domain: [f64; 3],
    /// Maximum number of attempts before declaring jamming.
    pub max_attempts: usize,
    /// Minimum gap between spheres (m). Usually 0.
    pub min_gap: f64,
}
/// Polydisperse sphere packing with various size distributions.
pub struct PolyDispersePacking {
    /// Distribution type.
    pub distribution: SizeDistribution,
    /// Distribution parameters.
    pub params: SizeDistributionParams,
    /// Domain dimensions.
    pub domain: [f64; 3],
    /// Placed spheres.
    pub spheres: Vec<PackingSphere>,
    /// Maximum placement attempts.
    pub max_attempts: usize,
}
impl PolyDispersePacking {
    /// Create a new polydisperse packing.
    pub fn new(
        distribution: SizeDistribution,
        params: SizeDistributionParams,
        domain: [f64; 3],
        max_attempts: usize,
    ) -> Self {
        Self {
            distribution,
            params,
            domain,
            spheres: Vec::new(),
            max_attempts,
        }
    }
    /// Sample a radius from the configured distribution.
    pub fn sample_radius(&self, rng: &mut impl Rng) -> f64 {
        match self.distribution {
            SizeDistribution::Monodisperse => self.params.mean_radius,
            SizeDistribution::Bidisperse => {
                if rng.random::<f64>() < self.params.large_fraction {
                    self.params.mean_radius * self.params.size_ratio
                } else {
                    self.params.mean_radius
                }
            }
            SizeDistribution::LogNormal => {
                let mu = self.params.mean_radius.ln();
                let sigma = self.params.std_radius;
                let z = self.normal_sample(rng);
                let r = (mu + sigma * z).exp();
                r.clamp(self.params.min_radius, self.params.max_radius)
            }
            SizeDistribution::Gaussian => {
                let r = self.params.mean_radius + self.params.std_radius * self.normal_sample(rng);
                r.clamp(self.params.min_radius, self.params.max_radius)
            }
            SizeDistribution::PowerLaw => {
                let alpha = self.params.power_law_exponent;
                let r_min = self.params.min_radius;
                let r_max = self.params.max_radius;
                let u: f64 = rng.random();
                if (alpha + 1.0).abs() < 1e-12 {
                    r_min * (r_max / r_min).powf(u)
                } else {
                    let p = alpha + 1.0;
                    (r_min.powf(p) + u * (r_max.powf(p) - r_min.powf(p))).powf(1.0 / p)
                }
            }
        }
    }
    /// Box-Muller normal sample.
    fn normal_sample(&self, rng: &mut impl Rng) -> f64 {
        let u1: f64 = rng.random::<f64>().max(1e-300);
        let u2: f64 = rng.random();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
    /// Check if a sphere with given center and radius can be placed.
    pub fn is_valid_placement(&self, center: [f64; 3], radius: f64) -> bool {
        let [lx, ly, lz] = self.domain;
        if center[0] - radius < 0.0
            || center[0] + radius > lx
            || center[1] - radius < 0.0
            || center[1] + radius > ly
            || center[2] - radius < 0.0
            || center[2] + radius > lz
        {
            return false;
        }
        let candidate = PackingSphere::new(center, radius, 0);
        for s in &self.spheres {
            if s.overlaps(&candidate) {
                return false;
            }
        }
        true
    }
    /// Run RSA with polydisperse radii.
    ///
    /// Returns the number of spheres placed.
    pub fn run(&mut self, n_spheres: usize) -> usize {
        let mut rng = rand::rng();
        let mut attempts = 0;
        while self.spheres.len() < n_spheres && attempts < self.max_attempts {
            let r = self.sample_radius(&mut rng);
            let [lx, ly, lz] = self.domain;
            let cx = rng.random_range(r..(lx - r).max(r + 1e-12));
            let cy = rng.random_range(r..(ly - r).max(r + 1e-12));
            let cz = rng.random_range(r..(lz - r).max(r + 1e-12));
            if self.is_valid_placement([cx, cy, cz], r) {
                let id = self.spheres.len();
                self.spheres.push(PackingSphere::new([cx, cy, cz], r, id));
            }
            attempts += 1;
        }
        self.spheres.len()
    }
    /// Compute the packing fraction.
    pub fn packing_fraction(&self) -> f64 {
        let [lx, ly, lz] = self.domain;
        let vol = lx * ly * lz;
        let sv: f64 = self.spheres.iter().map(|s| s.volume()).sum();
        sv / vol
    }
    /// Return radius statistics: (mean, std, min, max).
    pub fn radius_statistics(&self) -> (f64, f64, f64, f64) {
        if self.spheres.is_empty() {
            return (0.0, 0.0, 0.0, 0.0);
        }
        let radii: Vec<f64> = self.spheres.iter().map(|s| s.radius).collect();
        let n = radii.len() as f64;
        let mean = radii.iter().sum::<f64>() / n;
        let variance = radii.iter().map(|&r| (r - mean).powi(2)).sum::<f64>() / n;
        let std = variance.sqrt();
        let min = radii.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = radii.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (mean, std, min, max)
    }
    /// Compute the polydispersity index (PDI): std/mean.
    pub fn polydispersity_index(&self) -> f64 {
        let (mean, std, _, _) = self.radius_statistics();
        if mean > 1e-30 { std / mean } else { 0.0 }
    }
    /// Compute the size distribution histogram.
    ///
    /// Returns (bin_centers, counts).
    pub fn size_histogram(&self, n_bins: usize) -> (Vec<f64>, Vec<usize>) {
        if self.spheres.is_empty() || n_bins == 0 {
            return (Vec::new(), Vec::new());
        }
        let (_, _, min_r, max_r) = self.radius_statistics();
        let bin_width = (max_r - min_r) / n_bins as f64;
        if bin_width < 1e-30 {
            let centers = vec![min_r; n_bins];
            let mut counts = vec![0usize; n_bins];
            counts[0] = self.spheres.len();
            return (centers, counts);
        }
        let mut counts = vec![0usize; n_bins];
        for s in &self.spheres {
            let bin = ((s.radius - min_r) / bin_width).floor() as usize;
            let bin = bin.min(n_bins - 1);
            counts[bin] += 1;
        }
        let centers: Vec<f64> = (0..n_bins)
            .map(|i| min_r + (i as f64 + 0.5) * bin_width)
            .collect();
        (centers, counts)
    }
}
/// Void space analysis for sphere packings.
///
/// Analyses the pore space between spheres, computing pore size distributions,
/// percolation thresholds, and void connectivity.
pub struct VoidSpaceAnalysis {
    /// Spheres in the packing.
    pub spheres: Vec<PackingSphere>,
    /// Domain dimensions.
    pub domain: [f64; 3],
    /// Grid resolution for voxelization.
    pub grid_resolution: usize,
}
impl VoidSpaceAnalysis {
    /// Create a new void space analyzer.
    pub fn new(spheres: Vec<PackingSphere>, domain: [f64; 3], grid_resolution: usize) -> Self {
        Self {
            spheres,
            domain,
            grid_resolution,
        }
    }
    /// Check if a point is inside any sphere.
    pub fn point_in_sphere(&self, pt: [f64; 3]) -> bool {
        for s in &self.spheres {
            let dx = pt[0] - s.center[0];
            let dy = pt[1] - s.center[1];
            let dz = pt[2] - s.center[2];
            if dx * dx + dy * dy + dz * dz <= s.radius * s.radius {
                return true;
            }
        }
        false
    }
    /// Compute the void fraction by voxelization.
    pub fn void_fraction_by_voxelization(&self) -> f64 {
        let n = self.grid_resolution;
        let [lx, ly, lz] = self.domain;
        let dx = lx / n as f64;
        let dy = ly / n as f64;
        let dz = lz / n as f64;
        let mut void_count = 0usize;
        let total = n * n * n;
        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    let pt = [
                        (ix as f64 + 0.5) * dx,
                        (iy as f64 + 0.5) * dy,
                        (iz as f64 + 0.5) * dz,
                    ];
                    if !self.point_in_sphere(pt) {
                        void_count += 1;
                    }
                }
            }
        }
        void_count as f64 / total as f64
    }
    /// Compute the distance from point `pt` to the nearest sphere surface.
    ///
    /// Returns the maximum inscribed sphere radius at this point (pore radius).
    pub fn pore_radius_at(&self, pt: [f64; 3]) -> f64 {
        if self.point_in_sphere(pt) {
            return 0.0;
        }
        let mut min_dist_to_surface = f64::INFINITY;
        for s in &self.spheres {
            let dx = pt[0] - s.center[0];
            let dy = pt[1] - s.center[1];
            let dz = pt[2] - s.center[2];
            let dist_to_center = (dx * dx + dy * dy + dz * dz).sqrt();
            let dist_to_surface = (dist_to_center - s.radius).max(0.0);
            min_dist_to_surface = min_dist_to_surface.min(dist_to_surface);
        }
        let [lx, ly, lz] = self.domain;
        let wall_dist = [pt[0], lx - pt[0], pt[1], ly - pt[1], pt[2], lz - pt[2]]
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        min_dist_to_surface.min(wall_dist)
    }
    /// Compute the pore size distribution using random sampling.
    ///
    /// Samples `n_samples` random points in the void space and
    /// returns a histogram of pore radii.
    pub fn pore_size_distribution(&self, n_samples: usize, n_bins: usize) -> Vec<(f64, usize)> {
        let mut rng = rand::rng();
        let [lx, ly, lz] = self.domain;
        let mut pore_radii = Vec::new();
        for _ in 0..n_samples {
            let pt = [
                rng.random_range(0.0..lx),
                rng.random_range(0.0..ly),
                rng.random_range(0.0..lz),
            ];
            if !self.point_in_sphere(pt) {
                pore_radii.push(self.pore_radius_at(pt));
            }
        }
        if pore_radii.is_empty() || n_bins == 0 {
            return Vec::new();
        }
        let r_min = pore_radii.iter().cloned().fold(f64::INFINITY, f64::min);
        let r_max = pore_radii.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let bin_width = (r_max - r_min) / n_bins as f64;
        if bin_width < 1e-30 {
            return vec![(r_min, pore_radii.len())];
        }
        let mut counts = vec![0usize; n_bins];
        for &r in &pore_radii {
            let bin = ((r - r_min) / bin_width).floor() as usize;
            counts[bin.min(n_bins - 1)] += 1;
        }
        (0..n_bins)
            .map(|i| (r_min + (i as f64 + 0.5) * bin_width, counts[i]))
            .collect()
    }
    /// Estimate the percolation threshold.
    ///
    /// Uses the void fraction as a proxy: percolation typically occurs
    /// when void fraction > 0.32 in 3D.
    pub fn estimate_percolation_threshold(&self) -> f64 {
        0.3116
    }
    /// Check if the packing is above the percolation threshold.
    pub fn is_percolating(&self, void_fraction: f64) -> bool {
        void_fraction > self.estimate_percolation_threshold()
    }
    /// Compute the specific surface area (surface area per unit volume).
    pub fn specific_surface_area(&self) -> f64 {
        let [lx, ly, lz] = self.domain;
        let vol = lx * ly * lz;
        let sa: f64 = self
            .spheres
            .iter()
            .map(|s| 4.0 * std::f64::consts::PI * s.radius * s.radius)
            .sum();
        sa / vol
    }
    /// Perform a full analysis.
    pub fn analyze(&self, n_psd_samples: usize, n_bins: usize) -> VoidSpaceResult {
        let void_frac = self.void_fraction_by_voxelization();
        let [lx, ly, lz] = self.domain;
        let void_vol = void_frac * lx * ly * lz;
        let psd = self.pore_size_distribution(n_psd_samples, n_bins);
        VoidSpaceResult {
            pore_size_distribution: psd,
            percolation_threshold: self.estimate_percolation_threshold(),
            void_volume: void_vol,
            void_fraction: void_frac,
            connectivity: if void_frac > 0.32 { 1.0 } else { 0.0 },
        }
    }
}
/// Sphere packing inside various confining geometries.
///
/// Wall effects (reduced packing near walls) and boundary layers
/// are naturally captured by the placement algorithm.
pub struct ConfiningGeometry {
    /// Confinement configuration.
    pub config: ConfinementConfig,
    /// Placed spheres.
    pub spheres: Vec<PackingSphere>,
    /// Number of rejected attempts.
    pub n_rejected: usize,
}
impl ConfiningGeometry {
    /// Create a new confined packing.
    pub fn new(config: ConfinementConfig) -> Self {
        Self {
            config,
            spheres: Vec::new(),
            n_rejected: 0,
        }
    }
    /// Check if a point is inside the container (accounting for sphere radius).
    pub fn point_inside_container(&self, center: [f64; 3], radius: f64) -> bool {
        match self.config.shape {
            ConfinementShape::Box => {
                let [hx, hy, hz] = self.config.dimensions;
                center[0] - radius >= -hx
                    && center[0] + radius <= hx
                    && center[1] - radius >= -hy
                    && center[1] + radius <= hy
                    && center[2] - radius >= -hz
                    && center[2] + radius <= hz
            }
            ConfinementShape::Cylinder => {
                let r_cyl = self.config.dimensions[0];
                let h_cyl = self.config.dimensions[1];
                let r2 = center[0] * center[0] + center[1] * center[1];
                r2.sqrt() + radius <= r_cyl
                    && center[2] - radius >= -h_cyl
                    && center[2] + radius <= h_cyl
            }
            ConfinementShape::SphericalContainer => {
                let r_cont = self.config.dimensions[0];
                let dist =
                    (center[0] * center[0] + center[1] * center[1] + center[2] * center[2]).sqrt();
                dist + radius <= r_cont
            }
        }
    }
    /// Sample a random point inside the container.
    pub fn random_point_in_container(&self, rng: &mut impl Rng) -> [f64; 3] {
        match self.config.shape {
            ConfinementShape::Box => {
                let [hx, hy, hz] = self.config.dimensions;
                [
                    rng.random_range(-hx..hx),
                    rng.random_range(-hy..hy),
                    rng.random_range(-hz..hz),
                ]
            }
            ConfinementShape::Cylinder => {
                let r_cyl = self.config.dimensions[0];
                let h_cyl = self.config.dimensions[1];
                loop {
                    let x = rng.random_range(-r_cyl..r_cyl);
                    let y = rng.random_range(-r_cyl..r_cyl);
                    let z = rng.random_range(-h_cyl..h_cyl);
                    if (x * x + y * y).sqrt() <= r_cyl {
                        return [x, y, z];
                    }
                }
            }
            ConfinementShape::SphericalContainer => {
                let r_cont = self.config.dimensions[0];
                loop {
                    let x = rng.random_range(-r_cont..r_cont);
                    let y = rng.random_range(-r_cont..r_cont);
                    let z = rng.random_range(-r_cont..r_cont);
                    if (x * x + y * y + z * z).sqrt() <= r_cont {
                        return [x, y, z];
                    }
                }
            }
        }
    }
    /// Check if a candidate sphere overlaps any placed sphere.
    fn has_overlap(&self, center: [f64; 3], radius: f64) -> bool {
        let candidate = PackingSphere::new(center, radius, 0);
        for s in &self.spheres {
            if s.overlaps(&candidate) {
                return true;
            }
        }
        false
    }
    /// Try to place a sphere at a random location inside the container.
    pub fn try_place(&mut self, rng: &mut impl Rng) -> bool {
        let r = self.config.sphere_radius;
        let center = self.random_point_in_container(rng);
        if self.point_inside_container(center, r) && !self.has_overlap(center, r) {
            let id = self.spheres.len();
            self.spheres.push(PackingSphere::new(center, r, id));
            true
        } else {
            self.n_rejected += 1;
            false
        }
    }
    /// Run packing until jamming.
    pub fn run(&mut self) -> usize {
        let mut rng = rand::rng();
        let mut consecutive_failures = 0;
        loop {
            if self.try_place(&mut rng) {
                consecutive_failures = 0;
            } else {
                consecutive_failures += 1;
                if consecutive_failures >= self.config.max_attempts {
                    break;
                }
            }
        }
        self.spheres.len()
    }
    /// Compute the container volume.
    pub fn container_volume(&self) -> f64 {
        match self.config.shape {
            ConfinementShape::Box => {
                let [hx, hy, hz] = self.config.dimensions;
                8.0 * hx * hy * hz
            }
            ConfinementShape::Cylinder => {
                let r = self.config.dimensions[0];
                let h = self.config.dimensions[1];
                std::f64::consts::PI * r * r * 2.0 * h
            }
            ConfinementShape::SphericalContainer => {
                let r = self.config.dimensions[0];
                4.0 / 3.0 * std::f64::consts::PI * r * r * r
            }
        }
    }
    /// Compute the packing fraction.
    pub fn packing_fraction(&self) -> f64 {
        let vol = self.container_volume();
        let sv: f64 = self.spheres.iter().map(|s| s.volume()).sum();
        if vol > 1e-30 { sv / vol } else { 0.0 }
    }
    /// Compute the boundary layer fraction.
    ///
    /// Returns the fraction of spheres within one sphere diameter of the wall.
    pub fn boundary_layer_fraction(&self) -> f64 {
        if self.spheres.is_empty() {
            return 0.0;
        }
        let r = self.config.sphere_radius;
        let layer_thickness = 2.0 * r;
        let n_boundary = self
            .spheres
            .iter()
            .filter(|s| match self.config.shape {
                ConfinementShape::Box => {
                    let [hx, hy, hz] = self.config.dimensions;
                    let dist_x = hx - s.center[0].abs();
                    let dist_y = hy - s.center[1].abs();
                    let dist_z = hz - s.center[2].abs();
                    [dist_x, dist_y, dist_z]
                        .iter()
                        .any(|&d| d <= layer_thickness)
                }
                ConfinementShape::Cylinder => {
                    let r_cyl = self.config.dimensions[0];
                    let h_cyl = self.config.dimensions[1];
                    let dist_r =
                        r_cyl - (s.center[0] * s.center[0] + s.center[1] * s.center[1]).sqrt();
                    let dist_z = h_cyl - s.center[2].abs();
                    dist_r <= layer_thickness || dist_z <= layer_thickness
                }
                ConfinementShape::SphericalContainer => {
                    let r_cont = self.config.dimensions[0];
                    let dist_r = r_cont
                        - (s.center[0] * s.center[0]
                            + s.center[1] * s.center[1]
                            + s.center[2] * s.center[2])
                            .sqrt();
                    dist_r <= layer_thickness
                }
            })
            .count();
        n_boundary as f64 / self.spheres.len() as f64
    }
    /// Compute the radial distribution function g(r) for packing in a spherical container.
    ///
    /// Returns (r_values, g_r) histogram.
    pub fn radial_distribution_function(&self, n_bins: usize, r_max: f64) -> (Vec<f64>, Vec<f64>) {
        let n = self.spheres.len();
        if n < 2 || n_bins == 0 {
            return (Vec::new(), Vec::new());
        }
        let dr = r_max / n_bins as f64;
        let mut hist = vec![0usize; n_bins];
        let mut n_pairs = 0usize;
        for i in 0..n {
            for j in (i + 1)..n {
                let dist = self.spheres[i].distance_to(&self.spheres[j]);
                if dist < r_max {
                    let bin = (dist / dr).floor() as usize;
                    if bin < n_bins {
                        hist[bin] += 2;
                        n_pairs += 1;
                    }
                }
            }
        }
        let vol = self.container_volume();
        let rho = n as f64 / vol;
        let r_values: Vec<f64> = (0..n_bins).map(|i| (i as f64 + 0.5) * dr).collect();
        let g_r: Vec<f64> = (0..n_bins)
            .map(|i| {
                let r = r_values[i];
                let shell_vol = 4.0 * std::f64::consts::PI * r * r * dr;
                let expected = rho * shell_vol * n as f64;
                if expected > 1e-30 {
                    hist[i] as f64 / expected
                } else {
                    0.0
                }
            })
            .collect();
        let _ = n_pairs;
        (r_values, g_r)
    }
}
/// Result of a void space analysis.
#[derive(Debug, Clone)]
pub struct VoidSpaceResult {
    /// Pore size distribution: (radius, count) pairs.
    pub pore_size_distribution: Vec<(f64, usize)>,
    /// Estimated percolation threshold (porosity fraction).
    pub percolation_threshold: f64,
    /// Total void volume.
    pub void_volume: f64,
    /// Void fraction (porosity).
    pub void_fraction: f64,
    /// Connectivity (fraction of void grid points connected to inlet).
    pub connectivity: f64,
}
/// Parameters for a polydisperse size distribution.
#[derive(Debug, Clone)]
pub struct SizeDistributionParams {
    /// Mean radius (m).
    pub mean_radius: f64,
    /// Standard deviation or spread parameter.
    pub std_radius: f64,
    /// Minimum allowed radius (m).
    pub min_radius: f64,
    /// Maximum allowed radius (m).
    pub max_radius: f64,
    /// Size ratio for bidisperse (large/small).
    pub size_ratio: f64,
    /// Volume fraction of large spheres in bidisperse.
    pub large_fraction: f64,
    /// Power-law exponent.
    pub power_law_exponent: f64,
}
/// Configuration for packing in a confining geometry.
#[derive(Debug, Clone)]
pub struct ConfinementConfig {
    /// Shape of the container.
    pub shape: ConfinementShape,
    /// Container dimensions: \[half_x, half_y, half_z\] for box,
    /// \[radius, height, 0\] for cylinder, \[radius, 0, 0\] for sphere.
    pub dimensions: [f64; 3],
    /// Sphere radius for packing.
    pub sphere_radius: f64,
    /// Maximum placement attempts.
    pub max_attempts: usize,
}
/// Random packing via Random Sequential Adsorption (RSA).
///
/// In RSA, spheres are placed one at a time at random positions,
/// accepted only if they do not overlap with already-placed spheres.
/// The process terminates at the jamming limit (~64% for equal spheres in 3D).
pub struct RandomPacking {
    /// RSA configuration.
    pub config: RsaConfig,
    /// Placed spheres.
    pub spheres: Vec<PackingSphere>,
    /// Number of rejected attempts.
    pub n_rejected: usize,
    /// Whether the packing has reached the jamming limit.
    pub jammed: bool,
}
impl RandomPacking {
    /// Create a new RSA packing system.
    pub fn new(config: RsaConfig) -> Self {
        Self {
            config,
            spheres: Vec::new(),
            n_rejected: 0,
            jammed: false,
        }
    }
    /// Attempt to place one sphere at a random location.
    ///
    /// Returns `true` if placement succeeded.
    pub fn try_place_sphere(&mut self, rng: &mut impl Rng) -> bool {
        let r = self.config.radius + self.config.min_gap;
        let [lx, ly, lz] = self.config.domain;
        let cx = rng.random_range(r..(lx - r).max(r + 1e-12));
        let cy = rng.random_range(r..(ly - r).max(r + 1e-12));
        let cz = rng.random_range(r..(lz - r).max(r + 1e-12));
        let candidate = PackingSphere::new([cx, cy, cz], self.config.radius, self.spheres.len());
        if self.is_valid_placement(&candidate) {
            self.spheres.push(candidate);
            true
        } else {
            self.n_rejected += 1;
            false
        }
    }
    /// Check if a sphere can be placed without overlap.
    pub fn is_valid_placement(&self, candidate: &PackingSphere) -> bool {
        let gap = self.config.min_gap;
        for s in &self.spheres {
            let min_dist = s.radius + candidate.radius + gap;
            let dx = s.center[0] - candidate.center[0];
            let dy = s.center[1] - candidate.center[1];
            let dz = s.center[2] - candidate.center[2];
            if dx * dx + dy * dy + dz * dz < min_dist * min_dist {
                return false;
            }
        }
        true
    }
    /// Run the RSA algorithm until jammed.
    ///
    /// Returns the number of spheres placed.
    pub fn run(&mut self) -> usize {
        let mut rng = rand::rng();
        let mut consecutive_failures = 0;
        let jam_threshold = self.config.max_attempts;
        loop {
            if self.try_place_sphere(&mut rng) {
                consecutive_failures = 0;
            } else {
                consecutive_failures += 1;
                if consecutive_failures >= jam_threshold {
                    self.jammed = true;
                    break;
                }
            }
        }
        self.spheres.len()
    }
    /// Compute the current packing fraction (volume fraction).
    pub fn packing_fraction(&self) -> f64 {
        let [lx, ly, lz] = self.config.domain;
        let domain_volume = lx * ly * lz;
        let sphere_volume: f64 = self.spheres.iter().map(|s| s.volume()).sum();
        sphere_volume / domain_volume
    }
    /// Compute the theoretical jamming limit for monodisperse RSA in 3D.
    ///
    /// The RSA jamming limit in 3D is approximately 0.3841.
    pub fn theoretical_jamming_limit() -> f64 {
        0.3841
    }
    /// Return the number of contacts (overlapping pairs) in the packing.
    pub fn count_contacts(&self, tolerance: f64) -> usize {
        let mut count = 0;
        for i in 0..self.spheres.len() {
            for j in (i + 1)..self.spheres.len() {
                let gap = self.spheres[i].overlap_with(&self.spheres[j]);
                if gap.abs() < tolerance {
                    count += 1;
                }
            }
        }
        count
    }
    /// Place spheres at specific locations (for testing).
    pub fn place_at(&mut self, center: [f64; 3]) -> bool {
        let candidate = PackingSphere::new(center, self.config.radius, self.spheres.len());
        if self.is_valid_placement(&candidate) {
            self.spheres.push(candidate);
            true
        } else {
            false
        }
    }
}
/// Ordered sphere packing on crystal lattices.
pub struct OrderedPacking {
    /// Lattice type.
    pub lattice: LatticeType,
    /// Sphere radius.
    pub radius: f64,
    /// Lattice parameter (edge length of unit cell).
    pub lattice_parameter: f64,
    /// Number of unit cells in each direction.
    pub n_cells: [usize; 3],
    /// Generated sphere positions.
    pub spheres: Vec<PackingSphere>,
}
impl OrderedPacking {
    /// Create a new ordered packing.
    ///
    /// The lattice parameter is set so spheres are touching:
    /// - FCC: `a = 2*r*sqrt(2)`
    /// - BCC: `a = 4*r/sqrt(3)`
    /// - SC:  `a = 2*r`
    pub fn new(lattice: LatticeType, radius: f64, n_cells: [usize; 3]) -> Self {
        let a = match lattice {
            LatticeType::Fcc | LatticeType::Hcp => 2.0 * radius * 2.0_f64.sqrt(),
            LatticeType::Bcc => 4.0 * radius / 3.0_f64.sqrt(),
            LatticeType::SimpleCubic => 2.0 * radius,
        };
        Self {
            lattice,
            radius,
            lattice_parameter: a,
            n_cells,
            spheres: Vec::new(),
        }
    }
    /// Generate all sphere positions for the lattice.
    pub fn generate(&mut self) {
        self.spheres.clear();
        let a = self.lattice_parameter;
        let [nx, ny, nz] = self.n_cells;
        let basis: &[[f64; 3]] = match self.lattice {
            LatticeType::Fcc => &[
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.0],
                [0.5, 0.0, 0.5],
                [0.0, 0.5, 0.5],
            ],
            LatticeType::Hcp => &[
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.0],
                [0.0, 1.0 / 3.0, 0.5],
                [0.5, 5.0 / 6.0, 0.5],
            ],
            LatticeType::Bcc => &[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
            LatticeType::SimpleCubic => &[[0.0, 0.0, 0.0]],
        };
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    for b in basis {
                        let cx = (ix as f64 + b[0]) * a;
                        let cy = (iy as f64 + b[1]) * a;
                        let cz = (iz as f64 + b[2]) * a;
                        let id = self.spheres.len();
                        self.spheres
                            .push(PackingSphere::new([cx, cy, cz], self.radius, id));
                    }
                }
            }
        }
    }
    /// Theoretical packing efficiency for the lattice type.
    pub fn theoretical_packing_fraction(&self) -> f64 {
        match self.lattice {
            LatticeType::Fcc | LatticeType::Hcp => std::f64::consts::PI / (3.0 * 2.0_f64.sqrt()),
            LatticeType::Bcc => std::f64::consts::PI * 3.0_f64.sqrt() / 8.0,
            LatticeType::SimpleCubic => std::f64::consts::PI / 6.0,
        }
    }
    /// Coordination number (number of nearest neighbors) for the lattice.
    pub fn coordination_number(&self) -> usize {
        match self.lattice {
            LatticeType::Fcc | LatticeType::Hcp => 12,
            LatticeType::Bcc => 8,
            LatticeType::SimpleCubic => 6,
        }
    }
    /// Compute the packing fraction from generated spheres.
    pub fn actual_packing_fraction(&self) -> f64 {
        if self.spheres.is_empty() {
            return 0.0;
        }
        let [nx, ny, nz] = self.n_cells;
        let a = self.lattice_parameter;
        let volume = (nx as f64 * a) * (ny as f64 * a) * (nz as f64 * a);
        let sphere_vol: f64 = self.spheres.iter().map(|s| s.volume()).sum();
        sphere_vol / volume
    }
    /// Return the nearest-neighbor distance for the lattice.
    pub fn nearest_neighbor_distance(&self) -> f64 {
        2.0 * self.radius
    }
    /// Find neighbors within a given cutoff distance for a sphere index.
    pub fn find_neighbors(&self, idx: usize, cutoff: f64) -> Vec<usize> {
        let s = &self.spheres[idx];
        self.spheres
            .iter()
            .enumerate()
            .filter(|(i, other)| *i != idx && s.distance_to(other) <= cutoff)
            .map(|(i, _)| i)
            .collect()
    }
}
/// Particle size distribution type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SizeDistribution {
    /// Monodisperse: all spheres have the same radius.
    Monodisperse,
    /// Bidisperse: two sizes with a size ratio.
    Bidisperse,
    /// Log-normal distribution.
    LogNormal,
    /// Gaussian (truncated) distribution.
    Gaussian,
    /// Power-law distribution.
    PowerLaw,
}
/// Configuration for the force-based packing optimizer.
#[derive(Debug, Clone)]
pub struct OptimizerConfig {
    /// Number of optimization iterations.
    pub max_iterations: usize,
    /// Time step for displacement updates.
    pub time_step: f64,
    /// Convergence tolerance (max overlap).
    pub tolerance: f64,
    /// Damping coefficient to dissipate kinetic energy.
    pub damping: f64,
    /// Repulsion stiffness.
    pub stiffness: f64,
}
/// Crystal lattice type for ordered sphere packing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LatticeType {
    /// Face-centered cubic (FCC): packing fraction π/(3√2) ≈ 0.7405.
    Fcc,
    /// Hexagonally close-packed (HCP): same fraction as FCC ≈ 0.7405.
    Hcp,
    /// Body-centered cubic (BCC): packing fraction π√3/8 ≈ 0.6802.
    Bcc,
    /// Simple cubic (SC): packing fraction π/6 ≈ 0.5236.
    SimpleCubic,
}
/// Force-based packing optimizer.
///
/// Uses a soft-sphere repulsive potential to relax overlapping spheres
/// toward a mechanically stable configuration.
pub struct PackingOptimizer {
    /// Optimizer configuration.
    pub config: OptimizerConfig,
    /// Spheres being optimized.
    pub spheres: Vec<PackingSphere>,
    /// Domain bounds.
    pub domain: [f64; 3],
    /// Convergence history (max overlap at each iteration).
    pub convergence: Vec<f64>,
    /// Whether the optimizer has converged.
    pub converged: bool,
}
impl PackingOptimizer {
    /// Create a new optimizer.
    pub fn new(spheres: Vec<PackingSphere>, domain: [f64; 3], config: OptimizerConfig) -> Self {
        Self {
            config,
            spheres,
            domain,
            convergence: Vec::new(),
            converged: false,
        }
    }
    /// Compute the repulsive force on sphere `i` due to sphere `j`.
    ///
    /// Linear spring force: `f = k * overlap * normal`.
    fn pair_force(&self, i: usize, j: usize) -> [f64; 3] {
        let si = &self.spheres[i];
        let sj = &self.spheres[j];
        let dx = si.center[0] - sj.center[0];
        let dy = si.center[1] - sj.center[1];
        let dz = si.center[2] - sj.center[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        let min_dist = si.radius + sj.radius;
        if dist < min_dist && dist > 1e-30 {
            let overlap = min_dist - dist;
            let k = self.config.stiffness;
            let f = k * overlap / dist;
            [f * dx, f * dy, f * dz]
        } else {
            [0.0, 0.0, 0.0]
        }
    }
    /// Compute wall repulsion force to keep sphere inside domain.
    fn wall_force(&self, i: usize) -> [f64; 3] {
        let s = &self.spheres[i];
        let k = self.config.stiffness;
        let [lx, ly, lz] = self.domain;
        let mut f = [0.0f64; 3];
        let dims = [
            (s.center[0], 0.0, lx),
            (s.center[1], 0.0, ly),
            (s.center[2], 0.0, lz),
        ];
        for (d, (c, lo, hi)) in dims.iter().enumerate() {
            if c - s.radius < *lo {
                f[d] += k * (lo + s.radius - c);
            }
            if c + s.radius > *hi {
                f[d] -= k * (c + s.radius - hi);
            }
        }
        f
    }
    /// Perform one optimization step.
    ///
    /// Returns the maximum overlap after the step.
    pub fn step(&mut self) -> f64 {
        let n = self.spheres.len();
        let mut forces = vec![[0.0f64; 3]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let f = self.pair_force(i, j);
                for d in 0..3 {
                    forces[i][d] += f[d];
                    forces[j][d] -= f[d];
                }
            }
            let wf = self.wall_force(i);
            for d in 0..3 {
                forces[i][d] += wf[d];
            }
        }
        let dt = self.config.time_step;
        let damp = self.config.damping;
        for (i, s) in self.spheres.iter_mut().enumerate() {
            for d in 0..3 {
                s.center[d] += dt * (1.0 - damp) * forces[i][d];
            }
        }
        let mut max_overlap = 0.0f64;
        for i in 0..n {
            for j in (i + 1)..n {
                let ov = self.spheres[i].overlap_with(&self.spheres[j]);
                max_overlap = max_overlap.max(ov);
            }
        }
        max_overlap
    }
    /// Run the optimizer until convergence or max iterations.
    pub fn run(&mut self) {
        self.convergence.clear();
        for _ in 0..self.config.max_iterations {
            let max_ov = self.step();
            self.convergence.push(max_ov);
            if max_ov < self.config.tolerance {
                self.converged = true;
                break;
            }
        }
    }
    /// Compute the total potential energy of the packing.
    pub fn potential_energy(&self) -> f64 {
        let n = self.spheres.len();
        let k = self.config.stiffness;
        let mut energy = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let ov = self.spheres[i].overlap_with(&self.spheres[j]);
                if ov > 0.0 {
                    energy += 0.5 * k * ov * ov;
                }
            }
        }
        energy
    }
    /// Compute the packing fraction after optimization.
    pub fn packing_fraction(&self) -> f64 {
        let [lx, ly, lz] = self.domain;
        let vol = lx * ly * lz;
        let sv: f64 = self.spheres.iter().map(|s| s.volume()).sum();
        sv / vol
    }
    /// Return the maximum overlap in the current configuration.
    pub fn max_overlap(&self) -> f64 {
        let n = self.spheres.len();
        let mut max_ov = 0.0f64;
        for i in 0..n {
            for j in (i + 1)..n {
                let ov = self.spheres[i].overlap_with(&self.spheres[j]);
                max_ov = max_ov.max(ov);
            }
        }
        max_ov
    }
}
/// Shape of the confining container.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConfinementShape {
    /// Rectangular box with given half-extents.
    Box,
    /// Cylinder with axis along Z.
    Cylinder,
    /// Sphere container.
    SphericalContainer,
}
/// A sphere in 3D space, described by its center and radius.
#[derive(Debug, Clone, PartialEq)]
pub struct PackingSphere {
    /// Center of the sphere.
    pub center: [f64; 3],
    /// Radius of the sphere.
    pub radius: f64,
    /// Unique identifier.
    pub id: usize,
}
impl PackingSphere {
    /// Create a new sphere.
    pub fn new(center: [f64; 3], radius: f64, id: usize) -> Self {
        Self { center, radius, id }
    }
    /// Compute the volume of the sphere.
    pub fn volume(&self) -> f64 {
        4.0 / 3.0 * std::f64::consts::PI * self.radius.powi(3)
    }
    /// Check if this sphere overlaps with another sphere.
    pub fn overlaps(&self, other: &PackingSphere) -> bool {
        let min_dist = self.radius + other.radius;
        let dx = self.center[0] - other.center[0];
        let dy = self.center[1] - other.center[1];
        let dz = self.center[2] - other.center[2];
        dx * dx + dy * dy + dz * dz < min_dist * min_dist
    }
    /// Compute center-to-center distance.
    pub fn distance_to(&self, other: &PackingSphere) -> f64 {
        let dx = self.center[0] - other.center[0];
        let dy = self.center[1] - other.center[1];
        let dz = self.center[2] - other.center[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    /// Compute overlap (interpenetration) with another sphere.
    ///
    /// Returns the overlap distance (positive = overlap, negative = separation).
    pub fn overlap_with(&self, other: &PackingSphere) -> f64 {
        let min_dist = self.radius + other.radius;
        let dist = self.distance_to(other);
        min_dist - dist
    }
}
