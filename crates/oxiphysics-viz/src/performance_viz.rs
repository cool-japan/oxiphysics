// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Real-time performance visualization for the OxiPhysics engine.
//!
//! Provides CPU-side data structures and rendering helpers for:
//! - FPS counter and frame-time ring buffer
//! - Physics timing breakdown (collision / integration / constraints)
//! - Profiler flame-graph data
//! - Memory usage bar graphs
//! - Parallel thread-timeline strips
//! - GPU utilization graphs
//! - Benchmark comparison charts
//! - Spatial heat-maps of hot spots
//! - Iteration-convergence display
//! - Adaptive-timestep history

// ─────────────────────────────────────────────────────────────────────────────
// 1. Primitive rendering types (no GPU dependency)
// ─────────────────────────────────────────────────────────────────────────────

/// An RGBA colour packed into four `u8` channels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rgba {
    /// Red channel.
    pub r: u8,
    /// Green channel.
    pub g: u8,
    /// Blue channel.
    pub b: u8,
    /// Alpha channel (255 = fully opaque).
    pub a: u8,
}

impl Rgba {
    /// Construct from individual channels.
    #[must_use]
    pub const fn new(r: u8, g: u8, b: u8, a: u8) -> Self {
        Self { r, g, b, a }
    }

    /// Fully-opaque white.
    pub const WHITE: Self = Self::new(255, 255, 255, 255);
    /// Fully-opaque black.
    pub const BLACK: Self = Self::new(0, 0, 0, 255);
    /// Fully-opaque red.
    pub const RED: Self = Self::new(220, 50, 50, 255);
    /// Fully-opaque green.
    pub const GREEN: Self = Self::new(50, 200, 80, 255);
    /// Fully-opaque blue.
    pub const BLUE: Self = Self::new(50, 100, 220, 255);
    /// Fully-opaque yellow.
    pub const YELLOW: Self = Self::new(240, 200, 40, 255);
    /// Fully-opaque orange.
    pub const ORANGE: Self = Self::new(240, 140, 40, 255);
    /// Semi-transparent dark background.
    pub const PANEL_BG: Self = Self::new(20, 20, 30, 200);

    /// Linear interpolation between two colours.
    #[must_use]
    pub fn lerp(a: Self, b: Self, t: f32) -> Self {
        let t = t.clamp(0.0, 1.0);
        Self::new(
            lerp_u8(a.r, b.r, t),
            lerp_u8(a.g, b.g, t),
            lerp_u8(a.b, b.b, t),
            lerp_u8(a.a, b.a, t),
        )
    }
}

fn lerp_u8(a: u8, b: u8, t: f32) -> u8 {
    (a as f32 + (b as f32 - a as f32) * t) as u8
}

/// A 2-D integer rectangle used for widget layout.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rect {
    /// Left edge (pixels from left of screen).
    pub x: i32,
    /// Top edge (pixels from top of screen).
    pub y: i32,
    /// Width in pixels.
    pub w: i32,
    /// Height in pixels.
    pub h: i32,
}

impl Rect {
    /// Construct a rectangle.
    #[must_use]
    pub const fn new(x: i32, y: i32, w: i32, h: i32) -> Self {
        Self { x, y, w, h }
    }

    /// Return `true` if the pixel `(px, py)` lies inside the rectangle.
    #[must_use]
    pub fn contains_pixel(&self, px: i32, py: i32) -> bool {
        px >= self.x && px < self.x + self.w && py >= self.y && py < self.y + self.h
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 2. Ring buffer
// ─────────────────────────────────────────────────────────────────────────────

/// A fixed-capacity ring buffer for storing the last `N` `f64` samples.
#[derive(Debug, Clone)]
pub struct RingBuffer {
    data: Vec<f64>,
    head: usize,
    len: usize,
    capacity: usize,
}

impl RingBuffer {
    /// Create a new ring buffer with the given capacity (all zeros).
    #[must_use]
    pub fn new(capacity: usize) -> Self {
        assert!(capacity > 0, "capacity must be > 0");
        Self {
            data: vec![0.0; capacity],
            head: 0,
            len: 0,
            capacity,
        }
    }

    /// Push a new sample, overwriting the oldest if full.
    pub fn push(&mut self, value: f64) {
        self.data[self.head] = value;
        self.head = (self.head + 1) % self.capacity;
        if self.len < self.capacity {
            self.len += 1;
        }
    }

    /// Return the number of valid samples.
    #[must_use]
    pub fn len(&self) -> usize {
        self.len
    }

    /// Return `true` if no samples have been pushed yet.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Iterate samples oldest-first.
    #[must_use]
    pub fn iter_oldest_first(&self) -> Vec<f64> {
        if self.len == 0 {
            return Vec::new();
        }
        let start = if self.len == self.capacity {
            self.head
        } else {
            0
        };
        (0..self.len)
            .map(|i| self.data[(start + i) % self.capacity])
            .collect()
    }

    /// Return the arithmetic mean of stored samples.
    #[must_use]
    pub fn mean(&self) -> f64 {
        if self.len == 0 {
            return 0.0;
        }
        self.iter_oldest_first().iter().sum::<f64>() / self.len as f64
    }

    /// Return the maximum stored value.
    #[must_use]
    pub fn max(&self) -> f64 {
        self.iter_oldest_first()
            .into_iter()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Return the minimum stored value.
    #[must_use]
    pub fn min(&self) -> f64 {
        self.iter_oldest_first()
            .into_iter()
            .fold(f64::INFINITY, f64::min)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. FPS counter and frame-time graph
// ─────────────────────────────────────────────────────────────────────────────

/// Tracks frames-per-second and per-frame timings.
#[derive(Debug, Clone)]
pub struct FpsCounter {
    /// Ring buffer of recent frame durations [seconds].
    frame_times: RingBuffer,
    /// Exponential moving average of FPS.
    ema_fps: f64,
    /// EMA smoothing factor (0 < α ≤ 1).
    alpha: f64,
}

impl FpsCounter {
    /// Construct a new `FpsCounter` with a 120-sample history.
    #[must_use]
    pub fn new() -> Self {
        Self {
            frame_times: RingBuffer::new(120),
            ema_fps: 0.0,
            alpha: 0.1,
        }
    }

    /// Record a completed frame of duration `dt_s` seconds.
    pub fn record_frame(&mut self, dt_s: f64) {
        if dt_s > 0.0 {
            self.frame_times.push(dt_s);
            let instant_fps = 1.0 / dt_s;
            if self.ema_fps == 0.0 {
                self.ema_fps = instant_fps;
            } else {
                self.ema_fps = self.alpha * instant_fps + (1.0 - self.alpha) * self.ema_fps;
            }
        }
    }

    /// Smoothed FPS estimate.
    #[must_use]
    pub fn fps(&self) -> f64 {
        self.ema_fps
    }

    /// Mean frame time over the history window \[ms\].
    #[must_use]
    pub fn mean_frame_time_ms(&self) -> f64 {
        self.frame_times.mean() * 1000.0
    }

    /// Peak (worst) frame time over the history window \[ms\].
    #[must_use]
    pub fn peak_frame_time_ms(&self) -> f64 {
        self.frame_times.max() * 1000.0
    }

    /// Return the stored frame-time samples \[ms\], oldest first.
    #[must_use]
    pub fn frame_time_history_ms(&self) -> Vec<f64> {
        self.frame_times
            .iter_oldest_first()
            .into_iter()
            .map(|t| t * 1000.0)
            .collect()
    }

    /// Return a normalised sparkline of frame times suitable for a bar chart.
    ///
    /// Each value is in \[0, 1\] where 1 maps to the peak frame time.
    #[must_use]
    pub fn sparkline(&self) -> Vec<f64> {
        let samples = self.frame_time_history_ms();
        let peak = samples.iter().cloned().fold(0.0_f64, f64::max);
        if peak == 0.0 {
            return vec![0.0; samples.len()];
        }
        samples.into_iter().map(|v| v / peak).collect()
    }

    /// Format the current FPS and frame-time into a compact HUD string.
    #[must_use]
    pub fn hud_string(&self) -> String {
        format!(
            "FPS: {:.1}  Frame: {:.2}ms  Peak: {:.2}ms",
            self.fps(),
            self.mean_frame_time_ms(),
            self.peak_frame_time_ms()
        )
    }
}

impl Default for FpsCounter {
    fn default() -> Self {
        Self::new()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 4. Physics timing breakdown
// ─────────────────────────────────────────────────────────────────────────────

/// Timing breakdown for a single physics step.
#[derive(Debug, Clone, Copy, Default)]
pub struct PhysicsTimings {
    /// Time spent in broad-phase collision detection \[ms\].
    pub broadphase_ms: f64,
    /// Time spent in narrow-phase collision detection \[ms\].
    pub narrowphase_ms: f64,
    /// Time spent solving constraints \[ms\].
    pub constraints_ms: f64,
    /// Time spent integrating equations of motion \[ms\].
    pub integration_ms: f64,
    /// Time spent in all other sub-systems \[ms\].
    pub other_ms: f64,
}

impl PhysicsTimings {
    /// Construct from individual sub-system durations.
    #[must_use]
    pub fn new(
        broadphase_ms: f64,
        narrowphase_ms: f64,
        constraints_ms: f64,
        integration_ms: f64,
        other_ms: f64,
    ) -> Self {
        Self {
            broadphase_ms,
            narrowphase_ms,
            constraints_ms,
            integration_ms,
            other_ms,
        }
    }

    /// Total physics step time \[ms\].
    #[must_use]
    pub fn total_ms(&self) -> f64 {
        self.broadphase_ms
            + self.narrowphase_ms
            + self.constraints_ms
            + self.integration_ms
            + self.other_ms
    }

    /// Fractions of each stage as a proportion of the total.
    ///
    /// Returns `[broadphase, narrowphase, constraints, integration, other]`.
    #[must_use]
    pub fn fractions(&self) -> [f64; 5] {
        let t = self.total_ms();
        if t == 0.0 {
            return [0.0; 5];
        }
        [
            self.broadphase_ms / t,
            self.narrowphase_ms / t,
            self.constraints_ms / t,
            self.integration_ms / t,
            self.other_ms / t,
        ]
    }

    /// Return a coloured bar-chart segment description `(label, fraction, colour)`.
    #[must_use]
    pub fn bar_segments(&self) -> Vec<(&'static str, f64, Rgba)> {
        let f = self.fractions();
        vec![
            ("Broad", f[0], Rgba::BLUE),
            ("Narrow", f[1], Rgba::GREEN),
            ("Constr", f[2], Rgba::ORANGE),
            ("Integr", f[3], Rgba::RED),
            ("Other", f[4], Rgba::new(150, 150, 150, 255)),
        ]
    }
}

/// Accumulates per-step `PhysicsTimings` into ring buffers for display.
#[derive(Debug, Clone)]
pub struct PhysicsTimingHistory {
    /// Ring buffers for each sub-system (same order as `PhysicsTimings` fields).
    buffers: [RingBuffer; 5],
}

impl PhysicsTimingHistory {
    /// Create a new history with the given ring-buffer capacity.
    #[must_use]
    pub fn new(capacity: usize) -> Self {
        Self {
            buffers: [
                RingBuffer::new(capacity),
                RingBuffer::new(capacity),
                RingBuffer::new(capacity),
                RingBuffer::new(capacity),
                RingBuffer::new(capacity),
            ],
        }
    }

    /// Record a `PhysicsTimings` sample.
    pub fn record(&mut self, t: &PhysicsTimings) {
        self.buffers[0].push(t.broadphase_ms);
        self.buffers[1].push(t.narrowphase_ms);
        self.buffers[2].push(t.constraints_ms);
        self.buffers[3].push(t.integration_ms);
        self.buffers[4].push(t.other_ms);
    }

    /// Mean values \[ms\] over the stored history.
    #[must_use]
    pub fn means_ms(&self) -> [f64; 5] {
        [
            self.buffers[0].mean(),
            self.buffers[1].mean(),
            self.buffers[2].mean(),
            self.buffers[3].mean(),
            self.buffers[4].mean(),
        ]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 5. Profiler flame-graph data
// ─────────────────────────────────────────────────────────────────────────────

/// A single span in a flame graph.
#[derive(Debug, Clone)]
pub struct FlameSpan {
    /// Function or scope name.
    pub label: String,
    /// Start time within the frame \[ms\].
    pub start_ms: f64,
    /// Duration \[ms\].
    pub duration_ms: f64,
    /// Nesting depth (0 = top level).
    pub depth: u32,
    /// Display colour.
    pub colour: Rgba,
}

impl FlameSpan {
    /// Construct a new flame span.
    #[must_use]
    pub fn new(
        label: impl Into<String>,
        start_ms: f64,
        duration_ms: f64,
        depth: u32,
        colour: Rgba,
    ) -> Self {
        Self {
            label: label.into(),
            start_ms,
            duration_ms,
            depth,
            colour,
        }
    }

    /// End time of this span \[ms\].
    #[must_use]
    pub fn end_ms(&self) -> f64 {
        self.start_ms + self.duration_ms
    }

    /// Fraction of `frame_ms` this span occupies.
    #[must_use]
    pub fn width_fraction(&self, frame_ms: f64) -> f64 {
        if frame_ms == 0.0 {
            0.0
        } else {
            self.duration_ms / frame_ms
        }
    }
}

/// A complete flame-graph for one captured frame.
#[derive(Debug, Clone, Default)]
pub struct FlameGraph {
    /// All spans recorded in this frame.
    pub spans: Vec<FlameSpan>,
    /// Total frame duration \[ms\].
    pub frame_ms: f64,
}

impl FlameGraph {
    /// Create an empty flame graph.
    #[must_use]
    pub fn new(frame_ms: f64) -> Self {
        Self {
            spans: Vec::new(),
            frame_ms,
        }
    }

    /// Add a span.
    pub fn push(&mut self, span: FlameSpan) {
        self.spans.push(span);
    }

    /// Return spans at a specific depth level.
    #[must_use]
    pub fn spans_at_depth(&self, depth: u32) -> Vec<&FlameSpan> {
        self.spans.iter().filter(|s| s.depth == depth).collect()
    }

    /// Maximum nesting depth present in the graph.
    #[must_use]
    pub fn max_depth(&self) -> u32 {
        self.spans.iter().map(|s| s.depth).max().unwrap_or(0)
    }

    /// Return spans that contain time `t_ms` (i.e. t ∈ \[start, end\]).
    #[must_use]
    pub fn spans_at_time(&self, t_ms: f64) -> Vec<&FlameSpan> {
        self.spans
            .iter()
            .filter(|s| t_ms >= s.start_ms && t_ms <= s.end_ms())
            .collect()
    }

    /// Total time attributed to a label (sum of all matching spans).
    #[must_use]
    pub fn total_for_label(&self, label: &str) -> f64 {
        self.spans
            .iter()
            .filter(|s| s.label == label)
            .map(|s| s.duration_ms)
            .sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 6. Memory usage visualization
// ─────────────────────────────────────────────────────────────────────────────

/// A named memory allocation category.
#[derive(Debug, Clone)]
pub struct MemoryCategory {
    /// Human-readable category name.
    pub name: String,
    /// Bytes currently allocated.
    pub bytes: usize,
    /// Peak bytes allocated (high-water mark).
    pub peak_bytes: usize,
    /// Bar colour for the visualizer.
    pub colour: Rgba,
}

impl MemoryCategory {
    /// Construct a new category.
    #[must_use]
    pub fn new(name: impl Into<String>, bytes: usize, colour: Rgba) -> Self {
        Self {
            name: name.into(),
            bytes,
            peak_bytes: bytes,
            colour,
        }
    }

    /// Update the allocation count; raises the peak if necessary.
    pub fn update(&mut self, bytes: usize) {
        self.bytes = bytes;
        if bytes > self.peak_bytes {
            self.peak_bytes = bytes;
        }
    }

    /// Return the current allocation in MiB.
    #[must_use]
    pub fn mib(&self) -> f64 {
        self.bytes as f64 / (1024.0 * 1024.0)
    }

    /// Return the peak allocation in MiB.
    #[must_use]
    pub fn peak_mib(&self) -> f64 {
        self.peak_bytes as f64 / (1024.0 * 1024.0)
    }
}

/// A collection of memory categories with a total budget.
#[derive(Debug, Clone)]
pub struct MemoryUsagePanel {
    /// Registered categories.
    pub categories: Vec<MemoryCategory>,
    /// Total memory budget in bytes.
    pub budget_bytes: usize,
}

impl MemoryUsagePanel {
    /// Create a new panel with a given memory budget.
    #[must_use]
    pub fn new(budget_bytes: usize) -> Self {
        Self {
            categories: Vec::new(),
            budget_bytes,
        }
    }

    /// Register a new category and return its index.
    pub fn add_category(&mut self, cat: MemoryCategory) -> usize {
        let idx = self.categories.len();
        self.categories.push(cat);
        idx
    }

    /// Total bytes currently allocated across all categories.
    #[must_use]
    pub fn total_bytes(&self) -> usize {
        self.categories.iter().map(|c| c.bytes).sum()
    }

    /// Fraction of the budget consumed (0 … 1+).
    #[must_use]
    pub fn budget_fraction(&self) -> f64 {
        if self.budget_bytes == 0 {
            1.0
        } else {
            self.total_bytes() as f64 / self.budget_bytes as f64
        }
    }

    /// Format a compact status string.
    #[must_use]
    pub fn status_string(&self) -> String {
        let total_mib = self.total_bytes() as f64 / (1024.0 * 1024.0);
        let budget_mib = self.budget_bytes as f64 / (1024.0 * 1024.0);
        format!(
            "Mem: {:.1}/{:.1} MiB ({:.0}%)",
            total_mib,
            budget_mib,
            self.budget_fraction() * 100.0
        )
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 7. Parallel workload thread timeline
// ─────────────────────────────────────────────────────────────────────────────

/// A work item executed on one thread during one frame.
#[derive(Debug, Clone)]
pub struct ThreadWorkItem {
    /// Thread index (0 = main thread).
    pub thread_id: u32,
    /// Label describing the work.
    pub label: String,
    /// Start time within the frame \[ms\].
    pub start_ms: f64,
    /// Duration \[ms\].
    pub duration_ms: f64,
    /// Display colour.
    pub colour: Rgba,
}

impl ThreadWorkItem {
    /// Construct a work item.
    #[must_use]
    pub fn new(
        thread_id: u32,
        label: impl Into<String>,
        start_ms: f64,
        duration_ms: f64,
        colour: Rgba,
    ) -> Self {
        Self {
            thread_id,
            label: label.into(),
            start_ms,
            duration_ms,
            colour,
        }
    }

    /// End time of this work item \[ms\].
    #[must_use]
    pub fn end_ms(&self) -> f64 {
        self.start_ms + self.duration_ms
    }
}

/// A thread-timeline panel showing parallel workloads.
#[derive(Debug, Clone, Default)]
pub struct ThreadTimeline {
    /// All work items recorded for this frame.
    pub items: Vec<ThreadWorkItem>,
    /// Total number of worker threads.
    pub thread_count: u32,
    /// Frame duration \[ms\].
    pub frame_ms: f64,
}

impl ThreadTimeline {
    /// Construct a new timeline.
    #[must_use]
    pub fn new(thread_count: u32, frame_ms: f64) -> Self {
        Self {
            items: Vec::new(),
            thread_count,
            frame_ms,
        }
    }

    /// Add a work item.
    pub fn push(&mut self, item: ThreadWorkItem) {
        self.items.push(item);
    }

    /// Return work items for a specific thread, sorted by start time.
    #[must_use]
    pub fn items_for_thread(&self, thread_id: u32) -> Vec<&ThreadWorkItem> {
        let mut v: Vec<&ThreadWorkItem> = self
            .items
            .iter()
            .filter(|i| i.thread_id == thread_id)
            .collect();
        v.sort_by(|a, b| {
            a.start_ms
                .partial_cmp(&b.start_ms)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        v
    }

    /// Compute the utilisation fraction for each thread (work_time / frame_ms).
    #[must_use]
    pub fn thread_utilisation(&self) -> Vec<f64> {
        (0..self.thread_count)
            .map(|tid| {
                let work: f64 = self
                    .items_for_thread(tid)
                    .iter()
                    .map(|i| i.duration_ms)
                    .sum();
                if self.frame_ms == 0.0 {
                    0.0
                } else {
                    work / self.frame_ms
                }
            })
            .collect()
    }

    /// Average utilisation across all threads.
    #[must_use]
    pub fn mean_utilisation(&self) -> f64 {
        let u = self.thread_utilisation();
        if u.is_empty() {
            0.0
        } else {
            u.iter().sum::<f64>() / u.len() as f64
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 8. GPU utilization graph
// ─────────────────────────────────────────────────────────────────────────────

/// A single GPU utilization sample.
#[derive(Debug, Clone, Copy)]
pub struct GpuSample {
    /// Fractional GPU compute utilisation (0 … 1).
    pub compute: f64,
    /// Fractional GPU memory bandwidth utilisation (0 … 1).
    pub memory_bw: f64,
    /// GPU temperature in degrees Celsius.
    pub temp_c: f64,
    /// GPU power draw in watts.
    pub power_w: f64,
}

impl GpuSample {
    /// Construct a GPU sample.
    #[must_use]
    pub fn new(compute: f64, memory_bw: f64, temp_c: f64, power_w: f64) -> Self {
        Self {
            compute,
            memory_bw,
            temp_c,
            power_w,
        }
    }
}

/// Tracks GPU utilisation over time with a ring buffer per metric.
#[derive(Debug, Clone)]
pub struct GpuGraph {
    /// Compute utilisation history.
    pub compute: RingBuffer,
    /// Memory bandwidth utilisation history.
    pub memory_bw: RingBuffer,
    /// Temperature history \[°C\].
    pub temp: RingBuffer,
    /// Power history \[W\].
    pub power: RingBuffer,
}

impl GpuGraph {
    /// Create a new `GpuGraph` with the given history length.
    #[must_use]
    pub fn new(capacity: usize) -> Self {
        Self {
            compute: RingBuffer::new(capacity),
            memory_bw: RingBuffer::new(capacity),
            temp: RingBuffer::new(capacity),
            power: RingBuffer::new(capacity),
        }
    }

    /// Record a GPU sample.
    pub fn record(&mut self, s: GpuSample) {
        self.compute.push(s.compute);
        self.memory_bw.push(s.memory_bw);
        self.temp.push(s.temp_c);
        self.power.push(s.power_w);
    }

    /// Mean compute utilisation over the history window.
    #[must_use]
    pub fn mean_compute(&self) -> f64 {
        self.compute.mean()
    }

    /// Peak temperature \[°C\].
    #[must_use]
    pub fn peak_temp(&self) -> f64 {
        self.temp.max()
    }

    /// Format a compact status line.
    #[must_use]
    pub fn status_string(&self) -> String {
        format!(
            "GPU: {:.0}% compute  BW: {:.0}%  Temp: {:.0}°C  Power: {:.0}W",
            self.compute.mean() * 100.0,
            self.memory_bw.mean() * 100.0,
            self.temp.max(),
            self.power.mean()
        )
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 9. Benchmark comparison chart
// ─────────────────────────────────────────────────────────────────────────────

/// A single benchmark result.
#[derive(Debug, Clone)]
pub struct BenchmarkResult {
    /// Test case label.
    pub label: String,
    /// Mean execution time \[ms\].
    pub mean_ms: f64,
    /// Standard deviation \[ms\].
    pub std_ms: f64,
    /// Minimum observed time \[ms\].
    pub min_ms: f64,
    /// Maximum observed time \[ms\].
    pub max_ms: f64,
    /// Bar colour.
    pub colour: Rgba,
}

impl BenchmarkResult {
    /// Construct a benchmark result by computing statistics from raw samples.
    #[must_use]
    pub fn from_samples(label: impl Into<String>, samples: &[f64], colour: Rgba) -> Self {
        let n = samples.len() as f64;
        let mean_ms = if n > 0.0 {
            samples.iter().sum::<f64>() / n
        } else {
            0.0
        };
        let var = if n > 1.0 {
            samples.iter().map(|x| (x - mean_ms).powi(2)).sum::<f64>() / (n - 1.0)
        } else {
            0.0
        };
        let std_ms = var.sqrt();
        let min_ms = samples.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_ms = samples.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        Self {
            label: label.into(),
            mean_ms,
            std_ms,
            min_ms,
            max_ms,
            colour,
        }
    }

    /// Coefficient of variation (σ/μ).
    #[must_use]
    pub fn cv(&self) -> f64 {
        if self.mean_ms == 0.0 {
            0.0
        } else {
            self.std_ms / self.mean_ms
        }
    }

    /// Format a one-line summary.
    #[must_use]
    pub fn summary(&self) -> String {
        format!(
            "{}: {:.2}±{:.2}ms  [{:.2}, {:.2}]",
            self.label, self.mean_ms, self.std_ms, self.min_ms, self.max_ms
        )
    }
}

/// A comparison chart showing multiple benchmark results.
#[derive(Debug, Clone, Default)]
pub struct BenchmarkChart {
    /// Individual benchmark results.
    pub results: Vec<BenchmarkResult>,
    /// Chart title.
    pub title: String,
}

impl BenchmarkChart {
    /// Create a new chart with a title.
    #[must_use]
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            results: Vec::new(),
            title: title.into(),
        }
    }

    /// Add a result.
    pub fn add(&mut self, result: BenchmarkResult) {
        self.results.push(result);
    }

    /// Return the index of the fastest (lowest mean) result.
    #[must_use]
    pub fn fastest_idx(&self) -> Option<usize> {
        self.results
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                a.mean_ms
                    .partial_cmp(&b.mean_ms)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
    }

    /// Normalised mean times (relative to the fastest = 1.0).
    #[must_use]
    pub fn relative_means(&self) -> Vec<f64> {
        let fastest = self
            .results
            .iter()
            .map(|r| r.mean_ms)
            .fold(f64::INFINITY, f64::min);
        if fastest == 0.0 {
            return vec![0.0; self.results.len()];
        }
        self.results.iter().map(|r| r.mean_ms / fastest).collect()
    }

    /// Format the chart as a plain-text table.
    #[must_use]
    pub fn to_text_table(&self) -> String {
        let mut lines = vec![self.title.clone()];
        lines.push("-".repeat(60));
        for r in &self.results {
            lines.push(r.summary());
        }
        lines.join("\n")
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 10. Heat map of hot spots
// ─────────────────────────────────────────────────────────────────────────────

/// A 2-D heat-map grid where each cell stores an accumulated "heat" value.
#[derive(Debug, Clone)]
pub struct HeatMap {
    /// Number of columns.
    pub cols: usize,
    /// Number of rows.
    pub rows: usize,
    /// Row-major heat values (non-negative).
    pub data: Vec<f64>,
    /// World-space bounding box `[min_x, min_y, max_x, max_y]`.
    pub bounds: [f64; 4],
}

impl HeatMap {
    /// Construct a zero-filled heat map.
    #[must_use]
    pub fn new(cols: usize, rows: usize, bounds: [f64; 4]) -> Self {
        Self {
            cols,
            rows,
            data: vec![0.0; cols * rows],
            bounds,
        }
    }

    /// Splat `amount` of heat at world position `(wx, wy)` with a Gaussian
    /// kernel of radius `r` cells.
    pub fn splat(&mut self, wx: f64, wy: f64, amount: f64, r: f64) {
        let (col0, row0) = self.world_to_cell(wx, wy);
        let r_cells = r.ceil() as i64 + 1;
        for dr in -r_cells..=r_cells {
            for dc in -r_cells..=r_cells {
                let c = col0 as i64 + dc;
                let r_idx = row0 as i64 + dr;
                if c < 0 || r_idx < 0 || c >= self.cols as i64 || r_idx >= self.rows as i64 {
                    continue;
                }
                let dist2 = (dc as f64).powi(2) + (dr as f64).powi(2);
                let w = (-dist2 / (2.0 * r * r)).exp();
                self.data[r_idx as usize * self.cols + c as usize] += amount * w;
            }
        }
    }

    fn world_to_cell(&self, wx: f64, wy: f64) -> (usize, usize) {
        let fx = (wx - self.bounds[0]) / (self.bounds[2] - self.bounds[0]) * self.cols as f64;
        let fy = (wy - self.bounds[1]) / (self.bounds[3] - self.bounds[1]) * self.rows as f64;
        let col = (fx as usize).min(self.cols - 1);
        let row = (fy as usize).min(self.rows - 1);
        (col, row)
    }

    /// Return the maximum heat value.
    #[must_use]
    pub fn max_heat(&self) -> f64 {
        self.data.iter().cloned().fold(0.0_f64, f64::max)
    }

    /// Normalise all values to \[0, 1\] by dividing by the maximum.
    pub fn normalise(&mut self) {
        let m = self.max_heat();
        if m > 0.0 {
            self.data.iter_mut().for_each(|v| *v /= m);
        }
    }

    /// Map a normalised heat value to an RGBA colour using a hot-cold palette.
    #[must_use]
    pub fn heat_to_colour(v: f64) -> Rgba {
        let v = v.clamp(0.0, 1.0);
        if v < 0.5 {
            Rgba::lerp(Rgba::BLUE, Rgba::YELLOW, (v * 2.0) as f32)
        } else {
            Rgba::lerp(Rgba::YELLOW, Rgba::RED, ((v - 0.5) * 2.0) as f32)
        }
    }

    /// Convert the heat map to a flat RGBA image buffer (row-major).
    #[must_use]
    pub fn to_rgba_image(&self) -> Vec<Rgba> {
        let peak = self.max_heat();
        self.data
            .iter()
            .map(|&v| {
                let norm = if peak > 0.0 { v / peak } else { 0.0 };
                Self::heat_to_colour(norm)
            })
            .collect()
    }

    /// Zero out all heat values.
    pub fn clear(&mut self) {
        self.data.iter_mut().for_each(|v| *v = 0.0);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 11. Iteration convergence display
// ─────────────────────────────────────────────────────────────────────────────

/// Records the residual norm at each iteration of a solver.
#[derive(Debug, Clone, Default)]
pub struct ConvergencePlot {
    /// Residual norms (one per iteration).
    pub residuals: Vec<f64>,
    /// Convergence tolerance.
    pub tolerance: f64,
    /// Name of the solver.
    pub solver_name: String,
}

impl ConvergencePlot {
    /// Create a new convergence plot.
    #[must_use]
    pub fn new(solver_name: impl Into<String>, tolerance: f64) -> Self {
        Self {
            residuals: Vec::new(),
            tolerance,
            solver_name: solver_name.into(),
        }
    }

    /// Append a residual value.
    pub fn push(&mut self, residual: f64) {
        self.residuals.push(residual);
    }

    /// Return `true` if the last recorded residual is below the tolerance.
    #[must_use]
    pub fn has_converged(&self) -> bool {
        self.residuals
            .last()
            .map(|&r| r < self.tolerance)
            .unwrap_or(false)
    }

    /// Number of iterations recorded.
    #[must_use]
    pub fn iterations(&self) -> usize {
        self.residuals.len()
    }

    /// Reduction factor: ratio of first to last residual.
    #[must_use]
    pub fn reduction_factor(&self) -> f64 {
        match (self.residuals.first(), self.residuals.last()) {
            (Some(&r0), Some(&rn)) if rn > 0.0 => r0 / rn,
            _ => 1.0,
        }
    }

    /// Return log10 of each residual for a log-scale plot.
    #[must_use]
    pub fn log10_residuals(&self) -> Vec<f64> {
        self.residuals
            .iter()
            .map(|&r| {
                if r > 0.0 {
                    r.log10()
                } else {
                    f64::NEG_INFINITY
                }
            })
            .collect()
    }

    /// Format a one-line convergence summary.
    #[must_use]
    pub fn summary(&self) -> String {
        format!(
            "{}: {} iters, residual={:.2e}, converged={}",
            self.solver_name,
            self.iterations(),
            self.residuals.last().copied().unwrap_or(f64::NAN),
            self.has_converged()
        )
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 12. Adaptive timestep history
// ─────────────────────────────────────────────────────────────────────────────

/// A single adaptive-timestep decision.
#[derive(Debug, Clone, Copy)]
pub struct TimestepDecision {
    /// Simulation time at which the step was taken \[s\].
    pub sim_time_s: f64,
    /// Actual timestep used \[s\].
    pub dt_s: f64,
    /// Error estimate that triggered the step-size choice.
    pub error_estimate: f64,
    /// Whether the step was accepted (`true`) or rejected (`false`).
    pub accepted: bool,
}

impl TimestepDecision {
    /// Construct a timestep decision record.
    #[must_use]
    pub fn new(sim_time_s: f64, dt_s: f64, error_estimate: f64, accepted: bool) -> Self {
        Self {
            sim_time_s,
            dt_s,
            error_estimate,
            accepted,
        }
    }
}

/// Accumulates adaptive-timestep history for display.
#[derive(Debug, Clone, Default)]
pub struct AdaptiveTimestepHistory {
    /// All decisions in chronological order.
    pub decisions: Vec<TimestepDecision>,
    /// Maximum history length (oldest entries are dropped when exceeded).
    pub max_len: usize,
}

impl AdaptiveTimestepHistory {
    /// Create a history with a given maximum length.
    #[must_use]
    pub fn new(max_len: usize) -> Self {
        Self {
            decisions: Vec::new(),
            max_len,
        }
    }

    /// Record a timestep decision.
    pub fn push(&mut self, d: TimestepDecision) {
        if self.decisions.len() >= self.max_len {
            self.decisions.remove(0);
        }
        self.decisions.push(d);
    }

    /// Fraction of accepted steps.
    #[must_use]
    pub fn acceptance_rate(&self) -> f64 {
        if self.decisions.is_empty() {
            return 1.0;
        }
        let accepted = self.decisions.iter().filter(|d| d.accepted).count();
        accepted as f64 / self.decisions.len() as f64
    }

    /// Mean accepted timestep \[s\].
    #[must_use]
    pub fn mean_dt(&self) -> f64 {
        let v: Vec<f64> = self
            .decisions
            .iter()
            .filter(|d| d.accepted)
            .map(|d| d.dt_s)
            .collect();
        if v.is_empty() {
            0.0
        } else {
            v.iter().sum::<f64>() / v.len() as f64
        }
    }

    /// Return a `(sim_time, dt)` series for plotting accepted steps only.
    #[must_use]
    pub fn accepted_series(&self) -> Vec<(f64, f64)> {
        self.decisions
            .iter()
            .filter(|d| d.accepted)
            .map(|d| (d.sim_time_s, d.dt_s))
            .collect()
    }

    /// Format a compact status string.
    #[must_use]
    pub fn status_string(&self) -> String {
        format!(
            "dt: {:.2e}s  Accept: {:.0}%",
            self.mean_dt(),
            self.acceptance_rate() * 100.0
        )
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 13. Performance overlay compositor
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration for the performance overlay panel.
#[derive(Debug, Clone)]
pub struct OverlayConfig {
    /// Screen position and size of the overlay panel.
    pub panel_rect: Rect,
    /// Show the FPS counter.
    pub show_fps: bool,
    /// Show the physics timing breakdown bar.
    pub show_physics_timings: bool,
    /// Show the memory usage bar.
    pub show_memory: bool,
    /// Show the GPU utilisation graph.
    pub show_gpu: bool,
    /// Show the thread timeline.
    pub show_thread_timeline: bool,
    /// Show the adaptive timestep history.
    pub show_timestep_history: bool,
    /// Background colour.
    pub bg_colour: Rgba,
    /// Text colour.
    pub text_colour: Rgba,
}

impl OverlayConfig {
    /// Sensible defaults for a 300 × 400 pixel panel in the top-right corner.
    #[must_use]
    pub fn default_top_right(screen_w: i32) -> Self {
        Self {
            panel_rect: Rect::new(screen_w - 310, 10, 300, 400),
            show_fps: true,
            show_physics_timings: true,
            show_memory: true,
            show_gpu: true,
            show_thread_timeline: true,
            show_timestep_history: true,
            bg_colour: Rgba::PANEL_BG,
            text_colour: Rgba::WHITE,
        }
    }
}

/// Aggregates all performance sub-systems and produces text lines for an HUD.
#[derive(Debug, Clone)]
pub struct PerformanceOverlay {
    /// FPS counter.
    pub fps: FpsCounter,
    /// Physics timing history.
    pub physics: PhysicsTimingHistory,
    /// GPU utilisation graph.
    pub gpu: GpuGraph,
    /// Adaptive timestep history.
    pub timestep: AdaptiveTimestepHistory,
    /// Overlay display configuration.
    pub config: OverlayConfig,
}

impl PerformanceOverlay {
    /// Construct with default settings for a given screen width.
    #[must_use]
    pub fn new(screen_w: i32) -> Self {
        Self {
            fps: FpsCounter::new(),
            physics: PhysicsTimingHistory::new(120),
            gpu: GpuGraph::new(120),
            timestep: AdaptiveTimestepHistory::new(200),
            config: OverlayConfig::default_top_right(screen_w),
        }
    }

    /// Update all sub-systems with data from the most recently completed frame.
    pub fn update(
        &mut self,
        dt_s: f64,
        timings: &PhysicsTimings,
        gpu: GpuSample,
        ts_decision: TimestepDecision,
    ) {
        self.fps.record_frame(dt_s);
        self.physics.record(timings);
        self.gpu.record(gpu);
        self.timestep.push(ts_decision);
    }

    /// Generate a list of text lines to render on the HUD.
    #[must_use]
    pub fn hud_lines(&self) -> Vec<String> {
        let mut lines = Vec::new();
        if self.config.show_fps {
            lines.push(self.fps.hud_string());
        }
        if self.config.show_gpu {
            lines.push(self.gpu.status_string());
        }
        if self.config.show_physics_timings {
            let m = self.physics.means_ms();
            lines.push(format!(
                "Phys: bp={:.1} np={:.1} cs={:.1} int={:.1}ms",
                m[0], m[1], m[2], m[3]
            ));
        }
        if self.config.show_timestep_history {
            lines.push(self.timestep.status_string());
        }
        lines
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 14. Mini ASCII graph renderer
// ─────────────────────────────────────────────────────────────────────────────

/// Render a `Vec`f64` series as a fixed-width ASCII bar graph.
///
/// Each bar is one character wide; height is mapped to `height` rows.
/// Values are normalised to the maximum.
#[must_use]
pub fn ascii_bar_graph(data: &[f64], width: usize, height: usize) -> String {
    if data.is_empty() || width == 0 || height == 0 {
        return String::new();
    }
    // Downsample or pad to `width` bins
    let bins = resample(data, width);
    let max_val = bins.iter().cloned().fold(0.0_f64, f64::max);
    let bars: Vec<Vec<char>> = bins
        .iter()
        .map(|&v| {
            let filled = if max_val > 0.0 {
                ((v / max_val) * height as f64).round() as usize
            } else {
                0
            };
            (0..height)
                .map(|row| {
                    if (height - 1 - row) < filled {
                        '#'
                    } else {
                        ' '
                    }
                })
                .collect()
        })
        .collect();
    let mut out = String::new();
    for row in 0..height {
        for bar_col in bars.iter() {
            out.push(bar_col[row]);
        }
        out.push('\n');
    }
    out
}

fn resample(data: &[f64], target: usize) -> Vec<f64> {
    if data.len() == target {
        return data.to_vec();
    }
    (0..target)
        .map(|i| {
            let src = i as f64 / (target - 1).max(1) as f64 * (data.len() - 1) as f64;
            let lo = src.floor() as usize;
            let hi = (lo + 1).min(data.len() - 1);
            let t = src.fract();
            data[lo] * (1.0 - t) + data[hi] * t
        })
        .collect()
}

/// Render a convergence plot as a log10 ASCII line graph.
#[must_use]
pub fn ascii_convergence_graph(plot: &ConvergencePlot, width: usize, height: usize) -> String {
    let log_vals: Vec<f64> = plot
        .log10_residuals()
        .into_iter()
        .filter(|v| v.is_finite())
        .collect();
    ascii_bar_graph(&log_vals, width, height)
}

// ─────────────────────────────────────────────────────────────────────────────
// 15. Unit tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // --- Rgba -----------------------------------------------------------------

    #[test]
    fn test_rgba_lerp_half() {
        let c = Rgba::lerp(Rgba::BLACK, Rgba::WHITE, 0.5);
        // Allow ±1 due to f32 rounding
        assert!((c.r as i32 - 127).abs() <= 1);
    }

    #[test]
    fn test_rgba_lerp_extremes() {
        let zero = Rgba::lerp(Rgba::RED, Rgba::BLUE, 0.0);
        assert_eq!(zero, Rgba::RED);
        let one = Rgba::lerp(Rgba::RED, Rgba::BLUE, 1.0);
        assert_eq!(one, Rgba::BLUE);
    }

    #[test]
    fn test_rect_contains() {
        let r = Rect::new(10, 10, 100, 100);
        assert!(r.contains_pixel(50, 50));
        assert!(!r.contains_pixel(5, 50));
        assert!(!r.contains_pixel(50, 110));
    }

    // --- RingBuffer -----------------------------------------------------------

    #[test]
    fn test_ringbuffer_basic() {
        let mut rb = RingBuffer::new(5);
        for i in 1..=5 {
            rb.push(i as f64);
        }
        assert_eq!(rb.len(), 5);
        assert!((rb.mean() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_ringbuffer_overflow() {
        let mut rb = RingBuffer::new(3);
        for i in 0..10 {
            rb.push(i as f64);
        }
        assert_eq!(rb.len(), 3);
        // Newest three values: 7, 8, 9
        let v = rb.iter_oldest_first();
        assert_eq!(v, vec![7.0, 8.0, 9.0]);
    }

    #[test]
    fn test_ringbuffer_empty() {
        let rb = RingBuffer::new(10);
        assert!(rb.is_empty());
        assert_eq!(rb.mean(), 0.0);
    }

    #[test]
    fn test_ringbuffer_min_max() {
        let mut rb = RingBuffer::new(10);
        rb.push(3.0);
        rb.push(1.0);
        rb.push(4.0);
        assert!((rb.min() - 1.0).abs() < 1e-10);
        assert!((rb.max() - 4.0).abs() < 1e-10);
    }

    // --- FpsCounter -----------------------------------------------------------

    #[test]
    fn test_fps_counter_basic() {
        let mut fps = FpsCounter::new();
        fps.record_frame(1.0 / 60.0);
        assert!(fps.fps() > 0.0);
    }

    #[test]
    fn test_fps_mean_frame_time() {
        let mut fps = FpsCounter::new();
        fps.record_frame(0.016);
        assert!((fps.mean_frame_time_ms() - 16.0).abs() < 0.5);
    }

    #[test]
    fn test_fps_sparkline_len() {
        let mut fps = FpsCounter::new();
        for _ in 0..10 {
            fps.record_frame(0.016);
        }
        assert_eq!(fps.sparkline().len(), 10);
    }

    #[test]
    fn test_fps_hud_string() {
        let mut fps = FpsCounter::new();
        fps.record_frame(0.016_7);
        let s = fps.hud_string();
        assert!(s.contains("FPS:"));
    }

    // --- PhysicsTimings -------------------------------------------------------

    #[test]
    fn test_physics_timings_total() {
        let t = PhysicsTimings::new(1.0, 2.0, 3.0, 4.0, 0.5);
        assert!((t.total_ms() - 10.5).abs() < 1e-10);
    }

    #[test]
    fn test_physics_timings_fractions_sum() {
        let t = PhysicsTimings::new(1.0, 2.0, 3.0, 4.0, 0.5);
        let f: f64 = t.fractions().iter().sum();
        assert!((f - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_physics_timings_bar_segments_count() {
        let t = PhysicsTimings::default();
        assert_eq!(t.bar_segments().len(), 5);
    }

    // --- FlameGraph -----------------------------------------------------------

    #[test]
    fn test_flame_graph_max_depth() {
        let mut fg = FlameGraph::new(16.0);
        fg.push(FlameSpan::new("A", 0.0, 16.0, 0, Rgba::BLUE));
        fg.push(FlameSpan::new("B", 0.0, 8.0, 1, Rgba::GREEN));
        fg.push(FlameSpan::new("C", 0.0, 4.0, 2, Rgba::RED));
        assert_eq!(fg.max_depth(), 2);
    }

    #[test]
    fn test_flame_graph_spans_at_depth() {
        let mut fg = FlameGraph::new(16.0);
        fg.push(FlameSpan::new("A", 0.0, 16.0, 0, Rgba::BLUE));
        fg.push(FlameSpan::new("B", 0.0, 8.0, 1, Rgba::GREEN));
        assert_eq!(fg.spans_at_depth(0).len(), 1);
        assert_eq!(fg.spans_at_depth(1).len(), 1);
    }

    #[test]
    fn test_flame_graph_total_for_label() {
        let mut fg = FlameGraph::new(16.0);
        fg.push(FlameSpan::new("render", 0.0, 4.0, 0, Rgba::BLUE));
        fg.push(FlameSpan::new("render", 4.0, 3.0, 0, Rgba::BLUE));
        assert!((fg.total_for_label("render") - 7.0).abs() < 1e-10);
    }

    // --- MemoryUsagePanel -----------------------------------------------------

    #[test]
    fn test_memory_panel_total() {
        let mut panel = MemoryUsagePanel::new(1024 * 1024 * 512);
        panel.add_category(MemoryCategory::new("verts", 1024 * 1024, Rgba::BLUE));
        panel.add_category(MemoryCategory::new("textures", 1024 * 512, Rgba::GREEN));
        assert_eq!(panel.total_bytes(), 1024 * 1024 + 1024 * 512);
    }

    #[test]
    fn test_memory_category_peak() {
        let mut cat = MemoryCategory::new("buf", 100, Rgba::RED);
        cat.update(200);
        cat.update(150);
        assert_eq!(cat.peak_bytes, 200);
    }

    // --- ThreadTimeline -------------------------------------------------------

    #[test]
    fn test_thread_timeline_utilisation() {
        let mut tl = ThreadTimeline::new(2, 10.0);
        tl.push(ThreadWorkItem::new(0, "work", 0.0, 5.0, Rgba::BLUE));
        let u = tl.thread_utilisation();
        assert!((u[0] - 0.5).abs() < 1e-10);
        assert!((u[1]).abs() < 1e-10);
    }

    #[test]
    fn test_thread_timeline_mean_util() {
        let mut tl = ThreadTimeline::new(2, 10.0);
        tl.push(ThreadWorkItem::new(0, "A", 0.0, 10.0, Rgba::BLUE));
        tl.push(ThreadWorkItem::new(1, "B", 0.0, 10.0, Rgba::GREEN));
        assert!((tl.mean_utilisation() - 1.0).abs() < 1e-10);
    }

    // --- GpuGraph -------------------------------------------------------------

    #[test]
    fn test_gpu_graph_mean() {
        let mut g = GpuGraph::new(10);
        g.record(GpuSample::new(0.8, 0.5, 70.0, 200.0));
        g.record(GpuSample::new(0.6, 0.4, 72.0, 190.0));
        assert!((g.mean_compute() - 0.7).abs() < 1e-6);
    }

    #[test]
    fn test_gpu_peak_temp() {
        let mut g = GpuGraph::new(10);
        g.record(GpuSample::new(0.5, 0.5, 65.0, 150.0));
        g.record(GpuSample::new(0.9, 0.8, 90.0, 250.0));
        assert!((g.peak_temp() - 90.0).abs() < 1e-6);
    }

    // --- BenchmarkChart -------------------------------------------------------

    #[test]
    fn test_benchmark_from_samples_mean() {
        let samples = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let r = BenchmarkResult::from_samples("test", &samples, Rgba::BLUE);
        assert!((r.mean_ms - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_benchmark_chart_fastest() {
        let mut chart = BenchmarkChart::new("Test");
        chart.add(BenchmarkResult::from_samples(
            "slow",
            &[10.0, 12.0],
            Rgba::RED,
        ));
        chart.add(BenchmarkResult::from_samples(
            "fast",
            &[2.0, 3.0],
            Rgba::GREEN,
        ));
        assert_eq!(chart.fastest_idx(), Some(1));
    }

    #[test]
    fn test_benchmark_relative_means() {
        let mut chart = BenchmarkChart::new("Test");
        chart.add(BenchmarkResult::from_samples("A", &[4.0], Rgba::BLUE));
        chart.add(BenchmarkResult::from_samples("B", &[8.0], Rgba::RED));
        let rel = chart.relative_means();
        assert!((rel[0] - 1.0).abs() < 1e-10);
        assert!((rel[1] - 2.0).abs() < 1e-10);
    }

    // --- HeatMap --------------------------------------------------------------

    #[test]
    fn test_heatmap_splat_increases_value() {
        let mut hm = HeatMap::new(10, 10, [0.0, 0.0, 1.0, 1.0]);
        hm.splat(0.5, 0.5, 100.0, 1.0);
        assert!(hm.max_heat() > 0.0);
    }

    #[test]
    fn test_heatmap_normalise() {
        let mut hm = HeatMap::new(4, 4, [0.0, 0.0, 1.0, 1.0]);
        hm.splat(0.5, 0.5, 50.0, 2.0);
        hm.normalise();
        assert!((hm.max_heat() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_heatmap_to_rgba_image_length() {
        let hm = HeatMap::new(8, 8, [0.0, 0.0, 1.0, 1.0]);
        let img = hm.to_rgba_image();
        assert_eq!(img.len(), 64);
    }

    // --- ConvergencePlot ------------------------------------------------------

    #[test]
    fn test_convergence_not_converged() {
        let mut p = ConvergencePlot::new("cg", 1e-6);
        p.push(1.0);
        p.push(0.1);
        assert!(!p.has_converged());
    }

    #[test]
    fn test_convergence_converged() {
        let mut p = ConvergencePlot::new("cg", 1e-6);
        p.push(1e-7);
        assert!(p.has_converged());
    }

    #[test]
    fn test_convergence_reduction_factor() {
        let mut p = ConvergencePlot::new("test", 1e-6);
        p.push(1.0);
        p.push(0.01);
        assert!((p.reduction_factor() - 100.0).abs() < 1e-6);
    }

    #[test]
    fn test_convergence_log10() {
        let mut p = ConvergencePlot::new("test", 1e-6);
        p.push(100.0);
        p.push(10.0);
        let log = p.log10_residuals();
        assert!((log[0] - 2.0).abs() < 1e-6);
        assert!((log[1] - 1.0).abs() < 1e-6);
    }

    // --- AdaptiveTimestepHistory ----------------------------------------------

    #[test]
    fn test_timestep_acceptance_rate() {
        let mut h = AdaptiveTimestepHistory::new(100);
        h.push(TimestepDecision::new(0.0, 0.01, 1e-4, true));
        h.push(TimestepDecision::new(0.01, 0.005, 1e-3, false));
        assert!((h.acceptance_rate() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_timestep_mean_dt() {
        let mut h = AdaptiveTimestepHistory::new(100);
        h.push(TimestepDecision::new(0.0, 0.01, 1e-4, true));
        h.push(TimestepDecision::new(0.01, 0.02, 1e-4, true));
        assert!((h.mean_dt() - 0.015).abs() < 1e-10);
    }

    #[test]
    fn test_timestep_max_len_enforced() {
        let mut h = AdaptiveTimestepHistory::new(3);
        for i in 0..10 {
            h.push(TimestepDecision::new(i as f64, 0.01, 1e-4, true));
        }
        assert_eq!(h.decisions.len(), 3);
    }

    // --- PerformanceOverlay ---------------------------------------------------

    #[test]
    fn test_overlay_hud_lines_not_empty() {
        let mut overlay = PerformanceOverlay::new(1920);
        overlay.update(
            0.016,
            &PhysicsTimings::new(1.0, 2.0, 3.0, 4.0, 0.5),
            GpuSample::new(0.8, 0.5, 72.0, 200.0),
            TimestepDecision::new(0.0, 0.01, 1e-4, true),
        );
        let lines = overlay.hud_lines();
        assert!(!lines.is_empty());
    }

    // --- ASCII graph ----------------------------------------------------------

    #[test]
    fn test_ascii_bar_graph_dimensions() {
        let data: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let g = ascii_bar_graph(&data, 20, 8);
        let line_count = g.lines().count();
        assert_eq!(line_count, 8);
    }

    #[test]
    fn test_ascii_bar_graph_empty() {
        let g = ascii_bar_graph(&[], 10, 5);
        assert!(g.is_empty());
    }
}
