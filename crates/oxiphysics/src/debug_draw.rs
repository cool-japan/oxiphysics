// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Renderer-agnostic debug draw command buffer.
//!
//! Physics engines often need to visualise internal state (contact points,
//! AABB trees, constraint axes) without coupling to any specific renderer.
//! This module provides a lightweight, serialisable command buffer that a
//! renderer can consume each frame.
//!
//! ## Design
//!
//! - `DrawCommand` — an enum of primitive draw operations.
//! - `DrawList` — a growable buffer of commands with `Single` / `Persistent`
//!   lifetime semantics.
//!   `Single` commands are flushed by `DrawList::drain_single`; `Persistent`
//!   commands survive until `DrawList::clear` is called.
//! - `DebugDrawSession` — a thin wrapper that ties draw calls to simulation
//!   steps and produces a `DrawnFrame` per step.
//!
//! ## Quick Start
//!
//! ```rust,no_run
//! use oxiphysics::debug_draw::{DebugDrawSession, DrawDuration};
//!
//! let mut session = DebugDrawSession::new();
//!
//! // --- simulation loop ---
//! for step in 0..10 {
//!     session.begin_step(step as u64);
//!
//!     // Transient contact point visualisation
//!     session.sphere([0.0, 0.0, 0.0], 0.05, [1.0, 0.0, 0.0, 1.0]);
//!
//!     // Persistent gravity vector (rendered every frame until cleared)
//!     session.list.add_arrow(
//!         [0.0, 2.0, 0.0], [0.0, 0.0, 0.0],
//!         0.1, [0.5, 0.5, 1.0, 1.0], DrawDuration::Persistent,
//!     );
//!
//!     let _frame = session.end_step(); // send to renderer
//! }
//! ```

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

/// RGBA colour in linear [0, 1] range.
pub type DrawColor = [f32; 4];

/// How long a [`DrawCommand`] persists in a [`DrawList`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum DrawDuration {
    /// The command is removed from the list by the next call to
    /// [`DrawList::drain_single`].
    Single,
    /// The command remains in the list until [`DrawList::clear`] is called or
    /// the command is otherwise removed.
    Persistent,
}

/// A single renderable debug primitive.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DrawCommand {
    /// A line segment between two world-space points.
    Line {
        start: [f64; 3],
        end: [f64; 3],
        color: DrawColor,
        duration: DrawDuration,
    },
    /// A wire-frame sphere centred at a world-space position.
    Sphere {
        center: [f64; 3],
        radius: f64,
        color: DrawColor,
        duration: DrawDuration,
    },
    /// A wire-frame axis-aligned bounding box.
    Aabb {
        min: [f64; 3],
        max: [f64; 3],
        color: DrawColor,
        duration: DrawDuration,
    },
    /// An arrow from `from` to `to` with a cone head of the given size.
    Arrow {
        from: [f64; 3],
        to: [f64; 3],
        head_size: f64,
        color: DrawColor,
        duration: DrawDuration,
    },
    /// Three crossing line segments centred at a world-space point.
    Cross {
        center: [f64; 3],
        size: f64,
        color: DrawColor,
        duration: DrawDuration,
    },
    /// A world-space text label.
    Text {
        position: [f64; 3],
        label: String,
        color: DrawColor,
        duration: DrawDuration,
    },
}

impl DrawCommand {
    /// Return the duration of this command.
    pub fn duration(&self) -> &DrawDuration {
        match self {
            DrawCommand::Line { duration, .. } => duration,
            DrawCommand::Sphere { duration, .. } => duration,
            DrawCommand::Aabb { duration, .. } => duration,
            DrawCommand::Arrow { duration, .. } => duration,
            DrawCommand::Cross { duration, .. } => duration,
            DrawCommand::Text { duration, .. } => duration,
        }
    }

    /// Return `true` if this command has [`DrawDuration::Single`] lifetime.
    pub fn is_single(&self) -> bool {
        *self.duration() == DrawDuration::Single
    }
}

// ---------------------------------------------------------------------------
// DrawList
// ---------------------------------------------------------------------------

/// A growable buffer of [`DrawCommand`]s with `Single` / `Persistent` lifetime
/// semantics.
///
/// Renderers iterate over the list or call [`DrawList::drain_single`] to flush
/// the per-frame commands.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct DrawList {
    /// The underlying command buffer. Public to allow direct iteration.
    pub commands: Vec<DrawCommand>,
}

impl DrawList {
    /// Create an empty list.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a raw command.
    pub fn push(&mut self, cmd: DrawCommand) {
        self.commands.push(cmd);
    }

    /// Add a line segment.
    pub fn add_line(
        &mut self,
        start: [f64; 3],
        end: [f64; 3],
        color: DrawColor,
        duration: DrawDuration,
    ) {
        self.commands.push(DrawCommand::Line {
            start,
            end,
            color,
            duration,
        });
    }

    /// Add a wire-frame sphere.
    pub fn add_sphere(
        &mut self,
        center: [f64; 3],
        radius: f64,
        color: DrawColor,
        duration: DrawDuration,
    ) {
        self.commands.push(DrawCommand::Sphere {
            center,
            radius,
            color,
            duration,
        });
    }

    /// Add a wire-frame AABB.
    pub fn add_aabb(
        &mut self,
        min: [f64; 3],
        max: [f64; 3],
        color: DrawColor,
        duration: DrawDuration,
    ) {
        self.commands.push(DrawCommand::Aabb {
            min,
            max,
            color,
            duration,
        });
    }

    /// Add an arrow.
    pub fn add_arrow(
        &mut self,
        from: [f64; 3],
        to: [f64; 3],
        head_size: f64,
        color: DrawColor,
        duration: DrawDuration,
    ) {
        self.commands.push(DrawCommand::Arrow {
            from,
            to,
            head_size,
            color,
            duration,
        });
    }

    /// Add a cross (three perpendicular line segments).
    pub fn add_cross(
        &mut self,
        center: [f64; 3],
        size: f64,
        color: DrawColor,
        duration: DrawDuration,
    ) {
        self.commands.push(DrawCommand::Cross {
            center,
            size,
            color,
            duration,
        });
    }

    /// Add a text label.
    pub fn add_text(
        &mut self,
        position: [f64; 3],
        label: impl Into<String>,
        color: DrawColor,
        duration: DrawDuration,
    ) {
        self.commands.push(DrawCommand::Text {
            position,
            label: label.into(),
            color,
            duration,
        });
    }

    /// Remove all [`DrawDuration::Single`] commands from the list and return
    /// them.  [`DrawDuration::Persistent`] commands remain.
    pub fn drain_single(&mut self) -> Vec<DrawCommand> {
        let mut singles = Vec::new();
        let mut persistent = Vec::new();
        for cmd in self.commands.drain(..) {
            if cmd.is_single() {
                singles.push(cmd);
            } else {
                persistent.push(cmd);
            }
        }
        self.commands = persistent;
        singles
    }

    /// Iterate over all commands without consuming them.
    pub fn iter(&self) -> impl Iterator<Item = &DrawCommand> {
        self.commands.iter()
    }

    /// Remove all commands (both `Single` and `Persistent`).
    pub fn clear(&mut self) {
        self.commands.clear();
    }

    /// Number of commands currently in the list.
    pub fn len(&self) -> usize {
        self.commands.len()
    }

    /// Return `true` if the list contains no commands.
    pub fn is_empty(&self) -> bool {
        self.commands.is_empty()
    }
}

// ---------------------------------------------------------------------------
// DrawnFrame
// ---------------------------------------------------------------------------

/// A snapshot of [`DrawCommand`]s associated with a single simulation step.
///
/// Produced by [`DebugDrawSession::end_step`]. Contains only the `Single`-
/// duration commands that were active during the step; `Persistent` commands
/// remain in the [`DebugDrawSession::list`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DrawnFrame {
    /// Simulation step counter this frame was captured at.
    pub step: u64,
    /// Commands active during this step (`Single`-duration only).
    pub commands: Vec<DrawCommand>,
}

// ---------------------------------------------------------------------------
// DebugDrawSession
// ---------------------------------------------------------------------------

/// A drawing session that ties debug draw calls to simulation steps.
///
/// ```rust,no_run
/// use oxiphysics::debug_draw::{DebugDrawSession, DrawDuration};
///
/// let mut session = DebugDrawSession::new();
/// session.begin_step(42);
/// session.sphere([0.0, 1.0, 0.0], 0.5, [1.0, 0.0, 0.0, 1.0]);
/// let frame = session.end_step();
/// assert_eq!(frame.step, 42);
/// assert_eq!(frame.commands.len(), 1);
/// ```
pub struct DebugDrawSession {
    /// The underlying draw list.  Callers may add `Persistent` commands here
    /// directly via `list.add_*` methods.
    pub list: DrawList,
    step: u64,
}

impl DebugDrawSession {
    /// Create a new session at step 0.
    pub fn new() -> Self {
        Self {
            list: DrawList::new(),
            step: 0,
        }
    }

    /// Advance to the given simulation step.
    ///
    /// `Persistent` commands from prior steps continue to accumulate in
    /// `self.list`.  `Single` commands are consumed by [`end_step`][Self::end_step].
    pub fn begin_step(&mut self, step: u64) {
        self.step = step;
    }

    /// Flush `Single`-duration commands into a [`DrawnFrame`] and return it.
    ///
    /// `Persistent` commands remain in `self.list` for subsequent steps.
    pub fn end_step(&mut self) -> DrawnFrame {
        let single_cmds = self.list.drain_single();
        DrawnFrame {
            step: self.step,
            commands: single_cmds,
        }
    }

    /// Add a `Single`-duration line segment.
    pub fn line(&mut self, start: [f64; 3], end: [f64; 3], color: DrawColor) {
        self.list.add_line(start, end, color, DrawDuration::Single);
    }

    /// Add a `Single`-duration wire-frame sphere.
    pub fn sphere(&mut self, center: [f64; 3], radius: f64, color: DrawColor) {
        self.list
            .add_sphere(center, radius, color, DrawDuration::Single);
    }

    /// Add a `Single`-duration wire-frame AABB.
    pub fn aabb(&mut self, min: [f64; 3], max: [f64; 3], color: DrawColor) {
        self.list.add_aabb(min, max, color, DrawDuration::Single);
    }

    /// Add a `Single`-duration arrow.
    pub fn arrow(&mut self, from: [f64; 3], to: [f64; 3], head_size: f64, color: DrawColor) {
        self.list
            .add_arrow(from, to, head_size, color, DrawDuration::Single);
    }

    /// Add a `Single`-duration cross.
    pub fn cross(&mut self, center: [f64; 3], size: f64, color: DrawColor) {
        self.list
            .add_cross(center, size, color, DrawDuration::Single);
    }

    /// Add a `Single`-duration text label.
    pub fn text(&mut self, position: [f64; 3], label: impl Into<String>, color: DrawColor) {
        self.list
            .add_text(position, label, color, DrawDuration::Single);
    }

    /// Current simulation step counter.
    pub fn step(&self) -> u64 {
        self.step
    }
}

impl Default for DebugDrawSession {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Colour constants
// ---------------------------------------------------------------------------

/// Opaque red.
pub const RED: DrawColor = [1.0, 0.0, 0.0, 1.0];
/// Opaque green.
pub const GREEN: DrawColor = [0.0, 1.0, 0.0, 1.0];
/// Opaque blue.
pub const BLUE: DrawColor = [0.0, 0.0, 1.0, 1.0];
/// Opaque yellow.
pub const YELLOW: DrawColor = [1.0, 1.0, 0.0, 1.0];
/// Opaque cyan.
pub const CYAN: DrawColor = [0.0, 1.0, 1.0, 1.0];
/// Opaque magenta.
pub const MAGENTA: DrawColor = [1.0, 0.0, 1.0, 1.0];
/// Opaque white.
pub const WHITE: DrawColor = [1.0, 1.0, 1.0, 1.0];
/// Opaque orange.
pub const ORANGE: DrawColor = [1.0, 0.5, 0.0, 1.0];
