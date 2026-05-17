// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Physics event publish/subscribe system.
//!
//! The `EventBus` collects physics events emitted during a simulation step
//! and dispatches them to registered listeners.  Events are queued during the
//! step and either flushed (dispatching to all listeners) or drained (moved
//! into a `Vec` for manual inspection).
//!
//! ## Event types
//!
//! | Variant                  | Fired when …                                    |
//! |--------------------------|-------------------------------------------------|
//! | `CollisionBegin`       | Two bodies start overlapping                    |
//! | `CollisionEnd`         | Two bodies separate                             |
//! | `BodySleep`            | A body transitions to the sleeping state        |
//! | `BodyWake`             | A sleeping body is woken                        |
//! | `ConstraintBreak`      | A breakable constraint reaches its force limit  |
//! | `BodyEnterRegion`      | A body enters a tagged sensor region            |
//! | `BodyLeaveRegion`      | A body leaves a tagged sensor region            |
//! | `StepCompleted`        | A simulation step has finished                  |
//! | `UserDefined`          | Application-defined custom event                |
//!
//!
//!
//!
//!
//!
//!
//!
//!
//!
//!
//! ## Example
//!
//! ```rust
//! use oxiphysics::event_bus::{EventBus, PhysicsEvent};
//!
//! let mut bus = EventBus::new();
//!
//! // Register a listener that counts collisions
//! let id = bus.subscribe(|ev| {
//!     if matches!(ev, PhysicsEvent::CollisionBegin { .. }) {
//!         // handle collision
//!     }
//! });
//!
//! // Emit some events (normally done by the pipeline)
//! bus.publish(PhysicsEvent::CollisionBegin {
//!     body_a: 0, body_b: 1,
//!     contact_normal: [0.0, 1.0, 0.0],
//!     penetration_depth: 0.01,
//! });
//! bus.publish(PhysicsEvent::StepCompleted { step: 1, sim_time: 1.0 / 60.0, dt: 1.0 / 60.0 });
//!
//! // Flush: deliver all queued events to listeners
//! bus.flush();
//!
//! // Or drain: take ownership without calling listeners
//! let events = bus.drain();
//! assert!(events.is_empty()); // already flushed above
//!
//! bus.unsubscribe(id);
//! ```

#![allow(missing_docs)]
#![allow(dead_code)]

// ============================================================================
// PhysicsEvent
// ============================================================================

/// Contact information bundled with a [`CollisionBegin`] event.
///
/// [`CollisionBegin`]: PhysicsEvent::CollisionBegin
#[derive(Debug, Clone, PartialEq)]
pub struct ContactInfo {
    /// Index of the first colliding body.
    pub body_a: usize,
    /// Index of the second colliding body.
    pub body_b: usize,
    /// World-space contact normal pointing from A toward B.
    pub contact_normal: [f64; 3],
    /// Penetration depth at the contact point \[m\].
    pub penetration_depth: f64,
    /// World-space contact point.
    pub contact_point: [f64; 3],
}

/// All events that the physics engine may emit during a simulation step.
#[derive(Debug, Clone, PartialEq)]
pub enum PhysicsEvent {
    /// Two bodies have begun overlapping.
    CollisionBegin {
        body_a: usize,
        body_b: usize,
        contact_normal: [f64; 3],
        penetration_depth: f64,
    },

    /// Two previously overlapping bodies have separated.
    CollisionEnd { body_a: usize, body_b: usize },

    /// A dynamic body has entered the sleeping state.
    BodySleep { body: usize },

    /// A sleeping body has been woken (by collision, force, or API call).
    BodyWake { body: usize },

    /// A breakable constraint has exceeded its maximum force and been removed.
    ConstraintBreak {
        /// User-assigned constraint identifier.
        constraint_id: usize,
        /// The accumulated impulse magnitude that triggered the break \[N·s\].
        accumulated_impulse: f64,
    },

    /// A body has entered a named sensor region.
    BodyEnterRegion {
        body: usize,
        /// User-assigned region identifier.
        region_id: usize,
    },

    /// A body has left a named sensor region.
    BodyLeaveRegion { body: usize, region_id: usize },

    /// A simulation step has completed.
    StepCompleted {
        /// Step counter (0-based).
        step: u64,
        /// Accumulated simulation time at end of this step.
        sim_time: f64,
        /// Time step duration used.
        dt: f64,
    },

    /// Application-defined event with arbitrary payload.
    UserDefined {
        /// Application-assigned event type ID.
        id: u32,
        /// Arbitrary binary payload.
        data: Vec<u8>,
    },
}

impl PhysicsEvent {
    /// Returns a short human-readable name for this event type.
    pub fn kind_name(&self) -> &'static str {
        match self {
            PhysicsEvent::CollisionBegin { .. } => "CollisionBegin",
            PhysicsEvent::CollisionEnd { .. } => "CollisionEnd",
            PhysicsEvent::BodySleep { .. } => "BodySleep",
            PhysicsEvent::BodyWake { .. } => "BodyWake",
            PhysicsEvent::ConstraintBreak { .. } => "ConstraintBreak",
            PhysicsEvent::BodyEnterRegion { .. } => "BodyEnterRegion",
            PhysicsEvent::BodyLeaveRegion { .. } => "BodyLeaveRegion",
            PhysicsEvent::StepCompleted { .. } => "StepCompleted",
            PhysicsEvent::UserDefined { .. } => "UserDefined",
        }
    }

    /// Returns `true` if this is a collision event (begin or end).
    pub fn is_collision(&self) -> bool {
        matches!(
            self,
            PhysicsEvent::CollisionBegin { .. } | PhysicsEvent::CollisionEnd { .. }
        )
    }
}

impl std::fmt::Display for PhysicsEvent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PhysicsEvent::CollisionBegin {
                body_a,
                body_b,
                penetration_depth,
                ..
            } => write!(
                f,
                "CollisionBegin({body_a}↔{body_b}, depth={penetration_depth:.4}m)"
            ),
            PhysicsEvent::CollisionEnd { body_a, body_b } => {
                write!(f, "CollisionEnd({body_a}↔{body_b})")
            }
            PhysicsEvent::BodySleep { body } => write!(f, "BodySleep(body={body})"),
            PhysicsEvent::BodyWake { body } => write!(f, "BodyWake(body={body})"),
            PhysicsEvent::ConstraintBreak {
                constraint_id,
                accumulated_impulse,
            } => write!(
                f,
                "ConstraintBreak(id={constraint_id}, impulse={accumulated_impulse:.2}N·s)"
            ),
            PhysicsEvent::BodyEnterRegion { body, region_id } => {
                write!(f, "BodyEnterRegion(body={body}, region={region_id})")
            }
            PhysicsEvent::BodyLeaveRegion { body, region_id } => {
                write!(f, "BodyLeaveRegion(body={body}, region={region_id})")
            }
            PhysicsEvent::StepCompleted { step, sim_time, dt } => write!(
                f,
                "StepCompleted(step={step}, t={sim_time:.4}s, dt={dt:.6}s)"
            ),
            PhysicsEvent::UserDefined { id, data } => {
                write!(f, "UserDefined(id={id}, {n}B)", n = data.len())
            }
        }
    }
}

// ============================================================================
// ListenerId
// ============================================================================

/// Opaque handle returned by [`EventBus::subscribe`] for later unsubscription.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ListenerId(u64);

// ============================================================================
// EventBus
// ============================================================================

/// A boxed listener closure stored in the event bus.
type ListenerFn = Box<dyn Fn(&PhysicsEvent) + Send + 'static>;

/// A lightweight publish/subscribe event bus for physics events.
///
/// Events are published to an internal queue with [`publish`][EventBus::publish].
/// They are delivered to listeners either:
/// - On [`flush`][EventBus::flush] — calls every listener with each event.
/// - On [`drain`][EventBus::drain] — moves the events into the caller's `Vec`
///   without invoking any listeners.
///
/// Listeners are closures stored as `Box<dyn Fn(&PhysicsEvent)>` and are called
/// in registration order.
pub struct EventBus {
    /// Pending events not yet flushed.
    queue: Vec<PhysicsEvent>,
    /// Registered listener closures with their IDs.
    listeners: Vec<(ListenerId, ListenerFn)>,
    /// Monotonically increasing ID counter.
    next_id: u64,
    /// Total number of events published since creation.
    published_count: u64,
    /// Total number of listener invocations since creation.
    dispatch_count: u64,
}

impl Default for EventBus {
    fn default() -> Self {
        Self::new()
    }
}

impl EventBus {
    /// Create a new, empty event bus.
    pub fn new() -> Self {
        Self {
            queue: Vec::new(),
            listeners: Vec::new(),
            next_id: 0,
            published_count: 0,
            dispatch_count: 0,
        }
    }

    /// Publish an event to the queue.
    ///
    /// The event is not dispatched immediately — call [`flush`][Self::flush]
    /// or [`drain`][Self::drain] to process it.
    pub fn publish(&mut self, event: PhysicsEvent) {
        self.published_count += 1;
        self.queue.push(event);
    }

    /// Publish multiple events at once.
    pub fn publish_all(&mut self, events: impl IntoIterator<Item = PhysicsEvent>) {
        for e in events {
            self.publish(e);
        }
    }

    /// Register a listener closure.
    ///
    /// Returns a [`ListenerId`] that can be passed to [`unsubscribe`][Self::unsubscribe]
    /// to remove the listener.
    pub fn subscribe(&mut self, f: impl Fn(&PhysicsEvent) + Send + 'static) -> ListenerId {
        let id = ListenerId(self.next_id);
        self.next_id += 1;
        self.listeners.push((id, Box::new(f)));
        id
    }

    /// Register a listener that only receives events matching a predicate.
    pub fn subscribe_filtered(
        &mut self,
        predicate: impl Fn(&PhysicsEvent) -> bool + Send + 'static,
        handler: impl Fn(&PhysicsEvent) + Send + 'static,
    ) -> ListenerId {
        self.subscribe(move |ev| {
            if predicate(ev) {
                handler(ev);
            }
        })
    }

    /// Unsubscribe a listener by its ID.
    ///
    /// Returns `true` if the listener was found and removed.
    pub fn unsubscribe(&mut self, id: ListenerId) -> bool {
        if let Some(pos) = self.listeners.iter().position(|(lid, _)| *lid == id) {
            let _ = self.listeners.remove(pos);
            true
        } else {
            false
        }
    }

    /// Dispatch all queued events to every registered listener, then clear the queue.
    ///
    /// Listeners are called in registration order for each event in queue order.
    pub fn flush(&mut self) {
        for event in &self.queue {
            for (_, listener) in &self.listeners {
                listener(event);
                self.dispatch_count += 1;
            }
        }
        self.queue.clear();
    }

    /// Move all queued events out of the bus without calling listeners.
    ///
    /// Returns the events in FIFO order.  The queue is cleared after this call.
    pub fn drain(&mut self) -> Vec<PhysicsEvent> {
        std::mem::take(&mut self.queue)
    }

    /// Peek at the current queue without consuming or dispatching it.
    pub fn pending(&self) -> &[PhysicsEvent] {
        &self.queue
    }

    /// Number of events currently in the queue.
    pub fn pending_count(&self) -> usize {
        self.queue.len()
    }

    /// Clear the queue without dispatching or returning events.
    pub fn clear(&mut self) {
        self.queue.clear();
    }

    /// Number of currently registered listeners.
    pub fn listener_count(&self) -> usize {
        self.listeners.len()
    }

    /// Total events published since the bus was created.
    pub fn total_published(&self) -> u64 {
        self.published_count
    }

    /// Total listener invocations since the bus was created.
    pub fn total_dispatched(&self) -> u64 {
        self.dispatch_count
    }

    /// Count of events in the queue that satisfy `predicate`.
    pub fn count_pending<F>(&self, predicate: F) -> usize
    where
        F: Fn(&PhysicsEvent) -> bool,
    {
        self.queue.iter().filter(|e| predicate(e)).count()
    }
}

impl std::fmt::Debug for EventBus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("EventBus")
            .field("pending_count", &self.queue.len())
            .field("listener_count", &self.listeners.len())
            .field("published_count", &self.published_count)
            .field("dispatch_count", &self.dispatch_count)
            .finish()
    }
}

// ============================================================================
// Convenience: filter helpers
// ============================================================================

/// A helper that collects all [`CollisionBegin`] events into a `Vec`.
///
/// [`CollisionBegin`]: PhysicsEvent::CollisionBegin
pub fn collect_collision_begins(events: &[PhysicsEvent]) -> Vec<(usize, usize)> {
    events
        .iter()
        .filter_map(|e| {
            if let PhysicsEvent::CollisionBegin { body_a, body_b, .. } = e {
                Some((*body_a, *body_b))
            } else {
                None
            }
        })
        .collect()
}

/// A helper that extracts all bodies that entered or left any region in this batch.
pub fn region_crossing_events(events: &[PhysicsEvent]) -> Vec<(usize, usize, bool)> {
    // returns (body, region, is_enter)
    events
        .iter()
        .filter_map(|e| match e {
            PhysicsEvent::BodyEnterRegion { body, region_id } => Some((*body, *region_id, true)),
            PhysicsEvent::BodyLeaveRegion { body, region_id } => Some((*body, *region_id, false)),
            _ => None,
        })
        .collect()
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{Arc, Mutex};

    #[test]
    fn publish_and_drain() {
        let mut bus = EventBus::new();
        bus.publish(PhysicsEvent::BodySleep { body: 0 });
        bus.publish(PhysicsEvent::BodyWake { body: 1 });
        assert_eq!(bus.pending_count(), 2);
        let events = bus.drain();
        assert_eq!(events.len(), 2);
        assert_eq!(bus.pending_count(), 0);
    }

    #[test]
    fn flush_calls_listeners() {
        let mut bus = EventBus::new();
        let counter = Arc::new(Mutex::new(0u32));
        let c = Arc::clone(&counter);
        bus.subscribe(move |_| {
            *c.lock().unwrap_or_else(|e| e.into_inner()) += 1;
        });

        bus.publish(PhysicsEvent::BodySleep { body: 0 });
        bus.publish(PhysicsEvent::BodySleep { body: 1 });
        bus.flush();

        assert_eq!(*counter.lock().unwrap_or_else(|e| e.into_inner()), 2);
        assert_eq!(bus.pending_count(), 0);
    }

    #[test]
    fn unsubscribe_stops_delivery() {
        let mut bus = EventBus::new();
        let counter = Arc::new(Mutex::new(0u32));
        let c = Arc::clone(&counter);
        let id = bus.subscribe(move |_| {
            *c.lock().unwrap_or_else(|e| e.into_inner()) += 1;
        });
        bus.publish(PhysicsEvent::BodySleep { body: 0 });
        bus.flush();
        assert_eq!(*counter.lock().unwrap_or_else(|e| e.into_inner()), 1);

        bus.unsubscribe(id);
        bus.publish(PhysicsEvent::BodySleep { body: 0 });
        bus.flush();
        assert_eq!(*counter.lock().unwrap_or_else(|e| e.into_inner()), 1); // unchanged
    }

    #[test]
    fn filtered_subscription() {
        let mut bus = EventBus::new();
        let hits = Arc::new(Mutex::new(Vec::new()));
        let h = Arc::clone(&hits);
        bus.subscribe_filtered(
            |ev| matches!(ev, PhysicsEvent::CollisionBegin { .. }),
            move |ev| {
                h.lock()
                    .unwrap_or_else(|e| e.into_inner())
                    .push(ev.kind_name());
            },
        );
        bus.publish(PhysicsEvent::CollisionBegin {
            body_a: 0,
            body_b: 1,
            contact_normal: [0.0, 1.0, 0.0],
            penetration_depth: 0.1,
        });
        bus.publish(PhysicsEvent::BodySleep { body: 0 });
        bus.flush();
        let hits = hits.lock().unwrap_or_else(|e| e.into_inner());
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0], "CollisionBegin");
    }

    #[test]
    fn multiple_listeners_order() {
        let mut bus = EventBus::new();
        let log = Arc::new(Mutex::new(Vec::new()));
        for i in 0..3u32 {
            let l = Arc::clone(&log);
            bus.subscribe(move |_| {
                l.lock().unwrap_or_else(|e| e.into_inner()).push(i);
            });
        }
        bus.publish(PhysicsEvent::BodyWake { body: 0 });
        bus.flush();
        assert_eq!(
            *log.lock().unwrap_or_else(|e| e.into_inner()),
            vec![0, 1, 2]
        );
    }

    #[test]
    fn clear_discards_queue() {
        let mut bus = EventBus::new();
        bus.publish(PhysicsEvent::BodySleep { body: 0 });
        bus.clear();
        assert_eq!(bus.pending_count(), 0);
    }

    #[test]
    fn event_display() {
        let ev = PhysicsEvent::CollisionBegin {
            body_a: 3,
            body_b: 7,
            contact_normal: [0.0, 1.0, 0.0],
            penetration_depth: 0.005,
        };
        let s = format!("{}", ev);
        assert!(s.contains("CollisionBegin"));
        assert!(s.contains('3'));
        assert!(s.contains('7'));
    }

    #[test]
    fn collect_collision_begins_helper() {
        let events = vec![
            PhysicsEvent::CollisionBegin {
                body_a: 0,
                body_b: 1,
                contact_normal: [0., 1., 0.],
                penetration_depth: 0.01,
            },
            PhysicsEvent::BodySleep { body: 2 },
            PhysicsEvent::CollisionBegin {
                body_a: 3,
                body_b: 4,
                contact_normal: [1., 0., 0.],
                penetration_depth: 0.02,
            },
        ];
        let pairs = collect_collision_begins(&events);
        assert_eq!(pairs, vec![(0, 1), (3, 4)]);
    }

    #[test]
    fn total_stats() {
        let mut bus = EventBus::new();
        bus.subscribe(|_| {});
        bus.publish(PhysicsEvent::BodySleep { body: 0 });
        bus.publish(PhysicsEvent::BodyWake { body: 0 });
        bus.flush();
        assert_eq!(bus.total_published(), 2);
        assert_eq!(bus.total_dispatched(), 2);
    }
}
