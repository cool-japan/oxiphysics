// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ring-buffer physics logger with severity and category filtering.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::to_js_value;

/// Log level.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum LogLevel {
    /// Error — always shown.
    Error = 0,
    /// Warning.
    Warn = 1,
    /// Informational.
    Info = 2,
    /// Verbose debug output.
    Debug = 3,
}

/// A single log entry.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LogEntry {
    /// Severity.
    pub level: LogLevel,
    /// Category tag.
    #[wasm_bindgen(skip)]
    pub category: String,
    /// Log message.
    #[wasm_bindgen(skip)]
    pub message: String,
    /// Simulation time when logged.
    pub time: f64,
}

impl LogEntry {
    fn new(
        level: LogLevel,
        category: impl Into<String>,
        message: impl Into<String>,
        time: f64,
    ) -> Self {
        Self {
            level,
            category: category.into(),
            message: message.into(),
            time,
        }
    }
}

#[wasm_bindgen]
impl LogEntry {
    /// Category tag for the entry.
    #[wasm_bindgen(getter)]
    pub fn category(&self) -> String {
        self.category.clone()
    }

    /// Log message text.
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }

    /// Serialise to a `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

/// Ring-buffer logger with category filtering.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmPhysicsLogger {
    /// Ring buffer of entries.
    #[wasm_bindgen(skip)]
    pub entries: Vec<LogEntry>,
    /// Ring buffer capacity.
    pub capacity: usize,
    /// Minimum level to log.
    pub min_level: LogLevel,
    /// Optional category filter (`None` = all categories).
    #[wasm_bindgen(skip)]
    pub category_filter: Option<String>,
}

impl WasmPhysicsLogger {
    /// Create a logger.
    pub fn new(capacity: usize, min_level: LogLevel) -> Self {
        Self {
            entries: Vec::new(),
            capacity,
            min_level,
            category_filter: None,
        }
    }

    /// Filter to a specific category.
    pub fn set_category_filter(&mut self, category: impl Into<String>) {
        self.category_filter = Some(category.into());
    }

    /// Remove category filter.
    pub fn clear_filter(&mut self) {
        self.category_filter = None;
    }

    /// Log a message.
    pub fn log(
        &mut self,
        level: LogLevel,
        category: impl Into<String>,
        message: impl Into<String>,
        time: f64,
    ) {
        if level > self.min_level {
            return;
        }
        let cat: String = category.into();
        if let Some(ref f) = self.category_filter
            && &cat != f
        {
            return;
        }
        if self.capacity > 0 && self.entries.len() >= self.capacity {
            self.entries.remove(0);
        }
        self.entries.push(LogEntry::new(level, cat, message, time));
    }

    /// Shorthand error log.
    pub fn error(&mut self, cat: impl Into<String>, msg: impl Into<String>, time: f64) {
        self.log(LogLevel::Error, cat, msg, time);
    }

    /// Shorthand warn log.
    pub fn warn(&mut self, cat: impl Into<String>, msg: impl Into<String>, time: f64) {
        self.log(LogLevel::Warn, cat, msg, time);
    }

    /// Shorthand info log.
    pub fn info(&mut self, cat: impl Into<String>, msg: impl Into<String>, time: f64) {
        self.log(LogLevel::Info, cat, msg, time);
    }

    /// Shorthand debug log.
    pub fn debug(&mut self, cat: impl Into<String>, msg: impl Into<String>, time: f64) {
        self.log(LogLevel::Debug, cat, msg, time);
    }

    /// Entries at or above a given level.
    pub fn entries_at_level(&self, level: LogLevel) -> Vec<&LogEntry> {
        self.entries.iter().filter(|e| e.level <= level).collect()
    }

    /// Flush all entries.
    pub fn flush(&mut self) -> Vec<LogEntry> {
        std::mem::take(&mut self.entries)
    }
}

#[wasm_bindgen]
impl WasmPhysicsLogger {
    /// Create a logger (JS-compatible constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(capacity: usize, min_level: LogLevel) -> WasmPhysicsLogger {
        WasmPhysicsLogger::new(capacity, min_level)
    }

    /// Filter logging to a single category.
    #[wasm_bindgen(js_name = "set_category_filter")]
    pub fn set_category_filter_js(&mut self, category: String) {
        self.set_category_filter(category);
    }

    /// Remove the category filter.
    #[wasm_bindgen(js_name = "clear_filter")]
    pub fn clear_filter_js(&mut self) {
        self.clear_filter();
    }

    /// True if a category filter is currently set.
    #[wasm_bindgen(js_name = "has_category_filter")]
    pub fn has_category_filter_js(&self) -> bool {
        self.category_filter.is_some()
    }

    /// Get the current category filter, or an empty string when no filter is
    /// set. Use [`WasmPhysicsLogger::has_category_filter_js`] to disambiguate
    /// "filter set to empty" from "no filter".
    #[wasm_bindgen(js_name = "get_category_filter")]
    pub fn get_category_filter_js(&self) -> String {
        self.category_filter.clone().unwrap_or_default()
    }

    /// Log a message with the given severity, category, and simulation time.
    #[wasm_bindgen(js_name = "log")]
    pub fn log_js(&mut self, level: LogLevel, category: String, message: String, time: f64) {
        self.log(level, category, message, time);
    }

    /// Shorthand error log.
    #[wasm_bindgen(js_name = "error")]
    pub fn error_js(&mut self, category: String, message: String, time: f64) {
        self.error(category, message, time);
    }

    /// Shorthand warn log.
    #[wasm_bindgen(js_name = "warn")]
    pub fn warn_js(&mut self, category: String, message: String, time: f64) {
        self.warn(category, message, time);
    }

    /// Shorthand info log.
    #[wasm_bindgen(js_name = "info")]
    pub fn info_js(&mut self, category: String, message: String, time: f64) {
        self.info(category, message, time);
    }

    /// Shorthand debug log.
    #[wasm_bindgen(js_name = "debug")]
    pub fn debug_js(&mut self, category: String, message: String, time: f64) {
        self.debug(category, message, time);
    }

    /// Number of buffered entries.
    #[wasm_bindgen(js_name = "len")]
    pub fn len_js(&self) -> usize {
        self.entries.len()
    }

    /// Drain all entries from the buffer and return them as a `JsValue` array.
    #[wasm_bindgen(js_name = "flush")]
    pub fn flush_js(&mut self) -> Result<JsValue, JsValue> {
        let drained = self.flush();
        to_js_value(&drained)
    }

    /// Snapshot all current entries as a `JsValue` array, leaving the buffer
    /// untouched.
    #[wasm_bindgen(js_name = "entries_to_js_value")]
    pub fn entries_to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(&self.entries)
    }
}
