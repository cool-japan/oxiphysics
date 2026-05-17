// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Scene explorer — hierarchical body & joint enumeration.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::{err_to_jsvalue, to_js_value};

/// A node in the scene hierarchy.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SceneNode {
    /// Node (body) id.
    pub id: u32,
    /// Optional name.
    #[wasm_bindgen(skip)]
    pub name: Option<String>,
    /// Parent id (None = root).
    pub parent: Option<u32>,
    /// Child node ids.
    #[wasm_bindgen(skip)]
    pub children: Vec<u32>,
}

impl SceneNode {
    /// Create a new scene node.
    pub fn new(id: u32, name: Option<String>, parent: Option<u32>) -> Self {
        Self {
            id,
            name,
            parent,
            children: Vec::new(),
        }
    }
}

#[wasm_bindgen]
impl SceneNode {
    /// Optional name of the node, or an empty string if anonymous.
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.name.clone().unwrap_or_default()
    }

    /// True if the node has an explicit name.
    #[wasm_bindgen(js_name = "has_name")]
    pub fn has_name_js(&self) -> bool {
        self.name.is_some()
    }

    /// Number of child nodes.
    #[wasm_bindgen(js_name = "child_count")]
    pub fn child_count_js(&self) -> usize {
        self.children.len()
    }

    /// Child ids as a flat `Uint32Array`.
    #[wasm_bindgen(js_name = "children")]
    pub fn children_js(&self) -> Vec<u32> {
        self.children.clone()
    }

    /// Serialise the node to a `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

/// Enumerates and queries all scene objects.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct WasmSceneExplorer {
    /// Scene nodes.
    #[wasm_bindgen(skip)]
    pub nodes: Vec<SceneNode>,
    /// Joint records `(joint_id, body_a, body_b)`.
    #[wasm_bindgen(skip)]
    pub joints: Vec<(u32, u32, u32)>,
}

impl WasmSceneExplorer {
    /// Create an empty explorer.
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a body node.
    pub fn add_node(&mut self, id: u32, name: Option<String>, parent: Option<u32>) {
        // Link to parent
        if let Some(pid) = parent
            && let Some(p) = self.nodes.iter_mut().find(|n| n.id == pid)
        {
            p.children.push(id);
        }
        self.nodes.push(SceneNode::new(id, name, parent));
    }

    /// Register a joint.
    pub fn add_joint(&mut self, joint_id: u32, body_a: u32, body_b: u32) {
        self.joints.push((joint_id, body_a, body_b));
    }

    /// Enumerate all body ids.
    pub fn enumerate_bodies(&self) -> Vec<u32> {
        self.nodes.iter().map(|n| n.id).collect()
    }

    /// Enumerate all joint ids.
    pub fn enumerate_joints(&self) -> Vec<u32> {
        self.joints.iter().map(|&(jid, _, _)| jid).collect()
    }

    /// Get hierarchy as JSON.
    pub fn get_hierarchy(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(&self.nodes)
    }

    /// Find a body by name; returns id if found.
    pub fn find_by_name(&self, name: &str) -> Option<u32> {
        self.nodes.iter().find_map(|n| {
            if n.name.as_deref() == Some(name) {
                Some(n.id)
            } else {
                None
            }
        })
    }

    /// Root nodes (nodes with no parent).
    pub fn roots(&self) -> Vec<u32> {
        self.nodes
            .iter()
            .filter(|n| n.parent.is_none())
            .map(|n| n.id)
            .collect()
    }
}

#[wasm_bindgen]
impl WasmSceneExplorer {
    /// Create an empty explorer (JS-compatible constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmSceneExplorer {
        WasmSceneExplorer::new()
    }

    /// Register a body node. Pass an empty string for `name` to mean
    /// "anonymous", and a negative value for `parent` to mean "root".
    #[wasm_bindgen(js_name = "add_node")]
    pub fn add_node_js(&mut self, id: u32, name: String, parent: i64) {
        let name_opt = if name.is_empty() { None } else { Some(name) };
        let parent_opt = if parent < 0 {
            None
        } else {
            Some(parent as u32)
        };
        self.add_node(id, name_opt, parent_opt);
    }

    /// Register a joint between two bodies.
    #[wasm_bindgen(js_name = "add_joint")]
    pub fn add_joint_js(&mut self, joint_id: u32, body_a: u32, body_b: u32) {
        self.add_joint(joint_id, body_a, body_b);
    }

    /// Enumerate all body ids as a flat `Uint32Array`.
    #[wasm_bindgen(js_name = "enumerate_bodies")]
    pub fn enumerate_bodies_js(&self) -> Vec<u32> {
        self.enumerate_bodies()
    }

    /// Enumerate all joint ids as a flat `Uint32Array`.
    #[wasm_bindgen(js_name = "enumerate_joints")]
    pub fn enumerate_joints_js(&self) -> Vec<u32> {
        self.enumerate_joints()
    }

    /// Get the scene hierarchy as a JSON string.
    #[wasm_bindgen(js_name = "get_hierarchy")]
    pub fn get_hierarchy_js(&self) -> Result<String, JsValue> {
        self.get_hierarchy().map_err(err_to_jsvalue)
    }

    /// Find a body by name; returns `Some(id)` or `None` (undefined in JS).
    #[wasm_bindgen(js_name = "find_by_name")]
    pub fn find_by_name_js(&self, name: String) -> Option<u32> {
        self.find_by_name(&name)
    }

    /// Root nodes (those without a parent).
    #[wasm_bindgen(js_name = "roots")]
    pub fn roots_js(&self) -> Vec<u32> {
        self.roots()
    }

    /// Serialise the explorer to a `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}
