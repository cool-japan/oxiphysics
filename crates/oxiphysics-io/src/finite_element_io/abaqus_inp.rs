//! Reader for Abaqus `.inp` input files.
//!
//! The reader keeps what [`FeMesh`] can hold: nodes, elements, node and
//! element sets, isotropic elastic materials with their density, analysis
//! steps, displacement boundary conditions and concentrated loads.  Keywords it
//! does not know are skipped together with their data lines.
//!
//! # Data lines
//!
//! * `*NODE`: `node, x[, y[, z]]`; coordinates that are left out are zero, so
//!   two-dimensional models (`node, x, y`) are read.  `NSET=` collects the
//!   nodes in a set.
//! * `*ELEMENT`: `element, node...`.  The element type comes from `TYPE=`,
//!   `ELSET=` collects the elements in a set.  An element is one record that
//!   may span several lines: when the element type has a known node count the
//!   record is complete as soon as it has that many nodes; for other types a
//!   line that ends in a comma continues on the next line.
//! * `*NSET` / `*ELSET`: ids and names of earlier sets; with `GENERATE` each
//!   data line is `first, last[, increment]`.
//! * `*MATERIAL`, `*ELASTIC`, `*DENSITY`: `E, nu` of an isotropic material (the
//!   default `TYPE`; other elasticity types are reported and skipped) and its
//!   density, in either order.  Materials are numbered from 1 in file order.
//! * `*BOUNDARY`: `node|set, first[, last[, value]]` constrains the degrees of
//!   freedom `first..=last` (`last` defaults to `first`, `value` to zero);
//!   `node|set, ENCASTRE|PINNED|XSYMM|...` names a standard constraint.
//! * `*CLOAD`: `node|set, dof, magnitude`.
//!
//! Boundary conditions given in the model data, before the first `*STEP`
//! (or after an `*END STEP`), are model-level conditions and go to
//! [`FeMesh::global_bcs`]; those inside a step go to that step.  Sets are
//! expanded once the whole file is read, so a set may be defined after the
//! data line that uses it, and set names are case-insensitive.

use std::collections::HashMap;

use super::types::{
    AbaqusParseState, AnalysisStep, DirichletBc, FeElement, FeElementType, FeMesh, FeNode,
    LinearElasticMaterial, NodalForce,
};

/// Largest number of ids a single `GENERATE` range may expand to.
const MAX_GENERATED_IDS: usize = 1 << 26;

/// Parse an Abaqus input file.
///
/// Returns the mesh together with a message for every data line that could not
/// be interpreted and was skipped.
pub(super) fn parse_inp(source: &str) -> (FeMesh, Vec<String>) {
    Reader::default().run(source)
}

/// A keyword line: `*NAME, KEY=value, FLAG`.
struct Keyword {
    /// Upper-case keyword with single blanks, e.g. `*END STEP`.
    name: String,
    /// Parameters as `(UPPER-CASE KEY, value as written)`.
    params: Vec<(String, Option<String>)>,
}

impl Keyword {
    fn parse(line: &str) -> Self {
        let mut parts = line.split(',');
        let name = parts
            .next()
            .unwrap_or("")
            .split_whitespace()
            .collect::<Vec<_>>()
            .join(" ")
            .to_ascii_uppercase();
        let params = parts
            .map(str::trim)
            .filter(|part| !part.is_empty())
            .map(|part| match part.split_once('=') {
                Some((key, value)) => (
                    key.trim().to_ascii_uppercase(),
                    Some(value.trim().to_string()),
                ),
                None => (part.to_ascii_uppercase(), None),
            })
            .collect();
        Self { name, params }
    }

    /// The value of the parameter `key` (given in upper case).
    fn value(&self, key: &str) -> Option<&str> {
        self.params
            .iter()
            .find(|(name, _)| name == key)
            .and_then(|(_, value)| value.as_deref())
            .filter(|value| !value.is_empty())
    }

    /// Does the keyword carry the parameter `key` (with or without a value)?
    fn has(&self, key: &str) -> bool {
        self.params.iter().any(|(name, _)| name == key)
    }
}

/// A node, or a node set that is expanded at the end of the file.
#[derive(Debug, Clone, PartialEq)]
enum NodeTarget {
    Node(usize),
    Set(String),
}

/// One `*BOUNDARY` data line.
#[derive(Debug, Clone, PartialEq)]
struct BoundaryRecord {
    target: NodeTarget,
    dofs: Vec<u8>,
    value: f64,
}

/// One `*CLOAD` data line.
#[derive(Debug, Clone, PartialEq)]
struct LoadRecord {
    target: NodeTarget,
    dof: u8,
    magnitude: f64,
}

/// Where a boundary condition applies.
#[derive(Debug, Clone, Copy)]
enum Scope {
    /// Model data: every step.
    Model,
    /// The step with this index.
    Step(usize),
}

struct PendingBoundary {
    scope: Scope,
    line: usize,
    record: BoundaryRecord,
}

struct PendingLoad {
    step: usize,
    line: usize,
    record: LoadRecord,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SetKind {
    Node,
    Element,
}

impl SetKind {
    fn keyword(self) -> &'static str {
        match self {
            Self::Node => "*NSET",
            Self::Element => "*ELSET",
        }
    }
}

#[derive(Default)]
struct Reader {
    mesh: FeMesh,
    state: AbaqusParseState,
    boundaries: Vec<PendingBoundary>,
    loads: Vec<PendingLoad>,
    /// Density of each material by name; `*DENSITY` may come before or after
    /// the `*ELASTIC` of its material, so it is applied once the file is read.
    densities: HashMap<String, f64>,
    /// Inside the data of a `*DENSITY` keyword.
    in_density: bool,
    issues: Vec<String>,
}

impl Reader {
    fn run(mut self, source: &str) -> (FeMesh, Vec<String>) {
        let lines: Vec<&str> = source.lines().collect();
        let mut at = 0;
        while let Some(&raw) = lines.get(at) {
            let number = at + 1;
            at += 1;
            let text = raw.trim();
            if text.is_empty() || text.starts_with("**") {
                continue;
            }
            if text.starts_with('*') {
                // A keyword line that ends in a comma continues on the next
                // line.
                let mut keyword_line = text.to_string();
                while keyword_line.ends_with(',') {
                    match lines.get(at).map(|next| next.trim()) {
                        Some(next) if next.is_empty() || next.starts_with("**") => at += 1,
                        Some(next) if !next.starts_with('*') => {
                            keyword_line.push_str(next);
                            at += 1;
                        }
                        _ => break,
                    }
                }
                self.keyword(&keyword_line, number);
            } else {
                self.data_line(text, number);
            }
        }
        self.finish()
    }

    fn keyword(&mut self, line: &str, number: usize) {
        self.flush_pending_element();
        let keyword = Keyword::parse(line);
        self.in_density = false;
        let state = &mut self.state;
        state.node_section = false;
        state.element_section = false;
        state.node_set_section = None;
        state.elem_set_section = None;
        state.in_elastic = false;
        state.in_boundary = false;
        state.in_cload = false;
        state.generate = false;
        state.node_collect_set = None;
        state.element_collect_set = None;
        match keyword.name.as_str() {
            "*NODE" => {
                state.node_section = true;
                state.node_collect_set = keyword.value("NSET").map(str::to_string);
                if let Some(name) = state.node_collect_set.clone() {
                    self.set_entry(SetKind::Node, &name);
                }
            }
            "*ELEMENT" => {
                state.element_section = true;
                state.current_elem_type = element_type_from_name(keyword.value("TYPE"));
                state.element_collect_set = keyword.value("ELSET").map(str::to_string);
                if let Some(name) = state.element_collect_set.clone() {
                    self.set_entry(SetKind::Element, &name);
                }
            }
            "*NSET" | "*ELSET" => {
                let kind = if keyword.name == "*NSET" {
                    SetKind::Node
                } else {
                    SetKind::Element
                };
                let param = if kind == SetKind::Node {
                    "NSET"
                } else {
                    "ELSET"
                };
                match keyword.value(param) {
                    Some(name) => {
                        match kind {
                            SetKind::Node => state.node_set_section = Some(name.to_string()),
                            SetKind::Element => state.elem_set_section = Some(name.to_string()),
                        }
                        state.generate = keyword.has("GENERATE");
                        self.set_entry(kind, name);
                    }
                    None => self.issues.push(format!(
                        "line {number}: {} without a {param}= name",
                        kind.keyword()
                    )),
                }
            }
            "*MATERIAL" => {
                state.material_name = Some(keyword.value("NAME").unwrap_or("").to_string());
            }
            "*ELASTIC" => {
                let isotropic = keyword
                    .value("TYPE")
                    .is_none_or(|kind| kind.eq_ignore_ascii_case("ISOTROPIC"));
                state.in_elastic = isotropic;
                if !isotropic {
                    let kind = keyword.value("TYPE").unwrap_or("");
                    self.issues.push(format!(
                        "line {number}: *ELASTIC, TYPE={kind} is not isotropic elasticity and was ignored"
                    ));
                }
            }
            "*DENSITY" => self.in_density = true,
            "*STEP" => {
                state.in_step = true;
                self.mesh
                    .steps
                    .push(AnalysisStep::new(keyword.value("NAME").unwrap_or("")));
            }
            "*END STEP" => state.in_step = false,
            "*BOUNDARY" => {
                // Velocity and acceleration conditions, user subroutines and
                // external files are not displacement data.
                let displacement = keyword
                    .value("TYPE")
                    .is_none_or(|kind| kind.eq_ignore_ascii_case("DISPLACEMENT"));
                state.in_boundary =
                    displacement && !keyword.has("USER") && keyword.value("INPUT").is_none();
            }
            "*CLOAD" => state.in_cload = true,
            _ => {}
        }
    }

    fn data_line(&mut self, text: &str, number: usize) {
        if self.state.node_section {
            self.node_line(text, number);
        } else if self.state.element_section {
            self.element_line(text, number);
        } else if let Some(name) = self.state.node_set_section.clone() {
            self.set_line(SetKind::Node, &name, text, number);
        } else if let Some(name) = self.state.elem_set_section.clone() {
            self.set_line(SetKind::Element, &name, text, number);
        } else if self.state.in_boundary {
            self.boundary_line(text, number);
        } else if self.state.in_cload {
            self.cload_line(text, number);
        } else if self.state.in_elastic {
            self.elastic_line(text, number);
        } else if self.in_density {
            self.density_line(text, number);
        }
    }

    /// The set map for `kind`.
    fn sets_mut(&mut self, kind: SetKind) -> &mut HashMap<String, Vec<usize>> {
        match kind {
            SetKind::Node => &mut self.mesh.node_sets,
            SetKind::Element => &mut self.mesh.element_sets,
        }
    }

    /// Make sure the set `name` exists (set names are case-insensitive) and
    /// return the key it is stored under.
    fn set_entry(&mut self, kind: SetKind, name: &str) -> String {
        let sets = self.sets_mut(kind);
        let key = existing_key(sets, name).unwrap_or_else(|| name.to_string());
        sets.entry(key.clone()).or_default();
        key
    }

    fn node_line(&mut self, text: &str, number: usize) {
        let Some(node) = parse_node_line(text) else {
            self.issues.push(format!(
                "line {number}: *NODE: malformed data line `{text}`"
            ));
            return;
        };
        if let Some(name) = self.state.node_collect_set.clone() {
            let key = self.set_entry(SetKind::Node, &name);
            self.sets_mut(SetKind::Node)
                .entry(key)
                .or_default()
                .push(node.id);
        }
        self.mesh.nodes.push(node);
    }

    fn element_line(&mut self, text: &str, number: usize) {
        let (start, record) = match self.state.pending_element.take() {
            Some((start, mut previous)) => {
                if !previous.trim_end().ends_with(',') {
                    previous.push(',');
                }
                previous.push_str(text);
                (start, previous)
            }
            None => (number, text.to_string()),
        };
        if element_record_complete(&record, &self.state.current_elem_type) {
            self.finish_element(&record, start);
        } else {
            self.state.pending_element = Some((start, record));
        }
    }

    /// Store an element record that was cut off by the next keyword.
    fn flush_pending_element(&mut self) {
        if let Some((start, record)) = self.state.pending_element.take() {
            self.finish_element(&record, start);
        }
    }

    fn finish_element(&mut self, record: &str, number: usize) {
        let element_type = self.state.current_elem_type.clone();
        let Some(element) = parse_element_record(record, element_type.clone()) else {
            self.issues.push(format!(
                "line {number}: *ELEMENT: malformed data line `{record}`"
            ));
            return;
        };
        if let Some(expected) = element_type.node_count()
            && element.connectivity.len() != expected
        {
            self.issues.push(format!(
                "line {number}: *ELEMENT: element {} of type {:?} needs {expected} nodes, has {}",
                element.id,
                element_type,
                element.connectivity.len()
            ));
            return;
        }
        if let Some(name) = self.state.element_collect_set.clone() {
            let key = self.set_entry(SetKind::Element, &name);
            self.sets_mut(SetKind::Element)
                .entry(key)
                .or_default()
                .push(element.id);
        }
        self.mesh.elements.push(element);
    }

    fn set_line(&mut self, kind: SetKind, name: &str, text: &str, number: usize) {
        let ids = if self.state.generate {
            match generate_ids(text) {
                Ok(ids) => ids,
                Err(message) => {
                    self.issues.push(format!(
                        "line {number}: {} GENERATE: {message}",
                        kind.keyword()
                    ));
                    return;
                }
            }
        } else {
            self.explicit_ids(kind, text, number)
        };
        let key = self.set_entry(kind, name);
        self.sets_mut(kind).entry(key).or_default().extend(ids);
    }

    /// The ids on a set data line; a name stands for the members of an
    /// earlier set.
    fn explicit_ids(&mut self, kind: SetKind, text: &str, number: usize) -> Vec<usize> {
        let mut ids = Vec::new();
        for token in text.split(',').map(str::trim).filter(|t| !t.is_empty()) {
            if let Ok(id) = token.parse::<usize>() {
                ids.push(id);
                continue;
            }
            let members = {
                let sets = self.sets_mut(kind);
                existing_key(sets, token).and_then(|key| sets.get(&key).cloned())
            };
            match members {
                Some(members) => ids.extend(members),
                None => self.issues.push(format!(
                    "line {number}: {}: `{token}` is neither an id nor a defined set",
                    kind.keyword()
                )),
            }
        }
        ids
    }

    fn boundary_line(&mut self, text: &str, number: usize) {
        let Some(record) = parse_boundary_record(text) else {
            self.issues.push(format!(
                "line {number}: *BOUNDARY: malformed data line `{text}`"
            ));
            return;
        };
        let scope = match self.mesh.steps.len().checked_sub(1) {
            Some(index) if self.state.in_step => Scope::Step(index),
            _ => Scope::Model,
        };
        self.boundaries.push(PendingBoundary {
            scope,
            line: number,
            record,
        });
    }

    fn cload_line(&mut self, text: &str, number: usize) {
        let Some(record) = parse_load_record(text) else {
            self.issues.push(format!(
                "line {number}: *CLOAD: malformed data line `{text}`"
            ));
            return;
        };
        match self.mesh.steps.len().checked_sub(1) {
            Some(step) if self.state.in_step => self.loads.push(PendingLoad {
                step,
                line: number,
                record,
            }),
            _ => self.issues.push(format!(
                "line {number}: *CLOAD data line `{text}` is outside a *STEP and was ignored"
            )),
        }
    }

    fn elastic_line(&mut self, text: &str, number: usize) {
        let Some(name) = self.state.material_name.clone() else {
            return;
        };
        match parse_elastic_line(text, self.mesh.materials.len() + 1, &name) {
            Some(material) => self.mesh.materials.push(material),
            None => self.issues.push(format!(
                "line {number}: *ELASTIC: malformed data line `{text}`"
            )),
        }
        self.state.in_elastic = false;
    }

    /// The data line of `*DENSITY`: `density[, temperature]`.
    fn density_line(&mut self, text: &str, number: usize) {
        let Some(name) = self.state.material_name.clone() else {
            return;
        };
        match text.split(',').next().and_then(parse_finite) {
            Some(density) => {
                self.densities.insert(name, density);
            }
            None => self.issues.push(format!(
                "line {number}: *DENSITY: malformed data line `{text}`"
            )),
        }
        self.in_density = false;
    }

    /// Expand the boundary conditions and loads now that every set is known.
    fn finish(mut self) -> (FeMesh, Vec<String>) {
        self.flush_pending_element();
        for pending in std::mem::take(&mut self.boundaries) {
            let Some(nodes) = nodes_of(&pending.record.target, &self.mesh.node_sets) else {
                self.issues.push(undefined_set(
                    pending.line,
                    "*BOUNDARY",
                    &pending.record.target,
                ));
                continue;
            };
            let target = match pending.scope {
                Scope::Model => Some(&mut self.mesh.global_bcs),
                Scope::Step(index) => self.mesh.steps.get_mut(index).map(|step| &mut step.bcs),
            };
            if let Some(out) = target {
                for node_id in nodes {
                    for &dof in &pending.record.dofs {
                        out.push(DirichletBc {
                            node_id,
                            dof,
                            value: pending.record.value,
                        });
                    }
                }
            }
        }
        for pending in std::mem::take(&mut self.loads) {
            let Some(nodes) = nodes_of(&pending.record.target, &self.mesh.node_sets) else {
                self.issues.push(undefined_set(
                    pending.line,
                    "*CLOAD",
                    &pending.record.target,
                ));
                continue;
            };
            if let Some(step) = self.mesh.steps.get_mut(pending.step) {
                for node_id in nodes {
                    step.forces.push(NodalForce::new(
                        node_id,
                        pending.record.dof,
                        pending.record.magnitude,
                    ));
                }
            }
        }
        for material in &mut self.mesh.materials {
            if let Some(&density) = self.densities.get(&material.name) {
                material.density = density;
            }
        }
        (self.mesh, self.issues)
    }
}

fn undefined_set(line: usize, keyword: &str, target: &NodeTarget) -> String {
    let name = match target {
        NodeTarget::Set(name) => name.as_str(),
        NodeTarget::Node(_) => "",
    };
    format!("line {line}: {keyword} refers to the undefined node set `{name}`")
}

/// The key under which a set with this (case-insensitive) name is stored.
fn existing_key(sets: &HashMap<String, Vec<usize>>, name: &str) -> Option<String> {
    if sets.contains_key(name) {
        return Some(name.to_string());
    }
    sets.keys()
        .filter(|key| key.eq_ignore_ascii_case(name))
        .min()
        .cloned()
}

/// The nodes a `*BOUNDARY` / `*CLOAD` target designates; `None` for a set that
/// does not exist.
fn nodes_of(target: &NodeTarget, node_sets: &HashMap<String, Vec<usize>>) -> Option<Vec<usize>> {
    match target {
        NodeTarget::Node(id) => Some(vec![*id]),
        NodeTarget::Set(name) => {
            existing_key(node_sets, name).and_then(|key| node_sets.get(&key).cloned())
        }
    }
}

/// Parse a node data line: `id, x[, y[, z]]`, missing coordinates are zero.
fn parse_node_line(line: &str) -> Option<FeNode> {
    let parts: Vec<&str> = line.split(',').map(str::trim).collect();
    if parts.len() < 2 {
        return None;
    }
    let id: usize = parts[0].parse().ok()?;
    let mut coords = [0.0; 3];
    for (slot, text) in coords.iter_mut().zip(&parts[1..]) {
        if !text.is_empty() {
            *slot = text.parse::<f64>().ok().filter(|value| value.is_finite())?;
        }
    }
    Some(FeNode::new(id, coords))
}

/// The nodes in an element record: every field after the element number.
fn record_node_ids(record: &str) -> impl Iterator<Item = usize> + '_ {
    record
        .split(',')
        .skip(1)
        .filter_map(|text| text.trim().parse::<usize>().ok())
}

/// Is `record` the whole definition of an element of type `element_type`?
fn element_record_complete(record: &str, element_type: &FeElementType) -> bool {
    match element_type.node_count() {
        Some(expected) => record_node_ids(record).count() >= expected,
        None => !record.trim_end().ends_with(','),
    }
}

/// Parse an element record: `id, n1, n2, ...`.
fn parse_element_record(record: &str, element_type: FeElementType) -> Option<FeElement> {
    let id: usize = record.split(',').next()?.trim().parse().ok()?;
    let connectivity: Vec<usize> = record_node_ids(record).collect();
    if connectivity.is_empty() {
        return None;
    }
    Some(FeElement::new(id, element_type, connectivity))
}

/// Parse the data line of `*ELASTIC`: `E, nu`.
fn parse_elastic_line(line: &str, id: usize, name: &str) -> Option<LinearElasticMaterial> {
    let parts: Vec<&str> = line.split(',').collect();
    if parts.len() < 2 {
        return None;
    }
    let e = parse_finite(parts[0])?;
    let nu = parse_finite(parts[1])?;
    Some(LinearElasticMaterial::new(id, name, e, nu, 0.0))
}

/// Element type for the value of the `TYPE=` parameter.
fn element_type_from_name(name: Option<&str>) -> FeElementType {
    let name = name.unwrap_or("");
    let up = name.to_ascii_uppercase();
    if up.contains("C3D4") || up.contains("TET4") {
        FeElementType::Tet4
    } else if up.contains("C3D10") || up.contains("TET10") {
        FeElementType::Tet10
    } else if up.contains("C3D8") || up.contains("HEX8") {
        FeElementType::Hex8
    } else if up.contains("C3D20") || up.contains("HEX20") {
        FeElementType::Hex20
    } else if up.contains("S3") || up.contains("TRI3") {
        FeElementType::Tri3
    } else if up.contains("S4") || up.contains("QUAD4") {
        FeElementType::Quad4
    } else if up.contains("T3D2") {
        FeElementType::Line2
    } else {
        FeElementType::Unknown(name.to_string())
    }
}

/// Expand a `GENERATE` data line `first, last[, increment]`.
fn generate_ids(line: &str) -> Result<Vec<usize>, String> {
    let fields: Vec<&str> = line.split(',').map(str::trim).collect();
    let number = |index: usize| -> Result<Option<usize>, String> {
        match fields.get(index).copied().filter(|text| !text.is_empty()) {
            None => Ok(None),
            Some(text) => text
                .parse::<usize>()
                .map(Some)
                .map_err(|_| format!("`{text}` is not a valid id in `{line}`")),
        }
    };
    let (Some(first), Some(last)) = (number(0)?, number(1)?) else {
        return Err(format!("`{line}` needs `first, last[, increment]`"));
    };
    let step = number(2)?.unwrap_or(1);
    if step == 0 {
        return Err(format!("the increment of `{line}` must be at least 1"));
    }
    if first > last {
        return Err(format!("`{line}` has first > last"));
    }
    let count = (last - first) / step + 1;
    if count > MAX_GENERATED_IDS {
        return Err(format!(
            "`{line}` would generate {count} ids, more than the limit of {MAX_GENERATED_IDS}"
        ));
    }
    Ok((first..=last).step_by(step).collect())
}

fn parse_target(field: &str) -> Option<NodeTarget> {
    let field = field.trim();
    if field.is_empty() {
        return None;
    }
    Some(match field.parse::<usize>() {
        Ok(id) => NodeTarget::Node(id),
        Err(_) => NodeTarget::Set(field.to_string()),
    })
}

/// A degree of freedom number; they start at 1.
fn parse_dof(field: &str) -> Option<u8> {
    field.trim().parse::<u8>().ok().filter(|dof| *dof >= 1)
}

fn parse_finite(field: &str) -> Option<f64> {
    field.trim().parse::<f64>().ok().filter(|v| v.is_finite())
}

/// The degrees of freedom a named `*BOUNDARY` constraint fixes.
fn named_constraint_dofs(name: &str) -> Option<Vec<u8>> {
    let dofs: &[u8] = match name.trim().to_ascii_uppercase().as_str() {
        "ENCASTRE" => &[1, 2, 3, 4, 5, 6],
        "PINNED" => &[1, 2, 3],
        "XSYMM" => &[1, 5, 6],
        "YSYMM" => &[2, 4, 6],
        "ZSYMM" => &[3, 4, 5],
        "XASYMM" => &[2, 3, 4],
        "YASYMM" => &[1, 3, 5],
        "ZASYMM" => &[1, 2, 6],
        _ => return None,
    };
    Some(dofs.to_vec())
}

/// Parse a `*BOUNDARY` data line.
fn parse_boundary_record(line: &str) -> Option<BoundaryRecord> {
    let parts: Vec<&str> = line.split(',').map(str::trim).collect();
    if parts.len() < 2 {
        return None;
    }
    let target = parse_target(parts[0])?;
    if let Some(dofs) = named_constraint_dofs(parts[1]) {
        return Some(BoundaryRecord {
            target,
            dofs,
            value: 0.0,
        });
    }
    let first = parse_dof(parts[1])?;
    let given = |index: usize| parts.get(index).copied().filter(|field| !field.is_empty());
    let last = match given(2) {
        None => first,
        Some(field) => parse_dof(field)?,
    };
    let value = match given(3) {
        None => 0.0,
        Some(field) => parse_finite(field)?,
    };
    Some(BoundaryRecord {
        target,
        dofs: (first.min(last)..=first.max(last)).collect(),
        value,
    })
}

/// Parse a `*CLOAD` data line: `node|set, dof, magnitude`.
fn parse_load_record(line: &str) -> Option<LoadRecord> {
    let parts: Vec<&str> = line.split(',').map(str::trim).collect();
    if parts.len() < 3 {
        return None;
    }
    Some(LoadRecord {
        target: parse_target(parts[0])?,
        dof: parse_dof(parts[1])?,
        magnitude: parse_finite(parts[2])?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_element_io::types::AbaqusInpIo;

    fn parse(inp: &str) -> FeMesh {
        AbaqusInpIo::try_parse(inp).expect("input should parse without problems")
    }

    fn issues(inp: &str) -> String {
        AbaqusInpIo::try_parse(inp)
            .expect_err("input should report problems")
            .to_string()
    }

    const NODES: &str = "*NODE\n1, 0, 0, 0\n2, 1, 0, 0\n3, 1, 1, 0\n4, 0, 1, 0\n";

    #[test]
    fn two_dimensional_nodes_get_a_zero_z_coordinate() {
        let mesh = parse("*NODE\n1, 0.5, 2.5\n2, -1.0, 3.0, 4.0\n3, 7\n4, 1.0, , 2.0\n");
        assert_eq!(mesh.nodes.len(), 4);
        assert_eq!(mesh.nodes[0].coords, [0.5, 2.5, 0.0]);
        assert_eq!(mesh.nodes[1].coords, [-1.0, 3.0, 4.0]);
        assert_eq!(mesh.nodes[2].coords, [7.0, 0.0, 0.0]);
        assert_eq!(mesh.nodes[3].coords, [1.0, 0.0, 2.0]);
    }

    #[test]
    fn node_lines_with_direction_cosines_and_bad_numbers() {
        let mesh = AbaqusInpIo::parse("*NODE\n1, 0, 0, 0, 0, 0, 1\n2, 0, x, 0\n3\n");
        assert_eq!(mesh.nodes.len(), 1);
        let text = issues("*NODE\n1, 0, 0, 0\n2, 0, x, 0\n");
        assert!(text.contains("line 3") && text.contains("*NODE"), "{text}");
    }

    #[test]
    fn quadratic_hexahedron_spans_two_lines() {
        let inp = format!(
            "{NODES}*ELEMENT, TYPE=C3D20R\n\
             101, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,\n\
             16, 17, 18, 19, 20\n\
             102, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,\n\
             36, 37, 38, 39, 40\n"
        );
        let mesh = parse(&inp);
        assert_eq!(mesh.elements.len(), 2);
        assert_eq!(mesh.elements[0].id, 101);
        assert_eq!(mesh.elements[0].element_type, FeElementType::Hex20);
        assert_eq!(mesh.elements[0].connectivity, (1..=20).collect::<Vec<_>>());
        assert_eq!(mesh.elements[1].id, 102);
        assert_eq!(mesh.elements[1].connectivity, (21..=40).collect::<Vec<_>>());
    }

    #[test]
    fn a_continued_record_without_a_trailing_comma_is_joined_by_node_count() {
        let mesh = parse(
            "*ELEMENT, TYPE=C3D8\n1, 1, 2, 3, 4\n5, 6, 7, 8\n2, 11, 12, 13, 14, 15, 16, 17, 18\n",
        );
        assert_eq!(mesh.elements.len(), 2);
        assert_eq!(mesh.elements[0].connectivity, vec![1, 2, 3, 4, 5, 6, 7, 8]);
        assert_eq!(mesh.elements[1].id, 2);
    }

    #[test]
    fn unknown_element_types_use_the_trailing_comma_to_continue() {
        let mesh = parse("*ELEMENT, TYPE=WEIRD\n1, 1, 2, 3,\n4, 5\n2, 6, 7\n");
        assert_eq!(mesh.elements.len(), 2);
        assert_eq!(mesh.elements[0].connectivity, vec![1, 2, 3, 4, 5]);
        assert_eq!(
            mesh.elements[0].element_type,
            FeElementType::Unknown("WEIRD".to_string())
        );
        assert_eq!(mesh.elements[1].connectivity, vec![6, 7]);
    }

    #[test]
    fn element_with_the_wrong_node_count_is_dropped_and_reported() {
        let inp = "*ELEMENT, TYPE=C3D4\n1, 1, 2, 3, 4, 5\n2, 1, 2, 3, 4\n";
        let mesh = AbaqusInpIo::parse(inp);
        assert_eq!(mesh.elements.len(), 1);
        assert_eq!(mesh.elements[0].id, 2);
        let text = issues(inp);
        assert!(text.contains("needs 4 nodes, has 5"), "{text}");
        // An element cut off by the next keyword is reported, not kept.
        let cut = "*ELEMENT, TYPE=C3D8\n1, 1, 2, 3,\n*NSET, NSET=A\n1\n";
        assert!(AbaqusInpIo::parse(cut).elements.is_empty());
        assert!(issues(cut).contains("needs 8 nodes, has 3"));
    }

    #[test]
    fn element_type_comes_from_the_type_parameter_only() {
        // The set name contains "C3D4", which must not decide the type.
        let mesh = parse("*ELEMENT, ELSET=BLOCKS_C3D4, TYPE=C3D8R\n1, 1, 2, 3, 4, 5, 6, 7, 8\n");
        assert_eq!(mesh.elements[0].element_type, FeElementType::Hex8);
    }

    #[test]
    fn nset_and_elset_parameters_on_node_and_element_keywords_collect_ids() {
        let mesh = parse(
            "*NODE, NSET=ALLN\n1, 0, 0, 0\n2, 1, 0, 0\n*ELEMENT, TYPE=T3D2, ELSET=bars\n7, 1, 2\n8, 2, 1\n",
        );
        assert_eq!(mesh.node_sets.get("ALLN"), Some(&vec![1, 2]));
        assert_eq!(mesh.element_sets.get("bars"), Some(&vec![7, 8]));
    }

    #[test]
    fn nset_generate_expands_ranges_with_increments() {
        let mesh = parse(
            "*NSET, NSET=EVERY3, GENERATE\n1, 10, 3\n*NSET, NSET=RUN, generate\n20, 23\n*NSET, NSET=BOTH, GENERATE\n1, 3\n7, 9, 1\n",
        );
        assert_eq!(mesh.node_sets.get("EVERY3"), Some(&vec![1, 4, 7, 10]));
        assert_eq!(mesh.node_sets.get("RUN"), Some(&vec![20, 21, 22, 23]));
        assert_eq!(mesh.node_sets.get("BOTH"), Some(&vec![1, 2, 3, 7, 8, 9]));
    }

    #[test]
    fn elset_generate_and_explicit_ids_and_set_references() {
        let mesh = parse(
            "*ELSET, ELSET=A, GENERATE\n1, 5\n*ELSET, ELSET=B\n10, 11,\n12\n*ELSET, ELSET=C\nA, B, 99\n",
        );
        assert_eq!(mesh.element_sets.get("A"), Some(&vec![1, 2, 3, 4, 5]));
        assert_eq!(mesh.element_sets.get("B"), Some(&vec![10, 11, 12]));
        assert_eq!(
            mesh.element_sets.get("C"),
            Some(&vec![1, 2, 3, 4, 5, 10, 11, 12, 99])
        );
    }

    #[test]
    fn invalid_generate_ranges_are_reported_and_ignored() {
        for (data, fragment) in [
            ("1, 10, 0", "at least 1"),
            ("10, 1", "first > last"),
            ("1", "needs"),
            ("a, 3", "not a valid id"),
            ("0, 999999999999, 1", "limit"),
        ] {
            let inp = format!("*NSET, NSET=S, GENERATE\n{data}\n");
            let text = issues(&inp);
            assert!(text.contains(fragment), "{data}: {text}");
            let mesh = AbaqusInpIo::parse(&inp);
            assert_eq!(mesh.node_sets.get("S"), Some(&Vec::new()), "{data}");
        }
    }

    #[test]
    fn set_names_are_case_insensitive_and_sets_can_be_extended() {
        let mesh =
            parse("*NSET, NSET=Fixed\n1, 2\n*NSET, NSET=FIXED\n3\n*NSET, NSET=ALL\nfixed, 9\n");
        assert_eq!(mesh.node_sets.len(), 2);
        assert_eq!(mesh.node_sets.get("Fixed"), Some(&vec![1, 2, 3]));
        assert_eq!(mesh.node_sets.get("ALL"), Some(&vec![1, 2, 3, 9]));
    }

    #[test]
    fn undefined_set_references_in_set_data_are_reported() {
        let text = issues("*NSET, NSET=A\n1, MISSING\n");
        assert!(text.contains("MISSING"), "{text}");
        let text = issues("*NSET\n1\n");
        assert!(text.contains("without a NSET= name"), "{text}");
    }

    #[test]
    fn materials_get_sequential_ids_and_their_density_in_either_order() {
        let mesh = parse(
            "*MATERIAL, NAME=Steel\n*DENSITY\n7850.,\n*ELASTIC\n210.0E9, 0.3\n\
             *MATERIAL, NAME=Alu\n*ELASTIC, TYPE=ISOTROPIC\n70.0E9, 0.33\n*DENSITY\n2700., 20.\n",
        );
        assert_eq!(mesh.materials.len(), 2);
        assert_eq!(mesh.materials[0].id, 1);
        assert_eq!(mesh.materials[0].name, "Steel");
        assert_eq!(mesh.materials[0].density, 7850.0);
        assert_eq!(mesh.materials[0].young_modulus, 210.0e9);
        assert_eq!(mesh.materials[1].id, 2);
        assert_eq!(mesh.materials[1].name, "Alu");
        assert_eq!(mesh.materials[1].poisson_ratio, 0.33);
        assert_eq!(mesh.materials[1].density, 2700.0);
        let text = issues("*MATERIAL, NAME=M\n*DENSITY\nheavy\n");
        assert!(text.contains("*DENSITY: malformed"), "{text}");
    }

    #[test]
    fn anisotropic_elasticity_is_reported_not_read_as_isotropic() {
        let inp = "*MATERIAL, NAME=Ply\n*ELASTIC, TYPE=ORTHOTROPIC\n1.0E9, 2.0E9, 3.0E9, 0.3, 0.3, 0.3, 4.0E8, 4.0E8\n";
        assert!(AbaqusInpIo::parse(inp).materials.is_empty());
        let text = issues(inp);
        assert!(
            text.contains("TYPE=ORTHOTROPIC") && text.contains("line 2"),
            "{text}"
        );
    }

    #[test]
    fn boundary_in_a_step_constrains_each_dof_of_the_range() {
        let mesh = parse(&format!(
            "{NODES}*STEP, NAME=load\n*STATIC\n*BOUNDARY\n1, 1, 3\n2, 2\n3, 1, 2, 0.25\n4, 3, 3, -1.5e-3\n*END STEP\n"
        ));
        assert!(mesh.global_bcs.is_empty());
        assert_eq!(mesh.steps.len(), 1);
        assert_eq!(
            mesh.steps[0].bcs,
            vec![
                DirichletBc::fixed(1, 1),
                DirichletBc::fixed(1, 2),
                DirichletBc::fixed(1, 3),
                DirichletBc::fixed(2, 2),
                DirichletBc::prescribed(3, 1, 0.25),
                DirichletBc::prescribed(3, 2, 0.25),
                DirichletBc::prescribed(4, 3, -1.5e-3),
            ]
        );
    }

    #[test]
    fn model_level_boundary_goes_to_global_bcs_not_to_a_phantom_step() {
        let mesh = parse(&format!(
            "{NODES}*BOUNDARY\n1, 1, 3\n*STEP, NAME=push\n*BOUNDARY\n2, 2\n*END STEP\n"
        ));
        assert_eq!(
            mesh.global_bcs,
            vec![
                DirichletBc::fixed(1, 1),
                DirichletBc::fixed(1, 2),
                DirichletBc::fixed(1, 3)
            ]
        );
        assert_eq!(mesh.steps.len(), 1);
        assert_eq!(mesh.steps[0].name, "push");
        assert_eq!(mesh.steps[0].bcs, vec![DirichletBc::fixed(2, 2)]);
        // A model without any step still yields its supports.
        let only_model = parse(&format!("{NODES}*BOUNDARY\n4, 3\n"));
        assert!(only_model.steps.is_empty());
        assert_eq!(only_model.global_bcs, vec![DirichletBc::fixed(4, 3)]);
    }

    #[test]
    fn boundary_after_end_step_belongs_to_the_model_again() {
        let mesh = parse(&format!(
            "{NODES}*STEP\n*BOUNDARY\n1, 1\n*END STEP\n*BOUNDARY\n2, 2\n"
        ));
        assert_eq!(mesh.steps[0].bcs, vec![DirichletBc::fixed(1, 1)]);
        assert_eq!(mesh.global_bcs, vec![DirichletBc::fixed(2, 2)]);
    }

    #[test]
    fn boundary_on_a_node_set_expands_over_its_members() {
        let mesh = parse(&format!(
            "{NODES}*NSET, NSET=Left\n1, 4\n*STEP\n*BOUNDARY\nleft, 1, 2, 0.5\nLEFT, 3\n*END STEP\n"
        ));
        let bcs = &mesh.steps[0].bcs;
        assert_eq!(bcs.len(), 6);
        assert_eq!(bcs[0], DirichletBc::prescribed(1, 1, 0.5));
        assert_eq!(bcs[1], DirichletBc::prescribed(1, 2, 0.5));
        assert_eq!(bcs[2], DirichletBc::prescribed(4, 1, 0.5));
        assert_eq!(bcs[3], DirichletBc::prescribed(4, 2, 0.5));
        assert_eq!(bcs[4], DirichletBc::fixed(1, 3));
        assert_eq!(bcs[5], DirichletBc::fixed(4, 3));
    }

    #[test]
    fn a_set_defined_after_its_boundary_condition_is_still_expanded() {
        let mesh = parse(&format!(
            "{NODES}*BOUNDARY\nLATE, 2\n*NSET, NSET=late, GENERATE\n2, 4\n"
        ));
        assert_eq!(
            mesh.global_bcs,
            vec![
                DirichletBc::fixed(2, 2),
                DirichletBc::fixed(3, 2),
                DirichletBc::fixed(4, 2)
            ]
        );
    }

    #[test]
    fn named_constraints_fix_the_standard_degrees_of_freedom() {
        let expected: [(&str, &[u8]); 8] = [
            ("ENCASTRE", &[1, 2, 3, 4, 5, 6]),
            ("PINNED", &[1, 2, 3]),
            ("XSYMM", &[1, 5, 6]),
            ("YSYMM", &[2, 4, 6]),
            ("ZSYMM", &[3, 4, 5]),
            ("XASYMM", &[2, 3, 4]),
            ("YASYMM", &[1, 3, 5]),
            ("zasymm", &[1, 2, 6]),
        ];
        for (name, dofs) in expected {
            let mesh = parse(&format!("*NODE\n5, 0, 0, 0\n*BOUNDARY\n5, {name}\n"));
            let fixed: Vec<u8> = mesh.global_bcs.iter().map(|bc| bc.dof).collect();
            assert_eq!(fixed, dofs, "{name}");
            assert!(
                mesh.global_bcs
                    .iter()
                    .all(|bc| bc.node_id == 5 && bc.value == 0.0)
            );
        }
        let on_set = parse(
            "*NODE\n1, 0, 0, 0\n2, 1, 0, 0\n*NSET, NSET=FIX\n1, 2\n*BOUNDARY\nFIX, ENCASTRE\n",
        );
        assert_eq!(on_set.global_bcs.len(), 12);
    }

    #[test]
    fn velocity_boundaries_are_not_displacement_conditions() {
        let mesh = parse(&format!(
            "{NODES}*BOUNDARY, TYPE=VELOCITY\n1, 1, 1, 2.0\n*BOUNDARY, TYPE=DISPLACEMENT\n2, 1\n*BOUNDARY, OP=NEW\n3, 2\n"
        ));
        assert_eq!(
            mesh.global_bcs,
            vec![DirichletBc::fixed(2, 1), DirichletBc::fixed(3, 2)]
        );
    }

    #[test]
    fn malformed_and_unresolvable_boundary_lines_are_reported() {
        for (data, fragment) in [
            ("1, x", "malformed"),
            ("1", "malformed"),
            ("1, 0", "malformed"),
            ("1, 1, 2, abc", "malformed"),
            ("1, 1, 2, NaN", "malformed"),
            (", 1", "malformed"),
            ("MISSING, 1, 3", "undefined node set `MISSING`"),
        ] {
            let inp = format!("{NODES}*BOUNDARY\n{data}\n");
            let text = issues(&inp);
            assert!(
                text.contains("*BOUNDARY") && text.contains(fragment),
                "{data}: {text}"
            );
            assert!(AbaqusInpIo::parse(&inp).global_bcs.is_empty(), "{data}");
        }
    }

    #[test]
    fn cload_creates_one_force_per_node_and_expands_sets() {
        let mesh = parse(&format!(
            "{NODES}*NSET, NSET=top\n3, 4\n*STEP\n*CLOAD\n1, 2, -10.0\ntop, 3, 2.5\n*END STEP\n*STEP\n*CLOAD\n2, 1, 7\n*END STEP\n"
        ));
        assert_eq!(
            mesh.steps[0].forces,
            vec![
                NodalForce::new(1, 2, -10.0),
                NodalForce::new(3, 3, 2.5),
                NodalForce::new(4, 3, 2.5)
            ]
        );
        assert_eq!(mesh.steps[1].forces, vec![NodalForce::new(2, 1, 7.0)]);
    }

    #[test]
    fn cload_outside_a_step_and_malformed_cload_are_reported() {
        let outside = format!("{NODES}*CLOAD\n1, 2, 3.0\n");
        let text = issues(&outside);
        assert!(text.contains("outside a *STEP"), "{text}");
        let mesh = AbaqusInpIo::parse(&outside);
        assert!(mesh.steps.is_empty());
        let text = issues(&format!("{NODES}*STEP\n*CLOAD\n1, 2\n*END STEP\n"));
        assert!(text.contains("*CLOAD: malformed"), "{text}");
        let text = issues(&format!("{NODES}*STEP\n*CLOAD\nNOPE, 2, 1.0\n*END STEP\n"));
        assert!(text.contains("undefined node set `NOPE`"), "{text}");
    }

    #[test]
    fn keywords_are_matched_by_name_not_by_prefix() {
        // `*NODE OUTPUT` and `*ELEMENT OUTPUT` are output requests, not
        // definitions; their data lines must not become nodes or elements.
        let mesh = parse("*NODE OUTPUT\n1, 2, 3, 4\n*ELEMENT OUTPUT\n1, 2, 3, 4\n*NODE FILE\nU\n");
        assert!(mesh.nodes.is_empty());
        assert!(mesh.elements.is_empty());
    }

    #[test]
    fn keyword_lines_may_continue_on_the_next_line() {
        let mesh = parse(
            "*ELEMENT, TYPE=C3D4,\n ELSET=wrapped\n1, 1, 2, 3, 4\n*NSET,\n** a comment in between\n NSET=split\n9\n",
        );
        assert_eq!(mesh.elements[0].element_type, FeElementType::Tet4);
        assert_eq!(mesh.element_sets.get("wrapped"), Some(&vec![1]));
        assert_eq!(mesh.node_sets.get("split"), Some(&vec![9]));
    }

    #[test]
    fn a_complete_cantilever_model() {
        let inp = "\
** cantilever
*HEADING
cantilever beam
*NODE, NSET=ALL
1, 0.0, 0.0, 0.0
2, 1.0, 0.0, 0.0
3, 1.0, 1.0, 0.0
4, 0.0, 1.0, 0.0
5, 0.0, 0.0, 1.0
6, 1.0, 0.0, 1.0
7, 1.0, 1.0, 1.0
8, 0.0, 1.0, 1.0
*ELEMENT, TYPE=C3D8R, ELSET=BEAM
1, 1, 2, 3, 4, 5, 6, 7, 8
*NSET, NSET=FIXED
1, 4, 5, 8
*NSET, NSET=TIP
2, 3, 6, 7
*MATERIAL, NAME=Steel
*ELASTIC
200.0E9, 0.3
*SOLID SECTION, ELSET=BEAM, MATERIAL=Steel
*BOUNDARY
FIXED, ENCASTRE
*STEP, NAME=Load
*STATIC
0.1, 1.0
*CLOAD
TIP, 3, -250.0
*END STEP
";
        let mesh = parse(inp);
        assert_eq!(mesh.nodes.len(), 8);
        assert_eq!(mesh.elements.len(), 1);
        assert_eq!(mesh.elements[0].element_type, FeElementType::Hex8);
        assert_eq!(mesh.node_sets.get("ALL").map(Vec::len), Some(8));
        assert_eq!(mesh.element_sets.get("BEAM"), Some(&vec![1]));
        assert_eq!(mesh.materials.len(), 1);
        assert_eq!(mesh.materials[0].name, "Steel");
        assert_eq!(mesh.global_bcs.len(), 4 * 6);
        assert_eq!(mesh.steps.len(), 1);
        assert_eq!(mesh.steps[0].name, "Load");
        assert_eq!(mesh.steps[0].forces.len(), 4);
        assert!(
            mesh.steps[0]
                .forces
                .iter()
                .all(|f| f.dof == 3 && f.magnitude == -250.0)
        );
    }

    #[test]
    fn lenient_parse_keeps_everything_it_can_read() {
        let inp = format!(
            "{NODES}*BOUNDARY\nBROKEN\n1, 1\nMISSING, 2\n*STEP\n*CLOAD\n1, 2, 5.0\n*END STEP\n"
        );
        let mesh = AbaqusInpIo::parse(&inp);
        assert_eq!(mesh.global_bcs, vec![DirichletBc::fixed(1, 1)]);
        assert_eq!(mesh.steps[0].forces, vec![NodalForce::new(1, 2, 5.0)]);
        let text = issues(&inp);
        assert!(
            text.contains("BROKEN") && text.contains("MISSING"),
            "{text}"
        );
    }

    #[test]
    fn keyword_parameters_are_parsed() {
        let keyword = Keyword::parse("*Nset,  nset = Left Side , generate, INSTANCE=P-1");
        assert_eq!(keyword.name, "*NSET");
        assert_eq!(keyword.value("NSET"), Some("Left Side"));
        assert!(keyword.has("GENERATE"));
        assert_eq!(keyword.value("GENERATE"), None);
        assert_eq!(keyword.value("INSTANCE"), Some("P-1"));
        assert_eq!(Keyword::parse("*end   step").name, "*END STEP");
    }
}
