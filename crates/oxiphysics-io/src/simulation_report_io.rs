// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Simulation report generation: HTML, Markdown, and LaTeX output.
//!
//! Provides structured report data types ([`SimulationReport`], [`ReportSection`],
//! [`DataTable`]) and three writer implementations:
//! - [`HtmlReportWriter`] — styled HTML with embedded CSS
//! - [`MarkdownReportWriter`] — GitHub-flavored Markdown
//! - [`LatexReportWriter`] — LaTeX `article` class with tabular environments

#![allow(dead_code)]

use std::collections::HashMap;
use std::fmt::Write as FmtWrite;

// ---------------------------------------------------------------------------
// ColumnType
// ---------------------------------------------------------------------------

/// Type of data stored in a [`DataTable`] column.
#[derive(Debug, Clone, PartialEq)]
pub enum ColumnType {
    /// Integer column.
    Integer,
    /// Floating-point column.
    Float,
    /// Text column.
    Text,
    /// Boolean (true/false) column.
    Boolean,
}

impl ColumnType {
    /// Return a human-readable label for this column type.
    pub fn label(&self) -> &'static str {
        match self {
            ColumnType::Integer => "integer",
            ColumnType::Float => "float",
            ColumnType::Text => "text",
            ColumnType::Boolean => "boolean",
        }
    }
}

// ---------------------------------------------------------------------------
// CellValue
// ---------------------------------------------------------------------------

/// A single table cell value.
#[derive(Debug, Clone, PartialEq)]
pub enum CellValue {
    /// Integer cell.
    Int(i64),
    /// Float cell with optional precision (decimal places).
    Float(f64, Option<usize>),
    /// Text cell.
    Text(String),
    /// Boolean cell.
    Bool(bool),
    /// Empty / missing cell.
    Empty,
}

impl CellValue {
    /// Render the cell as a plain string.
    pub fn to_display(&self) -> String {
        match self {
            CellValue::Int(v) => v.to_string(),
            CellValue::Float(v, Some(prec)) => format!("{:.prec$}", v, prec = prec),
            CellValue::Float(v, None) => format!("{v:.6}"),
            CellValue::Text(s) => s.clone(),
            CellValue::Bool(b) => if *b { "true" } else { "false" }.to_string(),
            CellValue::Empty => String::new(),
        }
    }

    /// Render the cell as an HTML-safe string (escapes `<`, `>`, `&`).
    pub fn to_html(&self) -> String {
        html_escape(&self.to_display())
    }

    /// Render the cell as a LaTeX-safe string (escapes special characters).
    pub fn to_latex(&self) -> String {
        latex_escape(&self.to_display())
    }
}

// ---------------------------------------------------------------------------
// DataTable
// ---------------------------------------------------------------------------

/// A rectangular table of [`CellValue`]s with a header row and typed columns.
#[derive(Debug, Clone, Default)]
pub struct DataTable {
    /// Table caption / title.
    pub caption: String,
    /// Column header labels.
    pub headers: Vec<String>,
    /// Column types (must match `headers` length if non-empty).
    pub column_types: Vec<ColumnType>,
    /// Data rows; each row must have the same length as `headers`.
    pub rows: Vec<Vec<CellValue>>,
    /// Optional row highlighting: row index → CSS class or LaTeX command.
    pub row_highlights: HashMap<usize, String>,
}

impl DataTable {
    /// Create an empty table with the given headers.
    pub fn new(caption: impl Into<String>, headers: Vec<String>) -> Self {
        let ncols = headers.len();
        Self {
            caption: caption.into(),
            headers,
            column_types: vec![ColumnType::Text; ncols],
            rows: Vec::new(),
            row_highlights: HashMap::new(),
        }
    }

    /// Append a row. Panics in debug builds if the row width doesn't match headers.
    pub fn add_row(&mut self, row: Vec<CellValue>) {
        debug_assert!(
            self.headers.is_empty() || row.len() == self.headers.len(),
            "row width {} != header count {}",
            row.len(),
            self.headers.len()
        );
        self.rows.push(row);
    }

    /// Set the column type for column `idx`.
    pub fn set_column_type(&mut self, idx: usize, ct: ColumnType) {
        if idx < self.column_types.len() {
            self.column_types[idx] = ct;
        }
    }

    /// Highlight row `idx` with a label (e.g. `"highlight"`, `"\\rowcolor{yellow}"`).
    pub fn highlight_row(&mut self, idx: usize, label: impl Into<String>) {
        self.row_highlights.insert(idx, label.into());
    }

    /// Number of data rows.
    pub fn row_count(&self) -> usize {
        self.rows.len()
    }

    /// Number of columns.
    pub fn col_count(&self) -> usize {
        self.headers.len()
    }

    /// Returns `true` if the table has no data rows.
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// Compute the maximum display-string width of each column (including header).
    pub fn column_widths(&self) -> Vec<usize> {
        let n = self.headers.len();
        let mut widths: Vec<usize> = self.headers.iter().map(|h| h.len()).collect();
        for row in &self.rows {
            for (i, cell) in row.iter().enumerate() {
                if i < n {
                    widths[i] = widths[i].max(cell.to_display().len());
                }
            }
        }
        widths
    }
}

// ---------------------------------------------------------------------------
// FigureRef
// ---------------------------------------------------------------------------

/// A reference to an external figure file embedded in a report section.
#[derive(Debug, Clone)]
pub struct FigureRef {
    /// Unique label / id for the figure.
    pub label: String,
    /// File path or URL.
    pub path: String,
    /// Caption text.
    pub caption: String,
    /// Optional alt-text for HTML.
    pub alt: Option<String>,
    /// Optional width hint (e.g. `"80%"`, `"0.7\\linewidth"`).
    pub width: Option<String>,
}

impl FigureRef {
    /// Create a simple figure reference with just a path and caption.
    pub fn new(
        label: impl Into<String>,
        path: impl Into<String>,
        caption: impl Into<String>,
    ) -> Self {
        Self {
            label: label.into(),
            path: path.into(),
            caption: caption.into(),
            alt: None,
            width: None,
        }
    }

    /// Builder: set alt text.
    pub fn with_alt(mut self, alt: impl Into<String>) -> Self {
        self.alt = Some(alt.into());
        self
    }

    /// Builder: set width hint.
    pub fn with_width(mut self, width: impl Into<String>) -> Self {
        self.width = Some(width.into());
        self
    }
}

// ---------------------------------------------------------------------------
// ReportSection
// ---------------------------------------------------------------------------

/// A single section within a [`SimulationReport`].
#[derive(Debug, Clone, Default)]
pub struct ReportSection {
    /// Section title.
    pub title: String,
    /// Introductory / body text (may contain multiple paragraphs separated by `\n\n`).
    pub text: String,
    /// Data tables contained in this section.
    pub tables: Vec<DataTable>,
    /// Figure references in this section.
    pub figures: Vec<FigureRef>,
    /// Subsection level (0 = top-level section, 1 = subsection, …).
    pub level: u8,
}

impl ReportSection {
    /// Create a new top-level section.
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            level: 0,
            ..Default::default()
        }
    }

    /// Create a subsection (level 1).
    pub fn subsection(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            level: 1,
            ..Default::default()
        }
    }

    /// Append body text (separated from existing text by a blank line).
    pub fn append_text(&mut self, text: impl Into<String>) {
        let t = text.into();
        if self.text.is_empty() {
            self.text = t;
        } else {
            self.text.push_str("\n\n");
            self.text.push_str(&t);
        }
    }

    /// Add a data table to this section.
    pub fn add_table(&mut self, table: DataTable) {
        self.tables.push(table);
    }

    /// Add a figure reference to this section.
    pub fn add_figure(&mut self, figure: FigureRef) {
        self.figures.push(figure);
    }
}

// ---------------------------------------------------------------------------
// SimulationReport
// ---------------------------------------------------------------------------

/// A complete simulation report aggregating parameters, results, sections,
/// figures, and tables.
#[derive(Debug, Clone, Default)]
pub struct SimulationReport {
    /// Report title.
    pub title: String,
    /// Short description / abstract.
    pub description: String,
    /// Simulation parameters (key → value string).
    pub parameters: HashMap<String, String>,
    /// Key result values (key → value string).
    pub results: HashMap<String, String>,
    /// Ordered report sections.
    pub sections: Vec<ReportSection>,
    /// Top-level figure references (outside any section).
    pub figures: Vec<FigureRef>,
    /// Top-level data tables (outside any section).
    pub tables: Vec<DataTable>,
    /// Report author.
    pub author: String,
    /// Report date.
    pub date: String,
    /// Optional version string.
    pub version: Option<String>,
}

impl SimulationReport {
    /// Create a new empty report with the given title.
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            ..Default::default()
        }
    }

    /// Set the description / abstract.
    pub fn with_description(mut self, desc: impl Into<String>) -> Self {
        self.description = desc.into();
        self
    }

    /// Set the author.
    pub fn with_author(mut self, author: impl Into<String>) -> Self {
        self.author = author.into();
        self
    }

    /// Set the date string.
    pub fn with_date(mut self, date: impl Into<String>) -> Self {
        self.date = date.into();
        self
    }

    /// Set the version string.
    pub fn with_version(mut self, version: impl Into<String>) -> Self {
        self.version = Some(version.into());
        self
    }

    /// Insert a parameter key-value pair.
    pub fn add_parameter(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.parameters.insert(key.into(), value.into());
    }

    /// Insert a result key-value pair.
    pub fn add_result(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.results.insert(key.into(), value.into());
    }

    /// Append a report section.
    pub fn add_section(&mut self, section: ReportSection) {
        self.sections.push(section);
    }

    /// Append a top-level figure.
    pub fn add_figure(&mut self, figure: FigureRef) {
        self.figures.push(figure);
    }

    /// Append a top-level table.
    pub fn add_table(&mut self, table: DataTable) {
        self.tables.push(table);
    }

    /// Total number of tables across all sections and top-level.
    pub fn total_tables(&self) -> usize {
        self.tables.len() + self.sections.iter().map(|s| s.tables.len()).sum::<usize>()
    }

    /// Total number of figures across all sections and top-level.
    pub fn total_figures(&self) -> usize {
        self.figures.len() + self.sections.iter().map(|s| s.figures.len()).sum::<usize>()
    }
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Escape special HTML characters in `s`.
fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

/// Escape LaTeX special characters.
fn latex_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 8);
    for ch in s.chars() {
        match ch {
            '&' => out.push_str(r"\&"),
            '%' => out.push_str(r"\%"),
            '$' => out.push_str(r"\$"),
            '#' => out.push_str(r"\#"),
            '_' => out.push_str(r"\_"),
            '{' => out.push_str(r"\{"),
            '}' => out.push_str(r"\}"),
            '~' => out.push_str(r"\textasciitilde{}"),
            '^' => out.push_str(r"\textasciicircum{}"),
            '\\' => out.push_str(r"\textbackslash{}"),
            other => out.push(other),
        }
    }
    out
}

/// Build a Markdown table string from a [`DataTable`].
fn markdown_table(table: &DataTable) -> String {
    let mut out = String::new();
    // Caption
    if !table.caption.is_empty() {
        let _ = writeln!(out, "**{}**\n", table.caption);
    }
    if table.headers.is_empty() {
        return out;
    }
    // Header row
    let _ = write!(out, "|");
    for h in &table.headers {
        let _ = write!(out, " {} |", h);
    }
    out.push('\n');
    // Separator
    let _ = write!(out, "|");
    for _ in &table.headers {
        let _ = write!(out, " --- |");
    }
    out.push('\n');
    // Data rows
    for row in &table.rows {
        let _ = write!(out, "|");
        for cell in row {
            let _ = write!(out, " {} |", cell.to_display());
        }
        out.push('\n');
    }
    out
}

/// Build an HTML `<table>` string from a [`DataTable`].
fn html_table(table: &DataTable) -> String {
    let mut out = String::new();
    let _ = write!(out, "<table class=\"data-table\">");
    if !table.caption.is_empty() {
        let _ = write!(out, "<caption>{}</caption>", html_escape(&table.caption));
    }
    // thead
    if !table.headers.is_empty() {
        let _ = write!(out, "<thead><tr>");
        for h in &table.headers {
            let _ = write!(out, "<th>{}</th>", html_escape(h));
        }
        let _ = write!(out, "</tr></thead>");
    }
    // tbody
    let _ = write!(out, "<tbody>");
    for (i, row) in table.rows.iter().enumerate() {
        let class = if let Some(lbl) = table.row_highlights.get(&i) {
            format!(" class=\"{}\"", html_escape(lbl))
        } else {
            String::new()
        };
        let _ = write!(out, "<tr{}>", class);
        for cell in row {
            let _ = write!(out, "<td>{}</td>", cell.to_html());
        }
        let _ = write!(out, "</tr>");
    }
    let _ = write!(out, "</tbody></table>");
    out
}

/// Build a LaTeX `tabular` environment from a [`DataTable`].
fn latex_tabular(table: &DataTable) -> String {
    let ncols = table.headers.len();
    let col_spec: String = std::iter::repeat_n("l", ncols)
        .collect::<Vec<_>>()
        .join("|");
    let mut out = String::new();
    if !table.caption.is_empty() {
        let _ = writeln!(out, "\\begin{{table}}[htbp]");
        let _ = writeln!(out, "\\centering");
        let _ = writeln!(out, "\\caption{{{}}}", latex_escape(&table.caption));
    }
    let _ = writeln!(out, "\\begin{{tabular}}{{|{}|}}", col_spec);
    let _ = writeln!(out, "\\hline");
    // Header
    if !table.headers.is_empty() {
        let header_row = table
            .headers
            .iter()
            .map(|h| format!("\\textbf{{{}}}", latex_escape(h)))
            .collect::<Vec<_>>()
            .join(" & ");
        let _ = writeln!(out, "{} \\\\", header_row);
        let _ = writeln!(out, "\\hline");
    }
    // Data rows
    for row in &table.rows {
        let row_str = row
            .iter()
            .map(|c| c.to_latex())
            .collect::<Vec<_>>()
            .join(" & ");
        let _ = writeln!(out, "{} \\\\", row_str);
    }
    let _ = writeln!(out, "\\hline");
    let _ = writeln!(out, "\\end{{tabular}}");
    if !table.caption.is_empty() {
        let _ = writeln!(out, "\\end{{table}}");
    }
    out
}

// ---------------------------------------------------------------------------
// HtmlReportWriter
// ---------------------------------------------------------------------------

/// Writes a [`SimulationReport`] as a self-contained HTML document with
/// embedded CSS styling.
///
/// Usage:
/// ```no_run
/// use oxiphysics_io::simulation_report_io::{SimulationReport, HtmlReportWriter};
/// let report = SimulationReport::new("My Simulation");
/// let html = HtmlReportWriter::new().write(&report);
/// assert!(html.contains("<html"));
/// ```
#[derive(Debug, Clone, Default)]
pub struct HtmlReportWriter {
    /// Whether to embed the default CSS stylesheet.
    pub embed_css: bool,
    /// Custom CSS to append after the default stylesheet.
    pub extra_css: String,
}

impl HtmlReportWriter {
    /// Create a new writer with embedded CSS enabled.
    pub fn new() -> Self {
        Self {
            embed_css: true,
            extra_css: String::new(),
        }
    }

    /// Builder: disable embedded CSS.
    pub fn without_css(mut self) -> Self {
        self.embed_css = false;
        self
    }

    /// Builder: append extra CSS rules.
    pub fn with_extra_css(mut self, css: impl Into<String>) -> Self {
        self.extra_css = css.into();
        self
    }

    /// Render the report to an HTML string.
    pub fn write(&self, report: &SimulationReport) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "<!DOCTYPE html>");
        let _ = writeln!(out, "<html lang=\"en\">");
        let _ = writeln!(out, "<head>");
        let _ = writeln!(out, "<meta charset=\"UTF-8\">");
        let _ = writeln!(
            out,
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"
        );
        let _ = writeln!(out, "<title>{}</title>", html_escape(&report.title));
        if self.embed_css {
            let _ = writeln!(out, "<style>");
            out.push_str(DEFAULT_CSS);
            if !self.extra_css.is_empty() {
                out.push_str(&self.extra_css);
            }
            let _ = writeln!(out, "</style>");
        }
        let _ = writeln!(out, "</head>");
        let _ = writeln!(out, "<body>");
        let _ = writeln!(out, "<div class=\"container\">");

        // Title block
        let _ = writeln!(
            out,
            "<h1 class=\"report-title\">{}</h1>",
            html_escape(&report.title)
        );
        if !report.author.is_empty() || !report.date.is_empty() {
            let _ = write!(out, "<p class=\"report-meta\">");
            if !report.author.is_empty() {
                let _ = write!(
                    out,
                    "<strong>Author:</strong> {} ",
                    html_escape(&report.author)
                );
            }
            if !report.date.is_empty() {
                let _ = write!(out, "<strong>Date:</strong> {}", html_escape(&report.date));
            }
            if let Some(ref ver) = report.version {
                let _ = write!(out, " <strong>Version:</strong> {}", html_escape(ver));
            }
            let _ = writeln!(out, "</p>");
        }

        // Abstract
        if !report.description.is_empty() {
            let _ = writeln!(out, "<div class=\"abstract\">");
            let _ = writeln!(out, "<h2>Abstract</h2>");
            let _ = writeln!(out, "<p>{}</p>", html_escape(&report.description));
            let _ = writeln!(out, "</div>");
        }

        // Parameters table
        if !report.parameters.is_empty() {
            let _ = writeln!(out, "<h2>Parameters</h2>");
            let _ = writeln!(out, "<table class=\"data-table\">");
            let _ = writeln!(
                out,
                "<thead><tr><th>Parameter</th><th>Value</th></tr></thead>"
            );
            let _ = writeln!(out, "<tbody>");
            let mut params: Vec<(&String, &String)> = report.parameters.iter().collect();
            params.sort_by_key(|(k, _)| k.as_str());
            for (k, v) in params {
                let _ = writeln!(
                    out,
                    "<tr><td>{}</td><td>{}</td></tr>",
                    html_escape(k),
                    html_escape(v)
                );
            }
            let _ = writeln!(out, "</tbody></table>");
        }

        // Results table
        if !report.results.is_empty() {
            let _ = writeln!(out, "<h2>Results</h2>");
            let _ = writeln!(out, "<table class=\"data-table\">");
            let _ = writeln!(out, "<thead><tr><th>Metric</th><th>Value</th></tr></thead>");
            let _ = writeln!(out, "<tbody>");
            let mut res: Vec<(&String, &String)> = report.results.iter().collect();
            res.sort_by_key(|(k, _)| k.as_str());
            for (k, v) in res {
                let _ = writeln!(
                    out,
                    "<tr><td>{}</td><td>{}</td></tr>",
                    html_escape(k),
                    html_escape(v)
                );
            }
            let _ = writeln!(out, "</tbody></table>");
        }

        // Top-level tables
        for table in &report.tables {
            out.push_str(&html_table(table));
            out.push('\n');
        }

        // Top-level figures
        for fig in &report.figures {
            self.write_html_figure(&mut out, fig);
        }

        // Sections
        for section in &report.sections {
            self.write_html_section(&mut out, section);
        }

        let _ = writeln!(out, "</div>"); // container
        let _ = writeln!(out, "</body></html>");
        out
    }

    fn write_html_section(&self, out: &mut String, section: &ReportSection) {
        let tag = if section.level == 0 { "h2" } else { "h3" };
        let _ = writeln!(out, "<{tag}>{}</{tag}>", html_escape(&section.title));
        if !section.text.is_empty() {
            for para in section.text.split("\n\n") {
                let _ = writeln!(out, "<p>{}</p>", html_escape(para.trim()));
            }
        }
        for table in &section.tables {
            out.push_str(&html_table(table));
            out.push('\n');
        }
        for fig in &section.figures {
            self.write_html_figure(out, fig);
        }
    }

    fn write_html_figure(&self, out: &mut String, fig: &FigureRef) {
        let _ = writeln!(out, "<figure id=\"{}\">", html_escape(&fig.label));
        let alt = fig.alt.as_deref().unwrap_or(&fig.caption);
        let width_attr = if let Some(ref w) = fig.width {
            format!(" width=\"{}\"", html_escape(w))
        } else {
            String::new()
        };
        let _ = writeln!(
            out,
            "<img src=\"{}\" alt=\"{}\"{} />",
            html_escape(&fig.path),
            html_escape(alt),
            width_attr
        );
        let _ = writeln!(
            out,
            "<figcaption>{}</figcaption>",
            html_escape(&fig.caption)
        );
        let _ = writeln!(out, "</figure>");
    }
}

/// Default CSS embedded in HTML reports.
const DEFAULT_CSS: &str = r#"
body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
    line-height: 1.6;
    color: #333;
    background: #fafafa;
    margin: 0;
    padding: 0;
}
.container {
    max-width: 960px;
    margin: 2rem auto;
    background: #fff;
    padding: 2rem 3rem;
    box-shadow: 0 2px 8px rgba(0,0,0,.12);
    border-radius: 6px;
}
.report-title {
    font-size: 2rem;
    border-bottom: 2px solid #3b82f6;
    padding-bottom: .4rem;
    margin-bottom: .5rem;
}
.report-meta {
    color: #555;
    font-size: .9rem;
    margin-bottom: 1.5rem;
}
.abstract {
    background: #eff6ff;
    border-left: 4px solid #3b82f6;
    padding: .8rem 1rem;
    margin-bottom: 1.5rem;
    border-radius: 0 4px 4px 0;
}
h2 { font-size: 1.4rem; margin-top: 2rem; border-bottom: 1px solid #e5e7eb; padding-bottom: .3rem; }
h3 { font-size: 1.15rem; margin-top: 1.5rem; }
table.data-table {
    border-collapse: collapse;
    width: 100%;
    margin: 1rem 0;
    font-size: .9rem;
}
table.data-table th, table.data-table td {
    border: 1px solid #d1d5db;
    padding: .5rem .75rem;
    text-align: left;
}
table.data-table thead tr { background: #3b82f6; color: #fff; }
table.data-table tbody tr:nth-child(even) { background: #f3f4f6; }
table.data-table caption { caption-side: top; font-weight: bold; margin-bottom: .4rem; }
table.data-table tr.highlight { background: #fef9c3; }
figure { margin: 1.5rem 0; text-align: center; }
figure img { max-width: 100%; border-radius: 4px; }
figcaption { color: #555; font-size: .85rem; margin-top: .4rem; font-style: italic; }
"#;

// ---------------------------------------------------------------------------
// MarkdownReportWriter
// ---------------------------------------------------------------------------

/// Writes a [`SimulationReport`] as GitHub-flavored Markdown.
///
/// Usage:
/// ```no_run
/// use oxiphysics_io::simulation_report_io::{SimulationReport, MarkdownReportWriter};
/// let report = SimulationReport::new("My Simulation");
/// let md = MarkdownReportWriter::new().write(&report);
/// assert!(md.contains("# My Simulation"));
/// ```
#[derive(Debug, Clone, Default)]
pub struct MarkdownReportWriter {
    /// Whether to add a horizontal rule between sections.
    pub section_dividers: bool,
}

impl MarkdownReportWriter {
    /// Create a new writer with section dividers enabled.
    pub fn new() -> Self {
        Self {
            section_dividers: true,
        }
    }

    /// Render the report to a Markdown string.
    pub fn write(&self, report: &SimulationReport) -> String {
        let mut out = String::new();

        // Title
        let _ = writeln!(out, "# {}", report.title);
        out.push('\n');

        // Meta
        if !report.author.is_empty() {
            let _ = writeln!(out, "**Author:** {}  ", report.author);
        }
        if !report.date.is_empty() {
            let _ = writeln!(out, "**Date:** {}  ", report.date);
        }
        if let Some(ref ver) = report.version {
            let _ = writeln!(out, "**Version:** {}  ", ver);
        }
        if !report.author.is_empty() || !report.date.is_empty() {
            out.push('\n');
        }

        // Abstract
        if !report.description.is_empty() {
            let _ = writeln!(out, "## Abstract");
            out.push('\n');
            let _ = writeln!(out, "{}", report.description);
            out.push('\n');
        }

        // Parameters
        if !report.parameters.is_empty() {
            let _ = writeln!(out, "## Parameters");
            out.push('\n');
            let _ = writeln!(out, "| Parameter | Value |");
            let _ = writeln!(out, "| --- | --- |");
            let mut params: Vec<(&String, &String)> = report.parameters.iter().collect();
            params.sort_by_key(|(k, _)| k.as_str());
            for (k, v) in params {
                let _ = writeln!(out, "| {} | {} |", k, v);
            }
            out.push('\n');
        }

        // Results
        if !report.results.is_empty() {
            let _ = writeln!(out, "## Results");
            out.push('\n');
            let _ = writeln!(out, "| Metric | Value |");
            let _ = writeln!(out, "| --- | --- |");
            let mut res: Vec<(&String, &String)> = report.results.iter().collect();
            res.sort_by_key(|(k, _)| k.as_str());
            for (k, v) in res {
                let _ = writeln!(out, "| {} | {} |", k, v);
            }
            out.push('\n');
        }

        // Top-level tables
        for table in &report.tables {
            out.push_str(&markdown_table(table));
            out.push('\n');
        }

        // Top-level figures
        for fig in &report.figures {
            let alt = fig.alt.as_deref().unwrap_or(&fig.caption);
            let _ = writeln!(out, "![{}]({})", alt, fig.path);
            let _ = writeln!(out, "*Figure: {}*", fig.caption);
            out.push('\n');
        }

        // Sections
        for section in &report.sections {
            if self.section_dividers {
                let _ = writeln!(out, "---");
                out.push('\n');
            }
            self.write_md_section(&mut out, section);
        }

        out
    }

    fn write_md_section(&self, out: &mut String, section: &ReportSection) {
        let prefix = if section.level == 0 { "##" } else { "###" };
        let _ = writeln!(out, "{} {}", prefix, section.title);
        out.push('\n');

        if !section.text.is_empty() {
            let _ = writeln!(out, "{}", section.text);
            out.push('\n');
        }

        for table in &section.tables {
            out.push_str(&markdown_table(table));
            out.push('\n');
        }

        for fig in &section.figures {
            let alt = fig.alt.as_deref().unwrap_or(&fig.caption);
            let _ = writeln!(out, "![{}]({})", alt, fig.path);
            let _ = writeln!(out, "*Figure: {}*", fig.caption);
            out.push('\n');
        }
    }
}

// ---------------------------------------------------------------------------
// LatexReportWriter
// ---------------------------------------------------------------------------

/// Writes a [`SimulationReport`] as a compilable LaTeX `article` document.
///
/// Usage:
/// ```no_run
/// use oxiphysics_io::simulation_report_io::{SimulationReport, LatexReportWriter};
/// let report = SimulationReport::new("My Simulation");
/// let tex = LatexReportWriter::new().write(&report);
/// assert!(tex.contains("\\documentclass"));
/// ```
#[derive(Debug, Clone)]
pub struct LatexReportWriter {
    /// Paper size option for `geometry` package (e.g. `"a4paper"`).
    pub paper: String,
    /// Font size (e.g. `"11pt"`).
    pub font_size: String,
    /// Whether to include `\usepackage{booktabs}` and use `\toprule` etc.
    pub use_booktabs: bool,
    /// Whether to include the `graphicx` package.
    pub use_graphicx: bool,
}

impl Default for LatexReportWriter {
    fn default() -> Self {
        Self {
            paper: "a4paper".to_string(),
            font_size: "11pt".to_string(),
            use_booktabs: true,
            use_graphicx: true,
        }
    }
}

impl LatexReportWriter {
    /// Create a new writer with sensible defaults.
    pub fn new() -> Self {
        Self::default()
    }

    /// Render the report to a LaTeX string.
    pub fn write(&self, report: &SimulationReport) -> String {
        let mut out = String::new();

        // Preamble
        let _ = writeln!(
            out,
            "\\documentclass[{},{}]{{article}}",
            self.font_size, self.paper
        );
        let _ = writeln!(out, "\\usepackage[utf8]{{inputenc}}");
        let _ = writeln!(out, "\\usepackage[T1]{{fontenc}}");
        let _ = writeln!(out, "\\usepackage{{geometry}}");
        let _ = writeln!(out, "\\geometry{{margin=2.5cm}}");
        let _ = writeln!(out, "\\usepackage{{hyperref}}");
        if self.use_booktabs {
            let _ = writeln!(out, "\\usepackage{{booktabs}}");
        }
        if self.use_graphicx {
            let _ = writeln!(out, "\\usepackage{{graphicx}}");
        }
        let _ = writeln!(out, "\\usepackage{{array}}");
        let _ = writeln!(out, "\\usepackage{{colortbl}}");
        let _ = writeln!(out, "\\usepackage{{xcolor}}");
        out.push('\n');

        // Title, author, date
        let _ = writeln!(out, "\\title{{{}}}", latex_escape(&report.title));
        if !report.author.is_empty() {
            let _ = writeln!(out, "\\author{{{}}}", latex_escape(&report.author));
        }
        if !report.date.is_empty() {
            let _ = writeln!(out, "\\date{{{}}}", latex_escape(&report.date));
        } else {
            let _ = writeln!(out, "\\date{{\\today}}");
        }
        out.push('\n');

        // Begin document
        let _ = writeln!(out, "\\begin{{document}}");
        let _ = writeln!(out, "\\maketitle");
        if let Some(ref ver) = report.version {
            let _ = writeln!(
                out,
                "\\begin{{center}}\\small Version: {}\\end{{center}}",
                latex_escape(ver)
            );
        }
        out.push('\n');

        // Abstract
        if !report.description.is_empty() {
            let _ = writeln!(out, "\\begin{{abstract}}");
            let _ = writeln!(out, "{}", latex_escape(&report.description));
            let _ = writeln!(out, "\\end{{abstract}}");
            out.push('\n');
        }

        // Table of contents (optional, include for longer reports)
        let _ = writeln!(out, "\\tableofcontents");
        let _ = writeln!(out, "\\newpage");
        out.push('\n');

        // Parameters
        if !report.parameters.is_empty() {
            let _ = writeln!(out, "\\section{{Parameters}}");
            let mut params: Vec<(&String, &String)> = report.parameters.iter().collect();
            params.sort_by_key(|(k, _)| k.as_str());
            let _ = writeln!(out, "\\begin{{table}}[htbp]\\centering");
            let _ = writeln!(out, "\\caption{{Simulation Parameters}}");
            let _ = writeln!(out, "\\begin{{tabular}}{{|l|l|}}\\hline");
            let _ = writeln!(out, "\\textbf{{Parameter}} & \\textbf{{Value}} \\\\\\hline");
            for (k, v) in params {
                let _ = writeln!(out, "{} & {} \\\\", latex_escape(k), latex_escape(v));
            }
            let _ = writeln!(out, "\\hline\\end{{tabular}}\\end{{table}}");
            out.push('\n');
        }

        // Results
        if !report.results.is_empty() {
            let _ = writeln!(out, "\\section{{Results Summary}}");
            let mut res: Vec<(&String, &String)> = report.results.iter().collect();
            res.sort_by_key(|(k, _)| k.as_str());
            let _ = writeln!(out, "\\begin{{table}}[htbp]\\centering");
            let _ = writeln!(out, "\\caption{{Result Metrics}}");
            let _ = writeln!(out, "\\begin{{tabular}}{{|l|l|}}\\hline");
            let _ = writeln!(out, "\\textbf{{Metric}} & \\textbf{{Value}} \\\\\\hline");
            for (k, v) in res {
                let _ = writeln!(out, "{} & {} \\\\", latex_escape(k), latex_escape(v));
            }
            let _ = writeln!(out, "\\hline\\end{{tabular}}\\end{{table}}");
            out.push('\n');
        }

        // Top-level tables
        for table in &report.tables {
            out.push_str(&latex_tabular(table));
            out.push('\n');
        }

        // Top-level figures
        for fig in &report.figures {
            self.write_latex_figure(&mut out, fig);
        }

        // Sections
        for section in &report.sections {
            self.write_latex_section(&mut out, section);
        }

        let _ = writeln!(out, "\\end{{document}}");
        out
    }

    fn write_latex_section(&self, out: &mut String, section: &ReportSection) {
        let cmd = if section.level == 0 {
            "\\section"
        } else {
            "\\subsection"
        };
        let _ = writeln!(out, "{}{{{}}}", cmd, latex_escape(&section.title));
        out.push('\n');

        if !section.text.is_empty() {
            for para in section.text.split("\n\n") {
                let _ = writeln!(out, "{}\n", latex_escape(para.trim()));
            }
        }

        for table in &section.tables {
            out.push_str(&latex_tabular(table));
            out.push('\n');
        }

        for fig in &section.figures {
            self.write_latex_figure(out, fig);
        }
    }

    fn write_latex_figure(&self, out: &mut String, fig: &FigureRef) {
        let width = fig.width.as_deref().unwrap_or("0.8\\linewidth");
        let _ = writeln!(out, "\\begin{{figure}}[htbp]");
        let _ = writeln!(out, "\\centering");
        let _ = writeln!(
            out,
            "\\includegraphics[width={}]{{{}}}",
            width,
            latex_escape(&fig.path)
        );
        let _ = writeln!(out, "\\caption{{{}}}", latex_escape(&fig.caption));
        let _ = writeln!(out, "\\label{{fig:{}}}", latex_escape(&fig.label));
        let _ = writeln!(out, "\\end{{figure}}");
    }
}

// ---------------------------------------------------------------------------
// ReportBuilder — convenience builder
// ---------------------------------------------------------------------------

/// Fluent builder for constructing a [`SimulationReport`] and immediately
/// rendering it to a target format.
#[derive(Debug, Default)]
pub struct ReportBuilder {
    report: SimulationReport,
}

impl ReportBuilder {
    /// Create a new builder.
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            report: SimulationReport::new(title),
        }
    }

    /// Set the author.
    pub fn author(mut self, author: impl Into<String>) -> Self {
        self.report.author = author.into();
        self
    }

    /// Set the date.
    pub fn date(mut self, date: impl Into<String>) -> Self {
        self.report.date = date.into();
        self
    }

    /// Set the description.
    pub fn description(mut self, desc: impl Into<String>) -> Self {
        self.report.description = desc.into();
        self
    }

    /// Add a parameter.
    pub fn parameter(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.report.add_parameter(key, value);
        self
    }

    /// Add a result metric.
    pub fn result(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.report.add_result(key, value);
        self
    }

    /// Add a section.
    pub fn section(mut self, section: ReportSection) -> Self {
        self.report.add_section(section);
        self
    }

    /// Finalize and return the [`SimulationReport`].
    pub fn build(self) -> SimulationReport {
        self.report
    }

    /// Finalize and render as HTML.
    pub fn to_html(self) -> String {
        HtmlReportWriter::new().write(&self.report)
    }

    /// Finalize and render as Markdown.
    pub fn to_markdown(self) -> String {
        MarkdownReportWriter::new().write(&self.report)
    }

    /// Finalize and render as LaTeX.
    pub fn to_latex(self) -> String {
        LatexReportWriter::new().write(&self.report)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── CellValue ─────────────────────────────────────────────────────────────

    #[test]
    fn test_cell_int_display() {
        assert_eq!(CellValue::Int(42).to_display(), "42");
    }

    #[test]
    fn test_cell_float_default_precision() {
        let s = CellValue::Float(2.54321, None).to_display();
        assert!(s.starts_with("2.54321"), "got {s}");
    }

    #[test]
    fn test_cell_float_custom_precision() {
        let s = CellValue::Float(2.54321, Some(2)).to_display();
        assert_eq!(s, "2.54");
    }

    #[test]
    fn test_cell_text_display() {
        assert_eq!(CellValue::Text("hello".into()).to_display(), "hello");
    }

    #[test]
    fn test_cell_bool_true() {
        assert_eq!(CellValue::Bool(true).to_display(), "true");
    }

    #[test]
    fn test_cell_bool_false() {
        assert_eq!(CellValue::Bool(false).to_display(), "false");
    }

    #[test]
    fn test_cell_empty_display() {
        assert_eq!(CellValue::Empty.to_display(), "");
    }

    #[test]
    fn test_cell_html_escape() {
        let s = CellValue::Text("<b>bold & italic</b>".into()).to_html();
        assert!(s.contains("&lt;"), "got {s}");
        assert!(s.contains("&amp;"), "got {s}");
    }

    #[test]
    fn test_cell_latex_escape() {
        let s = CellValue::Text("x_1 & y_2".into()).to_latex();
        assert!(s.contains(r"\_"), "got {s}");
        assert!(s.contains(r"\&"), "got {s}");
    }

    // ── DataTable ──────────────────────────────────────────────────────────────

    #[test]
    fn test_data_table_add_row() {
        let mut t = DataTable::new("T", vec!["A".into(), "B".into()]);
        t.add_row(vec![CellValue::Int(1), CellValue::Int(2)]);
        assert_eq!(t.row_count(), 1);
    }

    #[test]
    fn test_data_table_col_count() {
        let t = DataTable::new("T", vec!["A".into(), "B".into(), "C".into()]);
        assert_eq!(t.col_count(), 3);
    }

    #[test]
    fn test_data_table_is_empty() {
        let t = DataTable::new("T", vec!["A".into()]);
        assert!(t.is_empty());
    }

    #[test]
    fn test_data_table_highlight_row() {
        let mut t = DataTable::new("T", vec!["A".into()]);
        t.add_row(vec![CellValue::Int(1)]);
        t.highlight_row(0, "highlight");
        assert_eq!(t.row_highlights.get(&0).unwrap(), "highlight");
    }

    #[test]
    fn test_data_table_column_widths() {
        let mut t = DataTable::new("T", vec!["Name".into(), "Val".into()]);
        t.add_row(vec![CellValue::Text("Alice".into()), CellValue::Int(100)]);
        let widths = t.column_widths();
        assert_eq!(widths[0], 5); // "Alice"
        assert!(widths[1] >= 3); // "Val" or "100"
    }

    #[test]
    fn test_data_table_set_column_type() {
        let mut t = DataTable::new("T", vec!["A".into(), "B".into()]);
        t.set_column_type(0, ColumnType::Integer);
        assert_eq!(t.column_types[0], ColumnType::Integer);
    }

    // ── FigureRef ─────────────────────────────────────────────────────────────

    #[test]
    fn test_figure_ref_new() {
        let f = FigureRef::new("fig1", "plot.png", "My Plot");
        assert_eq!(f.label, "fig1");
        assert_eq!(f.path, "plot.png");
        assert_eq!(f.caption, "My Plot");
        assert!(f.alt.is_none());
    }

    #[test]
    fn test_figure_ref_with_alt_and_width() {
        let f = FigureRef::new("f", "a.png", "cap")
            .with_alt("alt text")
            .with_width("80%");
        assert_eq!(f.alt.as_deref(), Some("alt text"));
        assert_eq!(f.width.as_deref(), Some("80%"));
    }

    // ── ReportSection ─────────────────────────────────────────────────────────

    #[test]
    fn test_report_section_new() {
        let s = ReportSection::new("Intro");
        assert_eq!(s.title, "Intro");
        assert_eq!(s.level, 0);
        assert!(s.tables.is_empty());
    }

    #[test]
    fn test_report_section_subsection_level() {
        let s = ReportSection::subsection("Sub");
        assert_eq!(s.level, 1);
    }

    #[test]
    fn test_report_section_append_text() {
        let mut s = ReportSection::new("S");
        s.append_text("Para 1.");
        s.append_text("Para 2.");
        assert!(s.text.contains("Para 1."));
        assert!(s.text.contains("Para 2."));
        assert!(s.text.contains("\n\n"));
    }

    #[test]
    fn test_report_section_add_table() {
        let mut s = ReportSection::new("S");
        s.add_table(DataTable::new("T", vec!["X".into()]));
        assert_eq!(s.tables.len(), 1);
    }

    #[test]
    fn test_report_section_add_figure() {
        let mut s = ReportSection::new("S");
        s.add_figure(FigureRef::new("f", "a.png", "cap"));
        assert_eq!(s.figures.len(), 1);
    }

    // ── SimulationReport ──────────────────────────────────────────────────────

    #[test]
    fn test_report_new() {
        let r = SimulationReport::new("Test Report");
        assert_eq!(r.title, "Test Report");
        assert!(r.parameters.is_empty());
        assert!(r.results.is_empty());
    }

    #[test]
    fn test_report_builder_chain() {
        let r = SimulationReport::new("T")
            .with_description("desc")
            .with_author("Alice")
            .with_date("2026-01-01")
            .with_version("1.0.0");
        assert_eq!(r.description, "desc");
        assert_eq!(r.author, "Alice");
        assert_eq!(r.version.as_deref(), Some("1.0.0"));
    }

    #[test]
    fn test_report_add_parameter() {
        let mut r = SimulationReport::new("T");
        r.add_parameter("dt", "0.001");
        assert_eq!(r.parameters.get("dt").map(|s| s.as_str()), Some("0.001"));
    }

    #[test]
    fn test_report_add_result() {
        let mut r = SimulationReport::new("T");
        r.add_result("energy", "42.5");
        assert_eq!(r.results.get("energy").map(|s| s.as_str()), Some("42.5"));
    }

    #[test]
    fn test_report_total_tables() {
        let mut r = SimulationReport::new("T");
        r.add_table(DataTable::new("T1", vec![]));
        let mut sec = ReportSection::new("S");
        sec.add_table(DataTable::new("T2", vec![]));
        r.add_section(sec);
        assert_eq!(r.total_tables(), 2);
    }

    #[test]
    fn test_report_total_figures() {
        let mut r = SimulationReport::new("T");
        r.add_figure(FigureRef::new("f1", "a.png", "A"));
        let mut sec = ReportSection::new("S");
        sec.add_figure(FigureRef::new("f2", "b.png", "B"));
        r.add_section(sec);
        assert_eq!(r.total_figures(), 2);
    }

    // ── HtmlReportWriter ──────────────────────────────────────────────────────

    #[test]
    fn test_html_writer_contains_doctype() {
        let r = SimulationReport::new("Test");
        let html = HtmlReportWriter::new().write(&r);
        assert!(html.contains("<!DOCTYPE html>"), "no doctype in:\n{html}");
    }

    #[test]
    fn test_html_writer_contains_title() {
        let r = SimulationReport::new("My Report");
        let html = HtmlReportWriter::new().write(&r);
        assert!(html.contains("My Report"), "title missing");
    }

    #[test]
    fn test_html_writer_contains_css() {
        let r = SimulationReport::new("T");
        let html = HtmlReportWriter::new().write(&r);
        assert!(html.contains("<style>"), "CSS missing");
    }

    #[test]
    fn test_html_writer_no_css_when_disabled() {
        let r = SimulationReport::new("T");
        let html = HtmlReportWriter::new().without_css().write(&r);
        assert!(!html.contains("<style>"), "CSS should be absent");
    }

    #[test]
    fn test_html_writer_parameters_table() {
        let mut r = SimulationReport::new("T");
        r.add_parameter("alpha", "0.5");
        let html = HtmlReportWriter::new().write(&r);
        assert!(html.contains("alpha"), "param key missing");
        assert!(html.contains("0.5"), "param value missing");
    }

    #[test]
    fn test_html_writer_results_table() {
        let mut r = SimulationReport::new("T");
        r.add_result("max_stress", "1.23e9");
        let html = HtmlReportWriter::new().write(&r);
        assert!(html.contains("max_stress"));
        assert!(html.contains("1.23e9"));
    }

    #[test]
    fn test_html_writer_section_heading() {
        let mut r = SimulationReport::new("T");
        r.add_section(ReportSection::new("Methodology"));
        let html = HtmlReportWriter::new().write(&r);
        assert!(html.contains("Methodology"));
    }

    #[test]
    fn test_html_writer_figure() {
        let mut r = SimulationReport::new("T");
        r.add_figure(FigureRef::new("f1", "energy.png", "Energy over time"));
        let html = HtmlReportWriter::new().write(&r);
        assert!(html.contains("energy.png"));
        assert!(html.contains("Energy over time"));
    }

    #[test]
    fn test_html_writer_data_table() {
        let mut r = SimulationReport::new("T");
        let mut tbl = DataTable::new("Results", vec!["Step".into(), "E".into()]);
        tbl.add_row(vec![CellValue::Int(1), CellValue::Float(1.5, Some(2))]);
        r.add_table(tbl);
        let html = HtmlReportWriter::new().write(&r);
        assert!(html.contains("<table"), "table tag missing");
        assert!(html.contains("Step"));
        assert!(html.contains("1.50"));
    }

    #[test]
    fn test_html_writer_html_escape_in_title() {
        let r = SimulationReport::new("Report <1> & <2>");
        let html = HtmlReportWriter::new().write(&r);
        assert!(html.contains("&lt;1&gt;"), "angle brackets not escaped");
        assert!(html.contains("&amp;"), "ampersand not escaped");
    }

    // ── MarkdownReportWriter ──────────────────────────────────────────────────

    #[test]
    fn test_md_writer_title_heading() {
        let r = SimulationReport::new("My Sim");
        let md = MarkdownReportWriter::new().write(&r);
        assert!(md.starts_with("# My Sim"), "title heading missing:\n{md}");
    }

    #[test]
    fn test_md_writer_parameters_table() {
        let mut r = SimulationReport::new("T");
        r.add_parameter("dt", "0.01");
        let md = MarkdownReportWriter::new().write(&r);
        assert!(md.contains("## Parameters"), "header missing");
        assert!(md.contains("| dt |"), "param row missing");
    }

    #[test]
    fn test_md_writer_results_table() {
        let mut r = SimulationReport::new("T");
        r.add_result("energy", "5.0");
        let md = MarkdownReportWriter::new().write(&r);
        assert!(md.contains("## Results"));
        assert!(md.contains("| energy |"));
    }

    #[test]
    fn test_md_writer_section_heading() {
        let mut r = SimulationReport::new("T");
        r.add_section(ReportSection::new("Discussion"));
        let md = MarkdownReportWriter::new().write(&r);
        assert!(md.contains("## Discussion"));
    }

    #[test]
    fn test_md_writer_subsection_heading() {
        let mut r = SimulationReport::new("T");
        r.add_section(ReportSection::subsection("Sub Details"));
        let md = MarkdownReportWriter::new().write(&r);
        assert!(md.contains("### Sub Details"));
    }

    #[test]
    fn test_md_writer_dividers() {
        let mut r = SimulationReport::new("T");
        r.add_section(ReportSection::new("A"));
        r.add_section(ReportSection::new("B"));
        let md = MarkdownReportWriter::new().write(&r);
        assert!(md.contains("---"), "divider missing");
    }

    #[test]
    fn test_md_writer_no_dividers_when_disabled() {
        let mut r = SimulationReport::new("T");
        r.add_section(ReportSection::new("A"));
        let w = MarkdownReportWriter {
            section_dividers: false,
        };
        let md = w.write(&r);
        assert!(!md.contains("\n---\n"), "divider should be absent");
    }

    #[test]
    fn test_md_writer_figure() {
        let mut r = SimulationReport::new("T");
        r.add_figure(FigureRef::new("f", "plot.png", "Energy").with_alt("alt"));
        let md = MarkdownReportWriter::new().write(&r);
        assert!(md.contains("![alt](plot.png)"));
    }

    #[test]
    fn test_md_writer_data_table() {
        let mut r = SimulationReport::new("T");
        let mut tbl = DataTable::new("T1", vec!["X".into(), "Y".into()]);
        tbl.add_row(vec![
            CellValue::Float(1.0, Some(1)),
            CellValue::Float(2.0, Some(1)),
        ]);
        r.add_table(tbl);
        let md = MarkdownReportWriter::new().write(&r);
        assert!(md.contains("| X |"), "header missing");
        assert!(md.contains("| 1.0 |"), "data missing");
    }

    #[test]
    fn test_md_writer_author_and_date() {
        let r = SimulationReport::new("T")
            .with_author("Bob")
            .with_date("2026-03-24");
        let md = MarkdownReportWriter::new().write(&r);
        assert!(md.contains("**Author:** Bob"));
        assert!(md.contains("**Date:** 2026-03-24"));
    }

    // ── LatexReportWriter ─────────────────────────────────────────────────────

    #[test]
    fn test_latex_writer_contains_documentclass() {
        let r = SimulationReport::new("T");
        let tex = LatexReportWriter::new().write(&r);
        assert!(tex.contains("\\documentclass"), "no documentclass");
    }

    #[test]
    fn test_latex_writer_begin_document() {
        let r = SimulationReport::new("T");
        let tex = LatexReportWriter::new().write(&r);
        assert!(tex.contains("\\begin{document}"));
        assert!(tex.contains("\\end{document}"));
    }

    #[test]
    fn test_latex_writer_title() {
        let r = SimulationReport::new("Simulation Report");
        let tex = LatexReportWriter::new().write(&r);
        assert!(tex.contains("\\title{Simulation Report}"));
    }

    #[test]
    fn test_latex_writer_parameters_section() {
        let mut r = SimulationReport::new("T");
        r.add_parameter("dt", "0.001");
        let tex = LatexReportWriter::new().write(&r);
        assert!(tex.contains("\\section{Parameters}"));
        assert!(tex.contains("dt"));
    }

    #[test]
    fn test_latex_writer_results_section() {
        let mut r = SimulationReport::new("T");
        r.add_result("E_max", "42.0");
        let tex = LatexReportWriter::new().write(&r);
        assert!(tex.contains("\\section{Results Summary}"));
        assert!(tex.contains("E"), "result key missing");
    }

    #[test]
    fn test_latex_writer_section() {
        let mut r = SimulationReport::new("T");
        r.add_section(ReportSection::new("Methodology"));
        let tex = LatexReportWriter::new().write(&r);
        assert!(tex.contains("\\section{Methodology}"));
    }

    #[test]
    fn test_latex_writer_subsection() {
        let mut r = SimulationReport::new("T");
        r.add_section(ReportSection::subsection("Numerical Setup"));
        let tex = LatexReportWriter::new().write(&r);
        assert!(tex.contains("\\subsection{Numerical Setup}"));
    }

    #[test]
    fn test_latex_writer_figure() {
        let mut r = SimulationReport::new("T");
        r.add_figure(FigureRef::new("fig1", "energy.pdf", "Energy Plot"));
        let tex = LatexReportWriter::new().write(&r);
        assert!(tex.contains("\\begin{figure}"));
        assert!(tex.contains("energy.pdf"));
        assert!(tex.contains("Energy Plot"));
    }

    #[test]
    fn test_latex_writer_data_table() {
        let mut r = SimulationReport::new("T");
        let mut tbl = DataTable::new("Results", vec!["A".into(), "B".into()]);
        tbl.add_row(vec![CellValue::Int(1), CellValue::Float(2.5, Some(1))]);
        r.add_table(tbl);
        let tex = LatexReportWriter::new().write(&r);
        assert!(tex.contains("\\begin{tabular}"), "tabular missing");
        assert!(tex.contains("2.5"));
    }

    #[test]
    fn test_latex_escape_special_chars() {
        assert_eq!(latex_escape("x_1 & y"), r"x\_1 \& y");
        assert_eq!(latex_escape("100%"), r"100\%");
        assert_eq!(latex_escape("$10"), r"\$10");
    }

    #[test]
    fn test_latex_writer_version() {
        let r = SimulationReport::new("T").with_version("2.0");
        let tex = LatexReportWriter::new().write(&r);
        assert!(tex.contains("Version: 2.0"));
    }

    // ── ReportBuilder ─────────────────────────────────────────────────────────

    #[test]
    fn test_report_builder_to_html() {
        let html = ReportBuilder::new("Builder Test")
            .author("Alice")
            .description("A test report")
            .parameter("n", "100")
            .result("max_err", "1e-5")
            .to_html();
        assert!(html.contains("Builder Test"));
        assert!(html.contains("Alice"));
        assert!(html.contains("max_err"));
    }

    #[test]
    fn test_report_builder_to_markdown() {
        let md = ReportBuilder::new("MD Report")
            .description("Desc")
            .parameter("dt", "0.01")
            .to_markdown();
        assert!(md.starts_with("# MD Report"));
        assert!(md.contains("dt"));
    }

    #[test]
    fn test_report_builder_to_latex() {
        let tex = ReportBuilder::new("LaTeX Report")
            .description("Abstract text")
            .to_latex();
        assert!(tex.contains("\\documentclass"));
        assert!(tex.contains("Abstract text"));
    }

    #[test]
    fn test_report_builder_section_in_output() {
        let mut sec = ReportSection::new("Methods");
        sec.append_text("We used finite differences.");
        let html = ReportBuilder::new("Report").section(sec).to_html();
        assert!(html.contains("Methods"));
        assert!(html.contains("finite differences"));
    }

    // ── html_escape helper ────────────────────────────────────────────────────

    #[test]
    fn test_html_escape_ampersand() {
        assert_eq!(html_escape("a & b"), "a &amp; b");
    }

    #[test]
    fn test_html_escape_angle_brackets() {
        assert_eq!(html_escape("<tag>"), "&lt;tag&gt;");
    }

    #[test]
    fn test_html_escape_quote() {
        assert_eq!(html_escape("\"quoted\""), "&quot;quoted&quot;");
    }

    #[test]
    fn test_html_escape_plain_text() {
        assert_eq!(html_escape("hello world"), "hello world");
    }

    // ── markdown_table helper ─────────────────────────────────────────────────

    #[test]
    fn test_markdown_table_empty_headers() {
        let t = DataTable::new("Cap", vec![]);
        let md = markdown_table(&t);
        // Should have at least the caption
        assert!(md.contains("Cap") || md.is_empty() || md.contains("**Cap**"));
    }

    #[test]
    fn test_markdown_table_with_data() {
        let mut t = DataTable::new("T", vec!["A".into(), "B".into()]);
        t.add_row(vec![CellValue::Text("x".into()), CellValue::Int(99)]);
        let md = markdown_table(&t);
        assert!(md.contains("| A |"), "header A missing");
        assert!(md.contains("| x |"), "cell x missing");
        assert!(md.contains("99"));
    }

    #[test]
    fn test_html_table_with_highlight() {
        let mut t = DataTable::new("T", vec!["A".into()]);
        t.add_row(vec![CellValue::Int(1)]);
        t.highlight_row(0, "highlight");
        let html = html_table(&t);
        assert!(
            html.contains("class=\"highlight\""),
            "highlight class missing: {html}"
        );
    }

    // ── ColumnType ────────────────────────────────────────────────────────────

    #[test]
    fn test_column_type_labels() {
        assert_eq!(ColumnType::Integer.label(), "integer");
        assert_eq!(ColumnType::Float.label(), "float");
        assert_eq!(ColumnType::Text.label(), "text");
        assert_eq!(ColumnType::Boolean.label(), "boolean");
    }
}
