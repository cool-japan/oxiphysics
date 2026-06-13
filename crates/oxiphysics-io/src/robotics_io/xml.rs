//! Minimal Pure-Rust XML parser for the URDF subset (no external XML crate).
//!
//! Supports: the `<?xml ?>` prolog, `<!-- -->` comments, self-closing tags,
//! open/close tags, both `"` and `'` attribute quoting, the five predefined
//! entity references, and nested elements. Whitespace-only text content is
//! discarded; non-whitespace text content is preserved on the element.

/// A parsed XML element node holding its name, attributes, children and text.
#[derive(Debug, Clone)]
pub struct XmlElement {
    /// The tag name of this element.
    pub name: String,
    /// The attributes of this element as `(key, value)` pairs, in document order.
    pub attrs: Vec<(String, String)>,
    /// The direct child elements of this element, in document order.
    pub children: Vec<XmlElement>,
    /// The concatenated, entity-decoded text content of this element.
    pub text: String,
}

impl XmlElement {
    /// Returns the value of the first attribute matching `name`, if any.
    pub fn attr(&self, name: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|(k, _)| k == name)
            .map(|(_, v)| v.as_str())
    }

    /// Returns an iterator over direct children whose tag name equals `name`.
    pub fn children_named<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a XmlElement> + 'a {
        self.children.iter().filter(move |c| c.name == name)
    }

    /// Returns the first direct child whose tag name equals `name`, if any.
    pub fn child(&self, name: &str) -> Option<&XmlElement> {
        self.children.iter().find(|c| c.name == name)
    }
}

/// Returns true if the slice starting at `i` matches `pat` character-by-character.
fn starts_with(chars: &[char], i: usize, pat: &str) -> bool {
    for (offset, pc) in pat.chars().enumerate() {
        match chars.get(i + offset) {
            Some(&c) if c == pc => {}
            _ => return false,
        }
    }
    true
}

/// Decodes the five predefined XML entities in a single left-to-right scan.
///
/// Unrecognised `&` sequences are passed through literally; this avoids the
/// double-decoding that chained `str::replace` calls would introduce.
fn decode_entities(s: &str) -> String {
    let chars: Vec<char> = s.chars().collect();
    let mut out = String::with_capacity(chars.len());
    let mut i = 0usize;
    while i < chars.len() {
        if chars[i] == '&' {
            if starts_with(&chars, i, "&lt;") {
                out.push('<');
                i += 4;
            } else if starts_with(&chars, i, "&gt;") {
                out.push('>');
                i += 4;
            } else if starts_with(&chars, i, "&amp;") {
                out.push('&');
                i += 5;
            } else if starts_with(&chars, i, "&quot;") {
                out.push('"');
                i += 6;
            } else if starts_with(&chars, i, "&apos;") {
                out.push('\'');
                i += 6;
            } else {
                out.push('&');
                i += 1;
            }
        } else {
            out.push(chars[i]);
            i += 1;
        }
    }
    out
}

/// Parses an XML document into a single root [`XmlElement`].
///
/// Returns `Err` with a human-readable message on any malformed input.
pub fn parse_xml(input: &str) -> Result<XmlElement, String> {
    let chars: Vec<char> = input.chars().collect();
    let mut i = 0usize;
    let mut stack: Vec<XmlElement> = Vec::new();
    let mut root: Option<XmlElement> = None;

    while i < chars.len() {
        // 1. Skip leading whitespace.
        while i < chars.len() && chars[i].is_whitespace() {
            i += 1;
        }
        // 2. Stop if we consumed the rest of the input.
        if i >= chars.len() {
            break;
        }

        // 3. Comment.
        if starts_with(&chars, i, "<!--") {
            let mut j = i + 4;
            let mut found = false;
            while j + 2 < chars.len() || j + 3 <= chars.len() {
                if starts_with(&chars, j, "-->") {
                    found = true;
                    break;
                }
                j += 1;
            }
            if !found {
                return Err("unterminated comment".into());
            }
            i = j + 3;
            continue;
        }

        // 4. Processing instruction.
        if starts_with(&chars, i, "<?") {
            let mut j = i + 2;
            let mut found = false;
            while j < chars.len() {
                if starts_with(&chars, j, "?>") {
                    found = true;
                    break;
                }
                j += 1;
            }
            if !found {
                return Err("unterminated processing instruction".into());
            }
            i = j + 2;
            continue;
        }

        // 5. Declaration (e.g. DOCTYPE).
        if starts_with(&chars, i, "<!") {
            let mut j = i + 2;
            let mut found = false;
            while j < chars.len() {
                if chars[j] == '>' {
                    found = true;
                    break;
                }
                j += 1;
            }
            if !found {
                return Err("unterminated declaration".into());
            }
            i = j + 1;
            continue;
        }

        // 6. Closing tag.
        if starts_with(&chars, i, "</") {
            i += 2;
            let start = i;
            while i < chars.len() && !chars[i].is_whitespace() && chars[i] != '>' {
                i += 1;
            }
            let name: String = chars[start..i].iter().collect();
            while i < chars.len() && chars[i].is_whitespace() {
                i += 1;
            }
            if i >= chars.len() || chars[i] != '>' {
                return Err("malformed closing tag".into());
            }
            i += 1;
            let elem = stack
                .pop()
                .ok_or_else(|| format!("unexpected closing tag </{}>", name))?;
            if elem.name != name {
                return Err(format!(
                    "mismatched closing tag: expected </{}>, found </{}>",
                    elem.name, name
                ));
            }
            if stack.is_empty() {
                if root.is_some() {
                    return Err("multiple root elements".into());
                }
                root = Some(elem);
            } else if let Some(parent) = stack.last_mut() {
                parent.children.push(elem);
            }
            continue;
        }

        // 7. Opening or self-closing tag.
        if chars[i] == '<' {
            i += 1;
            let start = i;
            while i < chars.len()
                && !chars[i].is_whitespace()
                && chars[i] != '/'
                && chars[i] != '>'
            {
                i += 1;
            }
            let name: String = chars[start..i].iter().collect();
            if name.is_empty() {
                return Err("empty tag name".into());
            }

            let mut attrs: Vec<(String, String)> = Vec::new();
            loop {
                while i < chars.len() && chars[i].is_whitespace() {
                    i += 1;
                }
                if i >= chars.len() {
                    return Err("unexpected end in tag".into());
                }
                if chars[i] == '/' || chars[i] == '>' {
                    break;
                }
                let astart = i;
                while i < chars.len()
                    && chars[i] != '='
                    && !chars[i].is_whitespace()
                    && chars[i] != '/'
                    && chars[i] != '>'
                {
                    i += 1;
                }
                let aname: String = chars[astart..i].iter().collect();
                if aname.is_empty() {
                    return Err("empty attribute name".into());
                }
                while i < chars.len() && chars[i].is_whitespace() {
                    i += 1;
                }
                if i >= chars.len() || chars[i] != '=' {
                    return Err("expected '=' in attribute".into());
                }
                i += 1;
                while i < chars.len() && chars[i].is_whitespace() {
                    i += 1;
                }
                if i >= chars.len() || (chars[i] != '"' && chars[i] != '\'') {
                    return Err("expected quote".into());
                }
                let quote = chars[i];
                i += 1;
                let vstart = i;
                while i < chars.len() && chars[i] != quote {
                    i += 1;
                }
                if i >= chars.len() {
                    return Err("unterminated attribute value".into());
                }
                let value: String = chars[vstart..i].iter().collect();
                i += 1;
                attrs.push((aname, decode_entities(&value)));
            }

            while i < chars.len() && chars[i].is_whitespace() {
                i += 1;
            }
            if i + 1 < chars.len() && chars[i] == '/' && chars[i + 1] == '>' {
                let elem = XmlElement {
                    name,
                    attrs,
                    children: Vec::new(),
                    text: String::new(),
                };
                i += 2;
                if stack.is_empty() {
                    if root.is_some() {
                        return Err("multiple root elements".into());
                    }
                    root = Some(elem);
                } else if let Some(parent) = stack.last_mut() {
                    parent.children.push(elem);
                }
                continue;
            } else if i < chars.len() && chars[i] == '>' {
                i += 1;
                stack.push(XmlElement {
                    name,
                    attrs,
                    children: Vec::new(),
                    text: String::new(),
                });
                continue;
            } else {
                return Err("malformed tag".into());
            }
        }

        // 8. Text content.
        let start = i;
        while i < chars.len() && chars[i] != '<' {
            i += 1;
        }
        let raw: String = chars[start..i].iter().collect();
        let decoded = decode_entities(&raw);
        if decoded.trim().is_empty() {
            continue;
        }
        let top = stack
            .last_mut()
            .ok_or_else(|| "text content outside root element".to_string())?;
        top.text.push_str(&decoded);
    }

    if !stack.is_empty() {
        return Err(format!(
            "unclosed tag <{}>",
            stack.last().map(|e| e.name.as_str()).unwrap_or("?")
        ));
    }
    root.ok_or_else(|| "no root element".to_string())
}

/// Parses a single `f64` from a string, trimming surrounding whitespace.
pub fn parse_float_str(s: &str) -> Result<f64, String> {
    s.trim()
        .parse::<f64>()
        .map_err(|e| format!("invalid float '{}': {}", s.trim(), e))
}

/// Parses exactly three whitespace-separated `f64` values into a `[f64; 3]`.
pub fn parse_vec3_str(s: &str) -> Result<[f64; 3], String> {
    let parts: Vec<&str> = s.split_whitespace().collect();
    if parts.len() != 3 {
        return Err(format!(
            "expected 3 floats, got {}: '{}'",
            parts.len(),
            s.trim()
        ));
    }
    let a = parse_float_str(parts[0])?;
    let b = parse_float_str(parts[1])?;
    let c = parse_float_str(parts[2])?;
    Ok([a, b, c])
}

/// Parses a single `f64` from the named attribute of `elem`.
pub fn parse_float_attr(elem: &XmlElement, attr: &str) -> Result<f64, String> {
    let v = elem
        .attr(attr)
        .ok_or_else(|| format!("missing attribute '{}' on <{}>", attr, elem.name))?;
    parse_float_str(v)
}

/// Parses three `f64` values from the named attribute of `elem`.
pub fn parse_vec3_attr(elem: &XmlElement, attr: &str) -> Result<[f64; 3], String> {
    let v = elem
        .attr(attr)
        .ok_or_else(|| format!("missing attribute '{}' on <{}>", attr, elem.name))?;
    parse_vec3_str(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn self_closing_tag() {
        let root = parse_xml("<a x=\"1\"/>").unwrap();
        assert_eq!(root.name, "a");
        assert_eq!(root.attr("x"), Some("1"));
    }

    #[test]
    fn nested_tags_with_text() {
        let root = parse_xml("<a><b>hi</b></a>").unwrap();
        assert_eq!(root.name, "a");
        let b = root.child("b").expect("child b");
        assert_eq!(b.text.trim(), "hi");
    }

    #[test]
    fn decode_amp_and_lt() {
        let root = parse_xml("<a v=\"x&amp;y\"/>").unwrap();
        assert_eq!(root.attr("v"), Some("x&y"));
        let root2 = parse_xml("<a v=\"x&lt;y\"/>").unwrap();
        assert_eq!(root2.attr("v"), Some("x<y"));
    }

    #[test]
    fn both_quote_styles() {
        let root = parse_xml("<a x='1' y=\"2\"/>").unwrap();
        assert_eq!(root.attr("x"), Some("1"));
        assert_eq!(root.attr("y"), Some("2"));
    }

    #[test]
    fn comment_skip() {
        let root = parse_xml("<!-- c --><a/>").unwrap();
        assert_eq!(root.name, "a");
    }

    #[test]
    fn unterminated_comment_errors() {
        assert!(parse_xml("<!-- unterminated").is_err());
    }
}
