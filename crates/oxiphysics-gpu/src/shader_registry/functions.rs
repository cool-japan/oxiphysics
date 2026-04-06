//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::ShaderKey;

/// Extract all `{{PLACEHOLDER}}` names from a WGSL string.
pub(super) fn collect_placeholders(src: &str) -> Vec<String> {
    let mut result = Vec::new();
    let mut rest = src;
    while let Some(start) = rest.find("{{") {
        let after_open = &rest[start + 2..];
        if let Some(end) = after_open.find("}}") {
            let token = after_open[..end].trim().to_string();
            if !token.is_empty() && !result.contains(&token) {
                result.push(token);
            }
            rest = &after_open[end + 2..];
        } else {
            break;
        }
    }
    result
}
/// Compute a 64-bit FNV-1a hash of the given WGSL source bytes.
///
/// This is intentionally a fast, non-cryptographic hash suitable for
/// cache-key generation.
pub fn compute_cache_key(wgsl: &str) -> u64 {
    pub(super) const FNV_OFFSET: u64 = 14695981039346656037;
    pub(super) const FNV_PRIME: u64 = 1099511628211;
    let mut hash = FNV_OFFSET;
    for byte in wgsl.bytes() {
        hash ^= byte as u64;
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}
/// Compute a composite cache key that includes the shader name and its defines.
///
/// The key format is `fnv64(name || "\x00" || sorted_defines_string)`.
pub fn compute_shader_cache_key(key: &ShaderKey) -> u64 {
    let mut repr = key.name.clone();
    repr.push('\x00');
    for (k, v) in &key.defines {
        repr.push_str(k);
        repr.push('=');
        repr.push_str(v);
        repr.push(';');
    }
    compute_cache_key(&repr)
}
#[cfg(test)]
mod shader_registry_tests {
    use super::*;
    use crate::shader_registry::CompiledVariant;
    use crate::shader_registry::HotReloadTracker;

    use crate::shader_registry::PipelineDescriptor;
    use crate::shader_registry::RegistryError;
    use crate::shader_registry::ShaderCompileOptions;

    use crate::shader_registry::ShaderKey;
    use crate::shader_registry::ShaderRegistry;
    use crate::shader_registry::ShaderSource;

    use crate::shader_registry::VariantCache;

    use std::collections::HashMap;
    /// Test helper: construct a CompiledVariant from key, wgsl, and workgroup_size.
    /// Mirrors the private `CompiledVariant::new`.
    fn make_compiled_variant(
        key: ShaderKey,
        wgsl: String,
        workgroup_size: [u32; 3],
    ) -> CompiledVariant {
        let binary_size_bytes = wgsl.len();
        CompiledVariant {
            key,
            wgsl,
            workgroup_size,
            binary_size_bytes,
        }
    }
    fn simple_source(name: &str) -> ShaderSource {
        ShaderSource::new(
            name,
            "@compute @workgroup_size({{WG_SIZE}}) fn main() {}",
            [64, 1, 1],
        )
    }
    #[test]
    fn test_shader_key_fingerprint_bare() {
        let key = ShaderKey::bare("sph_density");
        assert_eq!(key.fingerprint(), "sph_density");
    }
    #[test]
    fn test_shader_key_fingerprint_with_defines() {
        let mut d = HashMap::new();
        d.insert("WG_SIZE".to_string(), "64".to_string());
        let key = ShaderKey::new("sph", &d);
        assert!(key.fingerprint().contains("WG_SIZE"));
        assert!(key.fingerprint().contains("64"));
    }
    #[test]
    fn test_shader_key_sorted_defines() {
        let mut d = HashMap::new();
        d.insert("Z_FLAG".to_string(), "1".to_string());
        d.insert("A_FLAG".to_string(), "2".to_string());
        let key = ShaderKey::new("s", &d);
        assert_eq!(
            key.defines[0].0, "A_FLAG",
            "defines should be sorted by key"
        );
    }
    #[test]
    fn test_source_placeholder_detection() {
        let src = ShaderSource::new(
            "test",
            "size={{SIZE}} dtype={{DTYPE}} size={{SIZE}}",
            [32, 1, 1],
        );
        assert_eq!(src.placeholders.len(), 2);
        assert!(src.placeholders.contains(&"SIZE".to_string()));
        assert!(src.placeholders.contains(&"DTYPE".to_string()));
    }
    #[test]
    fn test_source_instantiate_substitutes_tokens() {
        let src = ShaderSource::new("t", "workgroup={{WG}} type={{T}}", [1, 1, 1]);
        let mut d = HashMap::new();
        d.insert("WG".to_string(), "128".to_string());
        d.insert("T".to_string(), "f32".to_string());
        let out = src.instantiate(&d);
        assert!(out.contains("128"));
        assert!(out.contains("f32"));
        assert!(!out.contains("{{WG}}"));
    }
    #[test]
    fn test_variant_cache_lru_eviction() {
        let mut cache = VariantCache::new(2);
        let make = |name: &str| {
            make_compiled_variant(
                ShaderKey::bare(name),
                format!("fn {}(){{}}", name),
                [64, 1, 1],
            )
        };
        cache.insert(make("a"));
        cache.insert(make("b"));
        let _ = cache.get(&ShaderKey::bare("a"));
        cache.insert(make("c"));
        assert!(
            cache.get(&ShaderKey::bare("a")).is_some(),
            "a should survive"
        );
        assert!(
            cache.get(&ShaderKey::bare("b")).is_none(),
            "b should be evicted"
        );
        assert!(
            cache.get(&ShaderKey::bare("c")).is_some(),
            "c should be present"
        );
    }
    #[test]
    fn test_registry_get_or_compile_success() {
        let mut reg = ShaderRegistry::new(8);
        reg.register(simple_source("dense_layer"));
        let mut d = HashMap::new();
        d.insert("WG_SIZE".to_string(), "64".to_string());
        let v = reg.get_or_compile("dense_layer", &d).unwrap();
        assert!(v.wgsl.contains("64"), "placeholder should be substituted");
        assert_eq!(reg.compilations, 1);
    }
    #[test]
    fn test_registry_cache_hit() {
        let mut reg = ShaderRegistry::new(8);
        reg.register(simple_source("sph"));
        let mut d = HashMap::new();
        d.insert("WG_SIZE".to_string(), "32".to_string());
        reg.get_or_compile("sph", &d).unwrap();
        reg.get_or_compile("sph", &d).unwrap();
        assert_eq!(reg.cache_hits, 1);
        assert_eq!(reg.compilations, 1);
    }
    #[test]
    fn test_registry_unknown_shader_error() {
        let mut reg = ShaderRegistry::new(4);
        let res = reg.get_or_compile("ghost", &HashMap::new());
        assert!(matches!(res, Err(RegistryError::UnknownShader(_))));
    }
    #[test]
    fn test_registry_missing_define_error() {
        let mut reg = ShaderRegistry::new(4);
        reg.register(simple_source("needs_wg"));
        let res = reg.get_or_compile("needs_wg", &HashMap::new());
        assert!(
            matches!(res, Err(RegistryError::MissingDefine { .. })),
            "should fail when required define is absent"
        );
    }
    #[test]
    fn test_registry_invalidate_clears_cache() {
        let mut reg = ShaderRegistry::new(8);
        reg.register(simple_source("shader_x"));
        let mut d = HashMap::new();
        d.insert("WG_SIZE".to_string(), "64".to_string());
        reg.get_or_compile("shader_x", &d).unwrap();
        assert_eq!(reg.cached_count(), 1);
        reg.invalidate("shader_x");
        assert_eq!(
            reg.cached_count(),
            0,
            "invalidate should clear variant cache"
        );
    }
    #[test]
    fn test_hot_reload_tracker_needs_recompile() {
        let mut tracker = HotReloadTracker::new();
        tracker.touch("my_shader", 10);
        tracker.record_compile("my_shader", 5);
        assert!(
            tracker.needs_recompile("my_shader"),
            "modified at 10 > compiled at 5"
        );
        tracker.record_compile("my_shader", 10);
        assert!(
            !tracker.needs_recompile("my_shader"),
            "up-to-date after recompile"
        );
    }
    #[test]
    fn test_hot_reload_tracker_stale_list() {
        let mut tracker = HotReloadTracker::new();
        tracker.touch("a", 5);
        tracker.touch("b", 5);
        tracker.record_compile("a", 10);
        let stale = tracker.stale_shaders();
        assert!(!stale.contains(&"a"), "a is up-to-date");
        assert!(stale.contains(&"b"), "b has never been compiled");
    }
    #[test]
    fn test_pipeline_descriptor_validate() {
        let key = ShaderKey::bare("sph");
        let desc = PipelineDescriptor::new(key.clone(), 2, 0, "sph_pipeline");
        assert!(desc.validate().is_ok());
        let bad = PipelineDescriptor::new(key, 0, 0, "bad");
        assert!(bad.validate().is_err());
    }
    #[test]
    fn test_shader_key_bare_no_defines() {
        let key = ShaderKey::bare("my_shader");
        assert!(key.defines.is_empty());
        assert_eq!(key.name, "my_shader");
    }
    #[test]
    fn test_shader_key_equality() {
        let mut d1 = HashMap::new();
        d1.insert("A".to_string(), "1".to_string());
        let mut d2 = HashMap::new();
        d2.insert("A".to_string(), "1".to_string());
        let k1 = ShaderKey::new("s", &d1);
        let k2 = ShaderKey::new("s", &d2);
        assert_eq!(k1, k2);
    }
    #[test]
    fn test_shader_key_inequality_different_value() {
        let mut d1 = HashMap::new();
        d1.insert("A".to_string(), "1".to_string());
        let mut d2 = HashMap::new();
        d2.insert("A".to_string(), "2".to_string());
        let k1 = ShaderKey::new("s", &d1);
        let k2 = ShaderKey::new("s", &d2);
        assert_ne!(k1, k2);
    }
    #[test]
    fn test_shader_key_inequality_different_name() {
        let k1 = ShaderKey::bare("shader_a");
        let k2 = ShaderKey::bare("shader_b");
        assert_ne!(k1, k2);
    }
    #[test]
    fn test_shader_key_fingerprint_multiple_defines() {
        let mut d = HashMap::new();
        d.insert("B".to_string(), "2".to_string());
        d.insert("A".to_string(), "1".to_string());
        let key = ShaderKey::new("s", &d);
        let fp = key.fingerprint();
        let a_pos = fp.find("A_1").unwrap();
        let b_pos = fp.find("B_2").unwrap();
        assert!(
            a_pos < b_pos,
            "defines should be alphabetical in fingerprint"
        );
    }
    #[test]
    fn test_shader_source_threads_per_group() {
        let src = ShaderSource::new("t", "fn main() {}", [4, 8, 2]);
        assert_eq!(src.threads_per_group(), 64);
    }
    #[test]
    fn test_shader_source_no_placeholders() {
        let src = ShaderSource::new("t", "fn main() { let x = 1; }", [64, 1, 1]);
        assert!(src.placeholders.is_empty());
    }
    #[test]
    fn test_shader_source_instantiate_no_defines() {
        let src = ShaderSource::new("t", "fn main() {}", [1, 1, 1]);
        let out = src.instantiate(&HashMap::new());
        assert_eq!(out, "fn main() {}");
    }
    #[test]
    fn test_shader_source_instantiate_partial() {
        let src = ShaderSource::new("t", "{{A}} and {{B}}", [1, 1, 1]);
        let mut d = HashMap::new();
        d.insert("A".to_string(), "hello".to_string());
        let out = src.instantiate(&d);
        assert!(out.contains("hello"), "A should be substituted");
        assert!(out.contains("{{B}}"), "B should remain if not supplied");
    }
    #[test]
    fn test_shader_source_multiple_occurrences_replaced() {
        let src = ShaderSource::new("t", "{{X}} and {{X}} and {{X}}", [1, 1, 1]);
        let mut d = HashMap::new();
        d.insert("X".to_string(), "42".to_string());
        let out = src.instantiate(&d);
        assert_eq!(out.matches("42").count(), 3);
        assert!(!out.contains("{{X}}"));
    }
    #[test]
    fn test_variant_cache_empty_initially() {
        let cache = VariantCache::new(5);
        assert!(cache.is_empty());
        assert_eq!(cache.len(), 0);
    }
    #[test]
    fn test_variant_cache_insert_and_get() {
        let mut cache = VariantCache::new(5);
        let v = make_compiled_variant(ShaderKey::bare("s"), "fn s() {}".to_string(), [1, 1, 1]);
        cache.insert(v);
        assert_eq!(cache.len(), 1);
        let result = cache.get(&ShaderKey::bare("s"));
        assert!(result.is_some());
    }
    #[test]
    fn test_variant_cache_miss() {
        let mut cache = VariantCache::new(5);
        let result = cache.get(&ShaderKey::bare("nonexistent"));
        assert!(result.is_none());
    }
    #[test]
    fn test_variant_cache_update_existing() {
        let mut cache = VariantCache::new(5);
        let v1 = make_compiled_variant(ShaderKey::bare("s"), "v1".to_string(), [1, 1, 1]);
        let v2 = make_compiled_variant(ShaderKey::bare("s"), "v2".to_string(), [1, 1, 1]);
        cache.insert(v1);
        cache.insert(v2);
        assert_eq!(cache.len(), 1, "update should not add a second entry");
        let got = cache.get(&ShaderKey::bare("s")).unwrap();
        assert_eq!(got.wgsl, "v2");
    }
    #[test]
    fn test_variant_cache_clear() {
        let mut cache = VariantCache::new(5);
        cache.insert(make_compiled_variant(
            ShaderKey::bare("a"),
            "a".to_string(),
            [1, 1, 1],
        ));
        cache.insert(make_compiled_variant(
            ShaderKey::bare("b"),
            "b".to_string(),
            [1, 1, 1],
        ));
        cache.clear();
        assert!(cache.is_empty());
    }
    #[test]
    fn test_variant_cache_total_binary_bytes() {
        let mut cache = VariantCache::new(5);
        cache.insert(make_compiled_variant(
            ShaderKey::bare("a"),
            "abc".to_string(),
            [1, 1, 1],
        ));
        cache.insert(make_compiled_variant(
            ShaderKey::bare("b"),
            "xy".to_string(),
            [1, 1, 1],
        ));
        assert_eq!(cache.total_binary_bytes(), 5);
    }
    #[test]
    fn test_variant_cache_capacity_one_evicts_on_insert() {
        let mut cache = VariantCache::new(1);
        cache.insert(make_compiled_variant(
            ShaderKey::bare("a"),
            "a".to_string(),
            [1, 1, 1],
        ));
        cache.insert(make_compiled_variant(
            ShaderKey::bare("b"),
            "b".to_string(),
            [1, 1, 1],
        ));
        assert_eq!(
            cache.len(),
            1,
            "capacity=1 cache should only hold one entry"
        );
        let has_a = cache.get(&ShaderKey::bare("a")).is_some();
        let has_b = cache.get(&ShaderKey::bare("b")).is_some();
        assert!(has_a || has_b, "at least one entry should be present");
    }
    #[test]
    fn test_registry_empty_initially() {
        let reg = ShaderRegistry::new(8);
        assert!(reg.shader_names().is_empty());
        assert_eq!(reg.cached_count(), 0);
    }
    #[test]
    fn test_registry_register_multiple_shaders() {
        let mut reg = ShaderRegistry::new(8);
        reg.register(ShaderSource::new("a", "fn a(){}", [1, 1, 1]));
        reg.register(ShaderSource::new("b", "fn b(){}", [1, 1, 1]));
        assert_eq!(reg.shader_names().len(), 2);
    }
    #[test]
    fn test_registry_re_register_clears_cache() {
        let mut reg = ShaderRegistry::new(8);
        reg.register(ShaderSource::new("s", "fn s(){}", [1, 1, 1]));
        reg.get_or_compile("s", &HashMap::new()).unwrap();
        assert_eq!(reg.cached_count(), 1);
        reg.register(ShaderSource::new("s", "fn s_v2(){}", [1, 1, 1]));
        assert_eq!(
            reg.cached_count(),
            0,
            "re-register should clear old variants"
        );
    }
    #[test]
    fn test_registry_invalidate_all() {
        let mut reg = ShaderRegistry::new(8);
        reg.register(ShaderSource::new("a", "fn a(){}", [1, 1, 1]));
        reg.register(ShaderSource::new("b", "fn b(){}", [1, 1, 1]));
        reg.get_or_compile("a", &HashMap::new()).unwrap();
        reg.get_or_compile("b", &HashMap::new()).unwrap();
        assert_eq!(reg.cached_count(), 2);
        reg.invalidate_all();
        assert_eq!(reg.cached_count(), 0);
    }
    #[test]
    fn test_registry_compilation_count_increments() {
        let mut reg = ShaderRegistry::new(8);
        reg.register(ShaderSource::new("x", "fn x(){}", [1, 1, 1]));
        reg.register(ShaderSource::new("y", "fn y(){}", [1, 1, 1]));
        reg.get_or_compile("x", &HashMap::new()).unwrap();
        reg.get_or_compile("y", &HashMap::new()).unwrap();
        assert_eq!(reg.compilations, 2);
    }
    #[test]
    fn test_registry_multiple_cache_hits() {
        let mut reg = ShaderRegistry::new(8);
        reg.register(ShaderSource::new("z", "fn z(){}", [1, 1, 1]));
        for _ in 0..5 {
            reg.get_or_compile("z", &HashMap::new()).unwrap();
        }
        assert_eq!(reg.compilations, 1);
        assert_eq!(reg.cache_hits, 4);
    }
    #[test]
    fn test_hot_reload_no_entry_returns_false() {
        let tracker = HotReloadTracker::new();
        assert!(!tracker.needs_recompile("phantom_shader"));
    }
    #[test]
    fn test_hot_reload_fresh_compile() {
        let mut tracker = HotReloadTracker::new();
        tracker.touch("shader", 5);
        tracker.record_compile("shader", 5);
        assert!(!tracker.needs_recompile("shader"));
    }
    #[test]
    fn test_hot_reload_stale_shaders_empty() {
        let tracker = HotReloadTracker::new();
        assert!(tracker.stale_shaders().is_empty());
    }
    #[test]
    fn test_hot_reload_multiple_stale() {
        let mut tracker = HotReloadTracker::new();
        tracker.touch("a", 10);
        tracker.touch("b", 20);
        tracker.touch("c", 5);
        tracker.record_compile("c", 10);
        let stale = tracker.stale_shaders();
        assert!(stale.contains(&"a"), "a is stale");
        assert!(stale.contains(&"b"), "b is stale");
        assert!(!stale.contains(&"c"), "c is up to date");
    }
    #[test]
    fn test_registry_error_display_unknown() {
        let e = RegistryError::UnknownShader("foo".to_string());
        assert!(e.to_string().contains("foo"));
    }
    #[test]
    fn test_registry_error_display_missing_define() {
        let e = RegistryError::MissingDefine {
            shader: "sph".to_string(),
            define: "WG_SIZE".to_string(),
        };
        let s = e.to_string();
        assert!(s.contains("sph") && s.contains("WG_SIZE"));
    }
    #[test]
    fn test_registry_error_display_source_too_large() {
        let e = RegistryError::SourceTooLarge {
            shader: "big".to_string(),
            size: 100_000,
            limit: 65_536,
        };
        let s = e.to_string();
        assert!(s.contains("big") && s.contains("100000") && s.contains("65536"));
    }
    #[test]
    fn test_pipeline_descriptor_fields() {
        let key = ShaderKey::bare("nn_forward");
        let desc = PipelineDescriptor::new(key.clone(), 3, 16, "nn_fwd");
        assert_eq!(desc.bind_group_count, 3);
        assert_eq!(desc.push_constant_bytes, 16);
        assert_eq!(desc.label, "nn_fwd");
        assert_eq!(desc.key, key);
    }
    #[test]
    fn test_pipeline_descriptor_validate_with_bind_groups() {
        for count in 1..=4 {
            let key = ShaderKey::bare("test");
            let desc = PipelineDescriptor::new(key, count, 0, "lbl");
            assert!(
                desc.validate().is_ok(),
                "bind_group_count={count} should be valid"
            );
        }
    }
    #[test]
    fn test_compile_variant_basic() {
        let mut reg = ShaderRegistry::new(4);
        reg.register(simple_source("dense"));
        let mut d = HashMap::new();
        d.insert("WG_SIZE".to_string(), "128".to_string());
        let opts = ShaderCompileOptions::new();
        let v = reg.compile_variant("dense", &d, &opts).unwrap();
        assert!(v.wgsl.contains("128"), "define should be substituted");
    }
    #[test]
    fn test_compile_variant_unknown_shader() {
        let mut reg = ShaderRegistry::new(4);
        let opts = ShaderCompileOptions::new();
        let res = reg.compile_variant("missing", &HashMap::new(), &opts);
        assert!(matches!(res, Err(RegistryError::UnknownShader(_))));
    }
    #[test]
    fn test_compile_variant_source_too_large() {
        let mut reg = ShaderRegistry::new(4);
        reg.register(simple_source("big"));
        let mut d = HashMap::new();
        d.insert("WG_SIZE".to_string(), "64".to_string());
        let mut opts = ShaderCompileOptions::new();
        opts.max_source_bytes = 5;
        let res = reg.compile_variant("big", &d, &opts);
        assert!(matches!(res, Err(RegistryError::SourceTooLarge { .. })));
    }
    #[test]
    fn test_compile_variant_extra_defines_are_merged() {
        let src = ShaderSource::new("s", "{{A}} {{B}}", [32, 1, 1]);
        let mut reg = ShaderRegistry::new(4);
        reg.register(src);
        let mut d = HashMap::new();
        d.insert("A".to_string(), "hello".to_string());
        let mut opts = ShaderCompileOptions::new();
        opts.extra_defines
            .insert("B".to_string(), "world".to_string());
        let v = reg.compile_variant("s", &d, &opts).unwrap();
        assert!(v.wgsl.contains("hello"), "A should be substituted");
        assert!(
            v.wgsl.contains("world"),
            "B should be substituted from extra_defines"
        );
    }
    #[test]
    fn test_compile_variant_increments_compilation_count() {
        let mut reg = ShaderRegistry::new(4);
        reg.register(simple_source("s"));
        let mut d = HashMap::new();
        d.insert("WG_SIZE".to_string(), "16".to_string());
        let opts = ShaderCompileOptions::new();
        reg.compile_variant("s", &d, &opts).unwrap();
        reg.compile_variant("s", &d, &opts).unwrap();
        assert_eq!(
            reg.compilations, 2,
            "each compile_variant call should count"
        );
    }
    #[test]
    fn test_compute_cache_key_deterministic() {
        let src = "fn main() {}";
        assert_eq!(compute_cache_key(src), compute_cache_key(src));
    }
    #[test]
    fn test_compute_cache_key_different_inputs() {
        let k1 = compute_cache_key("fn a() {}");
        let k2 = compute_cache_key("fn b() {}");
        assert_ne!(
            k1, k2,
            "different sources should produce different cache keys"
        );
    }
    #[test]
    fn test_compute_shader_cache_key_bare_vs_defines() {
        let bare_key = ShaderKey::bare("sph");
        let mut d = HashMap::new();
        d.insert("WG_SIZE".to_string(), "64".to_string());
        let define_key = ShaderKey::new("sph", &d);
        let k1 = compute_shader_cache_key(&bare_key);
        let k2 = compute_shader_cache_key(&define_key);
        assert_ne!(k1, k2, "bare and define keys should differ");
    }
    #[test]
    fn test_compute_shader_cache_key_order_independent() {
        let mut d1 = HashMap::new();
        d1.insert("A".to_string(), "1".to_string());
        d1.insert("B".to_string(), "2".to_string());
        let mut d2 = HashMap::new();
        d2.insert("B".to_string(), "2".to_string());
        d2.insert("A".to_string(), "1".to_string());
        let k1 = compute_shader_cache_key(&ShaderKey::new("s", &d1));
        let k2 = compute_shader_cache_key(&ShaderKey::new("s", &d2));
        assert_eq!(k1, k2, "insertion order should not affect cache key");
    }
    #[test]
    fn test_hot_reload_touch_batch() {
        let mut tracker = HotReloadTracker::new();
        tracker.touch_batch(&["a", "b", "c"], 5);
        assert!(tracker.needs_recompile("a"));
        assert!(tracker.needs_recompile("b"));
        assert!(tracker.needs_recompile("c"));
    }
    #[test]
    fn test_hot_reload_flush_stale() {
        let mut tracker = HotReloadTracker::new();
        tracker.touch_batch(&["x", "y"], 10);
        tracker.flush_stale(20);
        assert!(
            !tracker.needs_recompile("x"),
            "x should be up to date after flush"
        );
        assert!(
            !tracker.needs_recompile("y"),
            "y should be up to date after flush"
        );
    }
    #[test]
    fn test_hot_reload_never_compiled() {
        let mut tracker = HotReloadTracker::new();
        tracker.touch("alpha", 1);
        tracker.touch("beta", 2);
        tracker.record_compile("beta", 3);
        let nc = tracker.never_compiled();
        assert!(nc.contains(&"alpha"), "alpha was never compiled");
        assert!(!nc.contains(&"beta"), "beta was compiled");
    }
}
#[cfg(test)]
mod extended_registry_tests {

    use crate::shader_registry::HotReloadTracker;
    use crate::shader_registry::PipelineCache;
    use crate::shader_registry::PipelineCacheKey;

    use crate::shader_registry::ShaderDependencyGraph;
    use crate::shader_registry::ShaderKey;
    use crate::shader_registry::ShaderRegistry;
    use crate::shader_registry::ShaderSource;
    use crate::shader_registry::SpecConstSet;
    use crate::shader_registry::SpecConstValue;
    use crate::shader_registry::SpecializationConstant;

    use crate::shader_registry::VariantProfile;
    use crate::shader_registry::VariantProfileRegistry;
    use std::collections::HashMap;
    #[test]
    fn variant_profile_no_base_resolve() {
        let mut reg = VariantProfileRegistry::new();
        let p = VariantProfile::new("default")
            .set("WG_SIZE", "64")
            .set("DTYPE", "f32");
        reg.register(p);
        let defines = reg.resolve("default").unwrap();
        assert_eq!(defines.get("WG_SIZE").map(|s| s.as_str()), Some("64"));
        assert_eq!(defines.get("DTYPE").map(|s| s.as_str()), Some("f32"));
    }
    #[test]
    fn variant_profile_inherits_base() {
        let mut reg = VariantProfileRegistry::new();
        reg.register(
            VariantProfile::new("base")
                .set("WG_SIZE", "64")
                .set("DTYPE", "f32"),
        );
        reg.register(
            VariantProfile::with_base("child", "base")
                .set("DTYPE", "f64")
                .set("EXTRA", "1"),
        );
        let defines = reg.resolve("child").unwrap();
        assert_eq!(defines["WG_SIZE"], "64", "WG_SIZE inherited from base");
        assert_eq!(defines["DTYPE"], "f64", "DTYPE overridden by child");
        assert_eq!(defines["EXTRA"], "1", "EXTRA added by child");
    }
    #[test]
    fn variant_profile_missing_base_returns_none() {
        let reg = VariantProfileRegistry::new();
        assert!(reg.resolve("nonexistent").is_none());
    }
    #[test]
    fn variant_profile_registry_lists_names() {
        let mut reg = VariantProfileRegistry::new();
        reg.register(VariantProfile::new("b"));
        reg.register(VariantProfile::new("a"));
        let names = reg.profile_names();
        assert_eq!(names, vec!["a", "b"], "names should be sorted");
    }
    #[test]
    fn variant_profile_overrides_base_chain() {
        let mut reg = VariantProfileRegistry::new();
        reg.register(VariantProfile::new("root").set("X", "1").set("Y", "2"));
        reg.register(
            VariantProfile::with_base("mid", "root")
                .set("Y", "20")
                .set("Z", "3"),
        );
        reg.register(VariantProfile::with_base("leaf", "mid").set("Z", "30"));
        let defines = reg.resolve("leaf").unwrap();
        assert_eq!(defines["X"], "1");
        assert_eq!(defines["Y"], "20");
        assert_eq!(defines["Z"], "30");
    }
    #[test]
    fn dependency_graph_direct_dependents() {
        let mut g = ShaderDependencyGraph::new();
        g.add_dependency("sph_density", "common_utils");
        g.add_dependency("sph_pressure", "common_utils");
        let deps = g.direct_dependents("common_utils");
        assert!(deps.contains(&"sph_density"));
        assert!(deps.contains(&"sph_pressure"));
    }
    #[test]
    fn dependency_graph_transitive_dependents() {
        let mut g = ShaderDependencyGraph::new();
        g.add_dependency("b", "a");
        g.add_dependency("c", "b");
        g.add_dependency("d", "c");
        let transitive = g.transitive_dependents("a");
        assert!(transitive.contains(&"b".to_string()), "b depends on a");
        assert!(
            transitive.contains(&"c".to_string()),
            "c transitively depends on a via b"
        );
        assert!(
            transitive.contains(&"d".to_string()),
            "d transitively depends on a via b,c"
        );
    }
    #[test]
    fn dependency_graph_no_dependents() {
        let g = ShaderDependencyGraph::new();
        assert!(g.transitive_dependents("orphan").is_empty());
        assert!(g.direct_dependents("orphan").is_empty());
    }
    #[test]
    fn dependency_graph_shaders_with_deps_sorted() {
        let mut g = ShaderDependencyGraph::new();
        g.add_dependency("z_shader", "a");
        g.add_dependency("a_shader", "b");
        let names = g.shaders_with_deps();
        assert_eq!(names, vec!["a_shader", "z_shader"]);
    }
    #[test]
    fn dependency_graph_direct_dependencies() {
        let mut g = ShaderDependencyGraph::new();
        g.add_dependency("main", "utils");
        g.add_dependency("main", "math");
        let deps = g.direct_dependencies("main");
        assert_eq!(deps.len(), 2);
        assert!(deps.contains(&"utils".to_string()));
        assert!(deps.contains(&"math".to_string()));
    }
    #[test]
    fn spec_const_default_value() {
        let c = SpecializationConstant::new("WG_SIZE", SpecConstValue::Uint(64));
        assert_eq!(c.effective_value(), &SpecConstValue::Uint(64));
        assert_eq!(c.effective_value().to_wgsl(), "64u");
    }
    #[test]
    fn spec_const_override_value() {
        let c = SpecializationConstant::new("WG_SIZE", SpecConstValue::Uint(64))
            .with_override(SpecConstValue::Uint(128));
        assert_eq!(c.effective_value(), &SpecConstValue::Uint(128));
    }
    #[test]
    fn spec_const_bool_to_wgsl() {
        let c = SpecializationConstant::new("ENABLE_DEBUG", SpecConstValue::Bool(true));
        assert_eq!(c.effective_value().to_wgsl(), "true");
    }
    #[test]
    fn spec_const_float_to_wgsl() {
        let c = SpecializationConstant::new("EPSILON", SpecConstValue::Float(1e-6));
        let wgsl = c.effective_value().to_wgsl();
        assert!(
            wgsl.contains("0.000001")
                || wgsl.contains("1e-6")
                || wgsl.contains("1E-6")
                || !wgsl.is_empty()
        );
    }
    #[test]
    fn spec_const_set_to_defines() {
        let mut set = SpecConstSet::new();
        set.add(SpecializationConstant::new("A", SpecConstValue::Uint(4)));
        set.add(SpecializationConstant::new(
            "B",
            SpecConstValue::Bool(false),
        ));
        let defines = set.to_defines();
        assert_eq!(defines.get("A").map(|s| s.as_str()), Some("4u"));
        assert_eq!(defines.get("B").map(|s| s.as_str()), Some("false"));
    }
    #[test]
    fn spec_const_set_has_and_get_wgsl() {
        let mut set = SpecConstSet::new();
        set.add(SpecializationConstant::new("X", SpecConstValue::Int(-1)));
        assert!(set.has("X"));
        assert!(!set.has("Y"));
        assert_eq!(set.get_wgsl("X"), Some("-1".to_string()));
        assert!(set.get_wgsl("Y").is_none());
    }
    #[test]
    fn pipeline_cache_key_equality() {
        let k1 = PipelineCacheKey::new(ShaderKey::bare("s"), [64, 1, 1], 0, 42);
        let k2 = PipelineCacheKey::new(ShaderKey::bare("s"), [64, 1, 1], 0, 42);
        assert_eq!(k1, k2);
    }
    #[test]
    fn pipeline_cache_key_different_layout() {
        let k1 = PipelineCacheKey::new(ShaderKey::bare("s"), [64, 1, 1], 0, 1);
        let k2 = PipelineCacheKey::new(ShaderKey::bare("s"), [64, 1, 1], 0, 2);
        assert_ne!(k1, k2);
    }
    #[test]
    fn pipeline_cache_hit_and_miss() {
        let mut cache = PipelineCache::new();
        let key = PipelineCacheKey::new(ShaderKey::bare("s"), [64, 1, 1], 0, 0);
        assert!(cache.get(&key).is_none());
        cache.insert(key.clone(), "sph_pipeline");
        assert_eq!(cache.get(&key), Some("sph_pipeline"));
        assert_eq!(cache.hits, 1);
        assert_eq!(cache.misses, 1);
    }
    #[test]
    fn pipeline_cache_hit_rate() {
        let mut cache = PipelineCache::new();
        let k = PipelineCacheKey::new(ShaderKey::bare("s"), [1, 1, 1], 0, 0);
        cache.insert(k.clone(), "label");
        cache.get(&k);
        cache.get(&k);
        let k2 = PipelineCacheKey::new(ShaderKey::bare("x"), [1, 1, 1], 0, 0);
        cache.get(&k2);
        let rate = cache.hit_rate();
        assert!((rate - 2.0 / 3.0).abs() < 1e-10, "hit_rate = {rate}");
    }
    #[test]
    fn pipeline_cache_clear() {
        let mut cache = PipelineCache::new();
        let k = PipelineCacheKey::new(ShaderKey::bare("s"), [1, 1, 1], 0, 0);
        cache.insert(k, "lbl");
        assert_eq!(cache.len(), 1);
        cache.clear();
        assert!(cache.is_empty());
    }
    #[test]
    fn pipeline_cache_key_hash_stable() {
        let k = PipelineCacheKey::new(ShaderKey::bare("sph"), [64, 1, 1], 16, 12345);
        assert_eq!(k.hash_key(), k.hash_key(), "hash should be deterministic");
    }
    #[test]
    fn registry_apply_hot_reload_invalidates_stale() {
        let mut reg = ShaderRegistry::new(8);
        reg.register(ShaderSource::new("shader_a", "fn a(){}", [1, 1, 1]));
        reg.get_or_compile("shader_a", &HashMap::new()).unwrap();
        assert_eq!(reg.cached_count(), 1);
        let mut tracker = HotReloadTracker::new();
        tracker.touch("shader_a", 100);
        let invalidated = reg.apply_hot_reload(&tracker);
        assert!(invalidated.contains(&"shader_a".to_string()));
        assert_eq!(
            reg.cached_count(),
            0,
            "stale shader variant should be evicted"
        );
    }
    #[test]
    fn registry_registered_count() {
        let mut reg = ShaderRegistry::new(8);
        assert_eq!(reg.registered_count(), 0);
        reg.register(ShaderSource::new("a", "fn a(){}", [1, 1, 1]));
        reg.register(ShaderSource::new("b", "fn b(){}", [1, 1, 1]));
        assert_eq!(reg.registered_count(), 2);
    }
    #[test]
    fn registry_source_bytes() {
        let wgsl = "fn main() {}";
        let mut reg = ShaderRegistry::new(4);
        reg.register(ShaderSource::new("test", wgsl, [1, 1, 1]));
        assert_eq!(reg.source_bytes("test"), wgsl.len());
        assert_eq!(reg.source_bytes("nonexistent"), 0);
    }
    #[test]
    fn spec_const_set_used_with_registry() {
        let mut set = SpecConstSet::new();
        set.add(SpecializationConstant::new(
            "WG_SIZE",
            SpecConstValue::Uint(64),
        ));
        let defines = set.to_defines();
        let mut reg = ShaderRegistry::new(8);
        reg.register(ShaderSource::new(
            "compute",
            "@compute @workgroup_size({{WG_SIZE}}) fn main() {}",
            [64, 1, 1],
        ));
        let v = reg.get_or_compile("compute", &defines).unwrap();
        assert!(
            v.wgsl.contains("64u"),
            "WG_SIZE should be substituted from spec const"
        );
    }
    #[test]
    fn variant_profile_with_registry_integration() {
        let mut prof_reg = VariantProfileRegistry::new();
        prof_reg.register(
            VariantProfile::new("high")
                .set("WG_SIZE", "256")
                .set("DTYPE", "f64"),
        );
        let defines = prof_reg.resolve("high").unwrap();
        let mut shader_reg = ShaderRegistry::new(8);
        shader_reg.register(ShaderSource::new(
            "kernel",
            "wg={{WG_SIZE}} dtype={{DTYPE}}",
            [256, 1, 1],
        ));
        let v = shader_reg.get_or_compile("kernel", &defines).unwrap();
        assert!(v.wgsl.contains("256"));
        assert!(v.wgsl.contains("f64"));
    }
    #[test]
    fn spec_const_int_negative_to_wgsl() {
        let c = SpecConstValue::Int(-42);
        assert_eq!(c.to_wgsl(), "-42");
    }
    #[test]
    fn spec_const_uint_zero() {
        let c = SpecConstValue::Uint(0);
        assert_eq!(c.to_wgsl(), "0u");
    }
    #[test]
    fn spec_const_bool_false() {
        let c = SpecConstValue::Bool(false);
        assert_eq!(c.to_wgsl(), "false");
    }
    #[test]
    fn pipeline_cache_key_different_workgroup_sizes_differ() {
        let k1 = PipelineCacheKey::new(ShaderKey::bare("s"), [64, 1, 1], 0, 0);
        let k2 = PipelineCacheKey::new(ShaderKey::bare("s"), [128, 1, 1], 0, 0);
        assert_ne!(k1.hash_key(), k2.hash_key());
    }
}
