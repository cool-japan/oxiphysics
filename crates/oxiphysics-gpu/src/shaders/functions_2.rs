//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod shader_meta_tests {
    use super::super::functions::*;
    use crate::shaders::BytecodeShaderCache;
    use crate::shaders::ShaderMetaRegistry;
    use crate::shaders::ShaderMetadata;
    use crate::shaders::ShaderTemplateV2;
    use crate::shaders::ShaderVariant;
    use std::collections::HashMap;
    #[test]
    fn test_registry_lookup_after_register() {
        let mut reg = ShaderMetaRegistry::new();
        let meta = ShaderMetadata::new(ShaderVariant::Physics, "main", [64, 1, 1], 1);
        reg.register("physics_kernel", meta);
        let found = reg.lookup("physics_kernel");
        assert!(found.is_some());
        assert_eq!(found.unwrap().variant, ShaderVariant::Physics);
        assert_eq!(found.unwrap().entry_point, "main");
        assert!(reg.lookup("nonexistent").is_none());
    }
    #[test]
    fn test_registry_all_names() {
        let mut reg = ShaderMetaRegistry::new();
        reg.register(
            "a",
            ShaderMetadata::new(ShaderVariant::Sph, "main", [64, 1, 1], 2),
        );
        reg.register(
            "b",
            ShaderMetadata::new(ShaderVariant::Lbm, "cs_main", [32, 1, 1], 3),
        );
        let mut names = reg.all_names();
        names.sort();
        assert_eq!(names, vec!["a", "b"]);
    }
    #[test]
    fn test_bytecode_cache_insert_and_get() {
        let mut cache = BytecodeShaderCache::new(1024);
        cache.insert("kernel_a", vec![0u8; 100]);
        assert!(cache.get("kernel_a").is_some());
        assert_eq!(cache.get("kernel_a").unwrap().len(), 100);
        assert_eq!(cache.total_bytes(), 100);
    }
    #[test]
    fn test_bytecode_cache_eviction() {
        let mut cache = BytecodeShaderCache::new(150);
        cache.insert("a", vec![0u8; 100]);
        cache.insert("b", vec![0u8; 100]);
        assert!(cache.get("a").is_none(), "oldest entry should be evicted");
        assert!(cache.get("b").is_some(), "newest entry should remain");
        assert!(cache.total_bytes() <= 150);
    }
    #[test]
    fn test_template_instantiation() {
        let mut defines = HashMap::new();
        defines.insert("WORKGROUP_SIZE".to_string(), "128".to_string());
        defines.insert("DTYPE".to_string(), "f32".to_string());
        let tmpl = ShaderTemplateV2::new(
            "@compute @workgroup_size(WORKGROUP_SIZE)\nvar x: DTYPE;",
            defines,
        );
        let out = tmpl.instantiate();
        assert!(out.contains("128"), "WORKGROUP_SIZE should be replaced");
        assert!(out.contains("f32"), "DTYPE should be replaced");
        assert!(
            !out.contains("WORKGROUP_SIZE"),
            "placeholder should be gone"
        );
    }
    #[test]
    fn test_metadata_validation_valid() {
        let meta = ShaderMetadata::new(ShaderVariant::Collision, "main", [64, 1, 1], 2);
        assert!(validate_shader_metadata(&meta).is_ok());
    }
    #[test]
    fn test_metadata_validation_zero_workgroup() {
        let meta = ShaderMetadata::new(ShaderVariant::RigidBody, "main", [0, 1, 1], 1);
        assert!(validate_shader_metadata(&meta).is_err());
    }
    #[test]
    fn test_metadata_validation_too_many_threads() {
        let meta = ShaderMetadata::new(ShaderVariant::NeuralInference, "main", [32, 32, 2], 1);
        assert!(validate_shader_metadata(&meta).is_err());
    }
    #[test]
    fn test_metadata_validation_too_many_bind_groups() {
        let meta = ShaderMetadata::new(ShaderVariant::Sph, "main", [64, 1, 1], 5);
        assert!(validate_shader_metadata(&meta).is_err());
    }
    #[test]
    fn test_metadata_validation_empty_entry_point() {
        let meta = ShaderMetadata::new(ShaderVariant::Lbm, "", [64, 1, 1], 1);
        assert!(validate_shader_metadata(&meta).is_err());
    }
}
