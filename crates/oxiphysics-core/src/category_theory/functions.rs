//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Compose two functions `f: A -> B` and `g: B -> C` into `g . f: A -> C`.
///
/// The resulting closure applies `f` first, then `g`.
pub fn compose<A, B, C>(f: impl Fn(A) -> B, g: impl Fn(B) -> C) -> impl Fn(A) -> C {
    move |a| g(f(a))
}
/// Trait for objects in a category.
///
/// Every object has a unique label that identifies it within its category.
pub trait CatObject: Clone + PartialEq + std::fmt::Debug {
    /// Return the unique label for this object.
    fn label(&self) -> &str;
}
/// Trait for morphisms in a category.
///
/// A morphism connects a domain object to a codomain object.
pub trait CatMorphism: Clone + std::fmt::Debug {
    /// Return the name/label of this morphism.
    fn name(&self) -> &str;
    /// Return the domain object label.
    fn domain(&self) -> &str;
    /// Return the codomain object label.
    fn codomain(&self) -> &str;
}
/// Check the triangle identity for an adjoint pair at point `x`.
///
/// For an adjunction `F -| G`, one triangle identity states that
/// `counit(unit(x)) ~ x` (up to floating-point tolerance `1e-9`).
pub fn adjoint_pair_check(unit: impl Fn(f64) -> f64, counit: impl Fn(f64) -> f64, x: f64) -> bool {
    let result = counit(unit(x));
    (result - x).abs() < 1e-9
}
/// Compose two Kleisli arrows `f: A -> Option<B>` and `g: B -> Option`C`.
///
/// The result is `g . f` in the Kleisli category for `Option`.
pub fn kleisli_compose<A, B, C>(
    f: impl Fn(A) -> Option<B>,
    g: impl Fn(B) -> Option<C>,
) -> impl Fn(A) -> Option<C> {
    move |a| f(a).and_then(&g)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::category_theory::Adjunction;
    use crate::category_theory::Category;
    use crate::category_theory::Endofunctor;
    use crate::category_theory::FreeMonoid;
    use crate::category_theory::Functor;
    use crate::category_theory::Monoid;
    use crate::category_theory::MonoidalCategory;
    use crate::category_theory::Morphism;
    use crate::category_theory::NaturalTransformation;
    use crate::category_theory::OpticLens;
    use crate::category_theory::OptionMonad;
    use crate::category_theory::PhysicalSystem;
    use crate::category_theory::PhysicsCategory;
    use crate::category_theory::ProductCategory;
    use crate::category_theory::ProfunctorOptic;
    use crate::category_theory::SymmetryTransform;
    use crate::category_theory::VecMonad;
    use crate::category_theory::YonedaLemma;
    #[test]
    fn test_compose_applies_f_then_g() {
        let add1 = |x: i32| x + 1;
        let double = |x: i32| x * 2;
        let h = compose(add1, double);
        assert_eq!(h(3), 8);
    }
    #[test]
    fn test_compose_identity() {
        let id = |x: f64| x;
        let double = |x: f64| x * 2.0;
        let h = compose(id, double);
        assert!((h(5.0) - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_compose_string_pipeline() {
        let upper = |s: String| s.to_uppercase();
        let exclaim = |s: String| format!("{}!", s);
        let shout = compose(upper, exclaim);
        assert_eq!(shout("hello".to_string()), "HELLO!");
    }
    #[test]
    fn test_compose_three_stages() {
        let f = |x: i32| x + 1;
        let g = |x: i32| x * 2;
        let h = |x: i32| x - 3;
        let fgh = compose(compose(f, g), h);
        assert_eq!(fgh(5), 9);
    }
    #[test]
    fn test_compose_negation() {
        let negate = |x: f64| -x;
        let square = |x: f64| x * x;
        let neg_then_sq = compose(negate, square);
        assert!((neg_then_sq(3.0) - 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_morphism_name() {
        let m: Morphism<i32, f64> = Morphism::new("cast");
        assert_eq!(m.name(), "cast");
    }
    #[test]
    fn test_morphism_clone() {
        let m: Morphism<u8, u16> = Morphism::new("widen");
        let m2 = m.clone();
        assert_eq!(m2.name(), "widen");
    }
    #[test]
    fn test_category_empty() {
        let cat = Category::new();
        assert_eq!(cat.object_count(), 0);
        assert_eq!(cat.morphism_count(), 0);
    }
    #[test]
    fn test_category_add_objects_creates_identity() {
        let mut cat = Category::new();
        cat.add_object("A");
        cat.add_object("B");
        assert_eq!(cat.object_count(), 2);
        assert_eq!(cat.morphism_count(), 2);
    }
    #[test]
    fn test_category_add_morphisms() {
        let mut cat = Category::new();
        cat.add_object("A");
        cat.add_object("B");
        cat.add_morphism("f", "A", "B");
        assert_eq!(cat.morphism_count(), 3);
    }
    #[test]
    fn test_category_is_valid_morphism_true() {
        let mut cat = Category::new();
        cat.add_object("X");
        cat.add_object("Y");
        cat.add_morphism("phi", "X", "Y");
        assert!(cat.is_valid_morphism("X", "Y", "phi"));
    }
    #[test]
    fn test_category_is_valid_morphism_wrong_name() {
        let mut cat = Category::new();
        cat.add_object("X");
        cat.add_object("Y");
        cat.add_morphism("phi", "X", "Y");
        assert!(!cat.is_valid_morphism("X", "Y", "psi"));
    }
    #[test]
    fn test_category_is_valid_morphism_wrong_domain() {
        let mut cat = Category::new();
        cat.add_object("X");
        cat.add_object("Y");
        cat.add_object("Z");
        cat.add_morphism("phi", "X", "Y");
        assert!(!cat.is_valid_morphism("Z", "Y", "phi"));
    }
    #[test]
    fn test_category_compose_morphisms_valid() {
        let mut cat = Category::new();
        cat.add_object("A");
        cat.add_object("B");
        cat.add_object("C");
        cat.add_morphism("f", "A", "B");
        cat.add_morphism("g", "B", "C");
        let comp = cat.compose_morphisms("f", "g");
        assert!(comp.is_some());
        let (name, dom, cod) = comp.unwrap();
        assert_eq!(name, "f;g");
        assert_eq!(dom, "A");
        assert_eq!(cod, "C");
    }
    #[test]
    fn test_category_compose_morphisms_type_mismatch() {
        let mut cat = Category::new();
        cat.add_object("A");
        cat.add_object("B");
        cat.add_object("C");
        cat.add_morphism("f", "A", "B");
        cat.add_morphism("g", "C", "A");
        assert!(cat.compose_morphisms("f", "g").is_none());
    }
    #[test]
    fn test_category_compose_unknown_morphism() {
        let cat = Category::new();
        assert!(cat.compose_morphisms("unknown", "also_unknown").is_none());
    }
    #[test]
    fn test_category_left_identity_law() {
        let mut cat = Category::new();
        cat.add_object("A");
        cat.add_object("B");
        cat.add_morphism("f", "A", "B");
        assert!(cat.check_left_identity("f"));
    }
    #[test]
    fn test_category_right_identity_law() {
        let mut cat = Category::new();
        cat.add_object("A");
        cat.add_object("B");
        cat.add_morphism("f", "A", "B");
        assert!(cat.check_right_identity("f"));
    }
    #[test]
    fn test_category_composition_table() {
        let mut cat = Category::new();
        cat.add_object("A");
        cat.add_object("B");
        cat.add_object("C");
        cat.add_morphism("f", "A", "B");
        cat.add_morphism("g", "B", "C");
        cat.add_composition("f", "g", "h");
        let comp = cat.compose_morphisms("f", "g");
        assert!(comp.is_some());
        let (name, _, _) = comp.unwrap();
        assert_eq!(name, "h");
    }
    #[test]
    fn test_category_associativity() {
        let mut cat = Category::new();
        cat.add_object("A");
        cat.add_object("B");
        cat.add_object("C");
        cat.add_object("D");
        cat.add_morphism("f", "A", "B");
        cat.add_morphism("g", "B", "C");
        cat.add_morphism("h", "C", "D");
        assert!(cat.check_associativity("f", "g", "h"));
    }
    #[test]
    fn test_category_hom_set() {
        let mut cat = Category::new();
        cat.add_object("A");
        cat.add_object("B");
        cat.add_morphism("f", "A", "B");
        cat.add_morphism("g", "A", "B");
        let hom = cat.hom_set("A", "B");
        assert_eq!(hom.len(), 2);
        assert!(hom.contains(&"f".to_string()));
        assert!(hom.contains(&"g".to_string()));
    }
    #[test]
    fn test_functor_new() {
        let f = Functor::new("List", "Set", "Set");
        assert_eq!(f.name, "List");
        assert_eq!(f.source_category, "Set");
    }
    #[test]
    fn test_functor_apply_object() {
        let mut f = Functor::new("Maybe", "C", "C");
        f.map_object("Int", "Maybe<Int>");
        assert_eq!(f.apply_object("Int"), Some("Maybe<Int>"));
    }
    #[test]
    fn test_functor_apply_morphism() {
        let mut f = Functor::new("Maybe", "C", "C");
        f.map_morphism("succ", "fmap(succ)");
        assert_eq!(f.apply_morphism("succ"), Some("fmap(succ)"));
    }
    #[test]
    fn test_functor_apply_object_missing() {
        let f = Functor::new("F", "C", "D");
        assert_eq!(f.apply_object("X"), None);
    }
    #[test]
    fn test_functor_identity_law() {
        let mut f = Functor::new("F", "C", "D");
        f.map_object("A", "FA");
        f.map_morphism("id_A", "id_FA");
        assert!(f.check_identity_law("A"));
    }
    #[test]
    fn test_functor_composition_law() {
        let mut f = Functor::new("F", "C", "D");
        f.map_morphism("f", "Ff");
        f.map_morphism("g", "Fg");
        f.map_morphism("f;g", "Ff;Fg");
        assert!(f.check_composition_law("f", "g", "f;g"));
    }
    #[test]
    fn test_natural_transformation_new() {
        let nt = NaturalTransformation::new("F", "G");
        assert_eq!(nt.source, "F");
        assert_eq!(nt.target, "G");
        assert_eq!(nt.component_count(), 0);
    }
    #[test]
    fn test_natural_transformation_add_components() {
        let mut nt = NaturalTransformation::new("F", "G");
        nt.add_component("alpha_A");
        nt.add_component("alpha_B");
        assert_eq!(nt.component_count(), 2);
    }
    #[test]
    fn test_natural_transformation_component_for() {
        let mut nt = NaturalTransformation::new("F", "G");
        nt.add_component_for("A", "eta_A");
        nt.add_component_for("B", "eta_B");
        assert_eq!(nt.component_for("A"), Some("eta_A"));
        assert_eq!(nt.component_for("B"), Some("eta_B"));
        assert_eq!(nt.component_for("C"), None);
    }
    #[test]
    fn test_natural_transformation_naturality_check() {
        let nt = NaturalTransformation::new("F", "G");
        assert!(nt.check_naturality("f", "eta_B_after_Ff", "eta_B_after_Ff"));
        assert!(!nt.check_naturality("f", "eta_B_after_Ff", "Gf_after_eta_A"));
    }
    #[test]
    fn test_natural_transformation_vertical_compose() {
        let mut alpha = NaturalTransformation::new("F", "G");
        alpha.add_component_for("A", "alpha_A");
        let mut beta = NaturalTransformation::new("G", "H");
        beta.add_component_for("A", "beta_A");
        let gamma = alpha.vertical_compose(&beta);
        assert_eq!(gamma.source, "F");
        assert_eq!(gamma.target, "H");
        assert_eq!(gamma.component_for("A"), Some("alpha_A;beta_A"));
    }
    #[test]
    fn test_product_category_objects() {
        let mut c = Category::new();
        c.add_object("A");
        c.add_object("B");
        let mut d = Category::new();
        d.add_object("X");
        d.add_object("Y");
        let pc = ProductCategory::new(c, d);
        assert_eq!(pc.object_count(), 4);
    }
    #[test]
    fn test_product_category_morphisms() {
        let mut c = Category::new();
        c.add_object("A");
        c.add_object("B");
        c.add_morphism("f", "A", "B");
        let mut d = Category::new();
        d.add_object("X");
        let pc = ProductCategory::new(c, d);
        assert_eq!(pc.morphism_count(), 3);
    }
    #[test]
    fn test_product_category_projections() {
        let mut c = Category::new();
        c.add_object("A");
        let mut d = Category::new();
        d.add_object("X");
        let pc = ProductCategory::new(c, d);
        assert_eq!(pc.project_first(0), Some("A"));
        assert_eq!(pc.project_second(0), Some("X"));
    }
    #[test]
    fn test_product_category_bifunctor() {
        let mut c = Category::new();
        c.add_object("A");
        let mut d = Category::new();
        d.add_object("X");
        let pc = ProductCategory::new(c, d);
        let mut f = Functor::new("F", "C", "E");
        f.map_object("A", "FA");
        let mut g = Functor::new("G", "D", "E");
        g.map_object("X", "GX");
        let result = pc.apply_bifunctor(&f, &g);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], ("FA".to_string(), "GX".to_string()));
    }
    #[test]
    fn test_monoidal_category_tensor() {
        let mut cat = Category::new();
        cat.add_object("I");
        cat.add_object("A");
        cat.add_object("B");
        cat.add_object("A_tensor_B");
        let mut mc = MonoidalCategory::new(cat, "I");
        mc.add_tensor("A", "B", "A_tensor_B");
        assert_eq!(mc.tensor("A", "B"), Some("A_tensor_B"));
    }
    #[test]
    fn test_monoidal_category_unit() {
        let cat = Category::new();
        let mc = MonoidalCategory::new(cat, "I");
        assert_eq!(mc.unit(), "I");
    }
    #[test]
    fn test_monoidal_category_left_unit() {
        let mut cat = Category::new();
        cat.add_object("I");
        cat.add_object("A");
        let mut mc = MonoidalCategory::new(cat, "I");
        mc.add_tensor("I", "A", "A");
        assert!(mc.check_left_unit("A"));
    }
    #[test]
    fn test_monoidal_category_right_unit() {
        let mut cat = Category::new();
        cat.add_object("I");
        cat.add_object("A");
        let mut mc = MonoidalCategory::new(cat, "I");
        mc.add_tensor("A", "I", "A");
        assert!(mc.check_right_unit("A"));
    }
    #[test]
    fn test_monoidal_category_pentagon() {
        let cat = Category::new();
        let mut mc = MonoidalCategory::new(cat, "I");
        mc.add_associator("A", "B", "C", "alpha_ABC");
        mc.add_associator("B", "C", "D", "alpha_BCD");
        assert!(mc.check_pentagon("A", "B", "C", "D"));
    }
    #[test]
    fn test_monoidal_category_tensor_morphisms() {
        let cat = Category::new();
        let mc = MonoidalCategory::new(cat, "I");
        let name = mc.tensor_morphisms("f", "g");
        assert_eq!(name, "f_tensor_g");
    }
    #[test]
    fn test_adjunction_triangle_identity_1() {
        let mut left = Functor::new("F", "C", "D");
        left.map_object("A", "FA");
        let mut right = Functor::new("G", "D", "C");
        right.map_object("FA", "GFA");
        let mut unit_nt = NaturalTransformation::new("Id_C", "GF");
        unit_nt.add_component_for("A", "eta_A");
        let mut counit_nt = NaturalTransformation::new("FG", "Id_D");
        counit_nt.add_component_for("FA", "epsilon_FA");
        let adj = Adjunction::new(left, right, unit_nt, counit_nt);
        assert!(adj.check_triangle_identity_1("A"));
    }
    #[test]
    fn test_adjunction_triangle_identity_2() {
        let mut left = Functor::new("F", "C", "D");
        left.map_object("GB", "FGB");
        let mut right = Functor::new("G", "D", "C");
        right.map_object("B", "GB");
        let mut unit_nt = NaturalTransformation::new("Id_C", "GF");
        unit_nt.add_component_for("GB", "eta_GB");
        let mut counit_nt = NaturalTransformation::new("FG", "Id_D");
        counit_nt.add_component_for("B", "epsilon_B");
        let adj = Adjunction::new(left, right, unit_nt, counit_nt);
        assert!(adj.check_triangle_identity_2("B"));
    }
    #[test]
    fn test_adjunction_both_triangles() {
        let mut left = Functor::new("F", "C", "D");
        left.map_object("A", "FA");
        left.map_object("GB", "FGB");
        let mut right = Functor::new("G", "D", "C");
        right.map_object("FA", "GFA");
        right.map_object("B", "GB");
        let mut unit_nt = NaturalTransformation::new("Id_C", "GF");
        unit_nt.add_component_for("A", "eta_A");
        unit_nt.add_component_for("GB", "eta_GB");
        let mut counit_nt = NaturalTransformation::new("FG", "Id_D");
        counit_nt.add_component_for("FA", "epsilon_FA");
        counit_nt.add_component_for("B", "epsilon_B");
        let adj = Adjunction::new(left, right, unit_nt, counit_nt);
        assert!(adj.check_both_triangles("A", "B"));
    }
    #[test]
    fn test_adjunction_transpose() {
        let mut left = Functor::new("F", "C", "D");
        left.map_object("A", "FA");
        let mut right = Functor::new("G", "D", "C");
        right.map_object("B", "GB");
        let unit_nt = NaturalTransformation::new("Id_C", "GF");
        let counit_nt = NaturalTransformation::new("FG", "Id_D");
        let adj = Adjunction::new(left, right, unit_nt, counit_nt);
        let t = adj.transpose_left("A", "B", "f");
        assert_eq!(t, Some("transpose_left(f)".to_string()));
    }
    #[test]
    fn test_yoneda_hom_sets() {
        let mut cat = Category::new();
        cat.add_object("A");
        cat.add_object("B");
        cat.add_morphism("f", "A", "B");
        let yl = YonedaLemma::new("A", cat);
        let hom_ab = yl.hom("B").unwrap();
        assert!(hom_ab.contains(&"f".to_string()));
    }
    #[test]
    fn test_yoneda_identity_in_hom() {
        let mut cat = Category::new();
        cat.add_object("A");
        let yl = YonedaLemma::new("A", cat);
        assert!(yl.check_representable_identity());
    }
    #[test]
    fn test_yoneda_evaluate_at_identity() {
        let mut cat = Category::new();
        cat.add_object("A");
        let yl = YonedaLemma::new("A", cat);
        assert_eq!(yl.evaluate_at_identity(), Some("id_A".to_string()));
    }
    #[test]
    fn test_yoneda_embedding() {
        let mut cat = Category::new();
        cat.add_object("A");
        cat.add_object("B");
        cat.add_morphism("f", "A", "B");
        cat.add_morphism("g", "B", "B");
        let yl = YonedaLemma::new("A", cat);
        let emb = yl.yoneda_embedding("f", "B");
        let precomposed_bb = emb.get("B").unwrap();
        assert!(precomposed_bb.contains(&"f;g".to_string()));
    }
    #[test]
    fn test_physics_category_add_system() {
        let mut pc = PhysicsCategory::new();
        let sys = PhysicalSystem::new("Particle", 3, 10.0, "SO(3)");
        pc.add_system(sys);
        assert_eq!(pc.system_count(), 1);
        assert_eq!(pc.transform_count(), 1);
    }
    #[test]
    fn test_physics_category_identity_law() {
        let mut pc = PhysicsCategory::new();
        pc.add_system(PhysicalSystem::new("A", 3, 1.0, "SO(3)"));
        pc.add_system(PhysicalSystem::new("B", 3, 1.0, "SO(3)"));
        let rot = SymmetryTransform::new(
            "rot90",
            "A",
            "B",
            [0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            true,
        );
        pc.add_transform(rot);
        assert!(pc.check_identity_law("rot90"));
    }
    #[test]
    fn test_physics_category_associativity() {
        let mut pc = PhysicsCategory::new();
        pc.add_system(PhysicalSystem::new("A", 3, 1.0, "SO(3)"));
        pc.add_system(PhysicalSystem::new("B", 3, 1.0, "SO(3)"));
        pc.add_system(PhysicalSystem::new("C", 3, 1.0, "SO(3)"));
        pc.add_system(PhysicalSystem::new("D", 3, 1.0, "SO(3)"));
        let id_mat = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        pc.add_transform(SymmetryTransform::new("f", "A", "B", id_mat, true));
        pc.add_transform(SymmetryTransform::new("g", "B", "C", id_mat, true));
        pc.add_transform(SymmetryTransform::new("h", "C", "D", id_mat, true));
        assert!(pc.check_associativity("f", "g", "h"));
    }
    #[test]
    fn test_symmetry_transform_identity() {
        let id = SymmetryTransform::identity("Particle");
        assert!(id.is_identity(1e-12));
        assert!(id.is_rotation(1e-12));
    }
    #[test]
    fn test_symmetry_transform_compose() {
        let t1 = SymmetryTransform::new(
            "rot",
            "A",
            "B",
            [0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            true,
        );
        let t2 = SymmetryTransform::new(
            "rot2",
            "B",
            "C",
            [0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            true,
        );
        let composed = t1.compose(&t2).unwrap();
        assert_eq!(composed.source, "A");
        assert_eq!(composed.target, "C");
        assert!((composed.matrix[0] - (-1.0)).abs() < 1e-9);
    }
    #[test]
    fn test_symmetry_transform_determinant() {
        let rot = SymmetryTransform::new(
            "rot",
            "A",
            "A",
            [0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            true,
        );
        assert!((rot.determinant() - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_physics_category_energy_preserving() {
        let mut pc = PhysicsCategory::new();
        pc.add_system(PhysicalSystem::new("A", 3, 1.0, "SO(3)"));
        pc.add_system(PhysicalSystem::new("B", 3, 1.0, "SO(3)"));
        let id_mat = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        pc.add_transform(SymmetryTransform::new("f", "A", "B", id_mat, true));
        pc.add_transform(SymmetryTransform::new("g", "A", "B", id_mat, false));
        let ep = pc.energy_preserving_transforms();
        assert_eq!(ep.len(), 3);
    }
    #[test]
    fn test_physics_category_to_category() {
        let mut pc = PhysicsCategory::new();
        pc.add_system(PhysicalSystem::new("A", 3, 1.0, "SO(3)"));
        pc.add_system(PhysicalSystem::new("B", 3, 1.0, "SO(3)"));
        let id_mat = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        pc.add_transform(SymmetryTransform::new("f", "A", "B", id_mat, true));
        let cat = pc.to_category();
        assert_eq!(cat.object_count(), 2);
        assert!(cat.is_valid_morphism("A", "B", "f"));
    }
    #[test]
    fn test_option_monad_unit() {
        let m = OptionMonad::unit(42);
        assert_eq!(m, Some(42));
    }
    #[test]
    fn test_option_monad_bind_some() {
        let m = Some(10);
        let result = OptionMonad::bind(m, |x| Some(x * 2));
        assert_eq!(result, Some(20));
    }
    #[test]
    fn test_option_monad_bind_none() {
        let m: Option<i32> = None;
        let result = OptionMonad::bind(m, |x| Some(x * 2));
        assert_eq!(result, None);
    }
    #[test]
    fn test_option_monad_left_identity() {
        let x = 5;
        let f = |v: i32| Some(v + 1);
        let lhs = OptionMonad::bind(OptionMonad::unit(x), f);
        let rhs = f(x);
        assert_eq!(lhs, rhs);
    }
    #[test]
    fn test_option_monad_right_identity() {
        let m = Some(7);
        let result = OptionMonad::bind(m, OptionMonad::unit);
        assert_eq!(result, m);
    }
    #[test]
    fn test_option_monad_fmap() {
        let m = Some(3);
        let result = OptionMonad::fmap(m, |x| x * x);
        assert_eq!(result, Some(9));
    }
    #[test]
    fn test_vec_monad_unit() {
        assert_eq!(VecMonad::unit(42), vec![42]);
    }
    #[test]
    fn test_vec_monad_bind() {
        let xs = vec![1, 2, 3];
        let result = VecMonad::bind(xs, |x| vec![x, x * 10]);
        assert_eq!(result, vec![1, 10, 2, 20, 3, 30]);
    }
    #[test]
    fn test_vec_monad_bind_empty() {
        let xs: Vec<i32> = vec![];
        let result = VecMonad::bind(xs, |x| vec![x, x + 1]);
        assert!(result.is_empty());
    }
    #[test]
    fn test_vec_monad_fmap() {
        let xs = vec![1, 2, 3];
        let result = VecMonad::fmap(xs, |x| x * 2);
        assert_eq!(result, vec![2, 4, 6]);
    }
    #[test]
    fn test_vec_monad_left_identity() {
        let x = 3;
        let f = |v: i32| vec![v, v + 1, v + 2];
        let lhs = VecMonad::bind(VecMonad::unit(x), f);
        let rhs = f(x);
        assert_eq!(lhs, rhs);
    }
    #[test]
    fn test_free_monoid_empty_is_identity() {
        let e: FreeMonoid<char> = FreeMonoid::empty();
        assert!(e.is_empty());
        assert_eq!(e.len(), 0);
    }
    #[test]
    fn test_free_monoid_singleton() {
        let s = FreeMonoid::singleton('a');
        assert_eq!(s.len(), 1);
        assert_eq!(s.word, vec!['a']);
    }
    #[test]
    fn test_free_monoid_concat_identity() {
        let e: FreeMonoid<i32> = FreeMonoid::empty();
        let w = FreeMonoid::from_slice(&[1, 2, 3]);
        assert_eq!(e.concat(&w), w);
        assert_eq!(w.concat(&e), w);
    }
    #[test]
    fn test_free_monoid_concat_associativity() {
        let a = FreeMonoid::from_slice(&[1, 2]);
        let b = FreeMonoid::from_slice(&[3, 4]);
        let c = FreeMonoid::from_slice(&[5, 6]);
        assert_eq!(a.concat(&b).concat(&c), a.concat(&b.concat(&c)));
    }
    #[test]
    fn test_free_monoid_power() {
        let w = FreeMonoid::from_slice(&[1, 2]);
        let w3 = w.power(3);
        assert_eq!(w3.word, vec![1, 2, 1, 2, 1, 2]);
    }
    #[test]
    fn test_free_monoid_power_zero() {
        let w: FreeMonoid<u8> = FreeMonoid::from_slice(&[5]);
        assert!(w.power(0).is_empty());
    }
    #[test]
    fn test_free_monoid_reverse() {
        let w = FreeMonoid::from_slice(&[1, 2, 3]);
        let rev = w.reverse();
        assert_eq!(rev.word, vec![3, 2, 1]);
    }
    #[test]
    fn test_free_monoid_subword() {
        let w = FreeMonoid::from_slice(&[0, 1, 2, 3, 4]);
        let sub = w.subword(1, 4);
        assert_eq!(sub.word, vec![1, 2, 3]);
    }
    #[test]
    fn test_lens_get() {
        let lens = OpticLens::new(|s: &(i32, i32)| s.0, |s, a| (a, s.1));
        assert_eq!(lens.get(&(42, 7)), 42);
    }
    #[test]
    fn test_lens_set() {
        let lens = OpticLens::new(|s: &(i32, i32)| s.0, |s, a| (a, s.1));
        let result = lens.set((1, 2), 99);
        assert_eq!(result, (99, 2));
    }
    #[test]
    fn test_lens_put_get_law() {
        let lens = OpticLens::new(|s: &(i32, i32)| s.0, |s, a| (a, s.1));
        assert!(lens.check_put_get((10, 20), 99));
    }
    #[test]
    fn test_lens_get_put_law() {
        let lens = OpticLens::new(|s: &(i32, i32)| s.0, |s, a| (a, s.1));
        assert!(lens.check_get_put((42, 7)));
    }
    #[test]
    fn test_lens_modify() {
        let lens = OpticLens::new(|s: &(i32, i32)| s.0, |s, a| (a, s.1));
        let result = lens.modify((5, 3), |x| x * 2);
        assert_eq!(result, (10, 3));
    }
    #[test]
    fn test_lens_string_field() {
        let lens = OpticLens::<(String, i32), String>::new(
            |s: &(String, i32)| s.0.clone(),
            |s, a| (a, s.1),
        );
        let s = ("hello".to_string(), 42);
        assert_eq!(lens.get(&s), "hello");
        let s2 = lens.set(s, "world".to_string());
        assert_eq!(s2.0, "world");
    }
    #[test]
    fn test_profunctor_optic_apply() {
        let optic = ProfunctorOptic::<(i32, i32), (i32, i32), i32, i32>::new(
            |a: i32| a * 2,
            |fwd, s: (i32, i32)| (fwd(s.0), s.1),
        );
        let result = optic.apply((3, 7));
        assert_eq!(result, (6, 7));
    }
    #[test]
    fn test_profunctor_optic_string_transform() {
        let optic = ProfunctorOptic::<String, String, String, String>::new(
            |s: String| s.to_uppercase(),
            |fwd, s: String| fwd(s),
        );
        let result = optic.apply("hello".to_string());
        assert_eq!(result, "HELLO");
    }
    #[test]
    fn test_monoid_integer_addition_identity() {
        let m = Monoid::new(0i32, |a, b| a + b);
        assert!(m.is_identity(&0));
        assert!(!m.is_identity(&1));
    }
    #[test]
    fn test_monoid_combine() {
        let m = Monoid::new(0i32, |a, b| a + b);
        assert_eq!(m.combine(3, 4), 7);
    }
    #[test]
    fn test_monoid_left_identity() {
        let m = Monoid::new(0i32, |a, b| a + b);
        assert_eq!(m.combine(m.identity, 5), 5);
    }
    #[test]
    fn test_monoid_right_identity() {
        let m = Monoid::new(0i32, |a, b| a + b);
        assert_eq!(m.combine(7, m.identity), 7);
    }
    #[test]
    fn test_monoid_associativity() {
        let m = Monoid::new(0i32, |a, b| a + b);
        let lhs = m.combine(m.combine(1, 2), 3);
        let rhs = m.combine(1, m.combine(2, 3));
        assert_eq!(lhs, rhs);
    }
    #[test]
    fn test_monoid_power_zero() {
        let m = Monoid::new(0i32, |a, b| a + b);
        assert_eq!(m.power(5, 0), 0);
    }
    #[test]
    fn test_monoid_power_one() {
        let m = Monoid::new(0i32, |a, b| a + b);
        assert_eq!(m.power(5, 1), 5);
    }
    #[test]
    fn test_monoid_power_three() {
        let m = Monoid::new(0i32, |a, b| a + b);
        assert_eq!(m.power(5, 3), 15);
    }
    #[test]
    fn test_monoid_string_concatenation() {
        let m = Monoid::new(String::new(), |a: String, b: String| a + &b);
        assert_eq!(m.combine("foo".to_string(), "bar".to_string()), "foobar");
    }
    #[test]
    fn test_monoid_multiplication_identity() {
        let m = Monoid::new(1i32, |a, b| a * b);
        assert!(m.is_identity(&1));
        assert_eq!(m.power(2, 4), 16);
    }
    #[test]
    fn test_monoid_fold() {
        let m = Monoid::new(0i32, |a, b| a + b);
        let xs = vec![1, 2, 3, 4, 5];
        assert_eq!(m.fold(&xs), 15);
    }
    #[test]
    fn test_monoid_check_left_identity() {
        let m = Monoid::new(0i32, |a, b| a + b);
        assert!(m.check_left_identity(&7));
    }
    #[test]
    fn test_monoid_check_right_identity() {
        let m = Monoid::new(0i32, |a, b| a + b);
        assert!(m.check_right_identity(&7));
    }
    #[test]
    fn test_monoid_check_associativity() {
        let m = Monoid::new(0i32, |a, b| a + b);
        assert!(m.check_associativity(1, 2, 3));
    }
    #[test]
    fn test_kleisli_both_some() {
        let f = |x: i32| if x > 0 { Some(x * 2) } else { None };
        let g = |x: i32| if x < 100 { Some(x + 1) } else { None };
        let h = kleisli_compose(f, g);
        assert_eq!(h(3), Some(7));
    }
    #[test]
    fn test_kleisli_f_returns_none() {
        let f = |_x: i32| None::<i32>;
        let g = |x: i32| Some(x + 1);
        let h = kleisli_compose(f, g);
        assert_eq!(h(5), None);
    }
    #[test]
    fn test_kleisli_g_returns_none() {
        let f = |x: i32| Some(x * 2);
        let g = |_x: i32| None::<i32>;
        let h = kleisli_compose(f, g);
        assert_eq!(h(5), None);
    }
    #[test]
    fn test_kleisli_both_none() {
        let f = |_x: i32| None::<i32>;
        let g = |_x: i32| None::<i32>;
        let h = kleisli_compose(f, g);
        assert_eq!(h(5), None);
    }
    #[test]
    fn test_kleisli_string_pipeline() {
        let parse = |s: &str| s.parse::<i32>().ok();
        let double = |n: i32| if n < 1000 { Some(n * 2) } else { None };
        let h = kleisli_compose(parse, double);
        assert_eq!(h("21"), Some(42));
        assert_eq!(h("abc"), None);
    }
    #[test]
    fn test_adjoint_pair_identity_functions() {
        let ok = adjoint_pair_check(|x| x, |x| x, 3.125);
        assert!(ok);
    }
    #[test]
    fn test_adjoint_pair_scale_inverse() {
        let ok = adjoint_pair_check(|x| x * 2.0, |x| x / 2.0, 7.0);
        assert!(ok);
    }
    #[test]
    fn test_adjoint_pair_fails_when_not_inverse() {
        let ok = adjoint_pair_check(|x| x * 2.0, |x| x * 2.0, 1.0);
        assert!(!ok);
    }
    #[test]
    fn test_endofunctor_apply() {
        let ef = Endofunctor::new("List", |s| format!("List({})", s));
        assert_eq!(ef.apply("Int"), "List(Int)");
    }
    #[test]
    fn test_endofunctor_name() {
        let ef = Endofunctor::new("Maybe", |s| format!("Maybe({})", s));
        assert_eq!(ef.name(), "Maybe");
    }
    #[test]
    fn test_endofunctor_identity() {
        let ef = Endofunctor::new("Id", |s| s.to_string());
        assert_eq!(ef.apply("X"), "X");
    }
}
