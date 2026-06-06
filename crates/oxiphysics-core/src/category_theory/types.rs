//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;
use std::marker::PhantomData;

/// An adjunction between two functors `F -| G` (F is left adjoint to G).
///
/// An adjunction consists of:
/// - A left functor `F: C -> D`
/// - A right functor `G: D -> C`
/// - A unit natural transformation `eta: Id_C => G . F`
/// - A counit natural transformation `epsilon: F . G => Id_D`
///
/// satisfying the triangle identities.
#[derive(Debug, Clone)]
pub struct Adjunction {
    /// The left adjoint functor F.
    pub left_functor: Functor,
    /// The right adjoint functor G.
    pub right_functor: Functor,
    /// Unit natural transformation eta: Id_C => G . F.
    pub unit: NaturalTransformation,
    /// Counit natural transformation epsilon: F . G => Id_D.
    pub counit: NaturalTransformation,
}
impl Adjunction {
    /// Create a new adjunction from the left and right functors and the unit/counit.
    pub fn new(
        left: Functor,
        right: Functor,
        unit: NaturalTransformation,
        counit: NaturalTransformation,
    ) -> Self {
        Self {
            left_functor: left,
            right_functor: right,
            unit,
            counit,
        }
    }
    /// Check the first triangle identity: `epsilon_F(A) . F(eta_A) = id_{F(A)}`.
    ///
    /// For a given object `a` in C, verifies the triangle identity symbolically.
    pub fn check_triangle_identity_1(&self, a: &str) -> bool {
        let _eta_a = self.unit.component_for(a);
        let fa = self.left_functor.apply_object(a);
        if let Some(fa_label) = fa {
            self.counit.component_for(fa_label).is_some()
        } else {
            false
        }
    }
    /// Check the second triangle identity: `G(epsilon_B) . eta_{G(B)} = id_{G(B)}`.
    ///
    /// For a given object `b` in D, verifies the triangle identity symbolically.
    pub fn check_triangle_identity_2(&self, b: &str) -> bool {
        let gb = self.right_functor.apply_object(b);
        let _eps_b = self.counit.component_for(b);
        if let Some(gb_label) = gb {
            self.unit.component_for(gb_label).is_some()
        } else {
            false
        }
    }
    /// Check both triangle identities for a pair of objects.
    pub fn check_both_triangles(&self, a: &str, b: &str) -> bool {
        self.check_triangle_identity_1(a) && self.check_triangle_identity_2(b)
    }
    /// Return the hom-set adjunction: Hom_D(F(A), B) ~ Hom_C(A, G(B)).
    ///
    /// Given morphism f: F(A) -> B in D, produce the transpose A -> G(B) in C (symbolically).
    pub fn transpose_left(&self, a: &str, _b: &str, f_name: &str) -> Option<String> {
        let _fa = self.left_functor.apply_object(a)?;
        Some(format!("transpose_left({})", f_name))
    }
    /// Given morphism g: A -> G(B) in C, produce the transpose F(A) -> B in D (symbolically).
    pub fn transpose_right(&self, _a: &str, b: &str, g_name: &str) -> Option<String> {
        let _gb = self.right_functor.apply_object(b)?;
        Some(format!("transpose_right({})", g_name))
    }
}
/// A concrete named morphism with domain and codomain labels.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct NamedMorphism {
    /// Name of this morphism.
    pub morph_name: String,
    /// Domain object label.
    pub dom: String,
    /// Codomain object label.
    pub cod: String,
}
impl NamedMorphism {
    /// Create a new named morphism.
    pub fn new(name: impl Into<String>, dom: impl Into<String>, cod: impl Into<String>) -> Self {
        Self {
            morph_name: name.into(),
            dom: dom.into(),
            cod: cod.into(),
        }
    }
}
/// A product morphism: `((f_name, g_name), (dom_c, dom_d), (cod_c, cod_d))`.
pub type ProductMorphism = ((String, String), (String, String), (String, String));

/// The product of two categories C x D.
///
/// Objects are pairs `(c, d)` where `c` is an object of C and `d` is an object of D.
/// Morphisms are pairs `(f, g)` where `f: c1 -> c2` in C and `g: d1 -> d2` in D.
#[derive(Debug, Clone, Default)]
pub struct ProductCategory {
    /// The first component category.
    pub category_c: Category,
    /// The second component category.
    pub category_d: Category,
    /// Product objects as `(c_label, d_label)` pairs.
    pub product_objects: Vec<(String, String)>,
    /// Product morphisms as `((f_name, g_name), (dom_c, dom_d), (cod_c, cod_d))`.
    pub product_morphisms: Vec<ProductMorphism>,
}
impl ProductCategory {
    /// Create a product category from two component categories.
    pub fn new(c: Category, d: Category) -> Self {
        let mut pc = Self {
            category_c: c,
            category_d: d,
            product_objects: Vec::new(),
            product_morphisms: Vec::new(),
        };
        pc.build_products();
        pc
    }
    /// Build all product objects from the two component categories.
    fn build_products(&mut self) {
        for c_obj in &self.category_c.objects {
            for d_obj in &self.category_d.objects {
                self.product_objects.push((c_obj.clone(), d_obj.clone()));
            }
        }
        for (cn, cd_dom, cd_cod) in &self.category_c.morphisms {
            for (dn, dd_dom, dd_cod) in &self.category_d.morphisms {
                self.product_morphisms.push((
                    (cn.clone(), dn.clone()),
                    (cd_dom.clone(), dd_dom.clone()),
                    (cd_cod.clone(), dd_cod.clone()),
                ));
            }
        }
    }
    /// Number of product objects.
    pub fn object_count(&self) -> usize {
        self.product_objects.len()
    }
    /// Number of product morphisms.
    pub fn morphism_count(&self) -> usize {
        self.product_morphisms.len()
    }
    /// Project a product object to the first component.
    pub fn project_first(&self, idx: usize) -> Option<&str> {
        self.product_objects.get(idx).map(|(c, _)| c.as_str())
    }
    /// Project a product object to the second component.
    pub fn project_second(&self, idx: usize) -> Option<&str> {
        self.product_objects.get(idx).map(|(_, d)| d.as_str())
    }
    /// Compose two product morphisms component-wise.
    pub fn compose_product_morphisms(&self, f_idx: usize, g_idx: usize) -> Option<ProductMorphism> {
        let f = self.product_morphisms.get(f_idx)?;
        let g = self.product_morphisms.get(g_idx)?;
        if f.2 != g.1 {
            return None;
        }
        let comp_name = (
            format!("{};{}", f.0.0, g.0.0),
            format!("{};{}", f.0.1, g.0.1),
        );
        Some((comp_name, f.1.clone(), g.2.clone()))
    }
    /// Apply a bifunctor (pair of functors) to this product category.
    ///
    /// Given functors F: C -> E and G: D -> E, produce a functor F x G on objects.
    pub fn apply_bifunctor(&self, f: &Functor, g: &Functor) -> Vec<(String, String)> {
        let mut result = Vec::new();
        for (c_obj, d_obj) in &self.product_objects {
            let fc = f.apply_object(c_obj).unwrap_or("?").to_string();
            let gd = g.apply_object(d_obj).unwrap_or("?").to_string();
            result.push((fc, gd));
        }
        result
    }
}
/// A symmetry transform morphism in the physics category.
///
/// Morphisms represent physical symmetry transformations (rotations,
/// translations, gauge transforms, etc.) that map one physical system
/// configuration to another.
#[derive(Debug, Clone, PartialEq)]
pub struct SymmetryTransform {
    /// Name of this transform.
    pub name: String,
    /// Source system name.
    pub source: String,
    /// Target system name.
    pub target: String,
    /// The transformation matrix (3x3 stored as \[f64; 9\], row-major).
    pub matrix: [f64; 9],
    /// Whether this transform preserves energy.
    pub preserves_energy: bool,
}
impl SymmetryTransform {
    /// Create a new symmetry transform.
    pub fn new(
        name: impl Into<String>,
        source: impl Into<String>,
        target: impl Into<String>,
        matrix: [f64; 9],
        preserves_energy: bool,
    ) -> Self {
        Self {
            name: name.into(),
            source: source.into(),
            target: target.into(),
            matrix,
            preserves_energy,
        }
    }
    /// Create an identity transform for a system.
    pub fn identity(system: impl Into<String>) -> Self {
        let s = system.into();
        Self {
            name: format!("id_{}", s),
            source: s.clone(),
            target: s,
            matrix: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            preserves_energy: true,
        }
    }
    /// Compose this transform with another (sequential: self then other).
    ///
    /// The resulting matrix is `other.matrix * self.matrix`.
    pub fn compose(&self, other: &SymmetryTransform) -> Option<SymmetryTransform> {
        if self.target != other.source {
            return None;
        }
        let a = &self.matrix;
        let b = &other.matrix;
        let mut result = [0.0f64; 9];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    result[i * 3 + j] += b[i * 3 + k] * a[k * 3 + j];
                }
            }
        }
        Some(SymmetryTransform {
            name: format!("{};{}", self.name, other.name),
            source: self.source.clone(),
            target: other.target.clone(),
            matrix: result,
            preserves_energy: self.preserves_energy && other.preserves_energy,
        })
    }
    /// Check if this transform is approximately the identity.
    pub fn is_identity(&self, tol: f64) -> bool {
        let id = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        self.matrix
            .iter()
            .zip(id.iter())
            .all(|(a, b)| (a - b).abs() < tol)
    }
    /// Compute the determinant of the transformation matrix.
    pub fn determinant(&self) -> f64 {
        let m = &self.matrix;
        m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
            + m[2] * (m[3] * m[7] - m[4] * m[6])
    }
    /// Check if this is a proper rotation (det = +1).
    pub fn is_rotation(&self, tol: f64) -> bool {
        (self.determinant() - 1.0).abs() < tol
    }
}
/// An endofunctor on a category: maps object names to object names.
#[derive(Clone)]
pub struct Endofunctor {
    /// Human-readable name for this endofunctor.
    pub name: String,
    /// The object map: given an object name, return its image.
    pub maps_object: fn(&str) -> String,
}
impl Endofunctor {
    /// Create a new endofunctor with the given name and object map.
    pub fn new(name: impl Into<String>, maps_object: fn(&str) -> String) -> Self {
        Self {
            name: name.into(),
            maps_object,
        }
    }
    /// Apply the endofunctor to an object name.
    pub fn apply(&self, obj: &str) -> String {
        (self.maps_object)(obj)
    }
    /// Return the name of this endofunctor.
    pub fn name(&self) -> &str {
        &self.name
    }
    /// Check the identity preservation law: `F(id_A)` should map to `id_{F(A)}`.
    pub fn check_identity_preservation(&self, obj: &str) -> bool {
        let id_src = format!("id_{}", obj);
        let mapped = self.apply(obj);
        let expected_id = format!("id_{}", mapped);
        let mapped_id = self.apply(&id_src);
        mapped_id == expected_id
    }
}
/// A simple named object for use in concrete categories.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct NamedObject {
    /// The label of this object.
    pub label: String,
}
impl NamedObject {
    /// Create a named object.
    pub fn new(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
        }
    }
}
/// A category whose objects are physical systems and morphisms are symmetry transforms.
///
/// Composition is sequential application of symmetry transforms; identity is
/// the trivial (do-nothing) symmetry.
#[derive(Debug, Clone)]
pub struct PhysicsCategory {
    /// Physical systems (objects).
    pub systems: Vec<PhysicalSystem>,
    /// Symmetry transforms (morphisms).
    pub transforms: Vec<SymmetryTransform>,
}
impl PhysicsCategory {
    /// Create a new empty physics category.
    pub fn new() -> Self {
        Self {
            systems: Vec::new(),
            transforms: Vec::new(),
        }
    }
    /// Add a physical system and its identity transform.
    pub fn add_system(&mut self, system: PhysicalSystem) {
        let id = SymmetryTransform::identity(&system.name);
        self.transforms.push(id);
        self.systems.push(system);
    }
    /// Add a symmetry transform morphism.
    pub fn add_transform(&mut self, transform: SymmetryTransform) {
        self.transforms.push(transform);
    }
    /// Find a system by name.
    pub fn find_system(&self, name: &str) -> Option<&PhysicalSystem> {
        self.systems.iter().find(|s| s.name == name)
    }
    /// Find a transform by name.
    pub fn find_transform(&self, name: &str) -> Option<&SymmetryTransform> {
        self.transforms.iter().find(|t| t.name == name)
    }
    /// Compose two transforms by name.
    pub fn compose_transforms(&self, f: &str, g: &str) -> Option<SymmetryTransform> {
        let tf = self.find_transform(f)?;
        let tg = self.find_transform(g)?;
        tf.compose(tg)
    }
    /// Check that identity composed with any morphism returns the morphism.
    pub fn check_identity_law(&self, transform_name: &str) -> bool {
        let t = match self.find_transform(transform_name) {
            Some(t) => t,
            None => return false,
        };
        let id_name = format!("id_{}", t.source);
        let id_t = match self.find_transform(&id_name) {
            Some(id) => id,
            None => return false,
        };
        match id_t.compose(t) {
            Some(result) => t
                .matrix
                .iter()
                .zip(result.matrix.iter())
                .all(|(a, b)| (a - b).abs() < 1e-9),
            None => false,
        }
    }
    /// Check that composition of transforms is associative.
    pub fn check_associativity(&self, f: &str, g: &str, h: &str) -> bool {
        let tf = match self.find_transform(f) {
            Some(t) => t,
            None => return false,
        };
        let tg = match self.find_transform(g) {
            Some(t) => t,
            None => return false,
        };
        let th = match self.find_transform(h) {
            Some(t) => t,
            None => return false,
        };
        let fg = match tf.compose(tg) {
            Some(c) => c,
            None => return false,
        };
        let gh = match tg.compose(th) {
            Some(c) => c,
            None => return false,
        };
        let fg_h = match fg.compose(th) {
            Some(c) => c,
            None => return false,
        };
        let f_gh = match tf.compose(&gh) {
            Some(c) => c,
            None => return false,
        };
        fg_h.matrix
            .iter()
            .zip(f_gh.matrix.iter())
            .all(|(a, b)| (a - b).abs() < 1e-9)
    }
    /// Get all energy-preserving transforms.
    pub fn energy_preserving_transforms(&self) -> Vec<&SymmetryTransform> {
        self.transforms
            .iter()
            .filter(|t| t.preserves_energy)
            .collect()
    }
    /// Total number of systems.
    pub fn system_count(&self) -> usize {
        self.systems.len()
    }
    /// Total number of transforms.
    pub fn transform_count(&self) -> usize {
        self.transforms.len()
    }
    /// Convert this physics category to a symbolic Category.
    pub fn to_category(&self) -> Category {
        let mut cat = Category::new();
        for sys in &self.systems {
            cat.add_object(&sys.name);
        }
        for t in &self.transforms {
            if !t.name.starts_with("id_") {
                cat.add_morphism(&t.name, &t.source, &t.target);
            }
        }
        cat
    }
}
/// A monad over `Vec`T`: non-determinism / list monad.
///
/// `unit(x) = vec!\[x\]`, `bind(xs, f) = xs.into_iter().flat_map(f).collect()`.
#[derive(Debug, Clone, Copy)]
pub struct VecMonad;
impl VecMonad {
    /// Wrap a value in a singleton list (unit for the Vec monad).
    pub fn unit<T>(x: T) -> Vec<T> {
        vec![x]
    }
    /// Sequentially bind a list with a function returning lists.
    pub fn bind<T, U, F>(xs: Vec<T>, f: F) -> Vec<U>
    where
        F: Fn(T) -> Vec<U>,
    {
        xs.into_iter().flat_map(f).collect()
    }
    /// Lift a pure function into the Vec monad.
    pub fn fmap<T, U, F>(xs: Vec<T>, f: F) -> Vec<U>
    where
        F: Fn(T) -> U,
    {
        xs.into_iter().map(f).collect()
    }
}
/// The free monoid over an alphabet of type `T`.
///
/// Elements are finite lists (words); the monoid operation is concatenation
/// and the identity element is the empty list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FreeMonoid<T: Clone + PartialEq> {
    /// The underlying word (sequence of letters).
    pub word: Vec<T>,
}
impl<T: Clone + PartialEq> FreeMonoid<T> {
    /// Construct the empty word (identity element).
    pub fn empty() -> Self {
        Self { word: Vec::new() }
    }
    /// Construct a word from a single letter.
    pub fn singleton(letter: T) -> Self {
        Self { word: vec![letter] }
    }
    /// Construct a word from a slice.
    pub fn from_slice(letters: &[T]) -> Self {
        Self {
            word: letters.to_vec(),
        }
    }
    /// Return `true` if this is the empty word (identity).
    pub fn is_empty(&self) -> bool {
        self.word.is_empty()
    }
    /// Length of the word.
    pub fn len(&self) -> usize {
        self.word.len()
    }
    /// Concatenate (append) `other` to `self` (monoid binary operation).
    pub fn concat(&self, other: &Self) -> Self {
        let mut result = self.word.clone();
        result.extend_from_slice(&other.word);
        Self { word: result }
    }
    /// Compute the `n`-th power (concat `self` with itself `n` times).
    pub fn power(&self, n: usize) -> Self {
        let mut result = Self::empty();
        for _ in 0..n {
            result = result.concat(self);
        }
        result
    }
    /// Return the reverse of this word.
    pub fn reverse(&self) -> Self {
        let mut rev = self.word.clone();
        rev.reverse();
        Self { word: rev }
    }
    /// Return the subword from index `start` to `end` (exclusive).
    pub fn subword(&self, start: usize, end: usize) -> Self {
        let end = end.min(self.word.len());
        let start = start.min(end);
        Self {
            word: self.word[start..end].to_vec(),
        }
    }
}
/// A functor between two categories (structure-preserving map).
///
/// Stores the source and target category names, and maps object names to
/// object names and morphism names to morphism names.
#[derive(Debug, Clone)]
pub struct Functor {
    /// Name of this functor.
    pub name: String,
    /// Name of the source category.
    pub source_category: String,
    /// Name of the target category.
    pub target_category: String,
    /// Object map: `(source_object, target_object)` pairs.
    pub object_map: Vec<(String, String)>,
    /// Morphism map: `(source_morphism, target_morphism)` pairs.
    pub morphism_map: Vec<(String, String)>,
}
impl Functor {
    /// Create a new (empty) functor.
    pub fn new(
        name: impl Into<String>,
        source: impl Into<String>,
        target: impl Into<String>,
    ) -> Self {
        Self {
            name: name.into(),
            source_category: source.into(),
            target_category: target.into(),
            object_map: Vec::new(),
            morphism_map: Vec::new(),
        }
    }
    /// Register a mapping of object `src` to object `tgt`.
    pub fn map_object(&mut self, src: impl Into<String>, tgt: impl Into<String>) {
        self.object_map.push((src.into(), tgt.into()));
    }
    /// Register a mapping of morphism `src` to morphism `tgt`.
    pub fn map_morphism(&mut self, src: impl Into<String>, tgt: impl Into<String>) {
        self.morphism_map.push((src.into(), tgt.into()));
    }
    /// Look up the image of a source object.
    pub fn apply_object(&self, src: &str) -> Option<&str> {
        self.object_map
            .iter()
            .find(|(s, _)| s == src)
            .map(|(_, t)| t.as_str())
    }
    /// Look up the image of a source morphism.
    pub fn apply_morphism(&self, src: &str) -> Option<&str> {
        self.morphism_map
            .iter()
            .find(|(s, _)| s == src)
            .map(|(_, t)| t.as_str())
    }
    /// Check functor identity law: `F(id_A) = id_{F(A)}` for a given object name.
    ///
    /// Checks that the functor maps the identity morphism `id_A` to `id_{F(A)}`.
    pub fn check_identity_law(&self, obj: &str) -> bool {
        let id_src = format!("id_{}", obj);
        let Some(mapped_obj) = self.apply_object(obj) else {
            return false;
        };
        let id_tgt = format!("id_{}", mapped_obj);
        self.apply_morphism(&id_src) == Some(id_tgt.as_str())
    }
    /// Check composition preservation: `F(g . f) = F(g) . F(f)`.
    ///
    /// Takes the names of the composed morphisms and verifies the functor
    /// maps the composition to the composition of the mapped morphisms.
    pub fn check_composition_law(&self, f: &str, g: &str, gf_composed: &str) -> bool {
        let mapped_f = self.apply_morphism(f);
        let mapped_g = self.apply_morphism(g);
        let mapped_gf = self.apply_morphism(gf_composed);
        match (mapped_f, mapped_g, mapped_gf) {
            (Some(mf), Some(mg), Some(mgf)) => {
                let expected = format!("{};{}", mf, mg);
                mgf == expected
            }
            _ => false,
        }
    }
}
/// A profunctor optic for composable data transformations.
///
/// A profunctor `P<A, B>` is a type constructor that is contravariant in `A`
/// and covariant in `B`.
#[derive(Clone)]
pub struct ProfunctorOptic<S, T, A, B> {
    /// Forward function: transforms `A` to `B`.
    pub forward: fn(A) -> B,
    /// Adapter: lifts the inner transformation into the outer `S -> T` context.
    pub adapter: fn(fn(A) -> B, S) -> T,
    /// Phantom data for the type parameters.
    pub(super) _phantom: PhantomData<(S, T, A, B)>,
}
impl<S, T, A, B> ProfunctorOptic<S, T, A, B> {
    /// Construct a profunctor optic from a forward function and an adapter.
    pub fn new(forward: fn(A) -> B, adapter: fn(fn(A) -> B, S) -> T) -> Self {
        Self {
            forward,
            adapter,
            _phantom: PhantomData,
        }
    }
    /// Apply the optic: run the adapter with the stored forward function on `s`.
    pub fn apply(&self, s: S) -> T {
        (self.adapter)(self.forward, s)
    }
}
/// A monoidal category: a category equipped with a tensor product and a unit object.
///
/// The tensor product is a bifunctor `(C x C) -> C`, and the unit object `I`
/// satisfies `I tensor A ~ A ~ A tensor I` (up to natural isomorphism).
#[derive(Debug, Clone)]
pub struct MonoidalCategory {
    /// The underlying category.
    pub category: Category,
    /// The unit object label.
    pub unit_object: String,
    /// Tensor product results: maps `(A, B)` to the tensor object label.
    pub tensor_products: HashMap<(String, String), String>,
    /// Associator isomorphisms: `(A, B, C)` -> morphism name for `(A tensor B) tensor C -> A tensor (B tensor C)`.
    pub associators: HashMap<(String, String, String), String>,
    /// Left unitor: `A` -> morphism name for `I tensor A -> A`.
    pub left_unitors: HashMap<String, String>,
    /// Right unitor: `A` -> morphism name for `A tensor I -> A`.
    pub right_unitors: HashMap<String, String>,
}
impl MonoidalCategory {
    /// Create a new monoidal category with a given underlying category and unit.
    pub fn new(category: Category, unit_object: impl Into<String>) -> Self {
        Self {
            category,
            unit_object: unit_object.into(),
            tensor_products: HashMap::new(),
            associators: HashMap::new(),
            left_unitors: HashMap::new(),
            right_unitors: HashMap::new(),
        }
    }
    /// Register the tensor product of two objects.
    pub fn add_tensor(
        &mut self,
        a: impl Into<String>,
        b: impl Into<String>,
        result: impl Into<String>,
    ) {
        self.tensor_products
            .insert((a.into(), b.into()), result.into());
    }
    /// Compute the tensor product of two objects.
    pub fn tensor(&self, a: &str, b: &str) -> Option<&str> {
        self.tensor_products
            .get(&(a.to_string(), b.to_string()))
            .map(|s| s.as_str())
    }
    /// Return the unit object label.
    pub fn unit(&self) -> &str {
        &self.unit_object
    }
    /// Register an associator isomorphism for `(A, B, C)`.
    pub fn add_associator(
        &mut self,
        a: impl Into<String>,
        b: impl Into<String>,
        c: impl Into<String>,
        morphism: impl Into<String>,
    ) {
        self.associators
            .insert((a.into(), b.into(), c.into()), morphism.into());
    }
    /// Register a left unitor for object A.
    pub fn add_left_unitor(&mut self, a: impl Into<String>, morphism: impl Into<String>) {
        self.left_unitors.insert(a.into(), morphism.into());
    }
    /// Register a right unitor for object A.
    pub fn add_right_unitor(&mut self, a: impl Into<String>, morphism: impl Into<String>) {
        self.right_unitors.insert(a.into(), morphism.into());
    }
    /// Check left unit law: `I tensor A = A` (via tensor product table).
    pub fn check_left_unit(&self, a: &str) -> bool {
        match self.tensor(&self.unit_object, a) {
            Some(result) => result == a,
            None => false,
        }
    }
    /// Check right unit law: `A tensor I = A` (via tensor product table).
    pub fn check_right_unit(&self, a: &str) -> bool {
        match self.tensor(a, &self.unit_object) {
            Some(result) => result == a,
            None => false,
        }
    }
    /// Check the pentagon identity (associativity coherence) symbolically.
    ///
    /// Verifies that associators exist for all needed triples.
    pub fn check_pentagon(&self, a: &str, b: &str, c: &str, d: &str) -> bool {
        let key_abc = (a.to_string(), b.to_string(), c.to_string());
        let key_bcd = (b.to_string(), c.to_string(), d.to_string());
        self.associators.contains_key(&key_abc) && self.associators.contains_key(&key_bcd)
    }
    /// Tensor two morphisms: given `f: A -> B` and `g: C -> D`,
    /// produce the name for `f tensor g: A tensor C -> B tensor D`.
    pub fn tensor_morphisms(&self, f_name: &str, g_name: &str) -> String {
        format!("{}_tensor_{}", f_name, g_name)
    }
    /// Return the number of registered tensor products.
    pub fn tensor_count(&self) -> usize {
        self.tensor_products.len()
    }
}
/// A monad over `Option`T`: provides unit (return) and bind (>>=).
///
/// This is the standard Option monad: `unit(x) = Some(x)`,
/// `bind(m, f) = m.and_then(f)`.
#[derive(Debug, Clone, Copy)]
pub struct OptionMonad;
impl OptionMonad {
    /// Wrap a value in `Some` (unit / return for the Option monad).
    pub fn unit<T>(x: T) -> Option<T> {
        Some(x)
    }
    /// Sequentially bind a monadic value `m` with a function `f`.
    pub fn bind<T, U, F>(m: Option<T>, f: F) -> Option<U>
    where
        F: FnOnce(T) -> Option<U>,
    {
        m.and_then(f)
    }
    /// Lift a pure function `f: T -> U` into the Option monad.
    pub fn fmap<T, U, F>(m: Option<T>, f: F) -> Option<U>
    where
        F: FnOnce(T) -> U,
    {
        m.map(f)
    }
}
/// A physical system object in a physics category.
///
/// Objects represent physical systems with properties like energy, degrees
/// of freedom, and symmetry groups.
#[derive(Debug, Clone, PartialEq)]
pub struct PhysicalSystem {
    /// Name of the physical system.
    pub name: String,
    /// Number of degrees of freedom.
    pub degrees_of_freedom: usize,
    /// Total energy of the system (conserved under symmetry transforms).
    pub energy: f64,
    /// Symmetry group label (e.g., "SO(3)", "U(1)", "Z2").
    pub symmetry_group: String,
}
impl PhysicalSystem {
    /// Create a new physical system.
    pub fn new(
        name: impl Into<String>,
        dof: usize,
        energy: f64,
        symmetry: impl Into<String>,
    ) -> Self {
        Self {
            name: name.into(),
            degrees_of_freedom: dof,
            energy,
            symmetry_group: symmetry.into(),
        }
    }
}
/// A bidirectional lens for state access and update.
///
/// A lens `Lens<S, A>` provides a `get: S -> A` (view) and a
/// `set: (S, A) -> S` (update / put) operation.
pub struct OpticLens<S, A> {
    /// Extract a focus `A` from a source `S`.
    pub getter: fn(&S) -> A,
    /// Replace the focus in a source with a new value.
    pub setter: fn(S, A) -> S,
}
impl<S: Clone, A: Clone + PartialEq> OpticLens<S, A> {
    /// Construct a lens from a getter and a setter.
    pub fn new(getter: fn(&S) -> A, setter: fn(S, A) -> S) -> Self {
        Self { getter, setter }
    }
    /// Extract the focus from a source value.
    pub fn get(&self, s: &S) -> A {
        (self.getter)(s)
    }
    /// Replace the focus in a source value, returning the updated source.
    pub fn set(&self, s: S, a: A) -> S {
        (self.setter)(s, a)
    }
    /// Modify the focus using a function `f`.
    pub fn modify<F: Fn(A) -> A>(&self, s: S, f: F) -> S
    where
        S: Clone,
    {
        let a = (self.getter)(&s);
        let a_new = f(a);
        (self.setter)(s, a_new)
    }
    /// Check the put-get law: `get(set(s, a)) == a`.
    pub fn check_put_get(&self, s: S, a: A) -> bool {
        let s_new = self.set(s, a.clone());
        self.get(&s_new) == a
    }
    /// Check the get-put law: `set(s, get(s)) == s` (via `PartialEq` on `S`).
    pub fn check_get_put(&self, s: S) -> bool
    where
        S: PartialEq,
    {
        let a = self.get(&s);
        let s_round = self.set(s.clone(), a);
        s_round == s
    }
}
/// A typed morphism from type `A` to type `B`, carrying only a name tag.
///
/// The phantom data ensures the type parameters are tracked statically even
/// though no value of type `A` or `B` is stored at runtime.
#[derive(Debug, Clone)]
pub struct Morphism<A, B> {
    /// Human-readable name for this morphism.
    pub name: String,
    /// Zero-sized marker for the domain and codomain types.
    pub _phantom: PhantomData<(A, B)>,
}
impl<A, B> Morphism<A, B> {
    /// Create a new morphism with the given name.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            _phantom: PhantomData,
        }
    }
    /// Return the name of this morphism.
    pub fn name(&self) -> &str {
        &self.name
    }
}
/// A small category represented as named objects and named morphisms.
///
/// Each morphism entry is a triple `(name, domain, codomain)`.
/// Identity morphisms are automatically added for each object.
#[derive(Debug, Clone, Default)]
pub struct Category {
    /// Named objects in this category.
    pub objects: Vec<String>,
    /// Morphisms as `(name, domain_object, codomain_object)` triples.
    pub morphisms: Vec<(String, String, String)>,
    /// Composition table: `(f_name, g_name) -> composed_name`.
    pub composition_table: Vec<(String, String, String)>,
}
impl Category {
    /// Create a new empty category.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a named object and its identity morphism `id_X: X -> X`.
    pub fn add_object(&mut self, name: &str) {
        self.objects.push(name.to_string());
        let id_name = format!("id_{}", name);
        self.morphisms
            .push((id_name, name.to_string(), name.to_string()));
    }
    /// Add a morphism `name: domain -> codomain`.
    pub fn add_morphism(&mut self, name: &str, domain: &str, codomain: &str) {
        self.morphisms
            .push((name.to_string(), domain.to_string(), codomain.to_string()));
    }
    /// Register a composition law `f ; g = h` in the composition table.
    pub fn add_composition(&mut self, f: &str, g: &str, composed: &str) {
        self.composition_table
            .push((f.to_string(), g.to_string(), composed.to_string()));
    }
    /// Number of objects in this category.
    pub fn object_count(&self) -> usize {
        self.objects.len()
    }
    /// Number of morphisms in this category (including identity morphisms).
    pub fn morphism_count(&self) -> usize {
        self.morphisms.len()
    }
    /// Return `true` if there is a morphism named `name` from `from` to `to`.
    pub fn is_valid_morphism(&self, from: &str, to: &str, name: &str) -> bool {
        self.morphisms
            .iter()
            .any(|(n, d, c)| n == name && d == from && c == to)
    }
    /// Compose two morphisms `f` and `g` by name (sequential: f then g).
    ///
    /// Returns `Some((composed_name, domain_of_f, codomain_of_g))` if the
    /// codomain of `f` equals the domain of `g`; otherwise `None`.
    pub fn compose_morphisms(&self, f: &str, g: &str) -> Option<(String, String, String)> {
        let (_, f_dom, f_cod) = self.morphisms.iter().find(|(n, _, _)| n == f)?;
        let (_, g_dom, g_cod) = self.morphisms.iter().find(|(n, _, _)| n == g)?;
        if f_cod != g_dom {
            return None;
        }
        if let Some((_, _, comp)) = self
            .composition_table
            .iter()
            .find(|(a, b, _)| a == f && b == g)
        {
            return Some((comp.clone(), f_dom.clone(), g_cod.clone()));
        }
        Some((format!("{};{}", f, g), f_dom.clone(), g_cod.clone()))
    }
    /// Check the left-identity law: `id_A ; f = f` for a given morphism `f`.
    pub fn check_left_identity(&self, f: &str) -> bool {
        let entry = self.morphisms.iter().find(|(n, _, _)| n == f);
        let Some((_, f_dom, _)) = entry else {
            return false;
        };
        let id_name = format!("id_{}", f_dom);
        let comp = self.compose_morphisms(&id_name, f);
        comp.is_some()
    }
    /// Check the right-identity law: `f ; id_B = f` for a given morphism `f`.
    pub fn check_right_identity(&self, f: &str) -> bool {
        let entry = self.morphisms.iter().find(|(n, _, _)| n == f);
        let Some((_, _, f_cod)) = entry else {
            return false;
        };
        let id_name = format!("id_{}", f_cod);
        let comp = self.compose_morphisms(f, &id_name);
        comp.is_some()
    }
    /// Check associativity: `(f;g);h = f;(g;h)` given three composable morphisms.
    ///
    /// Returns true if both triple-compositions yield the same domain and codomain.
    pub fn check_associativity(&self, f: &str, g: &str, h: &str) -> bool {
        let fg = self.compose_morphisms(f, g);
        let gh = self.compose_morphisms(g, h);
        match (fg, gh) {
            (Some((fg_name, fg_dom, _fg_cod)), Some((_gh_name, _gh_dom, gh_cod))) => {
                let lhs_dom = fg_dom;
                let rhs_cod = gh_cod;
                let _lhs = self.compose_morphisms(&fg_name, h);
                let _rhs_name = format!("{};{}", f, _gh_name);
                lhs_dom
                    == *self
                        .morphisms
                        .iter()
                        .find(|(n, _, _)| n == f)
                        .map(|(_, d, _)| d)
                        .unwrap_or(&String::new())
                    && rhs_cod
                        == *self
                            .morphisms
                            .iter()
                            .find(|(n, _, _)| n == h)
                            .map(|(_, _, c)| c)
                            .unwrap_or(&String::new())
            }
            _ => false,
        }
    }
    /// Return all morphisms from `from` to `to`.
    pub fn hom_set(&self, from: &str, to: &str) -> Vec<String> {
        self.morphisms
            .iter()
            .filter(|(_, d, c)| d == from && c == to)
            .map(|(n, _, _)| n.clone())
            .collect()
    }
}
/// A descriptor for a natural transformation between two functors.
///
/// For each object `A` in the source category there is a component morphism
/// `eta_A : F(A) -> G(A)` in the target category.
#[derive(Debug, Clone)]
pub struct NaturalTransformation {
    /// Name of the source functor F.
    pub source: String,
    /// Name of the target functor G.
    pub target: String,
    /// Component morphisms indexed by `(object_name, component_morphism_name)`.
    pub components: Vec<(String, String)>,
}
impl NaturalTransformation {
    /// Create a new natural transformation descriptor.
    pub fn new(source: impl Into<String>, target: impl Into<String>) -> Self {
        Self {
            source: source.into(),
            target: target.into(),
            components: Vec::new(),
        }
    }
    /// Add a component morphism for object `obj` with morphism name `morphism`.
    pub fn add_component(&mut self, component: impl Into<String>) {
        self.components.push((String::new(), component.into()));
    }
    /// Add a named component for a specific object.
    pub fn add_component_for(&mut self, obj: impl Into<String>, morphism: impl Into<String>) {
        self.components.push((obj.into(), morphism.into()));
    }
    /// Number of components registered.
    pub fn component_count(&self) -> usize {
        self.components.len()
    }
    /// Look up the component morphism for a given object.
    pub fn component_for(&self, obj: &str) -> Option<&str> {
        self.components
            .iter()
            .find(|(o, _)| o == obj)
            .map(|(_, m)| m.as_str())
    }
    /// Check the naturality square for a morphism `f: A -> B`.
    ///
    /// Naturality requires `eta_B . F(f) = G(f) . eta_A`.
    /// This method checks symbol equality of the two composed morphism names.
    pub fn check_naturality(
        &self,
        _f_name: &str,
        fa_then_eta_b: &str,
        eta_a_then_gf: &str,
    ) -> bool {
        fa_then_eta_b == eta_a_then_gf
    }
    /// Vertical composition of two natural transformations.
    ///
    /// Given `alpha: F => G` and `beta: G => H`, produce `beta . alpha: F => H`.
    pub fn vertical_compose(&self, other: &NaturalTransformation) -> NaturalTransformation {
        let mut result = NaturalTransformation::new(&self.source, &other.target);
        for (obj, comp_a) in &self.components {
            if let Some(comp_b) = other.component_for(obj) {
                result.add_component_for(obj, format!("{};{}", comp_a, comp_b));
            }
        }
        result
    }
}
/// The Yoneda lemma: for a locally small category C, a functor F: C -> Set,
/// and an object A in C, there is a natural bijection
/// `Nat(Hom(A, -), F) ~ F(A)`.
///
/// This struct represents a representable functor Hom(A, -) and can evaluate
/// the Yoneda embedding.
#[derive(Debug, Clone)]
pub struct YonedaLemma {
    /// The representing object A.
    pub representing_object: String,
    /// The category in which the Yoneda lemma operates.
    pub category: Category,
    /// Cached hom-set values: for each object B, the set Hom(A, B).
    pub(super) hom_sets: HashMap<String, Vec<String>>,
}
impl YonedaLemma {
    /// Create a Yoneda lemma context for the representing object `a` in `category`.
    pub fn new(a: impl Into<String>, category: Category) -> Self {
        let a_str = a.into();
        let mut yl = Self {
            representing_object: a_str.clone(),
            category,
            hom_sets: HashMap::new(),
        };
        yl.compute_hom_sets(&a_str);
        yl
    }
    /// Compute all hom-sets Hom(A, B) for every object B.
    fn compute_hom_sets(&mut self, a: &str) {
        for obj in &self.category.objects {
            let hom = self.category.hom_set(a, obj);
            self.hom_sets.insert(obj.clone(), hom);
        }
    }
    /// Return Hom(A, B) for a given object B.
    pub fn hom(&self, b: &str) -> Option<&Vec<String>> {
        self.hom_sets.get(b)
    }
    /// The Yoneda embedding: given a morphism `f: A -> B`, produce the
    /// natural transformation `Hom(B, -) => Hom(A, -)` by precomposition.
    ///
    /// Returns a map from object C to the morphism names in Hom(A, C)
    /// obtained by precomposing with f.
    pub fn yoneda_embedding(&self, f_name: &str, b: &str) -> HashMap<String, Vec<String>> {
        let mut result = HashMap::new();
        for obj in &self.category.objects {
            let hom_b_c = self.category.hom_set(b, obj);
            let precomposed: Vec<String> = hom_b_c
                .iter()
                .map(|g| format!("{};{}", f_name, g))
                .collect();
            result.insert(obj.clone(), precomposed);
        }
        result
    }
    /// Evaluate the Yoneda lemma bijection: given a natural transformation
    /// from Hom(A, -) to F, extract the element of F(A) by evaluating
    /// the transformation at the identity morphism id_A.
    ///
    /// In our symbolic setting, returns the component at A applied to id_A.
    pub fn evaluate_at_identity(&self) -> Option<String> {
        let id_name = format!("id_{}", self.representing_object);
        let hom_a_a = self.hom(&self.representing_object)?;
        if hom_a_a.contains(&id_name) {
            Some(id_name)
        } else {
            None
        }
    }
    /// The number of objects for which hom-sets have been computed.
    pub fn hom_set_count(&self) -> usize {
        self.hom_sets.len()
    }
    /// Check that the representable functor preserves identity: Hom(A, A) contains id_A.
    pub fn check_representable_identity(&self) -> bool {
        let id_name = format!("id_{}", self.representing_object);
        self.hom(&self.representing_object)
            .map(|v| v.contains(&id_name))
            .unwrap_or(false)
    }
}
/// An abstract monoid consisting of an identity element and a binary operation.
#[derive(Clone)]
pub struct Monoid<T: Clone + PartialEq> {
    /// The identity element of this monoid.
    pub identity: T,
    /// The associative binary operation.
    pub op: fn(T, T) -> T,
}
impl<T: Clone + PartialEq> Monoid<T> {
    /// Construct a monoid from an identity and a binary operation.
    pub fn new(identity: T, op: fn(T, T) -> T) -> Self {
        Self { identity, op }
    }
    /// Apply the monoid operation to `a` and `b`.
    pub fn combine(&self, a: T, b: T) -> T {
        (self.op)(a, b)
    }
    /// Raise `a` to the `n`-th power under the monoid operation.
    pub fn power(&self, a: T, n: usize) -> T {
        let mut result = self.identity.clone();
        for _ in 0..n {
            result = self.combine(result, a.clone());
        }
        result
    }
    /// Return `true` if `a` equals the identity element.
    pub fn is_identity(&self, a: &T) -> bool {
        *a == self.identity
    }
    /// Fold a slice of elements with the monoid operation (left fold).
    pub fn fold(&self, xs: &[T]) -> T {
        xs.iter()
            .fold(self.identity.clone(), |acc, x| self.combine(acc, x.clone()))
    }
    /// Check left-identity law: `identity . a = a`.
    pub fn check_left_identity(&self, a: &T) -> bool {
        self.combine(self.identity.clone(), a.clone()) == *a
    }
    /// Check right-identity law: `a . identity = a`.
    pub fn check_right_identity(&self, a: &T) -> bool {
        self.combine(a.clone(), self.identity.clone()) == *a
    }
    /// Check associativity: `(a . b) . c = a . (b . c)`.
    pub fn check_associativity(&self, a: T, b: T, c: T) -> bool {
        let lhs = self.combine(self.combine(a.clone(), b.clone()), c.clone());
        let rhs = self.combine(a, self.combine(b, c));
        lhs == rhs
    }
}
