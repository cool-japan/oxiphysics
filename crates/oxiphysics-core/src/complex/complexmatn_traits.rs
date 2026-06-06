//! # ComplexMatN - Trait Implementations
//!
//! This module contains trait implementations for `ComplexMatN`.
//!
//! ## Implemented Traits
//!
//! - `IndexMut`
//! - `Index`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{Complex, ComplexMatN};

impl std::ops::IndexMut<(usize, usize)> for ComplexMatN {
    fn index_mut(&mut self, (r, c): (usize, usize)) -> &mut Complex {
        &mut self.data[r * self.n + c]
    }
}

impl std::ops::Index<(usize, usize)> for ComplexMatN {
    type Output = Complex;
    fn index(&self, (r, c): (usize, usize)) -> &Complex {
        &self.data[r * self.n + c]
    }
}
