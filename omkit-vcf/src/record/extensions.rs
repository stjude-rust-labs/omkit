//! A type-erased map for passing state between transformers.

use std::any::Any;
use std::any::TypeId;
use std::collections::HashMap;
use std::fmt;

use crate::Error;
use crate::Result;

/// A type-erased extension map for inter-transformer state.
///
/// Each transformer inserts its own state struct by type, and downstream
/// transformers retrieve it. This decouples transformers from a shared
/// state definition.
///
/// # Examples
///
/// ```
/// use omkit_vcf::record::extensions::Extension;
/// use omkit_vcf::record::extensions::Extensions;
///
/// #[derive(Clone, Debug, PartialEq)]
/// struct MyState {
///     value: u32,
/// }
///
/// impl Extension for MyState {}
///
/// let mut ext = Extensions::default();
/// ext.insert(MyState { value: 42 });
///
/// assert_eq!(ext.get::<MyState>().unwrap().value, 42);
/// ```
#[derive(Clone, Default)]
pub struct Extensions {
    /// The inner map of type-erased values.
    map: HashMap<TypeId, Box<dyn CloneableAny>>,
}

/// A marker trait for types that can be stored in [`Extensions`].
///
/// Implementors declare which transformer produces them so that error
/// messages can guide the user when a required extension is missing.
/// The default [`producer()`](Extension::producer) implementation uses
/// the type name.
pub trait Extension: Any + Clone + Send + Sync {
    /// The name of the transformer that produces this extension.
    ///
    /// The default implementation returns `std::any::type_name::<Self>()`.
    fn producer() -> &'static str {
        std::any::type_name::<Self>()
    }
}

/// An object-safe trait combining [`Any`] with [`Clone`].
trait CloneableAny: Any + Send + Sync {
    /// Clones this value into a new boxed trait object.
    fn clone_box(&self) -> Box<dyn CloneableAny>;

    /// Returns this value as an `&dyn Any` for downcasting.
    fn as_any(&self) -> &dyn Any;

    /// Returns the type name for error messages.
    fn type_name(&self) -> &'static str;
}

impl<T: Any + Clone + Send + Sync> CloneableAny for T {
    fn clone_box(&self) -> Box<dyn CloneableAny> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn type_name(&self) -> &'static str {
        std::any::type_name::<T>()
    }
}

impl Clone for Box<dyn CloneableAny> {
    fn clone(&self) -> Self {
        (**self).clone_box()
    }
}

impl Extensions {
    /// Inserts a value into the extensions, keyed by its type.
    ///
    /// If a value of the same type was already present, it is returned.
    pub fn insert<T: Extension>(&mut self, val: T) -> Option<T> {
        self.map
            .insert(TypeId::of::<T>(), Box::new(val))
            .and_then(|b| (*b).as_any().downcast_ref::<T>().cloned())
    }

    /// Returns a reference to the value of type `T`, if present.
    pub fn get<T: Extension>(&self) -> Option<&T> {
        self.map
            .get(&TypeId::of::<T>())
            .and_then(|b| (**b).as_any().downcast_ref())
    }

    /// Returns a reference to the value of type `T`, or an error if it
    /// has not been inserted by a prior transformer.
    ///
    /// The `consumer` parameter names the transformer calling this
    /// method so the error message can identify both ends of the
    /// dependency.
    pub fn require<T: Extension>(&self, consumer: &'static str) -> Result<&T> {
        self.get::<T>().ok_or_else(|| Error::MissingExtension {
            consumer,
            producer: T::producer(),
        })
    }
}

impl fmt::Debug for Extensions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let type_names: Vec<&'static str> = self.map.values().map(|v| (**v).type_name()).collect();
        f.debug_struct("Extensions")
            .field("types", &type_names)
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone, Debug, PartialEq)]
    struct Alpha {
        value: u32,
    }

    impl Extension for Alpha {}

    #[derive(Clone, Debug, PartialEq)]
    struct Beta {
        flag: bool,
    }

    impl Extension for Beta {}

    #[test]
    fn insert_and_get() {
        let mut ext = Extensions::default();
        ext.insert(Alpha { value: 42 });

        let alpha = ext.get::<Alpha>().unwrap();
        assert_eq!(alpha.value, 42);
    }

    #[test]
    fn get_returns_none_when_absent() {
        let ext = Extensions::default();
        assert!(ext.get::<Alpha>().is_none());
    }

    #[test]
    fn require_returns_value_when_present() {
        let mut ext = Extensions::default();
        ext.insert(Alpha { value: 7 });

        let alpha = ext.require::<Alpha>("test").unwrap();
        assert_eq!(alpha.value, 7);
    }

    #[test]
    fn require_returns_error_when_absent() {
        let ext = Extensions::default();
        let err = ext.require::<Alpha>("TestConsumer").unwrap_err();

        let msg = err.to_string();
        assert!(
            msg.contains("TestConsumer"),
            "error should name the consumer: {msg}"
        );
        assert!(
            msg.contains("Alpha"),
            "error should name the producer: {msg}"
        );
    }

    #[test]
    fn multiple_types_coexist() {
        let mut ext = Extensions::default();
        ext.insert(Alpha { value: 1 });
        ext.insert(Beta { flag: true });

        assert_eq!(ext.get::<Alpha>().unwrap().value, 1);
        assert!(ext.get::<Beta>().unwrap().flag);
    }

    #[test]
    fn insert_returns_none_on_first_insertion() {
        let mut ext = Extensions::default();
        assert!(ext.insert(Alpha { value: 1 }).is_none());
    }

    #[test]
    fn insert_returns_previous_value() {
        let mut ext = Extensions::default();
        ext.insert(Alpha { value: 1 });

        let previous = ext.insert(Alpha { value: 2 });
        assert_eq!(previous.unwrap().value, 1);
        assert_eq!(ext.get::<Alpha>().unwrap().value, 2);
    }

    #[test]
    fn clone_produces_independent_copy() {
        let mut ext = Extensions::default();
        ext.insert(Alpha { value: 10 });

        let mut cloned = ext.clone();
        cloned.insert(Alpha { value: 20 });

        assert_eq!(ext.get::<Alpha>().unwrap().value, 10);
        assert_eq!(cloned.get::<Alpha>().unwrap().value, 20);
    }

    #[test]
    fn debug_lists_type_names() {
        let mut ext = Extensions::default();
        ext.insert(Alpha { value: 1 });

        let debug = format!("{ext:?}");
        assert!(
            debug.contains("Alpha"),
            "debug output should contain type name: {debug}"
        );
    }

    #[test]
    fn default_is_empty() {
        let ext = Extensions::default();
        assert!(ext.get::<Alpha>().is_none());
    }
}
