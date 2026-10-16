//! Shared helpers of the validation suite: the shipped opacity table and the Julia reference values.
#![allow(dead_code)]
use qs_core::opacity::KapTable;
use serde_json::Value;
use std::sync::OnceLock;

pub fn table() -> &'static KapTable {
    static T: OnceLock<KapTable> = OnceLock::new();
    T.get_or_init(|| {
        let b = std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/data/kaptable.bin")).expect("data/kaptable.bin");
        KapTable::from_bytes(&b).expect("valid table")
    })
}

pub fn reference() -> &'static Value {
    static R: OnceLock<Value> = OnceLock::new();
    R.get_or_init(|| {
        let s = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/data/reference.json")).expect("data/reference.json");
        serde_json::from_str(&s).expect("valid json")
    })
}

/// Field `k` of a JSON object as f64 (`null` → NaN).
pub fn f(v: &Value, k: &str) -> f64 {
    v[k].as_f64().unwrap_or(f64::NAN)
}

pub fn rel(a: f64, b: f64) -> f64 {
    if a == b {
        0.0
    } else {
        (a - b).abs() / b.abs()
    }
}

/// Track and report the worst relative deviation of a quantity.
pub struct Worst {
    pub name: &'static str,
    pub tol: f64,
    pub max: f64,
}
impl Worst {
    pub fn new(name: &'static str, tol: f64) -> Self {
        Worst { name, tol, max: 0.0 }
    }
    pub fn check(&mut self, got: f64, want: f64, ctx: &str) {
        let r = rel(got, want);
        assert!(r <= self.tol, "{} at {}: got {:e}, Julia {:e}, rel {:.3e} > {:e}", self.name, ctx, got, want, r, self.tol);
        self.max = self.max.max(r);
    }
}
impl Drop for Worst {
    fn drop(&mut self) {
        if !std::thread::panicking() {
            println!("  {:<12} max rel deviation {:.2e} (tolerance {:.0e})", self.name, self.max, self.tol);
        }
    }
}
