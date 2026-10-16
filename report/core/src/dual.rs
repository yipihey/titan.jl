//! Minimal forward-mode dual numbers (value plus `N` first partial derivatives).
//!
//! The equation of state is written once, generically over [`Num`], and evaluated either in
//! `f64` or in `Dual<N>`. Seeding `ln ρ` and `ln T` as the two partials gives the exact first
//! derivatives ∂/∂ln ρ|_T and ∂/∂ln T|_ρ needed by Equation (eq-thermo) of the text.
use std::ops::{Add, Div, Mul, Neg, Sub};

/// A number that the generic equation-of-state code can be evaluated with.
pub(crate) trait Num:
    Copy
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + Add<f64, Output = Self>
    + Sub<f64, Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
{
    /// True for dual numbers (derivatives carried), false for plain `f64`.
    const IS_DUAL: bool;
    fn cst(x: f64) -> Self;
    fn val(self) -> f64;
    fn exp(self) -> Self;
    fn ln(self) -> Self;
    fn powf(self, p: f64) -> Self;
    fn recip(self) -> Self {
        Self::cst(1.0) / self
    }
}

impl Num for f64 {
    const IS_DUAL: bool = false;
    #[inline]
    fn cst(x: f64) -> Self {
        x
    }
    #[inline]
    fn val(self) -> f64 {
        self
    }
    #[inline]
    fn exp(self) -> Self {
        f64::exp(self)
    }
    #[inline]
    fn ln(self) -> Self {
        f64::ln(self)
    }
    #[inline]
    fn powf(self, p: f64) -> Self {
        f64::powf(self, p)
    }
}

/// Value `v` and partial derivatives `d`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Dual<const N: usize> {
    pub v: f64,
    pub d: [f64; N],
}

impl<const N: usize> Dual<N> {
    #[inline]
    pub fn new(v: f64, d: [f64; N]) -> Self {
        Dual { v, d }
    }
    #[inline]
    fn chain(self, v: f64, dv: f64) -> Self {
        let mut d = self.d;
        for x in d.iter_mut() {
            *x *= dv;
        }
        Dual { v, d }
    }
}

impl<const N: usize> Add for Dual<N> {
    type Output = Self;
    #[inline]
    fn add(self, o: Self) -> Self {
        let mut d = self.d;
        for i in 0..N {
            d[i] += o.d[i];
        }
        Dual { v: self.v + o.v, d }
    }
}
impl<const N: usize> Sub for Dual<N> {
    type Output = Self;
    #[inline]
    fn sub(self, o: Self) -> Self {
        let mut d = self.d;
        for i in 0..N {
            d[i] -= o.d[i];
        }
        Dual { v: self.v - o.v, d }
    }
}
impl<const N: usize> Mul for Dual<N> {
    type Output = Self;
    #[inline]
    fn mul(self, o: Self) -> Self {
        let mut d = [0.0; N];
        for i in 0..N {
            d[i] = self.d[i] * o.v + self.v * o.d[i];
        }
        Dual { v: self.v * o.v, d }
    }
}
impl<const N: usize> Div for Dual<N> {
    type Output = Self;
    #[inline]
    fn div(self, o: Self) -> Self {
        let v = self.v / o.v;
        let mut d = [0.0; N];
        for i in 0..N {
            d[i] = (self.d[i] - v * o.d[i]) / o.v;
        }
        Dual { v, d }
    }
}
impl<const N: usize> Neg for Dual<N> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        let mut d = self.d;
        for x in d.iter_mut() {
            *x = -*x;
        }
        Dual { v: -self.v, d }
    }
}
impl<const N: usize> Add<f64> for Dual<N> {
    type Output = Self;
    #[inline]
    fn add(self, o: f64) -> Self {
        Dual { v: self.v + o, d: self.d }
    }
}
impl<const N: usize> Sub<f64> for Dual<N> {
    type Output = Self;
    #[inline]
    fn sub(self, o: f64) -> Self {
        Dual { v: self.v - o, d: self.d }
    }
}
impl<const N: usize> Mul<f64> for Dual<N> {
    type Output = Self;
    #[inline]
    fn mul(self, o: f64) -> Self {
        let mut d = self.d;
        for x in d.iter_mut() {
            *x *= o;
        }
        Dual { v: self.v * o, d }
    }
}
impl<const N: usize> Div<f64> for Dual<N> {
    type Output = Self;
    #[inline]
    fn div(self, o: f64) -> Self {
        let mut d = self.d;
        for x in d.iter_mut() {
            *x /= o;
        }
        Dual { v: self.v / o, d }
    }
}

impl<const N: usize> Num for Dual<N> {
    const IS_DUAL: bool = true;
    #[inline]
    fn cst(x: f64) -> Self {
        Dual { v: x, d: [0.0; N] }
    }
    #[inline]
    fn val(self) -> f64 {
        self.v
    }
    #[inline]
    fn exp(self) -> Self {
        let e = self.v.exp();
        self.chain(e, e)
    }
    #[inline]
    fn ln(self) -> Self {
        self.chain(self.v.ln(), 1.0 / self.v)
    }
    #[inline]
    fn powf(self, p: f64) -> Self {
        let y = self.v.powf(p);
        self.chain(y, p * y / self.v)
    }
}
