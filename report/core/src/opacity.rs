//! Tabulated Rosseland and Planck mean opacities of the Z = 0 gas (port of the table part of
//! `opacity.jl`). The table itself (`kaptable.bin`) is computed by the Julia code from the
//! monochromatic continuum opacity (Thomson, H Rayleigh, H I and He II bound-free n ≤ 12 / 8,
//! He I bound-free, free-free, H⁻ bound-free and free-free; no lines, no molecules) and shipped as
//! data; it is not recomputed here.
//!
//! Layout of `kaptable.bin` (little endian): magic `"QSKAP1\0\0"`, `u32 nT`, `u32 nR`,
//! `f64 lT0, dlT, lR0, dlR`, then `f32 log10 κ_R[nT][nR]` (R index fastest), then
//! `f32 log10 κ_P` in the same layout. Axes: log T = lT0 + i dlT, log R = lR0 + j dlR,
//! R = ρ / T₆³.

/// The opacity table on a uniform (log T, log R) grid, R = ρ/T₆³, T₆ = T/10⁶ K.
#[derive(Clone, Debug)]
pub struct KapTable {
    /// Number of log T nodes.
    pub nt: usize,
    /// Number of log R nodes.
    pub nr: usize,
    /// First log T node.
    pub lt0: f64,
    /// log T spacing.
    pub dlt: f64,
    /// First log R node.
    pub lr0: f64,
    /// log R spacing.
    pub dlr: f64,
    lk: Vec<f32>,
    lkp: Vec<f32>,
}

/// Catmull–Rom cubic weights (interpolating, C¹ across cell edges).
#[inline]
fn crw(t: f64) -> [f64; 4] {
    let t2 = t * t;
    let t3 = t2 * t;
    [
        (-t3 + 2.0 * t2 - t) / 2.0,
        (3.0 * t3 - 5.0 * t2 + 2.0) / 2.0,
        (-3.0 * t3 + 4.0 * t2 + t) / 2.0,
        (t3 - t2) / 2.0,
    ]
}

impl KapTable {
    /// Parse `kaptable.bin` (layout in the module documentation).
    pub fn from_bytes(b: &[u8]) -> Result<KapTable, String> {
        if b.len() < 48 || &b[0..6] != b"QSKAP1" {
            return Err("not a QSKAP1 opacity table".into());
        }
        let u32_at = |o: usize| u32::from_le_bytes(b[o..o + 4].try_into().unwrap()) as usize;
        let f64_at = |o: usize| f64::from_le_bytes(b[o..o + 8].try_into().unwrap());
        let (nt, nr) = (u32_at(8), u32_at(12));
        let (lt0, dlt, lr0, dlr) = (f64_at(16), f64_at(24), f64_at(32), f64_at(40));
        let n = nt * nr;
        if b.len() != 48 + 8 * n {
            return Err(format!("table size {} does not match {}x{}", b.len(), nt, nr));
        }
        let read = |o: usize| -> Vec<f32> {
            b[o..o + 4 * n]
                .chunks_exact(4)
                .map(|c| f32::from_le_bytes(c.try_into().unwrap()))
                .collect()
        };
        Ok(KapTable { nt, nr, lt0, dlt, lr0, dlr, lk: read(48), lkp: read(48 + 4 * n) })
    }

    /// Last log T node.
    pub fn lt_max(&self) -> f64 {
        self.lt0 + (self.nt - 1) as f64 * self.dlt
    }
    /// Last log R node.
    pub fn lr_max(&self) -> f64 {
        self.lr0 + (self.nr - 1) as f64 * self.dlr
    }

    /// Tabulated log10 κ_R at node (i, j) (0-based; R index fastest).
    pub fn node_log_kappa_r(&self, i: usize, j: usize) -> f64 {
        self.lk[i * self.nr + j] as f64
    }

    /// C¹ bicubic (Catmull–Rom) interpolation of a log10 table in (log T, log R); outside the table
    /// the arguments are clamped to its edges.
    fn interp(&self, m: &[f32], rho: f64, t: f64) -> f64 {
        let lt = t.log10();
        let lr = rho.log10() - 3.0 * (lt - 6.0);
        let (nt, nr) = (self.nt, self.nr);
        let x = ((lt - self.lt0) / self.dlt).clamp(0.0, (nt - 1) as f64 - 1e-9);
        let y = ((lr - self.lr0) / self.dlr).clamp(0.0, (nr - 1) as f64 - 1e-9);
        let i = x.floor() as isize;
        let j = y.floor() as isize;
        let wx = crw(x - i as f64);
        let wy = crw(y - j as f64);
        let mut v = 0.0;
        for a in 0..4 {
            let ii = (i - 1 + a as isize).clamp(0, nt as isize - 1) as usize;
            let row = &m[ii * nr..(ii + 1) * nr];
            for b in 0..4 {
                let jj = (j - 1 + b as isize).clamp(0, nr as isize - 1) as usize;
                v += wx[a] * wy[b] * row[jj] as f64;
            }
        }
        10f64.powf(v)
    }

    /// Equation (eq-rosseland) of the text: the Rosseland mean
    /// `1/κ_R = ∫ (κ_ν^abs + κ_ν^sc)⁻¹ ∂B_ν/∂T dν / ∫ ∂B_ν/∂T dν` (cm² g⁻¹), interpolated from the table
    /// by C¹ bicubic Catmull–Rom interpolation of log10 κ_R in (log T, log R), R = ρ/T₆³;
    /// arguments outside the table are clamped to its edge (see [`KapTable::in_table`]).
    pub fn kappa_r(&self, rho: f64, t: f64) -> f64 {
        self.interp(&self.lk, rho, t)
    }

    /// Planck mean absorption opacity `κ_P = ∫ κ_ν^abs B_ν dν / ∫ B_ν dν` (cm² g⁻¹; scattering
    /// excluded), the companion of Equation (eq-rosseland) of the text, interpolated like
    /// [`KapTable::kappa_r`].
    pub fn kappa_p(&self, rho: f64, t: f64) -> f64 {
        self.interp(&self.lkp, rho, t)
    }

    /// Whether (ρ, T) lies inside the tabulated (log T, log R) rectangle (no clamping applied).
    pub fn in_table(&self, rho: f64, t: f64) -> bool {
        let lt = t.log10();
        let lr = rho.log10() - 3.0 * (lt - 6.0);
        self.lt0 <= lt && lt <= self.lt_max() && self.lr0 <= lr && lr <= self.lr_max()
    }
}
