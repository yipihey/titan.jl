#!/usr/bin/env python3
"""Read a TITAN 1.2 direct-access dump file (*.dmp) and write the final
(or a chosen) archived model to CSV at full double precision.

Record layout (from archive.f / start##.f; gfortran, recl in bytes):
  record 1 : header + problem-specific scalars
  record 2 : jdump, jstep of the latest archived model   (int32, int32)
  record jdump+2 : jdump, jstep (int32 x2), 12 real*8 scalars
      (cmax dmax smax timen dtime tfac tetot tee tes tew tel teq),
      7 arrays(mgr)  r0 xm0 d0 u0 egrt0 dvol0 tescr0,
      5 scalars      geddl geddr xmass xmext xmtot,
      38 arrays(mgr) (see ARRAYS below)
  ldump = 8*(45*mgr+19) bytes,  mgr = 330  (titan.par)

Staggering (TITAN convention, 1-based Fortran index k):
  interface (zone-edge) quantities: r, xme, u, fr, (xm, drdt, dmdt)
  zone-centred quantities: d, t, eg, pg, er, fedd, chif, xke, xkp, plf, ...
  zone k occupies r(k) <= r <= r(k+1).
  Real zones are k = ngrs .. ngre; interfaces ngrs .. ngre+1.
  Masses are stored divided by 4*pi (dvol = (r(k+1)^3 - r(k)^3)/3).

usage: read_dump.py FILE.dmp ngrs ncell out.csv [jdump]
"""
import sys
import numpy as np

MGR = 330
LDUMP = 8 * (45 * MGR + 19)
SCAL1 = ["cmax", "dmax", "smax", "timen", "dtime", "tfac",
         "tetot", "tee", "tes", "tew", "tel", "teq"]
ARR0 = ["r0", "xm0", "d0", "u0", "egrt0", "dvol0", "tescr0"]
SCAL2 = ["geddl", "geddr", "xmass", "xmext", "xmtot"]
ARRAYS = ["rn", "rmun", "rmup1n", "rmum1n", "avchin", "urel",
          "dvoln", "xmen", "dn", "dsn", "chifn", "drdt",
          "un", "usn", "tn", "asn", "xken", "dmdt",
          "egn", "egsn", "pgn", "xnen", "xkpn", "xnt",
          "ern", "ersn", "egrn", "egrsn", "plfn",
          "etn", "etsn", "egtn", "egtsn", "fedd",
          "frn", "frsn", "egrtn", "egrtsn"]
CAR = 7.5659141e-15  # titan.dat


def read_record(fname, irec):
    with open(fname, "rb") as f:
        f.seek((irec - 1) * LDUMP)
        return f.read(LDUMP)


def read_model(fname, jdump=None):
    if jdump is None:
        jdump, jstep = np.frombuffer(read_record(fname, 2)[:8], "<i4")
    buf = read_record(fname, int(jdump) + 2)
    jd, js = np.frombuffer(buf[:8], "<i4")
    x = np.frombuffer(buf[8:], "<f8")
    out = {"jdump": int(jd), "jstep": int(js)}
    p = 0
    for n in SCAL1:
        out[n] = x[p]; p += 1
    for n in ARR0:
        out[n] = x[p:p + MGR].copy(); p += MGR
    for n in SCAL2:
        out[n] = x[p]; p += 1
    for n in ARRAYS:
        out[n] = x[p:p + MGR].copy(); p += MGR
    return out


def write_csv(m, ngrs, ncell, outname):
    ngre = ngrs + ncell - 1
    k = np.arange(ngrs, ngre + 2)          # interfaces ngrs..ngre+1
    i = k - 1                              # 0-based
    nan = np.nan
    cols = {
        "k": k,
        # interface quantities
        "r": m["rn"][i], "xme_over_4pi": m["xmen"][i], "u": m["un"][i],
        "Fr": m["frn"][i],
        # zone-centred quantities (zone k between r(k) and r(k+1));
        # the last row (k = ngre+1) is the phantom zone and is set to NaN
    }
    zc = {"rho": "dn", "T": "tn", "eg": "egn", "pg": "pgn", "Er": "ern",
          "fedd": "fedd", "chiF": "chifn", "kapE": "xken", "kapP": "xkpn",
          "cs": "asn"}
    for c, a in zc.items():
        v = m[a][i].copy(); v[-1] = nan; cols[c] = v
    rc = 0.5 * (m["rn"][i] + m["rn"][i + 1]); rc[-1] = nan
    cols["r_center"] = rc
    tr = np.sign(cols["Er"]) * np.abs(cols["Er"] / CAR) ** 0.25
    cols["Trad"] = tr
    names = list(cols)
    with open(outname, "w") as f:
        f.write("# TITAN 1.2 dump %s jdump=%d jstep=%d timen=%.16e dtime=%.16e\n"
                % (outname, m["jdump"], m["jstep"], m["timen"], m["dtime"]))
        f.write("# geddl=%.16e geddr=%.16e xmtot=%.16e tetot=%.16e\n"
                % (m["geddl"], m["geddr"], m["xmtot"], m["tetot"]))
        f.write("# staggering: r,xme_over_4pi,u,Fr at interface k; "
                "rho,T,eg,pg,Er,fedd,chiF,kapE,kapP,cs,Trad at zone k "
                "[r(k),r(k+1)]; last row zone values NaN (phantom)\n")
        f.write("# units cgs; eg = specific gas energy (erg/g); Er = rad. energy "
                "density (erg/cm^3); Fr = rad. flux (comoving)\n")
        f.write(",".join(names) + "\n")
        for j in range(len(k)):
            f.write(",".join("%d" % cols[n][j] if n == "k" else "%.16e" % cols[n][j]
                             for n in names) + "\n")


if __name__ == "__main__":
    fn, ngrs, ncell, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    jd = int(sys.argv[5]) if len(sys.argv) > 5 else None
    m = read_model(fn, jd)
    write_csv(m, ngrs, ncell, out)
    print("%s: jdump=%d jstep=%d t=%.6e" % (out, m["jdump"], m["jstep"], m["timen"]))
