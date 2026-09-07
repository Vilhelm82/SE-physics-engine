"""Sweeps: parameter grids, incremental saving, and the decoder surface L_d(h, p) as an interpolable object.

The decoder sees only (herald h, unheralded Pauli rate p) per gate, so the map of the whole model factorises:
    figure of merit(config) = L_d( h(config), p(config) )
with L_d measured once on synthetic channels and (h, p) computed per configuration without any decoding.
"""
import itertools, json, time
import numpy as np
from scipy.interpolate import RegularGridInterpolator


def grid(**axes):
    """Cartesian product of named axes -> list of dicts, in a stable order."""
    names = list(axes); return [dict(zip(names, vals)) for vals in itertools.product(*(axes[n] for n in names))]


def progress(i, n, t0, every=25):
    if i % every == 0 or i == n:
        el = time.time() - t0; print(f"    {i}/{n}  {el:.0f}s elapsed, ~{el/max(i,1)*(n-i):.0f}s left", flush=True)


class DecoderSurface:
    """L_d(h, p) from a receipt of synthetic-channel runs; log-log bilinear interpolation with floors so h = 0 or p = 0 are usable."""
    def __init__(self, rows, h_floor=1e-6, p_floor=1e-7):
        self.hs = np.array(sorted({r['h'] for r in rows})); self.ps = np.array(sorted({r['p'] for r in rows}))
        self.h_floor, self.p_floor = h_floor, p_floor
        self.d = sorted({r['d'] for r in rows}); self.tables = {}
        for d in self.d:
            T = np.full((len(self.hs), len(self.ps)), np.nan)
            for r in rows:
                if r['d'] == d: T[np.searchsorted(self.hs, r['h']), np.searchsorted(self.ps, r['p'])] = r['rate']
            self.tables[d] = RegularGridInterpolator((np.log10(self.hs + h_floor), np.log10(self.ps + p_floor)), np.log10(np.maximum(T, 1e-9)),
                                                     bounds_error=False, fill_value=None)
    def __call__(self, d, h, p):
        v = self.tables[d](np.column_stack([np.log10(np.asarray(h, float) + self.h_floor), np.log10(np.asarray(p, float) + self.p_floor)]))
        return 10.0**v
    @classmethod
    def from_receipt(cls, path):
        rec = json.load(open(path)); return cls(rec['records']['surface'])
