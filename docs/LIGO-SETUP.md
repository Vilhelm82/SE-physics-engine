# LIGO Python environment

Prepared 7 September 2026 in `/home/williaml/seated-root/.venv-ligo`, using
Python 3.13.13. The environment is separate from system Python and the existing
research workloads.

## Start working

```bash
cd /home/williaml/seated-root
source .venv-ligo/bin/activate
jupyter lab
```

Select **Python (LIGO)** in Jupyter. This registered kernel uses the environment's
absolute Python path and is also available to other Jupyter installations.
For scripts, activation is optional:

```bash
/home/williaml/seated-root/.venv-ligo/bin/python your_analysis.py
```

## Installed tools

- **GWOSC 0.8.3 / GWpy 4.0.2:** discover and download public strain, read time
  series and quality channels, calculate spectra, filter and plot.
- **PyCBC 2.11.0 / LALSuite 7.26.15:** comparison waveforms, matched filtering,
  detector tools, and LIGO `.gwf` frame-file I/O through LALFrame.
- **Bilby 2.8.2 / Dynesty 3.1.0 / Corner 2.3.0:** likelihoods, Bayesian sampling
  and posterior plots.
- **NumPy, SciPy, Astropy, h5py, pandas, Matplotlib:** numerical operations,
  signal processing, units/GPS time, data storage and plotting.
- **JupyterLab / ipykernel:** interactive notebooks.

The exact resolved versions are in [`requirements-ligo.txt`](../requirements-ligo.txt);
[`requirements-ligo.in`](../requirements-ligo.in) records the direct dependencies.

## Data boundary

This setup targets the public **calibrated detector strain** `h(t)` distributed by
GWOSC. Instrumental readout and control channels before calibration are a distinct
data product; see [GWpy's data-type guide](https://gwpy.readthedocs.io/en/v4.0.1/timeseries/get/).
The downloaded archive is preserved as supplied. Reading a time interval does not
whiten, bandpass, or subtract a waveform from it.

The setup check uses H1 around GW150914: GPS `[1126259446, 1126259478)`, 32 seconds
at 4096 Hz. This is only an installation example; it does not select the event,
frequency range, noise assumptions or statistic for the proposed prediction test.
Data-quality and injection masks need to be applied when defining that analysis.

Download another interval with:

```python
from gwpy.timeseries import TimeSeries

strain = TimeSeries.fetch_open_data(
    "H1", 1126259446, 1126259478, sample_rate=4096, cache=True
)
```

GWOSC may serve a longer archive file for a short interval. Both 4096 Hz and
16384 Hz public products can be requested where available; choose the sampling
rate after specifying the prediction's frequency range.

Local setup artifacts are under `.ligo-data/`, including the original HDF5 archive,
its source URL and SHA-256 in `data-check.json`, `software-check.json`, and an ASD
plot. This directory and `.venv-ligo/` are excluded from Git.

Verification passed: dependency consistency and all main imports; a PyCBC/LAL
waveform and normalized self-match; a Bilby likelihood and short Dynesty run;
HDF5 and GWF round trips; a running Jupyter kernel; and the real H1 download,
131,072 finite strain samples, spectral calculation and plot export.
These checks establish software operation, not scientific validity of a prediction.

## Recreate the environment

From the repository root, using the installed `uv` package manager:

```bash
uv venv --python 3.13.13 --seed .venv-ligo
uv pip install --python .venv-ligo/bin/python -r requirements-ligo.txt
.venv-ligo/bin/python -m ipykernel install --user \
    --name seated-root-ligo --display-name 'Python (LIGO)'
.venv-ligo/bin/python -m pip check
```

The pins describe this Linux/Python 3.13 environment. To deliberately update
versions, re-resolve `requirements-ligo.in` and repeat the functional checks.

Package/source provenance is recorded in
[`EXTERNAL-MATHEMATICS-DEBT.md`](EXTERNAL-MATHEMATICS-DEBT.md), entry SR-015.
