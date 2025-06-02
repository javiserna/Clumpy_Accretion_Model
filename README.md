# Clumpy Accretion Toy Model

This repository contains a simple Python implementation of a *clumpy accretion* model designed to simulate the time-variable mass accretion rate onto a Classical T Tauri Star (CTTS). The model assumes that accretion occurs in discrete bursts associated with the stochastic infall of clumps of gas along magnetospheric field lines.

The resulting light curve is constructed by superimposing Gaussian-shaped bursts, each corresponding to an individual clump.

![Clumpy Accretion Illustration](clumpy_accretion_cartoon.png)

> *Illustration: Magnetospheric accretion from a protoplanetary disk with bursts modeled as Gaussian events, where the amplitude and width depend on the clump mass and time.*

This toy model was inspired by the ongoing work of **Tao Ji** and **Greg Herczeg**, who are developing observational and theoretical tools to understand burst-dominated accretion in TW Hya. The approach here allows others to generate synthetic light curves with physically motivated burst properties, suitable for testing fitting strategies or comparing with TESS/ground-based data.

## How to Use

Make sure you have Python 3 installed along with NumPy and Matplotlib.

```bash
pip install numpy matplotlib

Then, you can generate a synthetic clumpy accretion light curve with:

from clumpy_toy_model import generate_clumpy_lightcurve

time, Mdot = generate_clumpy_lightcurve(
    total_days=25,
    mean_bursts_per_day=2.5,
    mass_range=(1e-11, 1e-9),
    duration_range=(0.1, 0.6),
    seed=42
)

This returns:
	•	time: Array of time values (in days).
	•	Mdot: Mass accretion rate at each time step (in solar masses per year).

You can also plot the light curve:

import matplotlib.pyplot as plt

plt.plot(time, Mdot)
plt.xlabel("Time (days)")
plt.ylabel("Mass Accretion Rate (M$_\odot$/yr)")
plt.title("Synthetic Clumpy Accretion Light Curve")
plt.show()

Credits

This model is for exploratory and illustrative purposes only. It was inspired by the theoretical and observational framework developed by Tao Ji and Greg Herczeg.
