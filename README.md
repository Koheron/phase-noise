# Simulation of laser phase noise

## Requirements

* Python 3.12 or newer
* NumPy, PyQtGraph and PyQt6 (Qt 6)

Install dependencies:

```
python3 -m pip install -r requirements.txt
```

## Run simulation

```
python3 laser.py
```

## Physics

The laser field has constant amplitude and a diffusing phase:

$$E(t) = E_0 e^{i[2\pi\nu_0 t + \phi(t)]}, \qquad d\phi = \sqrt{2D_\phi\,dt}\,\xi, \quad \xi \sim \mathcal{N}(0,1).$$

This gives field coherence $|g^{(1)}(\tau)| = e^{-D_\phi |\tau|}$ and a
Lorentzian optical spectrum with FWHM $\Delta\nu = D_\phi/\pi$.

Two equal-power interferometer arms, separated by delay $\tau$, give the
normalized intensity

$$I(t) = \frac{1 + \cos[\theta + \phi(t) - \phi(t-\tau)]}{2}.$$

The carrier phase is absorbed into the fixed bias $\theta = \pi/4$.
The delay uses FFT interpolation with linear drift removed and restored.
The plot shows the two-sided PSD of intensity fluctuations after mean removal
and Hann windowing, in dB/Hz. Amplitude noise and detector noise are omitted.
