# 2DWD2_PhaseRetrieval

- The Matlab script `WD2_Demo.m` implements the 2DWD2: 2D Wigner Distribution Deconvolution algorithm for phase retrieval presented in:

&emsp; *Fast and Accurate 2D Phase Retrieval from Incomplete, Windowed Fourier Measurements* (link incoming...)\
&emsp; &emsp; Maya Hamka, Brandon Hutchison, Donald Liveoak

- The image files `image_magFourier.jpg` and `image_phaseRodenburg.png` are used by `WD2_Demo.m` to construct test signals.

- The folder `Numerics` contains the data which generated the timing plot below in `.mat` format.
<br/><br/>
## Background and Mathematical Model:
  
<p align="center">
  <img src="https://drive.google.com/uc?export=view&id=16z7quYEOxXqLlZKSVPjQ2hfzhmFzTwzW" />
</p>

  
  The algorithm recovers an unknown, 2D signal $X \in \mathbb{C}^{d_1 \times d_2}$ given phaseless (magnitude only) measurements of the form:
  
  $$Y_{\boldsymbol{\ell}} \ = \  \left| \mathcal{F}_D \left( X \circ S_{\ell_1, \ell_2} M \right) \right|^2 + N_{\boldsymbol{\ell}}$$
  
  where

| Symbol | Description |
| ----------- | ----------- |
| $\mathcal{F}_D$ | 2D Fourier Transform |
| $\circ$ | Componentwise Product |
| $S_{\ell_1, \ell_2}$ | 2D Shift Operator |
| $M \in \mathbb{C}^{d_1 \times d_2}$ | Known Mask, nonzero only in small $\delta_1 \times \delta_2$ block |
| $\vert \cdot \vert^2$ | Componentwise Squared Magnitude |
| $N_{\boldsymbol{\ell}} \in \mathbb{R}^{d_1 \times d_2}$ | Matrix of Arbitrary Noise |

We assume that we have access to measurements corresponding to all possible $d_1 \times d_2$ shifts of the *windowed* mask $M$, but only a $K_1 \times K_2$ subset of the possible Fourier frequencies.  Measurements of this form are found in **coherent diffraction imaging** (see above figure), where coherent radiation is diffracted within an unknown specimen, and the intensity of the resulting pattern is recorded as light and dark spots on a plane.  Particularly, by collecting multiple, localized diffration patterns, our model corresponds to the measurements which occur in **ptychography**.
<br/><br/>
## Available Test Signals:

Running `WD2_Demo.m` gives a choice of two test signals to recover:

- Entering \'a\' will generate a $228 \times 204$ random test signal drawn from the standard normal distribution. 
- Entering \'b\' will generate a $253 \times 200$ deterministic test signal by setting `image_magFourier.jpg` to the magnitude and `image_phaseRodenburg.png` to the phase:

<p align="center">
  <img src="https://drive.google.com/uc?export=view&id=1CT7FuqbOBERxvLDa1rL0V7xn8BNlibdk" />
</p>

- In either case, both the mask and noise are drawn from the standard normal distribution with the stregth of the additive noise governed by the Signal-to-Noise Ratio (SNR)


## Numerical Results:
<p align="center">
  <img src="https://drive.google.com/uc?export=view&id=1nQGbaH5B_RuxjUhkTQJvRIx6bZsTS_-z" />
</p>

- **Fast:** Empirically, execution timing (*i*) follows the known $\mathcal{O}(d_1 d_2 \log_2 (d_1 d_2))$ cost of the famous **Fast Fourier Transform**.
    - problem sizes ranged from $38 \times 34$ to $228 \times 204$ with a $10 \times 9$ mask window
    - 30 trials per problem size, median time recorded on log-log plot
    - the matrix `timeTrials_10x9.mat` contains the data for (*i*) and is found in the `Numerics` folder

- **Accurate:** Plot (*ii*) demonstrates the accuracy of 2DWD2 algorithm when **noise** is present.
    - window size of $\delta_1 \times \delta_2$
    - measurements with $(2\delta_1 -1) \times (2\delta_2 -1)$ Fourier frequencies
    - 30 trials per window size, mean relative error recorded in decibels
