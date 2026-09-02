# Rabi_fitting
Rabi_fitting

MATLAB fitting of Rabi oscillation data from a two-qubit NV⁻ centre quantum computer laboratory (University of Sydney, 2025). Rabi oscillations were measured at each of the two microwave resonance frequencies of the spin system; the fitted parameters were used to calibrate π and π/2 pulse widths for later pulse sequences, ending in an implementation of the Deutsch-Jozsa algorithm.

Model
Each dataset is fitted with a damped sine:

A(τ) = A₀ · sin(2π Ω_R τ + φ) · exp(−τ / T₂ᴿᵃᵇⁱ) + C

where τ is the microwave pulse duration (ns), Ω_R is the Rabi frequency, φ is the phase, T₂ᴿᵃᵇⁱ is the Rabi decay time, and C is a constant offset.

Fitting uses lsqcurvefit. The returned Jacobian is used to estimate the parameter covariance matrix, and 1σ uncertainties are taken from its diagonal. Fitted Rabi frequencies were approximately 6.4 MHz and 8.1 MHz for the first and second resonance respectively, with Rabi decay times of roughly 280 ns and 216 ns.

Files
freq1_uncertainties.m — fit for the first resonance frequency (reads rabi_fixed_freq1.xlsx)
freq2_uncertainties.m — fit for the second resonance frequency (reads rabifreq_2_fixed.xlsx)
rabi_fixed_freq1.xlsx, rabifreq_2_fixed.xlsx — measured Rabi oscillation data: pulse duration τ (ns) in column 1, signal amplitude in column 2

Any two-column data in this format can be fitted by editing the filename and initial guesses b0 in either script.

Requirements
MATLAB with the Optimization Toolbox (lsqcurvefit).
