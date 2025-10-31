import numpy as np

def exact_intensity_and_detuning(R_sc, Delta_LS, Gamma, I_sat):
    """
    Analytic solution for intensity (I) and detuning (Δ)
    given scattering rate and light shift for a two-level atom.
    Uses x = Δ / Γ as the dimensionless detuning variable.
    """
    # Compute coefficients
    term = (2 * Delta_LS / Gamma) - (Delta_LS / R_sc)
    disc = term**2 - 1
    if disc < 0:
        raise ValueError("No real physical solution: discriminant < 0.")

    # Two possible roots
    x1 = ( (Delta_LS / R_sc) - (2 * Delta_LS / Gamma) + np.sqrt(disc) ) / 2
    x2 = ( (Delta_LS / R_sc) - (2 * Delta_LS / Gamma) - np.sqrt(disc) ) / 2

    results = []
    for x in [x1, x2]:
        Delta = x * Gamma
        s0 = 8 * x * Delta_LS / Gamma
        I = s0 * I_sat
        if s0 > 0:
            results.append((I, Delta, s0, x))

    return results

# --- Testing ---
Gamma = 2*np.pi*29.13*10**6   # Hz (Yb171 1S0->1P1)
I_sat = 59.97*10**(-3)        # W/cm^2
R_sc = 0.449                  # Hz
Delta_LS = 2601535.348        # Hz

solutions = exact_intensity_and_detuning(R_sc, Delta_LS, Gamma, I_sat)

for i, (I, Delta, s0, x) in enumerate(solutions, 1):
    print(f"Solution {i}:")
    print(f"x = Δ/Γ = {x:.4f}")
    print(f"Δ = {Delta/(2*np.pi):.3f} Hz")
    print(f"s0 = {s0:.4e}")
    print(f"Intensity = {I:.3e} W/cm^2")
    print(f"Wavelength = {((2*np.pi*3*10**8)/((4.7242*10**15)-(Delta)))*10**9:.3e} nm\n")
