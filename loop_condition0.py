import mpmath as mp
import numpy as np
from datetime import datetime

# ---------------------------------------------------------
# 1. Setup High Precision
# ---------------------------------------------------------
mp.dps = 50  # 50 decimal digits of precision
mp.pretty = True

# ---------------------------------------------------------
# 2. Define Symbols & Parameters (Fixed)
# ---------------------------------------------------------
hbar_val = mp.mpf('1.0')
a0_val = mp.mpf('1.1')
m_val = mp.mpf('1.2')
vR_val = mp.mpf('2.0')
e_val = mp.mpf('-1.0')
ma_val = mp.mpf('200.0')
theta_a_val = 0.1 * mp.pi

# We will define beta_val dynamically in the loop later
beta_val = mp.mpf('1000.0')

# ---------------------------------------------------------
# 3. Physics Calculations (Pre-computation - Independent of Beta)
# ---------------------------------------------------------

# Pre-calculate rho10 (Fermi momentum location)
rho10_val = mp.sqrt(m_val ** 2 * vR_val ** 4 - a0_val ** 2 * hbar_val ** 4) / (mp.fabs(vR_val) * hbar_val ** 2)

# Calculate Chemical Potential muF
muF = (hbar_val ** 2 / (2 * m_val)) * (rho10_val ** 2) - (m_val * vR_val ** 2) / (hbar_val ** 2) - \
      (a0_val * hbar_val ** 4) / (2 * ma_val * m_val * vR_val ** 2) * (rho10_val ** 2)

# Pre-calculate numpy constants for asymptotic formula
rho11_n0_val_np = float(a0_val / (2 * vR_val ** 3) * (m_val ** 2 * vR_val ** 4 + a0_val ** 2 * hbar_val ** 4) / \
                        mp.sqrt(m_val ** 2 * vR_val ** 4 - a0_val ** 2 * hbar_val ** 4) * mp.sign(vR_val))

# ---------------------------------------------------------
# 4. Define Component Integrands
# ---------------------------------------------------------

def get_common_terms(rho, gamma):
    """
    Calculates common terms used in all I1 components.
    Uses the global 'beta_val' which is updated in the loop.
    """
    # Map parameters
    hbar = hbar_val
    a0 = a0_val
    m = m_val
    vR = vR_val
    ma = ma_val
    theta_a = theta_a_val

    # Coordinates
    kx = rho * mp.cos(gamma)
    ky = rho * mp.sin(gamma)
    k2 = rho ** 2

    # Mass term
    aniso = (kx ** 2 - ky ** 2) * mp.cos(theta_a) + 2 * kx * ky * mp.sin(theta_a)
    M_k = a0 + (hbar ** 2 / (2 * ma)) * aniso

    # Energies
    Ek = mp.sqrt(vR ** 2 * k2 + M_k ** 2)
    e_k_base = (hbar ** 2 / (2 * m)) * k2
    eps_k1 = e_k_base - Ek  # Lower band energy

    # Fermi Function for Lower Band
    # Note: beta_val is accessed globally here
    arg1 = 0.5 * beta_val * (eps_k1 - muF)
    fk1 = 0.5 * (1 - mp.tanh(arg1))

    return fk1, Ek

# --- I12 Integrand: Proportional to cos(alpha_k) / Ek^2 ---
# Ref Eq (145): I12 = i * (vR^2 / 2hbar^2) * Integral [ fk1 * cos(alpha_k)/Ek^2 ]
# Simplified: cos(alpha_k) = M_k / Ek, so term is M_k / Ek^3.
# But original code used: term = a0_val / (Ek ** 3).
# We stick to the original code's definition of integrand_I12.
def integrand_I12(rho, gamma):
    fk1, Ek = get_common_terms(rho, gamma)
    term = a0_val / (Ek ** 3)
    return fk1 * term * rho

# ---------------------------------------------------------
# 5. Execution Loop
# ---------------------------------------------------------

# Define the range of Beta values to test
beta_list = [10,50,100, 500, 1000, 2000, 5000, 10000]

print(f"{'Beta':<10} | {'Sigma Num':<25} | {'Asymp Val':<25} | {'Rel Error':<25}")
print("-" * 90)

start_total = datetime.now()

for b_val in beta_list:
    # Update global beta variable for the integrand function
    beta_val = mp.mpf(str(b_val))

    # --- Compute I12 ---
    # We only compute I12 because sigma_Hqco_numerical in your code
    # only depended on val_I12. I10 and I11 are skipped for speed.

    raw_I12 = mp.quad(integrand_I12, [0, rho10_val, mp.inf], [0, 2 * mp.pi], method='tanh-sinh')

    # Apply prefactor: i * (vR^2 / 2hbar^2)
    pref_I12 = (vR_val ** 2 / (2 * hbar_val ** 2))
    val_I12 = pref_I12 * raw_I12

    # Calculate Numerical Sigma
    sigma_Hqco_numerical = e_val**2 * hbar_val / (4 * mp.pi**2) * val_I12

    # Calculate Asymptotic Value (Numpy float precision as per original code)
    # Note: asymp_val depends on 1/beta
    asymp_val_np = np.sqrt(2) * float(e_val) ** 2 * np.log(2) / (4 * np.pi * float(beta_val)) * \
                   float(a0_val) ** (1.5) * float(m_val) ** 2 * float(ma_val) ** 0.5 * float(rho10_val) ** (-2) * \
                   (1 / (float(m_val) ** 3 * float(vR_val) ** 2) * float(rho10_val) -
                    3 * float(hbar_val) ** 4 / (float(ma_val) * float(m_val) ** 5 * float(vR_val) ** 4) * float(
                               rho10_val) ** 2 * rho11_n0_val_np +
                    1 / (float(ma_val) * float(m_val) ** 3 * float(vR_val) ** 2) * rho11_n0_val_np)

    # Calculate Relative Error
    rel_err = mp.fabs(sigma_Hqco_numerical - asymp_val_np) / mp.fabs(asymp_val_np)

    # Print row
    print(f"{int(b_val):<10} | {float(sigma_Hqco_numerical):<25.5e} | {asymp_val_np:<25.5e} | {float(rel_err):<25.5e}")

end_total = datetime.now()
print("-" * 90)
print(f"Total Time elapsed: {end_total - start_total}")