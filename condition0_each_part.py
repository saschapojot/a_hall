import mpmath as mp
import numpy as np
from datetime import datetime

# ---------------------------------------------------------
# 1. Setup High Precision
# ---------------------------------------------------------
mp.dps = 50  # 50 decimal digits of precision
mp.pretty = True

# ---------------------------------------------------------
# 2. Define Symbols & Parameters
# ---------------------------------------------------------
hbar_val = mp.mpf('1.0')
a0_val = mp.mpf('1.1')
m_val = mp.mpf('1.2')
vR_val = mp.mpf('2.0')
e_val = mp.mpf('-1.0')
ma_val = mp.mpf('200.0')
beta_val = mp.mpf('1000.0')
theta_a_val = 0.1 * mp.pi

# ---------------------------------------------------------
# 3. Physics Calculations (Pre-computation)
# ---------------------------------------------------------

# Pre-calculate rho10 (Fermi momentum location)
rho10_val = mp.sqrt(m_val ** 2 * vR_val ** 4 - a0_val ** 2 * hbar_val ** 4) / (mp.fabs(vR_val) * hbar_val ** 2)

# Calculate Chemical Potential muF
muF = (hbar_val ** 2 / (2 * m_val)) * (rho10_val ** 2) - (m_val * vR_val ** 2) / (hbar_val ** 2) - \
      (a0_val * hbar_val ** 4) / (2 * ma_val * m_val * vR_val ** 2) * (rho10_val ** 2)


# ---------------------------------------------------------
# 4. Define Component Integrands for Band 1
# ---------------------------------------------------------

def get_common_terms(rho, gamma):
    """
    Calculates common terms used in all I1 components.
    Returns: (fk1, Ek, kx, ky, cos_alpha_k)
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
    arg1 = 0.5 * beta_val * (eps_k1 - muF)
    fk1 = 0.5 * (1 - mp.tanh(arg1))

    # Angle term for I12
    cos_alpha_k = M_k / Ek

    return fk1, Ek, kx, ky, cos_alpha_k


# --- I10 Integrand: Proportional to kx*ky / Ek^3 ---
# Ref Eq (143): I10 = -i * (vR^2/ma) * sin(theta_a) * Integral [ fk1 * (kx*ky / Ek^3) ]
def integrand_I10(rho, gamma):
    fk1, Ek, kx, ky, _ = get_common_terms(rho, gamma)

    # Prefactor from Eq (143) (excluding 'i' for numerical integration of the magnitude)
    # We will compute the coefficient separately
    term = (kx * ky) / (Ek ** 3)
    return fk1 * term * rho  # Jacobian rho


# --- I11 Integrand: Proportional to (ky^2 - kx^2) / Ek^3 ---
# Ref Eq (144): I11 = i * (vR^2 / 2ma) * cos(theta_a) * Integral [ fk1 * (ky^2 - kx^2)/Ek^3 ]
def integrand_I11(rho, gamma):
    fk1, Ek, kx, ky, _ = get_common_terms(rho, gamma)

    term = (ky ** 2 - kx ** 2) / (Ek ** 3)
    return fk1 * term * rho


# --- I12 Integrand: Proportional to cos(alpha_k) / Ek^2 ---
# Ref Eq (145): I12 = i * (vR^2 / 2hbar^2) * Integral [ fk1 * cos(alpha_k)/Ek^2 ]
def integrand_I12(rho, gamma):
    fk1, Ek, _, _, _ = get_common_terms(rho, gamma)

    term = a0_val / (Ek ** 3)
    return fk1 * term * rho


# ---------------------------------------------------------
# 5. Execution
# ---------------------------------------------------------
print("Starting Component Integration for Band 1 (mpmath)...")
print(f"Beta: {beta_val}")
print(f"Integration range: rho [0, {rho10_val}, inf]")

start_t = datetime.now()

# 1. Compute I10 (Raw Integral part)
print("\nComputing I10 raw integral (kx*ky/Ek^3)...")
raw_I10 = mp.quad(integrand_I10, [0, rho10_val, mp.inf], [0, 2 * mp.pi], method='tanh-sinh')
# Apply prefactor: -i * (vR^2/ma) * sin(theta_a)
# Note: The result in the document includes 'i'. We print the imaginary magnitude.
pref_I10 = -1 * (vR_val ** 2 / ma_val) * mp.sin(theta_a_val)
val_I10 = pref_I10 * raw_I10
print(f"Raw Integral: {raw_I10}")
print(f"Value I10 (coeff of i): {val_I10}")

# 2. Compute I11 (Raw Integral part)
print("\nComputing I11 raw integral ((ky^2 - kx^2)/Ek^3)...")
raw_I11 = mp.quad(integrand_I11, [0, rho10_val, mp.inf], [0, 2 * mp.pi], method='tanh-sinh')
# Apply prefactor: i * (vR^2 / 2ma) * cos(theta_a)
pref_I11 = (vR_val ** 2 / (2 * ma_val)) * mp.cos(theta_a_val)
val_I11 = pref_I11 * raw_I11
print(f"Raw Integral: {raw_I11}")
print(f"Value I11 (coeff of i): {val_I11}")

# 3. Compute I12 (Raw Integral part)
print("\nComputing I12 raw integral (a0/Ek^3)...")
raw_I12 = mp.quad(integrand_I12, [0, rho10_val, mp.inf], [0, 2 * mp.pi], method='tanh-sinh')
# Apply prefactor: i * (vR^2 / 2hbar^2)
pref_I12 = (vR_val ** 2 / (2 * hbar_val ** 2))
val_I12 = pref_I12 * raw_I12
print(f"Raw Integral: {raw_I12}")
print(f"Value I12 (coeff of i): {val_I12}")
sigma_Hqco_numerical=e_val**2*hbar_val/(4*mp.pi**2)*val_I12
# 4. Total I1 (Sum of components)
# I1 = I10 + I11 + I12
# All terms have an 'i' factor in the document definition.
total_I1_mag = val_I10 + val_I11 + val_I12

end_t = datetime.now()
rho11_n0_val_np = float(a0_val / (2 * vR_val ** 3) * (m_val ** 2 * vR_val ** 4 + a0_val ** 2 * hbar_val ** 4) / \
                        mp.sqrt(m_val ** 2 * vR_val ** 4 - a0_val ** 2 * hbar_val ** 4) * mp.sign(vR_val))

asymp_val_np = np.sqrt(2) * float(e_val) ** 2 * np.log(2) / (4 * np.pi * float(beta_val)) * \
               float(a0_val) ** (1.5) * float(m_val) ** 2 * float(ma_val) ** 0.5 * float(rho10_val) ** (-2) * \
               (1 / (float(m_val) ** 3 * float(vR_val) ** 2) * float(rho10_val) -
                3 * float(hbar_val) ** 4 / (float(ma_val) * float(m_val) ** 5 * float(vR_val) ** 4) * float(
                           rho10_val) ** 2 * rho11_n0_val_np +
                1 / (float(ma_val) * float(m_val) ** 3 * float(vR_val) ** 2) * rho11_n0_val_np)

print("-" * 40)
print("FINAL COMPONENT BREAKDOWN (Coefficients of 'i'):")
print("-" * 40)
print(f"I10 (kx*ky term):      {val_I10}")
print(f"I11 (ky^2-kx^2 term):  {val_I11}")
print(f"I12 (cos alpha term):  {val_I12}")
print("-" * 40)
print(f"Total I1 Magnitude:    {total_I1_mag}")
print(f"sigma_Hqco_numerical={sigma_Hqco_numerical}")
print(f"asymp_val_np={asymp_val_np}")
print(f"rel err: {mp.fabs(sigma_Hqco_numerical-asymp_val_np)/mp.fabs(asymp_val_np)}")
print("-" * 40)
print(f"Time elapsed: {end_t - start_t}")