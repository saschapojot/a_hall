import mpmath as mp
import numpy as np
from datetime import datetime

# ---------------------------------------------------------
# 1. Setup High Precision (INCREASED)
# ---------------------------------------------------------
# Increased to 80 to prevent cancellation errors in the
# 1 - tanh(...) operation at high Beta.
mp.dps = 80
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
theta_a_val = 0.1 * mp.pi

# Placeholder (will be set in loop)
beta_val = mp.mpf('1000.0')

# ---------------------------------------------------------
# 3. Physics Calculations
# ---------------------------------------------------------
rho10_val = mp.sqrt(m_val ** 2 * vR_val ** 4 - a0_val ** 2 * hbar_val ** 4) / (mp.fabs(vR_val) * hbar_val ** 2)

muF = (hbar_val ** 2 / (2 * m_val)) * (rho10_val ** 2) - (m_val * vR_val ** 2) / (hbar_val ** 2) - \
      (a0_val * hbar_val ** 4) / (2 * ma_val * m_val * vR_val ** 2) * (rho10_val ** 2)

rho11_n0_val_np = float(a0_val / (2 * vR_val ** 3) * (m_val ** 2 * vR_val ** 4 + a0_val ** 2 * hbar_val ** 4) / \
                        mp.sqrt(m_val ** 2 * vR_val ** 4 - a0_val ** 2 * hbar_val ** 4) * mp.sign(vR_val))

# ---------------------------------------------------------
# 4. Helper Functions
# ---------------------------------------------------------

def stable_fermi(energy_diff, beta):
    """
    Computes 1 / (exp(beta * diff) + 1) numerically stably.
    Equivalent to 0.5 * (1 - tanh(0.5*beta*diff)) but safer.
    """
    x = beta * energy_diff

    # If x is huge positive, f -> 0
    if x > 200:
        return mp.mpf('0.0')
    # If x is huge negative, f -> 1
    elif x < -200:
        return mp.mpf('1.0')
    else:
        # Standard logistic function
        return 1 / (mp.exp(x) + 1)

def get_common_terms(rho, gamma):
    kx = rho * mp.cos(gamma)
    ky = rho * mp.sin(gamma)
    k2 = rho ** 2

    aniso = (kx ** 2 - ky ** 2) * mp.cos(theta_a_val) + 2 * kx * ky * mp.sin(theta_a_val)
    M_k = a0_val + (hbar_val ** 2 / (2 * ma_val)) * aniso

    Ek = mp.sqrt(vR_val ** 2 * k2 + M_k ** 2)
    e_k_base = (hbar_val ** 2 / (2 * m_val)) * k2
    eps_k1 = e_k_base - Ek

    # Use stable fermi function
    # Note: The original code used arg1 = 0.5 * beta * (eps - mu).
    # tanh(0.5*x) is related to logistic(x).
    # 0.5 * (1 - tanh(x/2)) == 1 / (e^x + 1)
    # So we pass (eps_k1 - muF) directly to the logistic function.
    fk1 = stable_fermi(eps_k1 - muF, beta_val)

    return fk1, Ek

def integrand_I12(rho, gamma):
    fk1, Ek = get_common_terms(rho, gamma)
    term = a0_val / (Ek ** 3)
    return fk1 * term * rho

# ---------------------------------------------------------
# 5. Execution Loop
# ---------------------------------------------------------

beta_list = [10,50,100, 500, 1000, 2000, 5000, 10000]

print(f"{'Beta':<10} | {'Sigma Num':<25} | {'Asymp Val':<25} | {'Rel Error':<25}")
print("-" * 90)

start_total = datetime.now()

for b_val in beta_list:
    beta_val = mp.mpf(str(b_val))

    # --- KEY CHANGE HERE ---
    # maxdegree=15 (default is usually 6-8).
    # This allows the integrator to double the grid density more times
    # to resolve the sharp step at high Beta.
    raw_I12 = mp.quad(integrand_I12,
                      [0, rho10_val, mp.inf],
                      [0, 2 * mp.pi],
                      method='tanh-sinh',
                      maxdegree=15)

    pref_I12 = (vR_val ** 2 / (2 * hbar_val ** 2))
    val_I12 = pref_I12 * raw_I12

    sigma_Hqco_numerical = e_val**2 * hbar_val / (4 * mp.pi**2) * val_I12

    # Asymptotic calculation (same as before)
    asymp_val_np = np.sqrt(2) * float(e_val) ** 2 * np.log(2) / (4 * np.pi * float(beta_val)) * \
                   float(a0_val) ** (1.5) * float(m_val) ** 2 * float(ma_val) ** 0.5 * float(rho10_val) ** (-2) * \
                   (1 / (float(m_val) ** 3 * float(vR_val) ** 2) * float(rho10_val) -
                    3 * float(hbar_val) ** 4 / (float(ma_val) * float(m_val) ** 5 * float(vR_val) ** 4) * float(
                               rho10_val) ** 2 * rho11_n0_val_np +
                    1 / (float(ma_val) * float(m_val) ** 3 * float(vR_val) ** 2) * rho11_n0_val_np)

    rel_err = mp.fabs(sigma_Hqco_numerical - asymp_val_np) / mp.fabs(asymp_val_np)

    print(f"{int(b_val):<10} | {float(sigma_Hqco_numerical):<25.5e} | {asymp_val_np:<25.5e} | {float(rel_err):<25.5e}")

end_total = datetime.now()
print("-" * 90)
print(f"Total Time elapsed: {end_total - start_total}")