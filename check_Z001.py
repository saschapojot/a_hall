import numpy as np
from scipy.integrate import quad, dblquad

hbar = 1
m = 1.0
ma = 10.0
vR = 2.0
a0 = 0.5
theta_a = 0.1 * np.pi

# --- Physics Functions ---

def Ek(rho, gamma):
    part1 = vR**2 * rho**2
    val2 = a0 + hbar**2 / (2 * ma) * rho**2 * np.cos(2 * gamma - theta_a)
    Ek2 = part1 + val2**2
    return np.sqrt(Ek2)

def G0(rho, gamma):
    return hbar**2 / (2 * m) * rho**2 + Ek(rho, gamma)

def g0(rho):
    part1 = hbar**2 / (2 * m) * rho**2
    part2 = np.sqrt(vR**2 * rho**2 + a0**2)
    return part1 + part2
def drho_g0(rho):
    leading = vR**2 * rho**2 + a0**2
    part1=hbar**2/m*rho
    part2=vR**2*rho*leading**(-1/2)
    return part1+part2
# --- Numerically Stable Fermi Functions ---

def fermi_dist(val):
    """
    Computes 1 / (exp(x) + 1) safely.
    """
    if val > 100: return 0.0
    if val < -100: return 1.0
    return 1.0 / (np.exp(val) + 1.0)

def fermi_deriv(val):
    """
    Computes exp(x) / (exp(x) + 1)^2 safely.
    Equivalent to: 1 / (4 * cosh^2(x/2))
    """
    if abs(val) > 100:
        return 0.0
    # Using cosh avoids the overflow of exp(val) entirely
    return 0.25 / (np.cosh(val / 2.0)**2)

def fk0(rho, gamma, beta, muF):
    val = beta * (G0(rho, gamma) - muF)
    return fermi_dist(val)

def c0(rho):
    leading = vR**2 * rho**2 + a0**2
    val = 0.5 * a0 * hbar**2 * rho**2 * leading**(-0.5)
    return val

def w00(rho, beta, muF):
    val = beta * (g0(rho) - muF)
    return fermi_dist(val)

def w01(rho, beta, muF):
    val = beta * (g0(rho) - muF)
    return fermi_deriv(val)


# --- Integrands and Asymptotics ---

def Z0010_approx_integrand(rho, beta, muF):
    """Integrand for Z0010 using the exponential approximation (Eq. 209)"""
    leading = vR**2 * rho**2 + a0**2
    # For rho < rho0, g0 < muF, so g0 - muF is negative
    val = rho**3 * leading**(-3/2) * c0(rho) * np.exp(beta * (g0(rho) - muF))
    val *= beta
    return val

def Z0011_approx_integrand(rho, beta, muF):
    """Integrand for Z0011 using the exponential approximation (Eq. 213, 214)"""
    leading = vR**2 * rho**2 + a0**2
    # For rho > rho0, g0 > muF, so muF - g0 is negative
    val = rho**3 * leading**(-3/2) * c0(rho) * np.exp(beta * (muF - g0(rho)))
    val *= beta
    return val

def Z001_exact_integrand(rho, beta, muF):
    """Exact integrand using the exact Fermi derivative w01 (Eq. 205 & 206)"""
    leading = vR**2 * rho**2 + a0**2
    val = beta * rho**3 * leading**(-3/2) * w01(rho, beta, muF) * c0(rho)
    return val

def get_rho0(muF):
    """Calculates rho_0 where g0(rho) = muF (Eq. 203)"""
    term1 = (muF * hbar**2 / m) + vR**2
    term2 = np.sqrt(vR**4 + (2 * muF * hbar**2 / m) * vR**2 + (a0**2 * hbar**4 / m**2))
    rho0_sq = (2 * m**2 / hbar**4) * (term1 - term2)
    return np.sqrt(rho0_sq)

def Z001_asymptotic(rho0):
    """Analytical asymptotic value for both Z0010 and Z0011 (Eq. 211 & 214)"""
    leading = vR**2 * rho0**2 + a0**2
    val = 2*rho0**3 * leading**(-3/2) * c0(rho0) * (1.0 / drho_g0(rho0))
    return val


# --- Comparison Execution ---

muF = 1.5
beta = 100.0  # Large beta (low temperature)

rho_0 = get_rho0(muF)



# Total Z001: Integration from 0 to infinity
num_exact_total, _ = quad(Z001_exact_integrand, 0, np.inf, args=(beta, muF))
Z001_asymp_val=Z001_asymptotic(rho_0)
print(f"--- Parameters ---")
print(f"muF:  {muF}")
print(f"beta: {beta}")
print(f"rho0: {rho_0:.6f}\n")
print(num_exact_total)
print(Z001_asymp_val)
print(num_exact_total/Z001_asymp_val)







