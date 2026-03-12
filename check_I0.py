import numpy as np
from scipy.integrate import quad, dblquad

hbar = 1
m = 1.0
ma = 20.0
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

# --- Integrands ---

def I00_integrand_complex(rho, gamma, beta, muF):
    # Calculate the complex value
    val = rho**3 * fk0(rho, gamma, beta, muF) * np.sin(2 * gamma) / Ek(rho, gamma)**3
    val *= -1j * vR**2 / (2 * ma) * np.sin(theta_a)
    return val

def I00_asymp_integrand_complex(rho, beta, muF):
    leading = vR**2 * rho**2 + a0**2

    # Term 1
    val1 = rho**5 * leading**(-2.5) * w00(rho, beta, muF)
    val1 *= (1/ma**2) * 1j * 0.75 * np.pi * vR**2 * np.sin(theta_a)**2 * a0 * hbar**2

    # Term 2
    val2 = rho**3 * leading**(-1.5) * w01(rho, beta, muF) * c0(rho)
    val2 *= (1/ma**2) * 1j * 0.5 * np.pi * vR**2 * beta * np.sin(theta_a)**2

    return val1 + val2

def I00_asymp_part1_integrand(rho, beta, muF):
    leading = vR**2 * rho**2 + a0**2
    # Term 1
    val1 = rho**5 * leading**(-2.5) * w00(rho, beta, muF)
    val1 *= (1/ma**2) * 1j * 0.75 * np.pi * vR**2 * np.sin(theta_a)**2 * a0 * hbar**2
    return val1

def get_rho0(muF):
    """Calculates rho_0 where g0(rho) = muF (Eq. 203)"""
    term1 = (muF * hbar**2 / m) + vR**2
    term2 = np.sqrt(vR**4 + (2 * muF * hbar**2 / m) * vR**2 + (a0**2 * hbar**4 / m**2))
    rho0_sq = (2 * m**2 / hbar**4) * (term1 - term2)
    return np.sqrt(rho0_sq)

def drho_g0(rho):
    leading = vR**2 * rho**2 + a0**2
    part1=hbar**2/m*rho
    part2=vR**2*rho*leading**(-1/2)
    return part1+part2

def I00_asymp_part2(rho, beta, muF):
    rho0=get_rho0(muF)
    leading = vR**2 * rho0**2 + a0**2
    val=rho0**3*leading**(-3/2)*c0(rho0)*1/drho_g0(rho0)
    val*=1/ma**2*1j*1/2*np.pi*vR**2*np.sin(theta_a)**2
    return val

# --- Comparison Execution ---

# Parameters for the test
beta = 500.0  # Large beta for the asymptotic limit (T -> 0)
muF = 2.0      # muF must be > a0 (0.5) to have a Fermi surface

print(f"Comparing I00 at beta={beta}, muF={muF}\n")

# 1. Exact Numerical Integration (2D)
# dblquad expects func(y, x). Let x = rho, y = gamma
# We integrate the imaginary part since the result is purely imaginary.
def exact_imag(gamma, rho):
    return np.imag(I00_integrand_complex(rho, gamma, beta, muF))

# Integrate rho from 0 to inf, gamma from 0 to 2*pi
exact_res_imag, exact_err = dblquad(exact_imag, 0, np.inf, lambda r: 0, lambda r: 2*np.pi)
exact_val = 1j * exact_res_imag

# 2. Asymptotic Calculation
# Part 1: Numerical integration over rho
def asymp_part1_imag(rho):
    return np.imag(I00_asymp_part1_integrand(rho, beta, muF))

part1_res_imag, part1_err = quad(asymp_part1_imag, 0, np.inf)
part1_val = 1j * part1_res_imag

# Part 2: Analytic boundary term evaluated at rho0
# (The `rho` argument is unused in the function body, so we pass 0)
part2_val = I00_asymp_part2(0, beta, muF)

asymp_total = part1_val + part2_val

# --- Output Results ---
print(f"Exact I00 (Numerical dblquad) : {exact_val}")
print(f"Asymptotic I00 (Total)        : {asymp_total}")
print(f"  -> Part 1 (Integral term)   : {part1_val}")
print(f"  -> Part 2 (Analytic term)   : {part2_val}")
print("-" * 50)
print(f"ratio         : {abs(exact_val/asymp_total):.5e}")