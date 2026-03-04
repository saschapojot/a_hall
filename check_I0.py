import numpy as np
from scipy.integrate import quad, dblquad

hbar = 1
m = 1.0
ma = 10
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

def G1(rho, gamma):
    return hbar**2 / (2 * m) * rho**2 - Ek(rho, gamma)

def g0(rho):
    part1 = hbar**2 / (2 * m) * rho**2
    part2 = np.sqrt(vR**2 * rho**2 + a0**2)
    return part1 + part2

def g1(rho):
    part1 = hbar**2 / (2 * m) * rho**2
    part2 = np.sqrt(vR**2 * rho**2 + a0**2)
    return part1 - part2

# Helper for numerical stability with large Beta
def fermi_dist(val):
    if val > 100: return 0.0
    if val < -100: return 1.0
    return 1 / (np.exp(val) + 1)

def fk0(rho, gamma, beta, muF):
    val = beta * G0(rho, gamma) - beta * muF
    return fermi_dist(val)

def fk1(rho, gamma, beta, muF):
    val = beta * G1(rho, gamma) - beta * muF
    return fermi_dist(val)

def c0(rho):
    leading = vR**2 * rho**2 + a0**2
    val = 1/2 * a0 * hbar**2 * rho**2 * leading**(-1/2)
    return val

def d0(rho):
    leading = vR**2 * rho**2 + a0**2
    val = 1/8 * hbar**4 * (
            rho**4 * leading**(-1/2) - a0**2 * rho**4 * leading**(-3/2)
    )
    return val

def c1(rho):
    leading = vR**2 * rho**2 + a0**2
    val = -1/2 * a0 * hbar**2 * rho**2 * leading**(-1/2)
    return val

def d1(rho):
    leading = vR**2 * rho**2 + a0**2
    val = -1/8 * hbar**4 * (
            rho**4 * leading**(-1/2) - a0**2 * rho**4 * leading**(-3/2)
    )
    return val

def w00(rho, beta, muF):
    val = beta * g0(rho) - beta * muF
    return fermi_dist(val)

def w10(rho, beta, muF):
    val = beta * g1(rho) - beta * muF
    return fermi_dist(val)

def I00_integrand(rho, gamma, beta, muF):
    val = rho**3 * fk0(rho, gamma, beta, muF) * np.sin(2 * gamma) / Ek(rho, gamma)**3
    val *= -1j * vR**2 / (2 * ma) * np.sin(theta_a)
    return val

def I00_asymp_integrand(rho, beta, muF):
    leading = vR**2 * rho**2 + a0**2
    val = rho**5 * leading**(-5/2) * w00(rho, beta, muF)
    val *= 1/ma**2 * 1j * 3/4 * np.pi * vR**2 * np.sin(theta_a)**2 * a0 * hbar**2
    return val  # Fixed: Added return statement

# --- Numerical Integration Logic ---

def compute_integrals_I00(beta_list, muF):
    print(f"{'Beta':<10} | {'ma^2 * Exact (Imag)':<25} | {'ma^2 * Asymp (Imag)':<25} | {'Ratio':<10}")
    print("-" * 80)

    # Integration limit for rho (infinity approximation)
    rho_limit = 50.0

    for beta in beta_list:

        # --- 1. Exact 2D Integral ---
        # dblquad integrates func(y, x). We map x->rho, y->gamma.

        def exact_real_integrand(gamma, rho):
            return np.real(I00_integrand(rho, gamma, beta, muF))

        def exact_imag_integrand(gamma, rho):
            return np.imag(I00_integrand(rho, gamma, beta, muF))

        # Integrate Real part
        res_exact_real, _ = dblquad(exact_real_integrand, 0, rho_limit, lambda x: 0, lambda x: 2*np.pi)
        # Integrate Imaginary part
        res_exact_imag, _ = dblquad(exact_imag_integrand, 0, rho_limit, lambda x: 0, lambda x: 2*np.pi)

        I00_exact = res_exact_real + 1j * res_exact_imag

        # --- 2. Asymptotic 1D Integral ---

        def asymp_real_integrand(rho):
            return np.real(I00_asymp_integrand(rho, beta, muF))

        def asymp_imag_integrand(rho):
            return np.imag(I00_asymp_integrand(rho, beta, muF))

        res_asymp_real, _ = quad(asymp_real_integrand, 0, rho_limit)
        res_asymp_imag, _ = quad(asymp_imag_integrand, 0, rho_limit)

        I00_asymp = res_asymp_real + 1j * res_asymp_imag

        # --- 3. Scaling and Output ---

        # Scale both results by ma^2
        scaled_exact_imag = I00_exact.imag * (ma**2)
        scaled_asymp_imag = I00_asymp.imag * (ma**2)

        # Ratio check
        if abs(scaled_asymp_imag) > 1e-12:
            ratio = scaled_exact_imag / scaled_asymp_imag
        else:
            ratio = 0.0

        print(f"{beta:<10.2f} | {scaled_exact_imag:<25.5e} | {scaled_asymp_imag:<25.5e} | {ratio:<10.4f}")


muF_val = 2 * np.abs(a0)
# Define increasing Beta values
betas = [5.0, 10.0, 20.0, 50.0,100,150]

compute_integrals_I00(betas, muF_val)