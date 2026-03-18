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

def dg0(rho):
    """Derivative of g0 with respect to rho, Eq. (164)."""
    return hbar**2 / m * rho + vR**2 * rho * (a0**2 + vR**2 * rho**2)**(-0.5)

# --- Numerically Stable Fermi Functions ---

def fermi_dist(val):
    if val > 100: return 0.0
    if val < -100: return 1.0
    return 1.0 / (np.exp(val) + 1.0)

def fermi_deriv(val):
    if abs(val) > 100:
        return 0.0
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

# --- Exact Integrand ---

def I01_integrand_exact(rho, gamma, beta, muF):
    """Exact integrand for I01 in polar coords, Eq. (207)."""
    val = rho**3 * fk0(rho, gamma, beta, muF) * np.cos(2 * gamma) / Ek(rho, gamma)**3
    val *= 1j * vR**2 / (2 * ma) * np.cos(theta_a)
    return -val

# --- Asymptotic Part 1: 1D integrand with w00, Eq. (209) first term / Eq. (210) ---

def I01_integrand_asymp_part1(rho, beta, muF):
    leading = vR**2 * rho**2 + a0**2
    part1 = rho**5 * leading**(-5 / 2) * w00(rho, beta, muF)
    part1 *= 1 / ma**2 * 1j * 3 / 4 * vR**2 * a0 * hbar**2 * np.pi * np.cos(theta_a)**2
    return part1

# --- Asymptotic Part 2: delta function evaluation, Eq. (211) ---

def get_rho0(muF):
    """Calculates rho_0 where g0(rho) = muF (Eq. 203)"""
    term1 = (muF * hbar**2 / m) + vR**2
    term2 = np.sqrt(vR**4 + (2 * muF * hbar**2 / m) * vR**2 + (a0**2 * hbar**4 / m**2))
    rho0_sq = (2 * m**2 / hbar**4) * (term1 - term2)
    return np.sqrt(max(rho0_sq, 0.0))

def drho_g0(rho):
    leading = vR**2 * rho**2 + a0**2
    part1 = hbar**2 / m * rho
    part2 = vR**2 * rho * leading**(-1 / 2)
    return part1 + part2

def I01_part2_asymp(beta, muF):
    """Second term of Eq. (209), evaluated via delta function (Eq. 211).
    Beta-independent (valid for beta >> 1)."""
    rho0 = get_rho0(muF)
    leading = vR**2 * rho0**2 + a0**2
    val = rho0**3 * leading**(-3 / 2) * c0(rho0) * 1 / drho_g0(rho0)
    val *= 1 / ma**2 * 1j * 1 / 2 * np.pi * vR**2 * np.cos(theta_a)**2
    return val

# --- Also keep part 2 as a 1D integral with beta*w01 for cross-check ---

def I01_integrand_asymp_part2_finite_beta(rho, beta, muF):
    """Second term of Eq. (209) as a 1D integrand (finite beta)."""
    leading = vR**2 * rho**2 + a0**2
    part2 = rho**3 * leading**(-3 / 2) * w01(rho, beta, muF) * c0(rho)
    part2 *= 1 / ma**2 * 1j * 1 / 2 * beta * np.pi * vR**2 * np.cos(theta_a)**2
    return part2

# --- Integration Routines ---

def compute_I01_numerical(beta, muF, rho_max=15.0):
    """Compute I01 via exact 2D numerical integration over (rho, gamma)."""
    def real_integrand(gamma, rho):
        return I01_integrand_exact(rho, gamma, beta, muF).real

    def imag_integrand(gamma, rho):
        return I01_integrand_exact(rho, gamma, beta, muF).imag

    real_part, _ = dblquad(real_integrand, 0, rho_max, 0, 2 * np.pi,
                           epsabs=1e-10, epsrel=1e-10)
    imag_part, _ = dblquad(imag_integrand, 0, rho_max, 0, 2 * np.pi,
                           epsabs=1e-10, epsrel=1e-10)
    return real_part + 1j * imag_part


def compute_I01_asymp_delta(beta, muF, rho_max=15.0):
    """Asymptotic I01: Part1 (1D integral w/ w00) + Part2 (delta func, Eq.211)."""
    def real_int(rho):
        return I01_integrand_asymp_part1(rho, beta, muF).real

    def imag_int(rho):
        return I01_integrand_asymp_part1(rho, beta, muF).imag

    rp, _ = quad(real_int, 0, rho_max, limit=200)
    ip, _ = quad(imag_int, 0, rho_max, limit=200)

    part1 = rp + 1j * ip
    part2 = I01_part2_asymp(beta, muF)

    return part1 + part2


def compute_I01_asymp_finite_beta(beta, muF, rho_max=15.0):
    """Asymptotic I01: Part1 + Part2, both as 1D integrals (finite beta)."""
    def real_int(rho):
        v = (I01_integrand_asymp_part1(rho, beta, muF)
             + I01_integrand_asymp_part2_finite_beta(rho, beta, muF))
        return v.real

    def imag_int(rho):
        v = (I01_integrand_asymp_part1(rho, beta, muF)
             + I01_integrand_asymp_part2_finite_beta(rho, beta, muF))
        return v.imag

    rp, _ = quad(real_int, 0, rho_max, limit=200)
    ip, _ = quad(imag_int, 0, rho_max, limit=200)

    return rp + 1j * ip


# --- Main Comparison ---

muF = 3.0

print("=" * 80)
print("Verification of I01: exact 2D vs asymptotic (Eq. 209-211)")
print("=" * 80)
print(f"  hbar={hbar}, m={m}, ma={ma}, vR={vR}, a0={a0}")
print(f"  theta_a = {theta_a / np.pi:.2f}*pi,  muF = {muF}")
print(f"  1/ma = {1 / ma:.4f}  (asymptotic condition: 1/ma << 1)")
print()

rho0 = get_rho0(muF)
print(f"  rho0 (Fermi momentum)  = {rho0:.6f}")
print(f"  g0(rho0)              = {g0(rho0):.6f}  (should = muF = {muF})")
print(f"  dg0/drho(rho0)        = {drho_g0(rho0):.6f}")
print(f"  c0(rho0)              = {c0(rho0):.6f}")
print()

part2_delta = I01_part2_asymp(100, muF)
print(f"  Part2 (Eq.211, delta func): Im = {part2_delta.imag:.12e}")
print()

betas = [5, 10, 20, 50, 100]

header = (f"{'beta':>6s} | {'Im(exact)':>20s} | {'Im(asy_delta)':>20s} "
          f"| {'Im(asy_finite)':>20s} | {'err_delta':>10s} | {'err_finite':>10s}")
print(header)
print("-" * len(header))

for beta in betas:
    I01_num = compute_I01_numerical(beta, muF)
    I01_delta = compute_I01_asymp_delta(beta, muF)
    I01_finite = compute_I01_asymp_finite_beta(beta, muF)

    denom = abs(I01_num) if abs(I01_num) > 1e-15 else 1.0
    err_d = abs(I01_num - I01_delta) / denom
    err_f = abs(I01_num - I01_finite) / denom

    print(f"I01_num/I01_delta={I01_num/I01_delta}")