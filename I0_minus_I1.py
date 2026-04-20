import numpy as np
from scipy.integrate import quad, dblquad

hbar = 1.2
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


def G1(rho, gamma):
    return hbar**2 / (2 * m) * rho**2 - Ek(rho, gamma)


def g1(rho):
    part1 = hbar**2 / (2 * m) * rho**2
    part2 = np.sqrt(vR**2 * rho**2 + a0**2)
    return part1 - part2




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
def d0(rho):
    leading = vR**2 * rho**2 + a0**2
    part1=rho**4*leading**(-1/2)
    part2=-a0**2*rho**4*leading**(-3/2)
    val=part1+part2
    val*=1/8*hbar**4
    return val


def fk1(rho, gamma, beta, muF):
    val = beta * (G1(rho, gamma) - muF)
    return fermi_dist(val)
def c1(rho):
    leading = vR**2 * rho**2 + a0**2
    val = -0.5 * a0 * hbar**2 * rho**2 * leading**(-0.5)
    return val

def d1(rho):
    leading = vR**2 * rho**2 + a0**2
    part1=rho**4*leading**(-1/2)
    part2=-a0**2*rho**4*leading**(-3/2)
    val=part1+part2
    val*=-1/8*hbar**4
    return val

def _stable_fermi(x):
    """Compute f(x) = 1/(exp(x) + 1) without overflow."""
    x = np.asarray(x, dtype=np.float64)
    out = np.empty_like(x)
    pos = x >= 0
    # x >= 0: use exp(-x) which is <= 1, no overflow
    em = np.exp(-x[pos])
    out[pos] = em / (1.0 + em)
    # x < 0: use exp(x) which is < 1, no overflow
    ep = np.exp(x[~pos])
    out[~pos] = 1.0 / (1.0 + ep)
    return out

def w00(rho, beta, muF):
    """Fermi function: 1/(e^x + 1)."""
    x=beta*g0(rho)-beta*muF
    return _stable_fermi(x)

def w01(rho, beta, muF):
    """w01(x) = e^x / (e^x+1)^2 = f(1-f)."""
    x=beta*g0(rho)-beta*muF
    f = _stable_fermi(x)
    return f * (1.0 - f)

def w02(rho, beta, muF):
    """
    Compute  w02(x) = (e^{2x} - e^x) / (e^x + 1)^3
    via the stable factorisation  f(1-f)(1-2f).

    Parameters
    ----------
    x : array_like
        Typically x = beta * (g0(rho) - mu_F).

    Returns
    -------
    ndarray
    """
    x=beta*g0(rho)-beta*muF
    f = _stable_fermi(x)
    return f * (1.0 - f) * (1.0 - 2.0 * f)

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

def I0_integrand_dominant_asymp(rho):
    leading = vR**2 * rho**2 + a0**2
    val=rho*leading**(-3/2)
    val*=1j*a0*vR**2/hbar**2*np.pi
    return val


def w10(rho, beta, muF):
    """Fermi function: 1/(e^x + 1)."""
    x=beta*g1(rho)-beta*muF
    return _stable_fermi(x)


def w11(rho, beta, muF):
    """w01(x) = e^x / (e^x+1)^2 = f(1-f)."""
    x=beta*g1(rho)-beta*muF
    f = _stable_fermi(x)
    return f * (1.0 - f)


def w12(rho, beta, muF):
    """
    Compute  w12(x) = (e^{2x} - e^x) / (e^x + 1)^3
    via the stable factorisation  f(1-f)(1-2f).

    Parameters
    ----------
    x : array_like
        Typically x = beta * (g0(rho) - mu_F).

    Returns
    -------
    ndarray
    """
    x=beta*g1(rho)-beta*muF
    f = _stable_fermi(x)
    return f * (1.0 - f) * (1.0 - 2.0 * f)


def get_rho1(muF):
    """Calculates rho_1 where g1(rho) = muF"""
    term1 = (muF * hbar**2 / m) + vR**2
    term2 = np.sqrt(vR**4 + (2 * muF * hbar**2 / m) * vR**2 + (a0**2 * hbar**4 / m**2))
    rho0_sq = (2 * m**2 / hbar**4) * (term1 + term2)
    return np.sqrt(max(rho0_sq, 0.0))


def drho_g1(rho):
    leading = vR**2 * rho**2 + a0**2
    part1 = hbar**2 / m * rho
    part2 = vR**2 * rho * leading**(-1 / 2)
    return part1 - part2


def I1_integrand_dominant_asymp(rho):
    leading = vR**2 * rho**2 + a0**2
    val=rho*leading**(-3/2)
    val*=1j*a0*vR**2/hbar**2*np.pi
    return val

muF = 1.5
rho0=get_rho0(muF)
rho1=get_rho1(muF)

def