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
def d0(rho):
    leading = vR**2 * rho**2 + a0**2
    part1=rho**4*leading**(-1/2)
    part2=-a0**2*rho**4*leading**(-3/2)
    val=part1+part2
    val*=1/8*hbar**4
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

def I0200_integrand_exact(rho, gamma, beta, muF):
    leading = vR**2 * rho**2 + a0**2
    val1=-3*rho**4*leading**(-5/2)
    val2=15*a0**2*rho**4*leading**(-7/2)
    inside=val1+val2

    inside*=1/ma**2*1/8*hbar**4*w00(rho, beta, muF)
    inside*=np.cos(2*gamma-theta_a)**2

    inside*=rho
    inside*=1j*a0*vR**2/(2*hbar**2)
    return inside

def I0200_asymp_integrand(rho):
    leading = vR**2 * rho**2 + a0**2
    val1=-3*rho**4*leading**(-5/2)
    val2=15*a0**2*rho**4*leading**(-7/2)
    inside=val1+val2
    inside*=rho
    inside*=1j*1/ma**2*1/16*np.pi*a0*vR**2*hbar**2

    return inside

def I0201_integrand_exact(rho, gamma, beta, muF):
    leading = vR**2 * rho**2 + a0**2
    val=1/ma**2*3/2*a0*hbar**2*beta*w01(rho, beta, muF)\
        *c0(rho)*rho**2*leading**(-5/2)*np.cos(2*gamma-theta_a)**2
    val*=1j*a0*vR**2/(2*hbar**2)
    return val

def I0201_asymp(muF):
    rho0=get_rho0(muF)
    leading = vR**2 * rho0**2 + a0**2
    val=1j*1/ma**2*3/4*np.pi*a0**2*vR**2*rho0**3*c0(rho0)\
        *leading**(-5/2)*1/drho_g0(rho0)
    return val


def I0202_integrand_exact(rho, gamma, beta, muF):
    leading = vR**2 * rho**2 + a0**2
    val1=1/2*beta**2*w02(rho,beta,muF)*c0(rho)**2
    val2=-beta*w01(rho,beta,muF)*d0(rho)
    inside=val1+val2
    inside*=np.cos(2*gamma-theta_a)**2

    inside*=1/ma**2*leading**(-3/2)*rho
    inside*=1j*a0*vR**2/(2*hbar**2)
    return inside


def F_I02020(rho):
    """
    The function inside the derivative for Eq. 222:
    F(rho) = [ rho * (vR^2 * rho^2 + a0^2)^{-3/2} * c0(rho)^2 ] / [ drho_g0(rho) ]
    """
    leading = vR**2 * rho**2 + a0**2
    numerator = rho * (leading**(-1.5)) * (c0(rho)**2)
    denominator = drho_g0(rho)
    return numerator / denominator

def dF_I02020_drho_num(rho, h=1e-5):
    """
    Numerically computes the derivative of F_I02020 with respect to rho
    using the central difference method.
    """
    return (F_I02020(rho + h) - F_I02020(rho - h)) / (2.0 * h)

def I02020_asymp(muF,h=1e-5):
    rho0=get_rho0(muF)
    val=dF_I02020_drho_num(rho0,h)
    val*=1/drho_g0(rho0)
    val*=1j*1/ma**2*np.pi*a0*vR**2/(4*hbar**2)
    return val


def I02021_asymp(muF):
    rho0=get_rho0(muF)
    leading = vR**2 * rho0**2 + a0**2
    val=rho0*leading**(-3/2)*d0(rho0)
    val*=1/drho_g0(rho0)
    val*=-1j*1/ma**2*np.pi*a0*vR**2/(2*hbar**2)
    return val
def I0202_asymp(muF,h=1e-5):
    return I02020_asymp(muF,h)+I02021_asymp(muF)


def I020_integrand_exact(rho, gamma, beta, muF):
    leading = vR**2 * rho**2 + a0**2
    part0=-3*rho**4*leading**(-5/2)+15*a0**2*rho**4*leading**(-7/2)
    part0*=np.cos(2*gamma-theta_a)**2
    part0*=1/ma**2*1/8*hbar**4*w00(rho, beta, muF)

    part1=c0(rho)*rho**2*leading**(-5/2)*np.cos(2*gamma-theta_a)**2
    part1*=1/ma**2*3/2*a0*hbar**2*beta*w01(rho, beta, muF)

    part2=1/2*beta**2*w02(rho, beta, muF)*c0(rho)**2\
        -beta*w01(rho, beta, muF)*d0(rho)
    part2*=np.cos(2*gamma-theta_a)**2
    part2*=1/ma**2*leading**(-3/2)

    val=part0+part1+part2
    val*=rho
    val*=1j*a0*vR**2/(2*hbar**2)
    return val

def I020_asymp(muF, h=1e-5):
    """
    Computes the total asymptotic value of I020 by summing the
    asymptotic contributions of I0200, I0201, and I0202.
    """
    rho0 = get_rho0(muF)

    # Integrate the I0200 asymptotic integrand from 0 to rho0
    def i0200_imag(rho):
        return np.imag(I0200_asymp_integrand(rho))

    i0200_val_imag, _ = quad(i0200_imag, 0, rho0, epsabs=1e-10, epsrel=1e-10)
    i0200_val = 1j * i0200_val_imag  # Reconstruct the complex value

    # Add the delta-function evaluated parts
    return i0200_val + I0201_asymp(muF) + I0202_asymp(muF, h)

def I0210_integrand_exact(rho, gamma, beta, muF):
    leading = vR**2 * rho**2 + a0**2
    val=rho**2*leading**(-5/2)*np.cos(2*gamma-theta_a)
    val*=-1/ma*w00(rho, beta, muF)*3/2*a0*hbar**2
    val*=rho**2*np.cos(2*gamma-theta_a)
    val*=rho*1j*vR**2/(4*ma)
    return val

def I0210_asymp_integrand(rho):
    leading = vR**2 * rho**2 + a0**2
    val=rho**5*leading**(-5/2)
    val*=-1j*1/ma**2*3/8*np.pi*vR**2*a0*hbar**2
    return val
def I0211_integrand_exact(rho, gamma, beta, muF):
    leading = vR**2 * rho**2 + a0**2
    val=-1/ma*beta*leading**(-3/2)*w01(rho, beta, muF)*c0(rho)*np.cos(2*gamma-theta_a)
    val*=rho**2*np.cos(2*gamma-theta_a)
    val*=rho
    val*=1j*vR**2/(4*ma)
    return val

def I0211_asymp(muF):
    rho0=get_rho0(muF)
    leading = vR**2 * rho0**2 + a0**2
    val=rho0**3*leading**(-3/2)*c0(rho0)*1/drho_g0(rho0)
    val*=-1j*1/ma**2*1/4*np.pi*vR**2
    return val

# --- Comparison Logic ---

def compare_integrals_I0200():
    # Parameters for the low-temperature (large beta) limit
    beta = 50.0
    muF = 1.5  # Ensure muF > a0 for a valid Fermi surface

    rho0 = get_rho0(muF)
    print(f"--- Parameters ---")
    print(f"beta = {beta}")
    print(f"muF  = {muF}")
    print(f"rho0 = {rho0:.6f}\n")

    # 1. Exact Integral
    # dblquad integrates func(y, x) where x is the outer variable (rho) and y is the inner variable (gamma)
    def exact_imag(gamma, rho):
        return np.imag(I0200_integrand_exact(rho, gamma, beta, muF))

    # The Fermi function decays exponentially for rho > rho0.
    # We can safely truncate the infinite upper bound to rho0 + margin to speed up integration.
    upper_limit_rho = rho0 + 1.0

    exact_val, exact_err = dblquad(
        exact_imag,
        0, upper_limit_rho,
        lambda x: 0, lambda x: 2 * np.pi,
        epsabs=1e-10, epsrel=1e-10
    )

    # 2. Asymptotic Integral
    # The asymptotic form assumes beta -> inf, so w00 becomes a step function from 0 to rho0.
    def asymp_imag(rho):
        return np.imag(I0200_asymp_integrand(rho))

    asymp_val, asymp_err = quad(
        asymp_imag,
        0, rho0,
        epsabs=1e-10, epsrel=1e-10
    )

    # Print Results
    print(f"--- Results (Imaginary Parts) ---")
    print(f"Exact Integral:      {exact_val:.10e}")
    print(f"Asymptotic Integral: {asymp_val:.10e}")

    abs_diff = abs(exact_val - asymp_val)
    rel_diff = abs_diff / abs(asymp_val) if asymp_val != 0 else float('inf')

    print(f"Absolute Difference: {abs_diff:.10e}")
    print(f"Relative Difference: {rel_diff:.10e}")
    print(f"beta={beta}, ratio={exact_val/asymp_val}")


def compare_integrals_I0201():
    beta = 100.0
    muF = 1.5

    rho0 = get_rho0(muF)
    print(f"=== I0201 Comparison ===")
    print(f"beta = {beta}, muF = {muF}, rho0 = {rho0:.6f}")

    def exact_imag(gamma, rho):
        # Multiply by rho here to account for the polar coordinate area element (rho d_rho d_gamma)
        return np.imag(I0201_integrand_exact(rho, gamma, beta, muF) * rho)

    # w01 is sharply peaked at rho0, so integrating slightly past rho0 is sufficient
    upper_limit_rho = rho0 + 1.0

    exact_val, exact_err = dblquad(
        exact_imag,
        0, upper_limit_rho,
        lambda x: 0, lambda x: 2 * np.pi,
        epsabs=1e-10, epsrel=1e-10
    )

    # The asymptotic form is already fully evaluated (no integration needed)
    asymp_val = np.imag(I0201_asymp(muF))

    abs_diff = abs(exact_val - asymp_val)
    rel_diff = abs_diff / abs(asymp_val) if asymp_val != 0 else float('inf')

    print(f"Exact Integral:      {exact_val:.10e}")
    print(f"Asymptotic Integral: {asymp_val:.10e}")
    print(f"Absolute Difference: {abs_diff:.10e}")
    print(f"Relative Difference: {rel_diff:.10e}")
    print(f"Ratio:               {exact_val/asymp_val:.10e}\n")


def compare_integrals_I0202(h=1e-5):
    beta = 500.0
    muF = 1.5

    rho0 = get_rho0(muF)
    print(f"=== I0202 Comparison ===")
    print(f"beta = {beta}, muF = {muF}, rho0 = {rho0:.6f}")

    def exact_imag(gamma, rho):
        # rho is already included in I0202_integrand_exact
        return np.imag(I0202_integrand_exact(rho, gamma, beta, muF))

    # w01 and w02 are sharply peaked at rho0
    upper_limit_rho = rho0 + 0.1

    exact_val, exact_err = dblquad(
        exact_imag,
        0, upper_limit_rho,
        lambda x: 0, lambda x: 2 * np.pi,
        epsabs=1e-10, epsrel=1e-10
    )

    asymp_val = np.imag(I0202_asymp(muF,h))

    abs_diff = abs(exact_val - asymp_val)
    rel_diff = abs_diff / abs(asymp_val) if asymp_val != 0 else float('inf')

    print(f"Exact Integral:      {exact_val:.10e}")
    print(f"Asymptotic Integral: {asymp_val:.10e}")
    print(f"Absolute Difference: {abs_diff:.10e}")
    print(f"Relative Difference: {rel_diff:.10e}")
    print(f"Ratio:               {exact_val/asymp_val:.10e}\n")

def compare_integrals_I020(h=1e-5):
    beta = 300.0
    muF = 1.5
    rho0 = get_rho0(muF)
    print(f"=== TOTAL I020 Comparison ===")
    print(f"beta = {beta}, muF = {muF}, rho0 = {rho0:.6f}")

    def exact_imag(gamma, rho):
        # rho is already fully included inside I020_integrand_exact
        return np.imag(I020_integrand_exact(rho, gamma, beta, muF))

    # Using a tight upper limit since beta is large and w01/w02 decay very fast
    upper_limit_rho = rho0 + 0.05

    exact_val, exact_err = dblquad(
        exact_imag,
        0, upper_limit_rho,
        lambda x: 0, lambda x: 2 * np.pi,
        epsabs=1e-10, epsrel=1e-10
    )

    asymp_val = np.imag(I020_asymp(muF, h))

    abs_diff = abs(exact_val - asymp_val)
    rel_diff = abs_diff / abs(asymp_val) if asymp_val != 0 else float('inf')

    print(f"Exact Integral:      {exact_val:.10e}")
    print(f"Asymptotic Integral: {asymp_val:.10e}")
    print(f"Absolute Difference: {abs_diff:.10e}")
    print(f"Relative Difference: {rel_diff:.10e}")
    print(f"Ratio:               {exact_val/asymp_val:.10e}\n")

def compare_integrals_I0210():
    beta = 100.0
    muF = 1.5
    rho0 = get_rho0(muF)
    print(f"=== I0210 Comparison ===")
    print(f"beta = {beta}, muF = {muF}, rho0 = {rho0:.6f}")

    def exact_imag(gamma, rho):
        # I0210_integrand_exact already includes the rho factor for the area element
        return np.imag(I0210_integrand_exact(rho, gamma, beta, muF))

    # w00 decays exponentially for rho > rho0
    upper_limit_rho = rho0 + 1.0

    exact_val, _ = dblquad(
        exact_imag,
        0, upper_limit_rho,
        lambda x: 0, lambda x: 2 * np.pi,
        epsabs=1e-10, epsrel=1e-10
    )

    def asymp_imag(rho):
        return np.imag(I0210_asymp_integrand(rho))

    asymp_val, _ = quad(
        asymp_imag,
        0, rho0,
        epsabs=1e-10, epsrel=1e-10
    )

    abs_diff = abs(exact_val - asymp_val)
    rel_diff = abs_diff / abs(asymp_val) if asymp_val != 0 else float('inf')

    print(f"Exact Integral:      {exact_val:.10e}")
    print(f"Asymptotic Integral: {asymp_val:.10e}")
    print(f"Absolute Difference: {abs_diff:.10e}")
    print(f"Relative Difference: {rel_diff:.10e}")
    print(f"Ratio:               {exact_val/asymp_val:.10e}\n")
def compare_integrals_I0211():
    beta = 100.0
    muF = 1.5
    rho0 = get_rho0(muF)
    print(f"=== I0211 Comparison ===")
    print(f"beta = {beta}, muF = {muF}, rho0 = {rho0:.6f}")

    def exact_imag(gamma, rho):
        # I0211_integrand_exact already includes the rho factor for the area element
        return np.imag(I0211_integrand_exact(rho, gamma, beta, muF))

    # w01 is sharply peaked at rho0, so integrating slightly past rho0 is sufficient
    upper_limit_rho = rho0 + 0.1

    exact_val, _ = dblquad(
        exact_imag,
        0, upper_limit_rho,
        lambda x: 0, lambda x: 2 * np.pi,
        epsabs=1e-10, epsrel=1e-10
    )

    asymp_val = np.imag(I0211_asymp(muF))

    abs_diff = abs(exact_val - asymp_val)
    rel_diff = abs_diff / abs(asymp_val) if asymp_val != 0 else float('inf')

    print(f"Exact Integral:      {exact_val:.10e}")
    print(f"Asymptotic Integral: {asymp_val:.10e}")
    print(f"Absolute Difference: {abs_diff:.10e}")
    print(f"Relative Difference: {rel_diff:.10e}")
    print(f"Ratio:               {exact_val/asymp_val:.10e}\n")

compare_integrals_I0211()