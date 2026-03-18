import numpy as np
from scipy.integrate import quad, dblquad

hbar = 1
m = 1.0
ma = 30.0
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

def I0_integrand_dominant_asymp(rho):
    leading = vR**2 * rho**2 + a0**2
    val=rho*leading**(-3/2)
    val*=1j*a0*vR**2/hbar**2*np.pi
    return val

def I0_integrand_asymp_part0(rho):
    leading = vR**2 * rho**2 + a0**2
    val=rho**5*leading**(-5/2)
    val*=1/ma**2*1j*3/4*np.pi*vR**2*np.sin(theta_a)**2*a0*hbar**2
    return val

# def I0_integrand_asymp_part1(rho):
#     leading = vR**2 * rho**2 + a0**2
#     val = rho**5 * leading**(-5/2)
#     # FIXED: Changed np.pi**np.cos to np.pi*np.cos
#     val *= -1/ma**2 * 1j * 3/4 * vR**2 * a0 * hbar**2 * np.pi * np.cos(theta_a)**2
#     return val

def I0_integrand_asymp_part1(rho):
    leading = vR**2 * rho**2 + a0**2
    val = rho**5 * leading**(-5/2)
    # FIXED: Sign error in document Eq. (207). k_y^2 - k_x^2 = -rho^2 cos(2*gamma).
    # The negative sign from the document is removed (making the term positive).
    val *= 1/ma**2 * 1j * 3/4 * vR**2 * a0 * hbar**2 * np.pi * np.cos(theta_a)**2
    return val

def I0_integrand_asymp_part2(rho):
    leading = vR**2 * rho**2 + a0**2
    val1=-3*rho**4*leading**(-5/2)
    val2=15*a0**2*rho**4*leading**(-7/2)
    val=rho*(val1+val2)
    val*=1j*1/ma**2*1/16*np.pi*a0*vR**2*hbar**2
    return val

def I0_integrand_asymp_part3(rho):
    leading = vR**2 * rho**2 + a0**2
    val=rho**5*leading**(-5/2)
    val*=-1j*1/ma**2*3/8*np.pi*vR**2*a0*hbar**2
    return val

def I0_integrand_asymp(rho):
    return I0_integrand_dominant_asymp(rho)+ I0_integrand_asymp_part0(rho)+I0_integrand_asymp_part1(rho)\
            +I0_integrand_asymp_part2(rho)+I0_integrand_asymp_part3(rho)

def I0_aymp_val0(muF):
    rho0=get_rho0(muF)
    leading = vR**2 * rho0**2 + a0**2
    val=rho0**3*leading**(-3/2)*c0(rho0)*1/drho_g0(rho0)
    val*=1/ma**2*1j*1/2*np.pi*vR**2*np.sin(theta_a)**2
    return val

# def I0_asymp_val1(muF):
#     rho0=get_rho0(muF)
#     leading = vR**2 * rho0**2 + a0**2
#     val=rho0**3*leading**(-3/2)*c0(rho0)*1/drho_g0(rho0)
#     val*=-1/ma**2*1j*1/2*np.pi*vR**2*np.cos(theta_a)**2
#     return val

def I0_asymp_val1(muF):
    rho0 = get_rho0(muF)
    leading = vR**2 * rho0**2 + a0**2
    val = rho0**3 * leading**(-3/2) * c0(rho0) * 1/drho_g0(rho0)
    # FIXED: Sign error in document Eq. (207) propagates here.
    # The negative sign from the document is removed (making the term positive).
    val *= 1/ma**2 * 1j * 1/2 * np.pi * vR**2 * np.cos(theta_a)**2
    return val

def I0_asymp_val2(muF):
    rho0=get_rho0(muF)
    leading = vR**2 * rho0**2 + a0**2
    val=rho0**3*c0(rho0)*leading**(-5/2)*1/drho_g0(rho0)
    val*=1j*1/ma**2*3/4*np.pi*a0**2*vR**2
    return val

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

def I0_asymp_val3(muF,h=1e-5):
    rho0=get_rho0(muF)
    leading = vR**2 * rho0**2 + a0**2
    val1=dF_I02020_drho_num(rho0,h)
    val1*=1/drho_g0(rho0)
    val1*=1j*1/ma**2*np.pi*a0*vR**2/(4*hbar**2)

    val2=rho0*leading**(-3/2)*d0(rho0)
    val2*=-1j*1/ma**2*np.pi*a0*vR**2/(2*hbar**2)*1/drho_g0(rho0)
    return val1+val2

def I0_asymp_val4(muF):
    rho0=get_rho0(muF)
    leading = vR**2 * rho0**2 + a0**2
    val=rho0**3*leading**(-3/2)*c0(rho0)*1/drho_g0(rho0)
    val*=-1j*1/ma**2*1/4*np.pi*vR**2
    return val

def I0_asymp_direct(muF,h=1e-5):
    return I0_aymp_val0(muF)+I0_asymp_val1(muF)+I0_asymp_val2(muF)\
            +I0_asymp_val3(muF,h)+I0_asymp_val4(muF)


def I0_integrand_exact(rho, gamma, beta, muF):
    kx=rho*np.cos(gamma)
    ky=rho*np.sin(gamma)
    Ek_val=Ek(rho, gamma)

    cos_alpha_k=(a0+hbar**2/(2*ma)*(kx**2+ky**2)*np.cos(2*gamma-theta_a))/Ek_val
    part1=-1j*4*vR**2/ma*np.sin(theta_a)*kx*ky/Ek_val

    part2=1j*2*vR**2/ma*np.cos(theta_a)*(ky**2-kx**2)/Ek_val

    part3=1j*2*vR**2/hbar**2*cos_alpha_k

    val=fk0(rho, gamma, beta, muF)* (part1+part2+part3)/(4*Ek_val**2)
    return val
# --- Comparison Logic ---

def compare_integrals_I0(h=1e-5):
    beta = 100.0
    muF = 1.5
    rho0 = get_rho0(muF)

    print(f"=== TOTAL I0 Comparison ===")
    print(f"beta = {beta}, muF = {muF}, rho0 = {rho0:.6f}")

    # 1. Exact 2D Integral
    def exact_imag(gamma, rho):
        # Multiply by rho for the polar coordinate area element (rho d_rho d_gamma)
        return np.imag(I0_integrand_exact(rho, gamma, beta, muF) * rho)

    # The Fermi function fk0 decays exponentially for rho > rho0.
    upper_limit_rho = rho0 + 0.1

    exact_val, exact_err = dblquad(
        exact_imag,
        0, upper_limit_rho,
        lambda x: 0, lambda x: 2 * np.pi,
        epsabs=1e-10, epsrel=1e-10
    )

    # 2. Asymptotic 1D Integral
    def asymp_imag(rho):
        return np.imag(I0_integrand_asymp(rho))

    asymp_int_val, asymp_int_err = quad(
        asymp_imag,
        0, rho0,
        epsabs=1e-10, epsrel=1e-10
    )

    # 3. Direct Asymptotic Boundary Terms
    asymp_dir_val = np.imag(I0_asymp_direct(muF, h))

    # Total Asymptotic Value
    total_asymp_val = asymp_int_val + asymp_dir_val

    # Print Results
    abs_diff = abs(exact_val - total_asymp_val)
    rel_diff = abs_diff / abs(total_asymp_val) if total_asymp_val != 0 else float('inf')

    print(f"Exact Integral:           {exact_val:.10e}")
    print(f"Asymp Integral Part:      {asymp_int_val:.10e}")
    print(f"Asymp Direct Part:        {asymp_dir_val:.10e}")
    print(f"Total Asymptotic:         {total_asymp_val:.10e}")
    print(f"Absolute Difference:      {abs_diff:.10e}")
    print(f"Relative Difference:      {rel_diff:.10e}")
    print(f"Ratio (Exact/Asymp):      {exact_val/total_asymp_val:.10e}\n")



def compare_integrals_I0_subtracted(h=1e-5):
    """
    Compares the exact and asymptotic integrals after subtracting
    the dominant asymptotic term from both.
    """
    beta = 10000.0  # Changed back to 1000.0 to match the first function
    muF = 1.5
    rho0 = get_rho0(muF)

    print(f"=== SUBTRACTED I0 Comparison (Without Dominant Term) ===")
    print(f"beta = {beta}, muF = {muF}, rho0 = {rho0:.6f}")

    # 1. Exact 2D Integral (Requires upper_limit_rho due to Fermi function tail)
    def exact_imag(gamma, rho):
        return np.imag(I0_integrand_exact(rho, gamma, beta, muF) * rho)

    upper_limit_rho = rho0 + 0.1

    exact_val, exact_err = dblquad(
        exact_imag,
        0, upper_limit_rho,
        lambda x: 0, lambda x: 2 * np.pi,
        epsabs=1e-10, epsrel=1e-10
    )

    # 2. Integral of the Dominant Asymptotic Term (Must stop exactly at rho0)
    def dom_asymp_imag(rho):
        return np.imag(I0_integrand_dominant_asymp(rho))

    dom_val, dom_err = quad(
        dom_asymp_imag,
        0, rho0,  # FIXED: Changed from upper_limit_rho to rho0
        epsabs=1e-10, epsrel=1e-10
    )

    # 3. Asymptotic 1D Integral (Full) (Must stop exactly at rho0)
    def asymp_imag(rho):
        return np.imag(I0_integrand_asymp(rho))

    asymp_int_val, asymp_int_err = quad(
        asymp_imag,
        0, rho0,  # FIXED: Changed from upper_limit_rho to rho0
        epsabs=1e-10, epsrel=1e-10
    )

    # 4. Direct Asymptotic Boundary Terms
    asymp_dir_val = np.imag(I0_asymp_direct(muF, h))

    # Total Asymptotic Value
    total_asymp_val = asymp_int_val + asymp_dir_val
    print(f"total_asymp_val={total_asymp_val}")
    # --- Subtractions ---
    exact_subtracted = exact_val - dom_val
    asymp_subtracted = total_asymp_val - dom_val

    # Print Results
    abs_diff = abs(exact_subtracted - asymp_subtracted)
    rel_diff = abs_diff / abs(asymp_subtracted) if asymp_subtracted != 0 else float('inf')

    print(f"Original Exact Integral:  {exact_val:.10e}")
    print(f"Dominant Asymp Integral:  {dom_val:.10e}")
    print(f"Exact Subtracted:         {exact_subtracted:.10e}")
    print(f"Asymp Subtracted:         {asymp_subtracted:.10e}")
    print(f"Absolute Difference:      {abs_diff:.10e}")
    print(f"Relative Difference:      {rel_diff:.10e}")
    print(f"Ratio (Exact/Asymp):      {exact_subtracted/asymp_subtracted if asymp_subtracted != 0 else float('inf'):.10e}\n")


compare_integrals_I0()
compare_integrals_I0_subtracted()

