from mpmath import mp

# Set the desired precision (number of decimal places)
mp.dps = 100

hbar = mp.mpf('0.1')
m = mp.mpf('1.0')
ma = mp.mpf('20.0')
vR = mp.mpf('2.0')
a0 = mp.mpf('0.5')
theta_a = mp.mpf('0.1') * mp.pi


# --- Physics Functions ---

def Ek(rho, gamma):
    part1 = vR**2 * rho**2
    val2 = a0 + hbar**2 / (mp.mpf('2.0') * ma) * rho**2 * mp.cos(mp.mpf('2.0') * gamma - theta_a)
    Ek2 = part1 + val2**2
    return mp.sqrt(Ek2)

def G0(rho, gamma):
    return hbar**2 / (mp.mpf('2.0') * m) * rho**2 + Ek(rho, gamma)

def g0(rho):
    part1 = hbar**2 / (mp.mpf('2.0') * m) * rho**2
    part2 = mp.sqrt(vR**2 * rho**2 + a0**2)
    return part1 + part2

def dg0(rho):
    """Derivative of g0 with respect to rho, Eq. (164)."""
    return hbar**2 / m * rho + vR**2 * rho * (a0**2 + vR**2 * rho**2)**(-mp.mpf('0.5'))


def G1(rho, gamma):
    return hbar**2 / (mp.mpf('2.0') * m) * rho**2 - Ek(rho, gamma)


def g1(rho):
    part1 = hbar**2 / (mp.mpf('2.0') * m) * rho**2
    part2 = mp.sqrt(vR**2 * rho**2 + a0**2)
    return part1 - part2


# --- Numerically Stable Fermi Functions ---

def fermi_dist(val):
    if val > 100: return mp.mpf('0.0')
    if val < -100: return mp.mpf('1.0')
    return mp.mpf('1.0') / (mp.exp(val) + mp.mpf('1.0'))

def fermi_deriv(val):
    if abs(val) > 100:
        return mp.mpf('0.0')
    return mp.mpf('0.25') / (mp.cosh(val / mp.mpf('2.0'))**2)

def fk0(rho, gamma, beta, muF):
    val = beta * (G0(rho, gamma) - muF)
    return fermi_dist(val)

def c0(rho):
    leading = vR**2 * rho**2 + a0**2
    val = mp.mpf('0.5') * a0 * hbar**2 * rho**2 * leading**(-mp.mpf('0.5'))
    return val

def d0(rho):
    leading = vR**2 * rho**2 + a0**2
    part1 = rho**4 * leading**(-mp.mpf('0.5'))
    part2 = -a0**2 * rho**4 * leading**(-mp.mpf('1.5'))
    val = part1 + part2
    val *= mp.mpf('0.125') * hbar**4
    return val


def fk1(rho, gamma, beta, muF):
    val = beta * (G1(rho, gamma) - muF)
    return fermi_dist(val)

def c1(rho):
    leading = vR**2 * rho**2 + a0**2
    val = -mp.mpf('0.5') * a0 * hbar**2 * rho**2 * leading**(-mp.mpf('0.5'))
    return val

def d1(rho):
    leading = vR**2 * rho**2 + a0**2
    part1 = rho**4 * leading**(-mp.mpf('0.5'))
    part2 = -a0**2 * rho**4 * leading**(-mp.mpf('1.5'))
    val = part1 + part2
    val *= -mp.mpf('0.125') * hbar**4
    return val

def _stable_fermi(x):
    """Compute f(x) = 1/(exp(x) + 1) without overflow for mpmath scalars."""
    if x >= 0:
        em = mp.exp(-x)
        return em / (mp.mpf('1.0') + em)
    else:
        ep = mp.exp(x)
        return mp.mpf('1.0') / (mp.mpf('1.0') + ep)

def w00(rho, beta, muF):
    """Fermi function: 1/(e^x + 1)."""
    x = beta * g0(rho) - beta * muF
    return _stable_fermi(x)

def w01(rho, beta, muF):
    """w01(x) = e^x / (e^x+1)^2 = f(1-f)."""
    x = beta * g0(rho) - beta * muF
    f = _stable_fermi(x)
    return f * (mp.mpf('1.0') - f)

def w02(rho, beta, muF):
    """
    Compute  w02(x) = (e^{2x} - e^x) / (e^x + 1)^3
    via the stable factorisation  f(1-f)(1-2f).
    """
    x = beta * g0(rho) - beta * muF
    f = _stable_fermi(x)
    return f * (mp.mpf('1.0') - f) * (mp.mpf('1.0') - mp.mpf('2.0') * f)

def get_rho0(muF):
    """Calculates rho_0 where g0(rho) = muF (Eq. 203) without catastrophic cancellation."""
    term1 = (muF * hbar**2 / m) + vR**2
    term2 = mp.sqrt(vR**4 + (mp.mpf('2.0') * muF * hbar**2 / m) * vR**2 + (a0**2 * hbar**4 / m**2))

    # Use the conjugate to avoid catastrophic cancellation:
    rho0_sq = mp.mpf('2.0') * (muF**2 - a0**2) / (term1 + term2)

    return mp.sqrt(max(rho0_sq, mp.mpf('0.0')))

def drho_g0(rho):
    leading = vR**2 * rho**2 + a0**2
    part1 = hbar**2 / m * rho
    part2 = vR**2 * rho * leading**(-mp.mpf('0.5'))
    return part1 + part2

def I0_integrand_dominant_asymp(rho):
    leading = vR**2 * rho**2 + a0**2
    val = rho * leading**(-mp.mpf('1.5'))
    val *= mp.mpc(0, 1) * a0 * vR**2 / hbar**2 * mp.pi
    return val


def w10(rho, beta, muF):
    """Fermi function: 1/(e^x + 1)."""
    x = beta * g1(rho) - beta * muF
    return _stable_fermi(x)


def w11(rho, beta, muF):
    """w01(x) = e^x / (e^x+1)^2 = f(1-f)."""
    x = beta * g1(rho) - beta * muF
    f = _stable_fermi(x)
    return f * (mp.mpf('1.0') - f)


def w12(rho, beta, muF):
    """
    Compute  w12(x) = (e^{2x} - e^x) / (e^x + 1)^3
    via the stable factorisation  f(1-f)(1-2f).
    """
    x = beta * g1(rho) - beta * muF
    f = _stable_fermi(x)
    return f * (mp.mpf('1.0') - f) * (mp.mpf('1.0') - mp.mpf('2.0') * f)


def get_rho1(muF):
    """Calculates rho_1 where g1(rho) = muF"""
    term1 = (muF * hbar**2 / m) + vR**2
    term2 = mp.sqrt(vR**4 + (mp.mpf('2.0') * muF * hbar**2 / m) * vR**2 + (a0**2 * hbar**4 / m**2))
    rho0_sq = (mp.mpf('2.0') * m**2 / hbar**4) * (term1 + term2)
    return mp.sqrt(max(rho0_sq, mp.mpf('0.0')))


def drho_g1(rho):
    leading = vR**2 * rho**2 + a0**2
    part1 = hbar**2 / m * rho
    part2 = vR**2 * rho * leading**(-mp.mpf('0.5'))
    return part1 - part2


def I1_integrand_dominant_asymp(rho):
    leading = vR**2 * rho**2 + a0**2
    val = rho * leading**(-mp.mpf('1.5'))
    val *= mp.mpc(0, 1) * a0 * vR**2 / hbar**2 * mp.pi
    return val


def get_rho0_sq(muF):
    term1 = (muF * hbar**2 / m) + vR**2
    term2 = mp.sqrt(vR**4 + (mp.mpf('2.0') * muF * hbar**2 / m) * vR**2 + (a0**2 * hbar**4 / m**2))
    rho0_sq = mp.mpf('2.0') * (muF**2 - a0**2) / (term1 + term2)
    return max(rho0_sq, mp.mpf('0.0'))

def get_rho1_sq(muF):
    term1 = (muF * hbar**2 / m) + vR**2
    term2 = mp.sqrt(vR**4 + (mp.mpf('2.0') * muF * hbar**2 / m) * vR**2 + (a0**2 * hbar**4 / m**2))
    rho1_sq = (mp.mpf('2.0') * m**2 / hbar**4) * (term1 + term2)
    return max(rho1_sq, mp.mpf('0.0'))

muF = mp.mpf('1.5')
rho0_sq = get_rho0_sq(muF)
rho1_sq = get_rho1_sq(muF)


# Use mp.power for exponents to prevent ANY standard Python float casting
lhs_part1 = mp.mpf('0.5') * mp.power(hbar, 4) * rho0_sq - (mp.power(hbar, 2) * m * muF + mp.power(m, 2) * mp.power(vR, 2))
lhs_part2 = mp.mpf('0.5') * mp.power(hbar, 4) * rho1_sq - (mp.power(hbar, 2) * m * muF + mp.power(m, 2) * mp.power(vR, 2))
lhs_part3 = mp.power(mp.power(a0, 2) + mp.power(vR, 2) * rho0_sq, 3)
lhs_part4 = mp.power(mp.power(a0, 2) + mp.power(vR, 2) * rho1_sq, 3)

lhs = lhs_part1 * lhs_part2 * lhs_part3 * lhs_part4

rhs_part1 = mp.power(m, 2) * mp.power(hbar, 4) * (mp.power(muF, 2) - mp.power(a0, 2)) - mp.power(m, 2) * mp.power(muF * mp.power(hbar, 2) + m * mp.power(vR, 2), 2)
rhs_part2 = mp.mpf('4.0') * mp.power(m, 2) * mp.power(vR, 4) / mp.power(hbar, 4) * (mp.power(muF, 2) - mp.power(a0, 2))
rhs_part3 = mp.mpf('4.0') * m * mp.power(vR, 2) * mp.power(a0, 2) / mp.power(hbar, 4) * (muF * mp.power(hbar, 2) + m * mp.power(vR, 2))
rhs_part4 = mp.power(a0, 4)

rhs = rhs_part1 * mp.power(rhs_part2 + rhs_part3 + rhs_part4, 3)

df = lhs - rhs

print("LHS:", lhs)
print("Difference (df):", df)