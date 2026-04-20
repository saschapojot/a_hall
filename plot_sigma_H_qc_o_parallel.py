from sympy import *
from sympy import expand_complex
from sympy.simplify.fu import TR5, TR11, TR9, fu
from sympy.physics.units import hbar
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import multiprocessing
from datetime import datetime

alpha_k = symbols("alpha_k", cls=Symbol, real=True)
theta_k = symbols("theta_k", cls=Symbol, real=True)
kx = symbols("k_x", cls=Symbol, real=True)
ky = symbols("k_y", cls=Symbol, real=True)
m = symbols("m", cls=Symbol, positive=True)
theta_a = symbols("theta_a", cls=Symbol, positive=True)
vR = symbols("v_R", cls=Symbol, positive=True, nonzero=True)
a0 = symbols("a0", cls=Symbol, positive=True)
rho = symbols("rho", cls=Symbol, positive=True)
gamma = symbols("gamma", cls=Symbol, real=True)
muF = symbols("muF", cls=Symbol, real=True)
eps = symbols("epsilon", cls=Symbol, positive=True)
beta = symbols("beta", cls=Symbol, positive=True)
ma = 1 / eps

k_abs = sqrt(kx ** 2 + ky ** 2)
cos_theta_k = kx / k_abs * sign(vR)
sin_theta_k = ky / k_abs * sign(vR)

Ek2 = vR ** 2 * (kx ** 2 + ky ** 2) \
      + (a0 + hbar ** 2 / (2 * ma) * ((kx ** 2 - ky ** 2) * cos(theta_a) + 2 * kx * ky * sin(theta_a))) ** 2
Ek = sqrt(Ek2)
epsilon_k_0 = hbar ** 2 / (2 * m) * (kx ** 2 + ky ** 2) + Ek
epsilon_k_1 = hbar ** 2 / (2 * m) * (kx ** 2 + ky ** 2) - Ek
fk0 = 1 / (exp(beta * (epsilon_k_0 - muF)) + 1)
fk1 = 1 / (exp(beta * (epsilon_k_1 - muF)) + 1)

cos_theta_a = cos(theta_a)
sin_theta_a = sin(theta_a)

sin_alpha_k = abs(vR) / Ek * k_abs

epsilon_L_k_a = a0 + hbar ** 2 / (2 * ma) * ((kx ** 2 - ky ** 2) * cos(theta_a) + 2 * kx * ky * sin(theta_a))
cos_alpha_k = epsilon_L_k_a / Ek

half = Rational(1, 2)
one_8 = Rational(1, 8)
three_2 = Rational(3, 2)

v_xk_0e_00 = ((1 / m + cos_theta_a / ma) * hbar * kx + sin_theta_a / ma * hbar * ky) * (half + half * cos_alpha_k) \
             + 1 / hbar * vR * sin_alpha_k * cos_alpha_k \
             + ((1 / m - cos_theta_a / ma) * hbar * kx - sin_theta_a / ma * hbar * ky) * (half - half * cos_alpha_k)

v_xk_0e_01 = hbar / ma * sin_alpha_k * (cos(theta_a) * kx + sin(theta_a) * ky) \
             - I * 1 / hbar * vR * sin_theta_k - 1 / hbar * vR * cos_alpha_k * cos_theta_k

v_xk_0e_10 = hbar / ma * sin_alpha_k * (cos(theta_a) * kx + sin(theta_a) * ky) + I * 1 / hbar * vR * sin_theta_k \
             - 1 / hbar * vR * cos_alpha_k * cos_theta_k

v_xk_0e_11 = ((1 / m + cos(theta_a) / ma) * hbar * kx + sin(theta_a) / ma * hbar * ky) * (half - half * cos_alpha_k) \
             - 1 / hbar * vR * sin_alpha_k * cos_theta_k \
             + ((1 / m - cos(theta_a) / ma) * hbar * kx - sin(theta_a) / ma * hbar * ky) * (half + half * cos_alpha_k)

v_yk_0e_00 = (half + half * cos_alpha_k) * (sin(theta_a) / ma * hbar * kx + (1 / m - cos(theta_a) / ma) * hbar * ky) \
             + 1 / hbar * vR * sin_alpha_k * sin_theta_k \
             + (-sin(theta_a) / ma * hbar * kx + (1 / m + cos(theta_a) / ma) * hbar * ky) * (half - half * cos_alpha_k)

v_yk_0e_01 = hbar / ma * sin_alpha_k * (sin(theta_a) * kx - cos(
    theta_a) * ky) + I * 1 / hbar * vR * cos_theta_k - 1 / hbar * vR * cos_alpha_k * sin_theta_k

v_yk_0e_10 = hbar / ma * sin_alpha_k * (sin(theta_a) * kx - cos(
    theta_a) * ky) - I * 1 / hbar * vR * cos_theta_k - 1 / hbar * vR * cos_alpha_k * sin_theta_k

v_yk_0e_11 = (sin(theta_a) / ma * hbar * kx + (1 / m - cos(theta_a) / ma) * hbar * ky) * (
        half - half * cos_alpha_k) - 1 / hbar * vR * sin_alpha_k * sin_theta_k \
             + (-sin(theta_a) / ma * hbar * kx + (1 / m + cos(theta_a) / ma) * hbar * ky) * (half + half * cos_alpha_k)

integrand_sigma = (fk0 - fk1) * (v_yk_0e_10 * v_xk_0e_01 - v_yk_0e_01 * v_xk_0e_10) / (epsilon_k_0 - epsilon_k_1) ** 2

# Define the numerical values in a dictionary
eV = 1.602176634e-19
meV = eV * 0.001
h = 6.62607015e-34
a = 3.5e-10
multiple = 5
subs_dict = {
    hbar: 6.62607015e-34 / (2 * np.pi),
    m: 9.1093837e-31,
    eps: 1 / (multiple * 9.1093837e-31),
    vR: 12 * meV * a,
    a0: 0.1 * meV,
    muF: 3 * meV
}

e = -1.602176634e-19
unit = e ** 2 / h


def sigma_analytical(subs_dict):
    muF_val = subs_dict[muF]
    m_val = subs_dict[m]
    vR_val = subs_dict[vR]
    a0_val = subs_dict[a0]
    hbar_val = subs_dict[hbar]
    I0_minus_I1 = -1j * a0_val * np.pi / hbar_val ** 2 \
                  * 2 * m_val * vR_val ** 2 / (2 * m_val * vR_val ** 2 * muF_val + a0_val ** 2 * hbar_val ** 2)
    sigma = I0_minus_I1 * 1j * e ** 2 * subs_dict[hbar] / (4 * np.pi ** 2)
    return sigma


integrand_sigma_num = integrand_sigma.subs(subs_dict) * 1j * e ** 2 * subs_dict[hbar] / (4 * np.pi ** 2)
integrand_func = lambdify((kx, ky, theta_a, beta), integrand_sigma_num, modules='numpy')


def get_rho0(subs_dict):
    muF_val = subs_dict[muF]
    m_val = subs_dict[m]
    vR_val = subs_dict[vR]
    a0_val = subs_dict[a0]
    hbar_val = subs_dict[hbar]
    term1 = (muF_val * hbar_val ** 2 / m_val) + vR_val ** 2
    term2 = np.sqrt(
        vR_val ** 4 + (2 * muF_val * hbar_val ** 2 / m_val) * vR_val ** 2 + (a0_val ** 2 * hbar_val ** 4 / m_val ** 2))
    rho0_sq = (2 * m_val ** 2 / hbar_val ** 4) * (term1 - term2)
    return np.sqrt(max(rho0_sq, 0.0))


def get_rho1(subs_dict):
    muF_val = subs_dict[muF]
    m_val = subs_dict[m]
    vR_val = subs_dict[vR]
    a0_val = subs_dict[a0]
    hbar_val = subs_dict[hbar]
    term1 = (muF_val * hbar_val ** 2 / m_val) + vR_val ** 2
    term2 = np.sqrt(
        vR_val ** 4 + (2 * muF_val * hbar_val ** 2 / m_val) * vR_val ** 2 + (a0_val ** 2 * hbar_val ** 4 / m_val ** 2))
    rho1_sq = (2 * m_val ** 2 / hbar_val ** 4) * (term1 + term2)
    return np.sqrt(max(rho1_sq, 0.0))


rho0 = get_rho0(subs_dict)
rho1 = get_rho1(subs_dict)

kB = 1.380649e-23
T = 0.01
beta_val = 1 / (kB * T)

# Determine the maximum radius for the integration domain
rho_max = max(rho0, rho1)
theoretical_value = np.real(sigma_analytical(subs_dict))


# --- Worker Function for Multiprocessing ---
def compute_integral(ta_val):
    """
    Worker function to compute the double integral for a specific theta_a value.
    Returns a tuple of (theta_a, integral_result) to ensure order can be explicitly reconstructed.
    """

    def polar_integrand_real(r, phi):
        kx_val = r * np.cos(phi)
        ky_val = r * np.sin(phi)
        val = integrand_func(kx_val, ky_val, ta_val, beta_val)
        return np.real(val) * r

    integral_real, error_real = integrate.dblquad(
        polar_integrand_real,
        0, 2 * np.pi,
        lambda phi: 0, lambda phi: rho_max,
        epsabs=1e-12, epsrel=1e-12
    )

    print(f"Completed theta_a = {ta_val:.3f} rad | Real Integral = {integral_real:.6e}")
    return (ta_val, integral_real)


# --- Main Execution Block ---
if __name__ == '__main__':
    t_int_start = datetime.now()
    num_points = 300
    theta_a_array = np.linspace(0, 2 * np.pi, num_points)

    # Define the number of parallel processes
    # multiprocessing.cpu_count() uses all available cores.
    # You can change this to a specific integer, e.g., num_processes = 4
    num_processes = 24

    print(f"Starting parallel integration over a circle of radius rho_max = {rho_max:.4f}...")
    print(f"Calculating for {num_points} values of theta_a from 0 to 2*pi using {num_processes} processes.\n")

    # Use multiprocessing.Pool to compute integrals in parallel
    with multiprocessing.Pool(processes=num_processes) as pool:
        # We can use imap_unordered for better performance since we sort manually afterwards
        results = list(pool.imap_unordered(compute_integral, theta_a_array))

    print("\nIntegration complete! Sorting results and generating plot...")
    t_int_end = datetime.now()
    print("time: ", t_int_end - t_int_start)

    # Sort the results explicitly based on ta_val (the first element of the tuple)
    results.sort(key=lambda x: x[0])

    # Unpack the sorted results
    sorted_theta_a = np.array([res[0] for res in results])
    real_integrals = np.array([res[1] for res in results])

    # --- Plotting ---
    plt.figure(figsize=(8, 5))
    plt.scatter(sorted_theta_a, real_integrals / unit, marker='.', color='b', label=r'Numerical Re($\sigma$)')

    # Add the theoretical value as a horizontal line
    plt.axhline(y=theoretical_value / unit, color='r', linestyle='--', linewidth=2,
                label=rf'Theoretical Re($\sigma$) = {theoretical_value / unit:.4e}')

    plt.xlabel(r'$\theta_a$ (radians)', fontsize=12)
    plt.ylabel(r'Re($\sigma$)', fontsize=12)
    plt.title(r'Real part of the Integral vs $\theta_a$', fontsize=14)
    plt.xticks(
        [0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi],
        ['0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$']
    )
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"sigma_multiple{multiple}.png")
    plt.close()