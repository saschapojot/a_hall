from sympy import *
from sympy import expand_complex
from sympy.simplify.fu import TR5, TR11, TR9,fu
from sympy.physics.units import hbar
import numpy as np

alpha_k=symbols("alpha_k",cls=Symbol,real=True)
theta_k=symbols("theta_k",cls=Symbol,real=True)
kx=symbols("k_x",cls=Symbol,real=True)
ky=symbols("k_y",cls=Symbol,real=True)
m=symbols("m",cls=Symbol,positive=True)
# ma=symbols("m_a",cls=Symbol,positive=True)
theta_a=symbols("theta_a",cls=Symbol,positive=True)
vR=symbols("v_R",cls=Symbol,positive=True,nonzero=True)
# Ek=symbols("Ek",cls=Symbol,positive=True)
a0=symbols("a0",cls=Symbol,positive=True)
rho=symbols("rho",cls=Symbol,positive=True)
gamma=symbols("gamma",cls=Symbol,real=True)
muF=symbols("muF",cls=Symbol,real=True)
eps=symbols("epsilon",cls=Symbol,positive=True)
beta=symbols("beta",cls=Symbol,positive=True)
ma=1/eps

k_abs=sqrt(kx**2+ky**2)
cos_theta_k=kx/k_abs*sign(vR)
sin_theta_k=ky/k_abs*sign(vR)

Ek2=vR**2*(kx**2+ky**2)\
    +(a0+hbar**2/(2*ma)*((kx**2-ky**2)*cos(theta_a)+2*kx*ky*sin(theta_a) ) )**2
Ek=sqrt(Ek2)
epsilon_k_0=hbar**2/(2*m)*(kx**2+ky**2)+Ek
epsilon_k_1=hbar**2/(2*m)*(kx**2+ky**2)-Ek
fk0=1/(exp(beta*(epsilon_k_0-muF))+1)
fk1=1/(exp(beta*(epsilon_k_1-muF))+1)
#cos(theta_a)
cos_theta_a=cos(theta_a)
#sin(theta_a)
sin_theta_a=sin(theta_a)

sin_alpha_k=abs(vR)/Ek*k_abs

epsilon_L_k_a=a0+hbar**2/(2*ma)*((kx**2-ky**2)*cos(theta_a)+2*kx*ky*sin(theta_a))
cos_alpha_k=epsilon_L_k_a/Ek

half=Rational(1,2)
one_8=Rational(1,8)
three_2=Rational(3,2)

v_xk_0e_00=((1/m+cos_theta_a/ma)*hbar*kx+sin_theta_a/ma*hbar*ky)*(half+half*cos_alpha_k) \
           +1/hbar*vR*sin_alpha_k*cos_alpha_k \
           +((1/m-cos_theta_a/ma)*hbar*kx-sin_theta_a/ma*hbar*ky)*(half-half*cos_alpha_k)

v_xk_0e_01=hbar/ma*sin_alpha_k*(cos(theta_a)*kx+sin(theta_a)*ky) \
           -I*1/hbar*vR*sin_theta_k-1/hbar*vR*cos_alpha_k*cos_theta_k

v_xk_0e_10=hbar/ma*sin_alpha_k*(cos(theta_a)*kx+sin(theta_a)*ky)+I*1/hbar*vR*sin_theta_k \
           -1/hbar*vR*cos_alpha_k*cos_theta_k

v_xk_0e_11=((1/m+cos(theta_a)/ma)*hbar*kx+sin(theta_a)/ma*hbar*ky)*(half-half*cos_alpha_k) \
           -1/hbar*vR*sin_alpha_k*cos_theta_k \
           +((1/m-cos(theta_a)/ma)*hbar*kx-sin(theta_a)/ma*hbar*ky)*(half+half*cos_alpha_k)


v_yk_0e_00=(half+half*cos_alpha_k)*(sin(theta_a)/ma*hbar*kx+(1/m-cos(theta_a)/ma)*hbar*ky) \
           +1/hbar*vR*sin_alpha_k*sin_theta_k \
           +(-sin(theta_a)/ma*hbar*kx+(1/m+cos(theta_a)/ma)*hbar*ky)*(half-half*cos_alpha_k)

v_yk_0e_01=hbar/ma*sin_alpha_k*(sin(theta_a)*kx-cos(theta_a)*ky)+I*1/hbar*vR*cos_theta_k-1/hbar*vR*cos_alpha_k*sin_theta_k

v_yk_0e_10=hbar/ma*sin_alpha_k*(sin(theta_a)*kx-cos(theta_a)*ky)-I*1/hbar*vR*cos_theta_k-1/hbar*vR*cos_alpha_k*sin_theta_k

v_yk_0e_11=(sin(theta_a)/ma*hbar*kx+(1/m-cos(theta_a)/ma)*hbar*ky)*(half-half*cos_alpha_k)-1/hbar*vR*sin_alpha_k*sin_theta_k \
           +(-sin(theta_a)/ma*hbar*kx+(1/m+cos(theta_a)/ma)*hbar*ky)*(half+half*cos_alpha_k)


integrand_sigma=(fk0-fk1)*(v_yk_0e_10*v_xk_0e_01- v_yk_0e_01*v_xk_0e_10)/(epsilon_k_0-epsilon_k_1)**2



# Define the numerical values in a dictionary
# Note: Since ma = 1/eps, substituting eps will automatically evaluate ma
subs_dict = {
    hbar: 1.1,
    m: 0.9,
    eps: 1/20,
    vR: 2.0,
    a0: 0.5,
    muF: 1.5
}
def sigma_analytical(subs_dict):
    muF_val = subs_dict[muF]
    m_val = subs_dict[m]
    vR_val = subs_dict[vR]
    a0_val = subs_dict[a0]
    hbar_val = subs_dict[hbar]
    I0_minus_I1=-1j*a0_val*np.pi/hbar_val**2\
                *2*m_val*vR_val**2/(2*m_val*vR_val**2*muF_val+a0_val**2*hbar_val**2)
    sigma=I0_minus_I1*1j*e**2*subs_dict[hbar]/(4*np.pi**2)
    return sigma

# Substitute the numerical values into the integrand
e=-1
integrand_sigma_num = integrand_sigma.subs(subs_dict)*1j*e**2*subs_dict[hbar]/(4*np.pi**2)
# Simplify the result (optional, but recommended for cleaner expressions)
# integrand_sigma_num = simplify(integrand_sigma_num)

integrand_func = lambdify((kx, ky, theta_a, beta), integrand_sigma_num, modules='numpy')


def get_rho0(subs_dict):
    # Use the SymPy symbols (no quotes) to get the numerical values
    muF_val = subs_dict[muF]
    m_val = subs_dict[m]
    vR_val = subs_dict[vR]
    a0_val = subs_dict[a0]
    hbar_val = subs_dict[hbar]

    """Calculates rho_0 where g0(rho) = muF (Eq. 203)"""
    term1 = (muF_val * hbar_val**2 / m_val) + vR_val**2
    term2 = np.sqrt(vR_val**4 + (2 * muF_val * hbar_val**2 / m_val) * vR_val**2 + (a0_val**2 * hbar_val**4 / m_val**2))
    rho0_sq = (2 * m_val**2 / hbar_val**4) * (term1 - term2)
    return np.sqrt(max(rho0_sq, 0.0))


def get_rho1(subs_dict):
    # Extract numerical values here as well to avoid symbolic math errors
    muF_val = subs_dict[muF]
    m_val = subs_dict[m]
    vR_val = subs_dict[vR]
    a0_val = subs_dict[a0]
    hbar_val = subs_dict[hbar]

    """Calculates rho_1 where g1(rho) = muF"""
    term1 = (muF_val * hbar_val**2 / m_val) + vR_val**2
    term2 = np.sqrt(vR_val**4 + (2 * muF_val * hbar_val**2 / m_val) * vR_val**2 + (a0_val**2 * hbar_val**4 / m_val**2))
    rho1_sq = (2 * m_val**2 / hbar_val**4) * (term1 + term2)
    return np.sqrt(max(rho1_sq, 0.0))


rho0=get_rho0(subs_dict)
rho1=get_rho1(subs_dict)

beta_val = 100.0
theta_a_val = 0.5 * np.pi


# Determine the maximum radius for the integration domain
rho_max = max(rho0, rho1)
from scipy import integrate

# Define the integrands in polar coordinates (kx = r*cos(phi), ky = r*sin(phi))
# Note: We multiply by 'r' which is the Jacobian determinant for polar coordinates.
def polar_integrand_real(r, phi):
    kx_val = r * np.cos(phi)
    ky_val = r * np.sin(phi)
    val = integrand_func(kx_val, ky_val, theta_a_val, beta_val)
    return np.real(val) * r

def polar_integrand_imag(r, phi):
    kx_val = r * np.cos(phi)
    ky_val = r * np.sin(phi)
    val = integrand_func(kx_val, ky_val, theta_a_val, beta_val)
    return np.imag(val) * r


print(f"Starting integration over a circle of radius rho_max = {rho_max:.4f}...")

# scipy.integrate.dblquad signature: func(y, x), [x_min, x_max], [y_min(x), y_max(x)]
# Let x = phi (0 to 2*pi) and y = r (0 to rho_max)
integral_real, error_real = integrate.dblquad(
    polar_integrand_real,
    0, 2 * np.pi,                  # phi limits
    lambda phi: 0, lambda phi: rho_max  # r limits
)

integral_imag, error_imag = integrate.dblquad(
    polar_integrand_imag,
    0, 2 * np.pi,                  # phi limits
    lambda phi: 0, lambda phi: rho_max  # r limits
)

# Combine the results
total_integral = integral_real + 1j * integral_imag

print(f"Integration complete!")
print(f"Real part:      {integral_real:.8e} (Estimated error: {error_real:.2e})")
print(f"Imaginary part: {integral_imag:.8e} (Estimated error: {error_imag:.2e})")
print(f"Total Integral: {total_integral}")
print(sigma_analytical(subs_dict))