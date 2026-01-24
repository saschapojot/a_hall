import os
os.environ['TERM'] = 'xterm-256color'
from sympy import  *
from sympy.physics.units import hbar
from sympy.simplify.fu import TR5, TR11

eps_k=symbols("epsilon_k",cls=Symbol,real=True)
vR=symbols("v_R",cls=Symbol,real=True,nonzero=True)
eps_L=symbols("epsilon_L",cls=Symbol,real=True)
eps_k_a=symbols(r"epsilon_ka",cls=Symbol,real=True)
kx=symbols("k_x",cls=Symbol,real=True)
ky=symbols("k_y",cls=Symbol,real=True)
m=symbols("m",cls=Symbol,positive=True)
ma=symbols("m_a",cls=Symbol,positive=True)
theta_a=symbols("theta_a",cls=Symbol,positive=True)

sigma_0=Matrix(
 [[1,0],
  [0,1]]
)

sigma_x=Matrix(
    [[0,1],
     [1,0]]
)

sigma_y=Matrix([
    [0,-I],
    [I,0]
])

sigma_z=Matrix(
    [[1,0],
     [0,-1]]
)

#|k|
k_abs=sqrt(kx**2+ky**2)

#h0
h0=eps_k*sigma_0-vR*ky*sigma_x+vR*kx*sigma_y+(eps_L+eps_k_a)*sigma_z

#cos(theta_k)
cos_theta_k=kx/k_abs*sign(vR)

#sin(theta_k)
sin_theta_k=ky/k_abs*sign(vR)

#dcos(theta_k)/dkx
d_cos_theta_k_d_kx=diff(cos_theta_k,kx)

#dcos(theta_k)/dky
d_cos_theta_k_d_ky=diff(cos_theta_k,ky)

#dsin(theta_k)/dkx
d_sin_theta_k_d_kx=diff(sin_theta_k,kx)

#dsin(theta_k)/dky
d_sin_theta_k_d_ky=diff(sin_theta_k,ky)

#theta_k computed from arccos
theta_k=acos(cos_theta_k)# in [0,pi]

#cos(2theta_k)
cos_2_theta_k=2*cos_theta_k**2-1

#sin(2theta_k)
sin_2_theta_k=2*sin_theta_k*cos_theta_k

#dcos(2theta_k)/dkx
d_cos_2_theta_k_d_kx=diff(cos_2_theta_k,kx)

#dcos(2theta_k)/dky
d_cos_2_theta_k_d_ky=diff(cos_2_theta_k,ky)

#dsin(2theta_k)/dkx
d_sin_2_theta_k_d_kx=diff(sin_2_theta_k,kx)

#dsin(2theta_k)/dky
d_sin_2_theta_k_d_ky=diff(sin_2_theta_k,ky)

#d|k|/dkx
d_k_abs_d_kx=diff(k_abs,kx)

#d|k|/dky
d_k_abs_d_ky=diff(k_abs,ky)

#eps_k's value
eps_k_val=hbar**2/(2*m)*k_abs**2

#d eps_k/dkx
d_eps_k_d_kx=diff(eps_k_val,kx)

#d eps_k/dky
d_eps_k_d_ky=diff(eps_k_val,ky)

#cos(theta_a)
cos_theta_a=cos(theta_a)

#sin(theta_a)
sin_theta_a=sin(theta_a)

#cos(2theta_k-theta_a)
cos_2theta_k_minus_theta_a=cos_2_theta_k*cos_theta_a+sin_2_theta_k*sin_theta_a

#d cos(2theta_k-theta_a)/dkx
d_cos_2theta_k_minus_theta_a_dkx=diff(cos_2theta_k_minus_theta_a,kx)

#d cos(2theta_k-theta_a)/dky
d_cos_2theta_k_minus_theta_a_dky=diff(cos_2theta_k_minus_theta_a,ky)

#cos(2theta_k-theta_a)

cos_2theta_k_minus_theta_a=cos_2_theta_k*cos_theta_a+sin_2_theta_k*sin_theta_a

#value of epsilon_{k}^{a}
eps_k_a_val=hbar**2*k_abs**2/(2*ma)*cos_2theta_k_minus_theta_a

# d epsilon_ka/dkx
d_eps_k_a_d_kx=diff(eps_k_a_val,kx)

# d epsilon_ka/dky
d_eps_k_a_d_ky=diff(eps_k_a_val,ky)

# v_{x,k}^{0}
v_xk_0=1/hbar*diff(eps_k_val,kx)*sigma_0+1/hbar*vR*sigma_y+1/hbar*diff(eps_k_a_val,kx)*sigma_z

#v_{y,k}^{0}
v_yk_0=1/hbar*diff(eps_k_val,ky)*sigma_0-1/hbar*vR*sigma_x+1/hbar*diff(eps_k_a_val,ky)*sigma_z

# Ek
# 对根式内部进行通分
Ek_inner = vR**2*k_abs**2+(eps_L+eps_k_a_val)**2
Ek_inner_together = expand(together(Ek_inner))
Ek = sqrt(Ek_inner_together)

#cos(alpha_k)
cos_alpha_k=(eps_L+eps_k_a_val)/Ek
h0_val=eps_k_val*sigma_0-vR*ky*sigma_x+vR*kx*sigma_y+(eps_L+eps_k_a_val)*sigma_z

