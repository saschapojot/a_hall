import os
os.environ['TERM'] = 'xterm-256color'
from sympy import  *
from sympy.physics.units import hbar
eps_k=symbols("epsilon_k",cls=Symbol,real=True)
vR=symbols("v_R",cls=Symbol,real=True,nonzero=True)
eps_L=symbols("epsilon_L",cls=Symbol,real=True)
eps_k_a=symbols(r"epsilon_ka",cls=Symbol,real=True)
kx=symbols("k_x",cls=Symbol,real=True)
ky=symbols("k_y",cls=Symbol,real=True)
m=symbols("m",cls=Symbol,positive=True)
ma=symbols("m_a",cls=Symbol,positive=True)

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

k_abs=sqrt(kx**2+ky**2)

h0=eps_k*sigma_0-vR*ky*sigma_x+vR*kx*sigma_y+(eps_L+eps_k_a)*sigma_z

cos_theta_k=kx/k_abs*sign(vR)

sin_theta_k=ky/k_abs*sign(vR)

d_cos_theta_k_d_kx=diff(cos_theta_k,kx)

d_cos_theta_k_d_ky=diff(cos_theta_k,ky)

d_sin_theta_k_d_kx=diff(sin_theta_k,kx)

d_sin_theta_k_d_ky=diff(sin_theta_k,ky)

theta_k=acos(cos_theta_k)# in [0,pi]
cos_2_theta_k=2*cos_theta_k**2-1
# sin_2_theta_k=sin(2*theta_k)
sin_2_theta_k=2*sin_theta_k*cos_theta_k
d_cos_2_theta_k_d_kx=diff(cos_2_theta_k,kx)

d_cos_2_theta_k_d_ky=diff(cos_2_theta_k,ky)

d_sin_2_theta_k_d_kx=diff(sin_2_theta_k,kx)

d_sin_2_theta_k_d_ky=diff(sin_2_theta_k,ky)

d_k_abs_d_kx=diff(k_abs,kx)
d_k_abs_d_ky=diff(k_abs,ky)

eps_k_val=hbar**2/(2*m)*k_abs**2

d_eps_k_d_kx=diff(eps_k_val,kx)

d_eps_k_d_ky=diff(eps_k_val,ky)

rhs=hbar**2/(m)*ky
df=d_eps_k_d_ky-rhs

pprint(df)