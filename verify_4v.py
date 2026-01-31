from sympy import *
from sympy import expand_complex
from sympy.simplify.fu import TR5, TR11, TR9,fu
from sympy.physics.units import hbar


alpha_k=symbols("alpha_k",cls=Symbol,real=True)
theta_k=symbols("theta_k",cls=Symbol,real=True)
kx=symbols("k_x",cls=Symbol,real=True)
ky=symbols("k_y",cls=Symbol,real=True)
m=symbols("m",cls=Symbol,positive=True)
ma=symbols("m_a",cls=Symbol,positive=True)
theta_a=symbols("theta_a",cls=Symbol,positive=True)
vR=symbols("v_R",cls=Symbol,real=True,nonzero=True)

a00,a01,a10,a11=symbols("a00,a01,a10,a11",cls=Symbol,real=True)

A=Matrix([[a00,a01],[a10,a11]])

U=Matrix(
    [
        [cos(alpha_k/2)*exp(-I*theta_k),sin(alpha_k/2)*exp(-I*theta_k)],
        [I*sin(alpha_k/2) ,-I*cos(alpha_k/2)]
    ]
)
k_abs=sqrt(kx**2+ky**2)
cos_theta_k=kx/k_abs
sin_theta_k=ky/k_abs

half=Rational(1,2)




v_xk_0e_00=((1/m+cos(theta_a)/ma)*hbar*kx+sin(theta_a)/ma*hbar*ky)*(half+half*cos(alpha_k))\
    +1/hbar*vR*sin(alpha_k)*cos_theta_k\
    +((1/m-cos(theta_a)/ma)*hbar*kx-sin(theta_a)/ma*hbar*ky)*(half-half*cos(alpha_k))


v_xk_0e_01=hbar/ma*sin(alpha_k)*(cos(theta_a)*kx+sin(theta_a)*ky)\
           -I*1/hbar*vR*sin_theta_k-1/hbar*vR*cos(alpha_k)*cos_theta_k

v_xk_0e_10=hbar/ma*sin(alpha_k)*(cos(theta_a)*kx+sin(theta_a)*ky)+I*1/hbar*vR*sin_theta_k\
    -1/hbar*vR*cos(alpha_k)*cos_theta_k


v_xk_0e_11=((1/m+cos(theta_a)/ma)*hbar*kx+sin(theta_a)/ma*hbar*ky)*(half-half*cos(alpha_k))\
    -1/hbar*vR*sin(alpha_k)*cos_theta_k\
    +((1/m-cos(theta_a)/ma)*hbar*kx-sin(theta_a)/ma*hbar*ky)*(half+half*cos(alpha_k))

v_xk_0e_mat=Matrix([[v_xk_0e_00,v_xk_0e_01],[v_xk_0e_10,v_xk_0e_11]])



v_yk_0=Matrix([
    [sin(theta_a)/ma*hbar*kx+(1/m-cos(theta_a)/ma)*hbar*ky,-1/hbar*vR],
    [-1/hbar*vR,-sin(theta_a)/ma*hbar*kx+(1/m+cos(theta_a)/ma)*hbar*ky]
])


v_yk_0e_00=(half+half*cos(alpha_k))*(sin(theta_a)/ma*hbar*kx+(1/m-cos(theta_a)/ma)*hbar*ky)\
    +1/hbar*vR*sin(alpha_k)*sin_theta_k\
    +(-sin(theta_a)/ma*hbar*kx+(1/m+cos(theta_a)/ma)*hbar*ky)*(half-half*cos(alpha_k))

v_yk_0e_01=hbar/ma*sin(alpha_k)*(sin(theta_a)*kx-cos(theta_a)*ky)+I*1/hbar*vR*cos_theta_k-1/hbar*vR*cos(alpha_k)*sin_theta_k


v_yk_0e_10=hbar/ma*sin(alpha_k)*(sin(theta_a)*kx-cos(theta_a)*ky)-I*1/hbar*vR*cos_theta_k\
           -1/hbar*vR*cos(alpha_k)*sin_theta_k

v_yk_0e_11=(sin(theta_a)/ma*hbar*kx+(1/m-cos(theta_a)/ma)*hbar*ky)*(half-half*cos(alpha_k))-1/hbar*vR*sin(alpha_k)*sin_theta_k\
    +(-sin(theta_a)/ma*hbar*kx+(1/m+cos(theta_a)/ma)*hbar*ky)*(half+half*cos(alpha_k))

v_yk_0e_mat=Matrix([[v_yk_0e_00,v_yk_0e_01],[v_yk_0e_10,v_yk_0e_11]])




lhs=v_yk_0e_10*v_xk_0e_01-v_yk_0e_01*v_xk_0e_10

rhs=-I*4*vR/ma*sin(theta_a)*sin(alpha_k)*kx*ky/k_abs\
    +I*2*vR/ma*cos(theta_a)*sin(alpha_k)*(ky**2-kx**2)/k_abs\
    +I*2*vR**2/hbar**2*cos(alpha_k)



df=fu(lhs-rhs)

pprint(simplify(expand(df)))