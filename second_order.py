from sympy import *
from sympy.physics.units import hbar

# Define symbols
b0,b1,b2=symbols("b0,b1,b2",cls=Symbol,positive=True)
eps=symbols("epsilon",cls=Symbol,real=True)
a=symbols("a",cls=Symbol,real=True)
alpha_k=symbols("alpha_k",cls=Symbol,real=True)
theta_k=symbols("theta_k",cls=Symbol,real=True)
kx=symbols("k_x",cls=Symbol,real=True)
ky=symbols("k_y",cls=Symbol,real=True)
m=symbols("m",cls=Symbol,positive=True)
ma=symbols("m_a",cls=Symbol,positive=True)
theta_a=symbols("theta_a",cls=Symbol,positive=True)
vR=symbols("v_R",cls=Symbol,positive=True,nonzero=True)
Ek=symbols("Ek",cls=Symbol,positive=True)
a0=symbols("a0",cls=Symbol,positive=True)
rho=symbols("rho",cls=Symbol,positive=True)
gamma=symbols("gamma",cls=Symbol,real=True)
a00,a01,a10,a11=symbols("a00,a01,a10,a11",cls=Symbol,real=True)

# Note: This overwrites the previous definition of eps.
# The 'positive=True' assumption here is the root cause of the series error later.
eps=symbols("epsilon",cls=Symbol,positive=True)

beta=symbols("beta",cls=Symbol,positive=True)
e=symbols("e",cls=Symbol,real=True)

# Define ma in terms of eps
ma=1/eps

# hbar = 1.0
# m = 1.0
# ma = 20.0
# vR = 2.0
# a0 = 0.5
# theta_a = 0.1 * np.pi
# eps=1/ma
A=Matrix([[a00,a01],[a10,a11]])

U=Matrix(
    [
        [cos(alpha_k/2)*exp(-I*theta_k),sin(alpha_k/2)*exp(-I*theta_k)],
        [I*sin(alpha_k/2) ,-I*cos(alpha_k/2)]
    ]
)
k_abs=sqrt(kx**2+ky**2)
cos_theta_k=kx/k_abs*sign(vR)
sin_theta_k=ky/k_abs*sign(vR)

#cos(theta_a)
cos_theta_a=cos(theta_a)
#sin(theta_a)
sin_theta_a=sin(theta_a)
#sin(alpha_k)
sin_alpha_k=abs(vR)/Ek*k_abs
cos_alpha_k=cos(alpha_k)

half=Rational(1,2)
one_8=Rational(1,8)
three_2=Rational(3,2)


v_xk_0e_00=((1/m+cos_theta_a/ma)*hbar*kx+sin_theta_a/ma*hbar*ky)*(half+half*cos_alpha_k) \
           +1/hbar*sin_alpha_k*cos_alpha_k \
           +((1/m-cos_theta_a/ma)*hbar*kx-sin_theta_a/ma*hbar*ky)*(half-half*cos_alpha_k)


v_xk_0e_01=half*((1/m+cos_theta_a/ma)*hbar*kx+sin_theta_a/ma*hbar*ky)*sin_alpha_k \
           -half*((1/m-cos_theta_a/ma)*hbar*kx-sin_theta_a/ma*hbar*ky)*sin_alpha_k \
           -I*1/hbar*vR*sin_theta_k-1/hbar*vR*cos_alpha_k*cos_theta_k

v_xk_0e_10=half*((1/m+cos_theta_a/ma)*hbar*kx+sin_theta_a/ma*hbar*ky)*sin_alpha_k \
           -half*((1/m-cos_theta_a/ma)*hbar*kx-sin_theta_a/ma*hbar*ky)*sin_alpha_k \
           +I*1/hbar*vR*sin_theta_k-1/hbar*vR*cos_alpha_k*cos_theta_k

v_xk_0e_11=((1/m+cos_theta_a/ma)*hbar*kx+sin_theta_a/ma*hbar*ky)*(half-half*cos_alpha_k) \
           -1/hbar*vR*sin_alpha_k*cos_theta_k \
           +((1/m-cos_theta_a/ma)*hbar*kx-sin_theta_a/ma*hbar*ky)*(half+half*cos_alpha_k)


v_xk_0e_mat=Matrix([[v_xk_0e_00,v_xk_0e_01],[v_xk_0e_10,v_xk_0e_11]])


v_yk_0e_00=(half+half*cos_alpha_k)*(sin(theta_a)/ma*hbar*kx+(1/m-cos(theta_a)/ma)*hbar*ky) \
           +1/hbar*vR*sin_alpha_k*sin_theta_k \
           +(-sin(theta_a)/ma*hbar*kx+(1/m+cos(theta_a)/ma)*hbar*ky)*(half-half*cos_alpha_k)

v_yk_0e_01=hbar/ma*sin_alpha_k*(sin(theta_a)*kx-cos(theta_a)*ky)+I*1/hbar*vR*cos_theta_k-1/hbar*vR*cos_alpha_k*sin_theta_k


v_yk_0e_10=hbar/ma*sin_alpha_k*(sin(theta_a)*kx-cos(theta_a)*ky)-I*1/hbar*vR*cos_theta_k-1/hbar*vR*cos_alpha_k*sin_theta_k

v_yk_0e_11=(sin(theta_a)/ma*hbar*kx+(1/m-cos(theta_a)/ma)*hbar*ky)*(half-half*cos_alpha_k)-1/hbar*vR*sin_alpha_k*sin_theta_k \
           +(-sin(theta_a)/ma*hbar*kx+(1/m+cos(theta_a)/ma)*hbar*ky)*(half+half*cos_alpha_k)




half=Rational(1,2)


# Define the function g1
g1 = hbar**2/(2*m)*rho**2 - sqrt(a0**2 + vR**2*rho**2)



c1=-a0*hbar**2/(2*ma)*rho**2*(a0**2+vR**2*rho**2)**(-half)
n=0

G1=g1+c1*cos(2*gamma-theta_a)

drho2_G1=diff(G1,(rho,2))

drho_c1=diff(c1,rho)

drho2_c1=diff(c1,(rho,2))


rho10=sqrt(m**2*vR**4-a0**2*hbar**4)/(abs(vR)*hbar**2)
# rho10=symbols("rho10",cls=Symbol,real=True)

# rho11=symbols("rho11",cls=Symbol,real=True)
rho11=a0/(2*vR**3)*(m**2*vR**4+a0**2*hbar**4)/sqrt(m**2*vR**4-a0**2*hbar**4)*sign(vR)*(-1)**n
rho1=rho10+eps*rho11
# f=a0**2+vR**2*rho1**2
Gamma0=Matrix([
    [hbar**6/(m**3*vR**2)*rho10**2,0],
    [0,0]
])
# Gamma1=Matrix([
#     [3*hbar**6/(m**3*vR**2)*rho10*rho11-3*hbar**10/(m**5*vR**4)*rho10**3*rho11+(a0**3*hbar**8/(2*m**3*vR**6)-3*a0**5*hbar**12/(2*m**5*vR**10))*(-1)**n ,0 ],
#     [0,2*a0*hbar**4/(m*vR**2)*rho10**2*(-1)**n]
# ])

Gamma1=Matrix([
    [0 ,0 ],
    [0,2*a0*hbar**4/(m*vR**2)*rho10**2*(-1)**n]
])
Gamma=Gamma0+1/ma*Gamma1


leading_Ek=vR**2*rho**2+a0**2
g0=hbar**2/(2*m)*rho**2+(a0**2+vR**2*rho**2)**half

c0=a0*hbar**2/2*rho**2/sqrt(a0**2+vR**2*rho**2)

d0=one_8*hbar**4*(
        rho**4*leading_Ek**(-half)-a0**2*rho**4*leading_Ek**(-3*half)
)



Ek2=vR**2*rho**2+(a0+hbar**2/(2*ma)*rho**2*cos(2*gamma-theta_a) )**2

Ek=sqrt(Ek2)


G0=hbar**2/(2*m)*rho**2+Ek

G1=hbar**2/(2*m)*rho**2-Ek


g1=hbar**2/(2*m)*rho**2-leading_Ek**half
c1=-half*a0*hbar**2*rho**2*leading_Ek**(-half)
d1=-one_8*hbar**4*(
        rho**4*leading_Ek**(-half)-a0**2*rho**4*leading_Ek**(-3*half)
)

half=Rational(1,2)


y,muF=symbols("y,mu_F",cls=Symbol,positive=True)
y=1/ma*c0*cos(2*gamma-theta_a)+1/ma**2*d0*(cos(2*gamma-theta_a))**2
z=1/ma*c1*cos(2*gamma-theta_a)+1/ma**2*d1*(cos(2*gamma-theta_a))**2



leading=beta*g0-beta*muF

w00=1/(exp(leading)+1)

w01=exp(leading)/(exp(leading)+1)**2

w02=(exp(2*leading)-exp(leading))/(exp(leading)+1)**3

leading1=beta*g1-beta*muF
w10=1/(exp(leading1)+1)
w11=exp(leading)/(exp(leading)+1)**2
w12=(exp(2*leading)-exp(leading))/(exp(leading)+1)**3

f1=1/(exp(beta*g1))
rhs=w10-1/ma*beta*w11*c1*cos(2*gamma-theta_a)+1/ma**2*(half*beta**2*w12*c1**2-beta*w11*d1)*(cos(2*gamma-theta_a))**2



df=f1-rhs

# --- FIX START ---
# Use a dummy variable for series expansion to avoid NotImplementedError
# caused by positive=True assumption on eps in nested exponents.
eps_dummy = symbols("eps_dummy", real=True)

# Substitute eps with the dummy variable
df_dummy = df.subs(eps, eps_dummy)

# Perform the series expansion on the dummy variable
val_dummy = df_dummy.series(eps_dummy, 0, 2)

# Substitute the dummy variable back to eps for the final result
val = val_dummy.subs(eps_dummy, eps)
# --- FIX END ---

pprint(simplify(val))