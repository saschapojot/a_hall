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
Ek=symbols("Ek",cls=Symbol,positive=True)
a0=symbols("a0",cls=Symbol,positive=True)
rho=symbols("rho",cls=Symbol,positive=True)
gamma=symbols("gamma",cls=Symbol,real=True)
a00,a01,a10,a11=symbols("a00,a01,a10,a11",cls=Symbol,real=True)
eps=symbols("epsilon",cls=Symbol,positive=True)
beta=symbols("beta",cls=Symbol,positive=True)
ma=1/eps

muF = symbols("mu_F", cls=Symbol, real=True)
rho0 = symbols("rho_0", cls=Symbol, positive=True)
rho1 = symbols("rho_1", cls=Symbol, positive=True)

half=Rational(1,2)
one_2m=1/(2*m)

chem=2*m*vR**2*muF+a0**2*hbar**2

B0=-I*1/ma**2*Rational(3,8)*m*pi*a0/vR**2\
    +I*1/ma**2*Rational(3,4)*m*pi*a0**3*hbar**2/vR**2*1/chem\
    -I*1/ma**2*Rational(1,8)*m*pi*a0**5*hbar**2/vR**2\
   *(6*m*muF*vR**2*hbar**2+4*m**2*vR**4+3*a0**2*hbar**4)/chem**3



# --- Partial Fraction Decomposition ---
# 1. Define a dummy variable for the substitution
u = symbols("u", real=True)
ma_val=symbols("m_a", real=True)
B0_frac=B0.subs(eps,1/ma_val)

# 2. The composite factor you want in the denominator
factor_expr = a0**2 * hbar**2 + 2 * m * muF * vR**2
# 3. Solve factor_expr = u for muF
muF_sub = (u - a0**2 * hbar**2) / (2 * m * vR**2)
# 4. Substitute muF with the u-expression in your fraction
B0_frac_u = B0_frac.subs(muF, muF_sub)
# Cancel out common terms to ensure it's a clean rational function in u
B0_frac_u = cancel(B0_frac_u)


# 5. Perform partial fraction decomposition with respect to the dummy variable u
B0_partial_frac_u = apart(B0_frac_u, u)

B0_partial_frac=B0_partial_frac_u.subs(u, factor_expr)

with open("B0_partial_frac.txt","w") as fptr:
    fptr.write(latex(B0_partial_frac))
