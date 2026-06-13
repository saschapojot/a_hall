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

B2n=(hbar**2*one_2m*rho0**2-muF)**5+(hbar**2*one_2m*rho1**2-muF)**5

#symmetrize B2n
B2n_expanded=expand(B2n)
R0, R1 = symbols("R0,R1",cls=Symbol,positive=True)
B2n_expanded_R=expand(B2n_expanded.subs({rho0: sqrt(R0), rho1: sqrt(R1)}))
B2n_expanded_R_sym,_=symmetrize(B2n_expanded_R, R0, R1)
B2n_expanded_R_sym_together=together(B2n_expanded_R_sym)

val_sum=4*m/hbar**4*(muF*hbar**2+m*vR**2)
val_prod=4*m**2/hbar**4*(muF**2-a0**2)
B2n_expanded_R_sym_together_replaced=B2n_expanded_R_sym_together.subs({
    R0 + R1: val_sum,
    R0 * R1: val_prod
})



B2=-I*1/ma**2*Rational(15,8)*pi*a0**3*hbar**2*m/vR**2*1/chem\
    +I*1/ma**2*Rational(5,4)*pi*a0**5*hbar**2*m/vR**2*(6*m*muF*hbar**2*vR**2+4*m**2*vR**4+3*a0**2*hbar**4)/chem**3\
    -I*1/ma**2*Rational(3,16)*pi*a0**7*hbar**12/vR**4*B2n_expanded_R_sym_together_replaced/chem**5


# --- Partial Fraction Decomposition ---
# 1. Define a dummy variable for the substitution
u = symbols("u", real=True)
ma_val=symbols("m_a", real=True)
B2_frac=B2.subs(eps,1/ma_val)

# 2. The composite factor you want in the denominator
factor_expr = a0**2 * hbar**2 + 2 * m * muF * vR**2
# 3. Solve factor_expr = u for muF
muF_sub = (u - a0**2 * hbar**2) / (2 * m * vR**2)
# 4. Substitute muF with the u-expression in your fraction
B2_frac_u = B2_frac.subs(muF, muF_sub)
# Cancel out common terms to ensure it's a clean rational function in u
B2_frac_u = cancel(B2_frac_u)


# 5. Perform partial fraction decomposition with respect to the dummy variable u
B2_partial_frac_u = apart(B2_frac_u, u)

B2_partial_frac=B2_partial_frac_u.subs(u, factor_expr)

with open("B2_partial_frac.txt","w") as fptr:
    fptr.write(latex(B2_partial_frac))