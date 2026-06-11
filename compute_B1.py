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

term0_up=half*a0*hbar**2*rho0**4
chem0=muF-hbar**2*one_2m*rho0**2
term0_down=hbar**2/m*chem0**4+vR**2*chem0**3

term1_up=half*a0*hbar**2*rho1**4
chem1=hbar**2*one_2m*rho1**2-muF
term1_down=hbar**2/m*chem1**4-vR**2*chem1**3

B1_term_up=term0_up*term1_down+term1_up*term0_down
B1_term_down=term0_down*term1_down

#up
B1_term_up_expanded=expand(B1_term_up)
R0, R1 = symbols("R0,R1",cls=Symbol,positive=True)
B1_term_up_expanded_R=expand(B1_term_up_expanded.subs({rho0: sqrt(R0), rho1: sqrt(R1)}))
B1_term_up_sym,_=symmetrize(B1_term_up_expanded_R, R0, R1)
B1_term_up_sym_together=together(B1_term_up_sym)
val_sum=4*m/hbar**4*(muF*hbar**2+m*vR**2)
val_prod=4*m**2/hbar**4*(muF**2-a0**2)
B1_term_up_sym_together_replaced=B1_term_up_sym_together.subs({
    R0 + R1: val_sum,
    R0 * R1: val_prod
})
B1_term_up_sym_together_replaced_num,B1_term_up_sym_together_replaced_den=factor(B1_term_up_sym_together_replaced).as_numer_denom()
B1_term_up_num_factor_list=B1_term_up_sym_together_replaced_num.as_ordered_factors()
B1_term_up_num_prefactor=Mul(*B1_term_up_num_factor_list[0:3])
B1_term_up_num_poly=Mul(*B1_term_up_num_factor_list[3:])
with open("B1_term_up_num_poly.txt","w") as fptr:
    fptr.write(latex(B1_term_up_num_poly))

#down
B1_term_down_expanded=expand(B1_term_down)
B1_term_down_expanded_R=expand(B1_term_down_expanded.subs({rho0: sqrt(R0), rho1: sqrt(R1)}))
B1_term_down_sym,_=symmetrize(B1_term_down_expanded_R,R0,R1)
B1_term_down_together=together(B1_term_down_sym)
B1_term_down_together_replaced=B1_term_down_together.subs({
    R0 + R1: val_sum,
    R0 * R1: val_prod
})
B1_term_down_together_replaced_num,B1_term_down_together_replaced_den=factor(B1_term_down_together_replaced).as_numer_denom()
B1_term_down_factor_list=B1_term_down_together_replaced_num.as_ordered_factors()
B1_term_down_num_poly=Mul(*B1_term_down_factor_list[:])
with open("B1_term_down_num_poly.txt","w") as fptr:
    fptr.write(latex(B1_term_down_num_poly))


B1=I*1/ma**2*pi*vR**2*a0*m**3*B1_term_up_num_poly/B1_term_down_num_poly

# --- Partial Fraction Decomposition ---
# 1. Define a dummy variable for the substitution
u = symbols("u", real=True)
ma_val=symbols("m_a", real=True)
B1_frac=B1.subs(eps,1/ma_val)
# 2. The composite factor you want in the denominator
factor_expr = a0**2 * hbar**2 + 2 * m * muF * vR**2
# 3. Solve factor_expr = u for muF
muF_sub = (u - a0**2 * hbar**2) / (2 * m * vR**2)
# 4. Substitute muF with the u-expression in your fraction
B1_frac_u = B1_frac.subs(muF, muF_sub)
# Cancel out common terms to ensure it's a clean rational function in u
B1_frac_u = cancel(B1_frac_u)

# 5. Perform partial fraction decomposition with respect to the dummy variable u
B1_partial_frac_u = apart(B1_frac_u, u)

B1_partial_frac=B1_partial_frac_u.subs(u, factor_expr)

with open("B1_partial_frac.txt","w") as fptr:
    fptr.write(latex(B1_partial_frac))