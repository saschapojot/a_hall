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
# 构造第一项 (涉及 rho1)
B4b_term1_num=hbar**2/m*rho1**2*(hbar**2*one_2m*rho1**2-muF)**3\
    +vR**4*rho1**4-vR**2*rho1**2*(hbar**2*one_2m*rho1**2-muF)**2

B4b_term1_den=(hbar**2/m*(hbar**2*one_2m*rho1**2-muF)-vR**2 )**3\
    *(hbar**2*one_2m*rho1**2-muF)**5

# 构造第二项 (涉及 rho0)
B4b_term2_num=hbar**2/m*rho0**2*(muF-hbar**2*one_2m*rho0**2)**3\
    +vR**2*rho0**2*(muF-hbar**2*one_2m*rho0**2)**2\
    -vR**4*rho0**4
B4b_term2_den=(hbar**2/m*(muF-hbar**2*one_2m*rho0**2)+vR**2)**3\
    *(muF-hbar**2*one_2m*rho0**2)**5

B4b_up=B4b_term1_num*B4b_term2_den-B4b_term2_num*B4b_term1_den

B4b_down=B4b_term1_den*B4b_term2_den

## B4b up
B4b_up_expanded=expand(B4b_up)
R0, R1 = symbols("R0,R1",cls=Symbol,positive=True)
B4b_up_R=expand(B4b_up_expanded.subs({rho0: sqrt(R0), rho1: sqrt(R1)}))
B4b_up_sym, _ = symmetrize(B4b_up_R, R0, R1)
B4b_up_sym_together = together(B4b_up_sym)
val_sum=4*m/hbar**4*(muF*hbar**2+m*vR**2)
val_prod=4*m**2/hbar**4*(muF**2-a0**2)

B4b_up_sym_together_replaced=B4b_up_sym_together.subs({
    R0 + R1: val_sum,
    R0 * R1: val_prod
})

B4b_up_sym_together_replaced_num, B4b_up_sym_together_replaced_den=factor(B4b_up_sym_together_replaced).as_numer_denom()

B4b_up_num_factors_list=B4b_up_sym_together_replaced_num.as_ordered_factors()

B4b_up_num_prefactor=Mul(*B4b_up_num_factors_list[0:2])
B4b_up_num_poly = Mul(*B4b_up_num_factors_list[2:])
with open("B4b_up_num_poly.txt","w") as fptr:
    fptr.write(latex(B4b_up_num_poly))

#B4b down
B4b_down_expanded=expand(B4b_down)
B4b_down_R =expand(B4b_down_expanded.subs({rho0: sqrt(R0), rho1: sqrt(R1)}))
B4b_down_sym, _ = symmetrize(B4b_down_R, R0, R1)
B4b_down_sym_together = together(B4b_down_sym)
B4b_down_sym_together_replaced=B4b_down_sym_together.subs({
    R0 + R1: val_sum,
    R0 * R1: val_prod
})

B4b_down_sym_together_replaced_num, B4b_down_sym_together_replaced_den=factor(B4b_down_sym_together_replaced).as_numer_denom()
B4a_down_num_poly=B4b_down_sym_together_replaced_num
with open("B4a_down_num_poly.txt","w") as fptr:
    fptr.write(latex(B4a_down_num_poly))


B4b_up=1/B4b_up_sym_together_replaced_den*B4b_up_num_prefactor*B4b_up_num_poly


B4b=a0**4*m**3*hbar**4*B4b_up_num_poly/B4a_down_num_poly

# --- Partial Fraction Decomposition ---
# 1. Define a dummy variable for the substitution
u = symbols("u", real=True)
ma_val=symbols("m_a", real=True)
# 2. The composite factor you want in the denominator
factor_expr = a0**2 * hbar**2 + 2 * m * muF * vR**2

# 3. Solve factor_expr = u for muF
muF_sub = (u - a0**2 * hbar**2) / (2 * m * vR**2)
B4b_with_factor=I*1/ma**2*Rational(1,4)*pi*a0*vR**2/hbar**2*B4b
B4b_with_factor=B4b_with_factor.subs(ma,ma_val)
# 4. Substitute muF with the u-expression in your fraction
B4b_with_factor_u = B4b_with_factor.subs(muF, muF_sub)
# Cancel out common terms to ensure it's a clean rational function in u
B4b_with_factor_u = cancel(B4b_with_factor_u)
# 5. Perform partial fraction decomposition with respect to the dummy variable u.
# This will automatically separate terms with `u` and `u * hbar**2 + m**2 * vR**4`
partial_frac_u = apart(B4b_with_factor_u, u)

# 6. Substitute the original expression back in place of u
partial_frac_B4b_full = partial_frac_u.subs(u, factor_expr)
with open("partial_frac_B4b_full.txt","w") as fptr:
    fptr.write(latex(partial_frac_B4b_full))
