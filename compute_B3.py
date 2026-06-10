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

chem=hbar**2*m*muF+m**2*vR**2

term0_up=half*a0*hbar**2*m**2\
        * rho0**4*(hbar**2*one_2m*rho0**2-muF)

term0_down=(half*hbar**4*rho0**2-chem)*(a0**2+vR**2*rho0**2)**3

term1_up=half*a0*hbar**2*m**2\
        *rho1**4*(hbar**2*one_2m*rho1**2-muF)
term1_down=(half*hbar**4*rho1**2-chem)*(a0**2+vR**2*rho1**2)**3

B3_term_up=term0_up*term1_down+term1_up*term0_down
B3_term_down=term0_down*term1_down

#up
B3_term_up_expanded=expand(B3_term_up)
R0, R1 = symbols("R0,R1",cls=Symbol,positive=True)
B3_term_up_expanded_R=expand(B3_term_up_expanded.subs({rho0: sqrt(R0), rho1: sqrt(R1)}))
B3_term_up_sym,_=symmetrize(B3_term_up_expanded_R, R0, R1)
B3_term_up_sym_together=together(B3_term_up_sym)

val_sum=4*m/hbar**4*(muF*hbar**2+m*vR**2)
val_prod=4*m**2/hbar**4*(muF**2-a0**2)
B3_term_up_sym_together_replaced=B3_term_up_sym_together.subs({
    R0 + R1: val_sum,
    R0 * R1: val_prod
})
B3_term_up_sym_together_replaced_num,B3_term_up_sym_together_replaced_den=factor(B3_term_up_sym_together_replaced).as_numer_denom()

B3_term_up_num_factor_list=B3_term_up_sym_together_replaced_num.as_ordered_factors()
B3_term_up_num_prefactor=Mul(*B3_term_up_num_factor_list[0:4])
B3_term_up_num_poly=Mul(*B3_term_up_num_factor_list[4:])
with open("B3_term_up_num_poly.txt","w") as fptr:
    fptr.write(latex(B3_term_up_num_poly))

#down
B3_term_down_expanded=expand(B3_term_down)
B3_term_down_expanded_R=expand(B3_term_down_expanded.subs({rho0: sqrt(R0), rho1: sqrt(R1)}))
B3_term_down_sym,_=symmetrize(B3_term_down_expanded_R,R0,R1)
B3_term_down_together=together(B3_term_down_sym)
B3_term_down_together_replaced=B3_term_down_together.subs({
    R0 + R1: val_sum,
    R0 * R1: val_prod
})
B3_term_down_together_replaced_num,B3_term_down_together_replaced_den=factor(B3_term_down_together_replaced).as_numer_denom()
B3_term_down_factor_list=B3_term_down_together_replaced_num.as_ordered_factors()
B3_term_down_num_prefactor=Mul(*B3_term_down_factor_list[0:2])
B3_term_down_num_poly=Mul(*B3_term_down_factor_list[2:])
with open("B3_term_down_num_poly.txt","w") as fptr:
    fptr.write(latex(B3_term_down_num_poly))
