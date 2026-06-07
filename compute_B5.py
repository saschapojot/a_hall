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

chem0=muF-hbar**2*one_2m*rho0**2
x0up=rho0**4*chem0**2-a0**2*rho0**4
x0down=hbar**2/m*chem0**6+vR**2*chem0**5

chem1=hbar**2*one_2m*rho1**2-muF
x1up=rho1**4*chem1**2-a0**2*rho1**4
x1down=hbar**2/m*chem1**6-vR**2*chem1**5


term_up=x0up*x1down+x1up*x0down

term_down=x0down*x1down

#up
term_up_expanded=expand(term_up)
R0, R1 = symbols("R0,R1",cls=Symbol,positive=True)
term_up_R=expand(term_up_expanded.subs({rho0: sqrt(R0), rho1: sqrt(R1)}))

term_up_sym, _= symmetrize(term_up_R, R0, R1)
term_up_sym_together=together(term_up_sym)

val_sum=4*m/hbar**4*(muF*hbar**2+m*vR**2)
val_prod=4*m**2/hbar**4*(muF**2-a0**2)

term_up_sym_together_replaced=term_up_sym_together.subs({
    R0 + R1: val_sum,
    R0 * R1: val_prod
})

term_up_sym_together_replaced_num,term_up_sym_together_replaced_den=factor(term_up_sym_together_replaced).as_numer_denom()

term_up_num_factor_list=term_up_sym_together_replaced_num.as_ordered_factors()
term_up_num_prefactpr=Mul(*term_up_num_factor_list[0:4])
term_up_num_poly=Mul(*term_up_num_factor_list[4:])
B5_up_num_poly=term_up_num_poly
with open("B5_up_num_poly.txt","w") as fptr:
    fptr.write(latex(B5_up_num_poly))


#down
term_down_expanded=expand(term_down)
term_down_R=expand(term_down.subs({rho0: sqrt(R0), rho1: sqrt(R1)}))
term_down_sym,res=symmetrize(term_down_R,R0,R1)

term_down_sym_together=together(term_down_sym)
term_down_sym_together_replaced=term_down_sym_together.subs({
    R0 + R1: val_sum,
    R0 * R1: val_prod
})
term_down_sym_together_replaced_num,term_down_sym_together_replaced_den=factor(term_down_sym_together_replaced).as_numer_denom()

term_down_num_factor_list=term_down_sym_together_replaced_num.as_ordered_factors()
term_down_num_poly=Mul(*term_down_num_factor_list[:])
B5_down_num_poly=term_down_num_poly
with open("B5_down_num_poly.txt","w") as fptr:
    fptr.write(latex(B5_down_num_poly))


#compute partial fraction
ma_val=symbols("m_a", real=True)
frac=B5_up_num_poly/B5_down_num_poly*I*1/ma**2*pi*a0*m**4*vR**4
frac=frac.subs(eps,1/ma_val)
# 1. Define a dummy variable for the substitution
u = symbols("u", real=True)

# 2. The composite factor you want in the denominator
factor_expr = a0**2 * hbar**2 + 2 * m * muF * vR**2

# 3. Solve factor_expr = u for one of the variables, e.g., muF
# u = a0**2 * hbar**2 + 2 * m * muF * vR**2  =>  muF = (u - a0**2 * hbar**2) / (2 * m * vR**2)
muF_sub = (u - a0**2 * hbar**2) / (2 * m * vR**2)


# 4. Substitute muF with the u-expression in your fraction
frac_u = frac.subs(muF, muF_sub)
# Cancel out common terms to ensure it's a clean rational function in u
frac_u = cancel(frac_u)
# 5. Perform partial fraction decomposition with respect to the dummy variable u
partial_frac_u = apart(frac_u, u)

# 6. Substitute the original expression back in place of u
partial_frac_B5 = partial_frac_u.subs(u, factor_expr)
pprint(partial_frac_B5)
with open("partial_frac_B5.txt","w") as fptr:
    fptr.write(latex(partial_frac_B5))


