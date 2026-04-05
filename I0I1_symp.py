from sympy import *
from sympy.physics.units import hbar


m=symbols("m",cls=Symbol,positive=True)
vR=symbols("v_R",cls=Symbol,real=True,nonzero=True)

half=Rational(1,2)
muF=symbols("mu_F",cls=Symbol,real=True,nonzero=True)
a0=symbols("a0",cls=Symbol,real=True,nonzero=True)


rho0rho1_sum=4*m*muF/hbar**2+4*m**2/hbar**4*vR**2

rho0rho1_prod=4*m**2/hbar**4*(muF**2-a0**2)

rho0=sqrt(2*m**2/hbar**4*(muF*hbar**2/m+vR**2-sqrt(vR**4+2*muF*hbar**2/m*vR**2+a0**2*hbar**4/m**2 ) ) )

rho1=sqrt(2*m**2/hbar**4*(muF*hbar**2/m+vR**2+sqrt(vR**4+2*muF*hbar**2/m*vR**2+a0**2*hbar**4/m**2 ) ) )


lhs=(vR**2*rho0**2+a0**2)**(-3*half)-(vR**2*rho1**2+a0**2)**(-3*half)

rhs=2*m/hbar**2*vR**2*(6*m*muF/hbar**2*vR**2+4*m**2/hbar**4*vR**4+3*a0**2 )/(2*m/hbar**2*muF*vR**2+a0**2)**3
df=lhs-rhs
# Define the substitution dictionary
values = {
    hbar: 1.0,
    m: 1.0,
    muF: 1.5,
    vR: 2.0,
    a0: 0.5
}

# Substitute the values into df and evaluate to a floating-point number
df_evaluated = df.subs(values).evalf()



print("\nEvaluated df with substituted values:")
print(df_evaluated)