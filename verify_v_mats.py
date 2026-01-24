from sympy import *
from sympy.simplify.fu import TR5, TR11, TR9,fu

alpha_k=symbols("alpha_k",cls=Symbol,real=True)
theta_k=symbols("theta_k",cls=Symbol,real=True)

a00,a01,a10,a11=symbols("a00,a01,a10,a11",cls=Symbol,real=True)

A=Matrix([[a00,a01],[a10,a11]])

U=Matrix(
    [
        [cos(alpha_k/2)*exp(-I*theta_k),sin(alpha_k/2)*exp(-I*theta_k)],
        [I*sin(alpha_k/2) ,-I*cos(alpha_k/2)]
    ]
)

half=Rational(1,2)

Ae00=(half+half*cos(alpha_k))*a00-I*half*a10*sin(alpha_k)*exp(-I*theta_k)\
     +I*half*a01*sin(alpha_k)*exp(I*theta_k)+a11*(half-half*cos(alpha_k))


Ae01=half*a00*sin(alpha_k)-half*a11*sin(alpha_k)\
    -I*a10*(half-half*cos(alpha_k))*exp(-I*theta_k)-I*a01*(half+half*cos(alpha_k))*exp(I*theta_k)

Ae10=half*a00*sin(alpha_k)-half*a11*sin(alpha_k)\
    +I*a10*(half+half*cos(alpha_k))*exp(-I*theta_k)+I*a01*(half-half*cos(alpha_k))*exp(I*theta_k)

Ae11=a00*(half-half*cos(alpha_k))+I*half*a10*sin(alpha_k)*exp(-I*theta_k)\
    -I*half*a01*sin(alpha_k)*exp(I*theta_k)+a11*(half+half*cos(alpha_k))


Ae=Matrix([
    [Ae00,Ae01],
    [Ae10,Ae11]
])

df = fu(expand_complex(Ae - U.H@A@U))


pprint(expand(df))

