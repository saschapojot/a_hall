from sympy import  *

eps_k=symbols("epsilon_k",cls=Symbol,real=True)
vR=symbols("v_R",cls=Symbol,real=True)
eps_L=symbols("epsilon_L",cls=Symbol,real=True)
eps_k_a=symbols(r"epsilon_ka",cls=Symbol,real=True)
kx,ky=symbols("kx,ky",cls=Symbol,real=True)

sigma_0=Matrix(
 [[1,0],
  [0,1]]
)

sigma_x=Matrix(
    [[0,1],
     [1,0]]
)

sigma_y=Matrix([
    [0,-I],
    [I,0]
])

sigma_z=Matrix(
    [[1,0],
     [0,-1]]
)

k_abs=sqrt(kx**2+ky**2)

h0=eps_k
