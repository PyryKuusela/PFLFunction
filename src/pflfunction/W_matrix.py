from sympy import symbols, integrate, Rational, exp, diff, simplify, expand, Matrix
from fractions import Fraction as frac

"""
    **Description**

    `d_log_n`: Defining the logarithmic derivative for order 0, order 1 and higher orders. 

    **Arguments**

    `f`: function
    `x`: variable
    `n`: order of the derivative

    **Returns** 
    The logarithmic derivative for order 0 returns f, for order 1 returns the logarithmic derivative we defined before,
    and for higher orders up tu n-1.

"""
def d_log(f,x,n=1):
    if n==0:
        return f
    elif n==1:
        return x*diff(f,x)
    elif n>0:
        return d_log(x*diff(f,x),x,n-1)
    else:
        raise Exception(f"Cannot take a derivative of order {n}")
    
"""
    **Description**

    `W_matrix`: period matrix constructed by giving every entry in terms of differential equations in terms of just one variable:
                for the case of 3-folds, the relations are expressed in terms of W03
                for the case of 4-folds, the relations are expressed in terms of W04    

    **Arguments**

    `L`: differential operator
    `hodge_type`: Hodge type for 3-folds[1,1,1,1] and 4-folds [1,1,1,1,1]

    **Returns** 
    A matrix

"""

def W_matrix(L,hodge_type):
    if(hodge_type==[1,1,1]):

        theta, z = symbols('theta z')

        L = expand(L)

        S1 = L.coeff(theta, 1)
        S2 = L.coeff(theta, 2) 
        S3 = L.coeff(theta, 3)

        R1 = S1/S3
        R2 = S2/S3

        W11 = simplify(exp(integrate(-R2*2/(3*z), z))*(exp(-integrate(-R2*2/(3*z), z)).subs(z,0)))
        W02 = -W11
        W12 = simplify(-d_log(W11,z) - R2*W11)
        W22 = simplify(-d_log(W11,z,2) - d_log(R2,z)*W11 - frac(1,2)*R2*d_log(W11,z) + R1*W11)
        return -Matrix([[0,0,W02],[0,W11,W12],[W02,W12,W22]])    
    
    elif(hodge_type==[1,1,1,1]):
        theta, z = symbols('theta z')

        L = expand(L)

        S2 = L.coeff(theta, 2)
        S4 = L.coeff(theta, 4)  
        S3 = L.coeff(theta, 3) 

        R = S3/S4

        W03 = simplify(-exp(integrate(-R/(2*z), z))*(exp(-integrate(-R/(2*z), z)).subs(z,0)))
        W12 = -W03
        W13 = simplify(-d_log(W03,z))
        W23 = simplify(-d_log(d_log(W03,z),z) + R*W13 + S2/S4*W12)
        return Matrix([[0,0,0,W03],[0,0,W12,W13],[0,-W12,0,W23],[-W03,-W13,-W23,0]])

    elif(hodge_type==[1,1,1,1,1]):

        theta, z = symbols('theta z')

        L = expand(L)

        S1 = L.coeff(theta, 1)
        S2 = L.coeff(theta, 2)
        S3 = L.coeff(theta, 3)
        S4 = L.coeff(theta, 4)  
        S5 = L.coeff(theta, 5) 
       
        R = S4/S5
        R2 = S3/S5
        R3 = S2/S5
        R4 = S1/S5

        W04 = simplify(exp(integrate(-R*2/(5*z), z))*exp(-integrate(-R*2/(5*z), z)).subs(z,0))
        W13 = -W04
        W22 = W04
        W14 = simplify(-frac(3,2) * d_log(W04,z))
        W23 = simplify(frac(1,2)*d_log(W04,z))
        W24 = simplify(-frac(3,2) * d_log(W04,z,2) - frac(3,2) *R*d_log(W04,z) - R2*W04)
        W33 = simplify(2*d_log(W04,z,2) + frac(3,2) *R*d_log(W04,z) + R2*W04)
        W34 = simplify(d_log(W04,z,3) + frac(3,4)*R*d_log(W04,z,2) + frac(1,2)*R2*d_log(W04,z) + frac(3,4)*d_log(R,z)*d_log(W04,z)+ frac(1,2)*d_log(R2,z)*W04)
        W44 = simplify(d_log(W04,z,4) + frac(7,4)*R*d_log(W04,z,3) + d_log(W04,z,2)*(frac(3,2) * d_log(R,z) + frac(5,2)*R2 + frac(3,4)*R*R)+ d_log(W04,z)*(frac(3,4)*d_log(R,z,2) + d_log(R2,z) + frac(1,2)*R3 + 2*R*R2 + frac(3,4)*R*d_log(R,z)) + W04*(-R4 + R2*R2 + frac(1,2)*R*d_log(R2,z) + frac(1,2)*d_log(R2,z,2)))
        return Matrix([[0,0,0,0,W04],[0,0,0,W13,W14],[0,0,W22,W23,W24],[0,W13,W23,W33,W34],[W04,W14,W24,W34,W44]])

    else:
        raise Exception("The Hodge type currently not supported.")

