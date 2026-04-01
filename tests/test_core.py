import sympy as sp
import flint
import pflfunction
from pflfunction.one_parameter_zeta_functions import L_functions

def test_4_1_1(): # mirror quintic, only conifold sings
    theta = sp.symbols('theta')
    z = sp.symbols('z')
    assert (
        L_functions(theta**4-5*z*(5*theta+1)*(5*theta+2)*(5*theta+3)*(5*theta+4),1,7,1-3125*z,1,1,5,"mirror_quintic",pacc_init=20,nadd=100) 
            == [[7, 1, [5, 385]], [7, 2, [25, 350]], [7, 3, [10, 420]], [7, 4, [-5, -210]], [7, 5, [1, 301], 'C'], [7, 6, [-35, 805]]]
    )

def test_4_3_1(): # HV CY3
    theta = sp.symbols('theta')
    z = sp.symbols('z')
    assert (
        L_functions(theta**4-225*z**3*(theta+1)**2*(theta+2)**2+z**2*(theta+1)**2*(259*theta**2+518*theta+285)-z*(35*theta**4+70*theta**3+63*theta**2+28*theta+5),3,19,(1-z)*(1-9*z)*(1-25*z),1,1,flint.fmpq(2,3),"4.3.1",pacc_init=100,nadd=300)
            == [[19, 1, [-39, 7239], 'C'], [19, 2, [76, 38]], [19, 3, [-8, 4598]], [19, 4, [16, 9158]], [19, 5, [16, 9158]], [19, 6, [8, -6042]], [19, 7, [-44, -4522]], [19, 8, [-118, 16758]], [19, 9, [-84, 1558]], [19, 10, [12, 10678]], [19, 11, [-64, 3078]], [19, 12, [12, 1558]], [19, 13, [178, 20558]], [19, 14, [12, -3002]], [19, 15, [42, -722]], [19, 16, [57, 5415], 'C'], [19, 17, [-39, 7239], 'C'], [19, 18, [-54, 6118]]]
    )


def test_4_3_7(): # https://cycluster.mpim-bonn.mpg.de/operator.html?nn=4.3.7
    theta = sp.symbols('theta')
    z = sp.symbols('z')
    assert (
        L_functions(theta**4-944784*z**3*(theta+1)*(theta+2)*(2*theta+1)*(2*theta+5)-2916*z**2*(theta+1)**2*(20*theta**2+40*theta+17)-18*z*(6*theta**4+12*theta**3+3*theta**2-3*theta-1),5,7,1-324*z,1,1+108*z,flint.fmpq(-4,3),"4.3.7",pacc_init=50,nadd=300)
            == [[7, 1, [26, 567]], [7, 2, 0], [7, 3, [-28, 819]], [7, 4, [0, 294], 'C'], [7, 5, [16, 56]], [7, 6, [-20, 560]]]
    )

def test_HVK3():
    theta = sp.symbols('theta')
    z = sp.symbols('z')
    assert (
        L_functions(theta**3-2*z*(1+2*theta)*(2+5*theta*(1+theta))+64*z**2*(1+theta)**3,1.5,7,(-1+4*z)*(-1+16*z),1,1,0,"HVK3",pacc_init=20,nadd=50)
            == [[7, 1, [-7, -49]], [7, 2, [-2, 49], 'C'], [7, 3, [5, 35]], [7, 4, [-2, 49], 'C'], [7, 5, [5, 35]], [7, 6, [3, -21]]]
    )

