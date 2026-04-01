from fractions import Fraction as frac
import sympy as sp
import flint 
import numpy as np
import math

from . import multipoint_evaluation as mpe
from . import W_matrix as Wm

import os
import sys


from .PicardFuchs import CYnOperatorPeriods 
from .PicardFuchs import EtildeOneParam

from .pAdic import rational_multiparameter_series as rs
from .pAdic import p_adic_utilities as p_utils


#This class holds the "global" data of the computation
class computation_data:
    def __init__(self,L,hodge_type,nmax,p,pacc,pacc_init,conifold_locus,apparent_sing_locus,other_sing_locus,chiratio,ctx,nadd,label):
        self.L = L
        self.hodge_type = hodge_type
        self.nmax = nmax
        self.p = p
        self.pacc = pacc
        self.pacc_init = pacc_init
        self.ctx = ctx
        self.conifold_locus = sp.sympify(conifold_locus)
        self.apparent_sing_locus = sp.sympify(apparent_sing_locus)
        self.other_sing_locus = sp.sympify(other_sing_locus)
        self.chiratio = chiratio

        self.nadd = nadd #The number of terms that should vanish at the end of the series for the purposes of checking that the series terminates
        self.label = label #unique label for the computation (e.g. the new number of the corresponding operator). Should be a string.

def period_matrix(cd):
    padic_periods=CYnOperatorPeriods.padic_periods(cd.L,len(cd.hodge_type)-1,cd.nmax,cd.p,cd.pacc_init)

    mat = np.array(EtildeOneParam.period_matrix(padic_periods,cd.ctx),dtype=object).T  

    return mat 

#Returns the intersection matrix.
def sigma_mat(hodge_type):
    if(hodge_type == [1,1,1]):
        return np.array([[0,0,1],[0,-1,0],[1,0,0]])
    elif(hodge_type == [1,1,1,1]):
        return np.array([[0,0,0,-1],[0,0,1,0],[0,-1,0,0],[1,0,0,0]])
    else:
        raise Exception("The Hodge type currently not supported.")

def inverse_W_series_mat(L,hodge_type,nmax,ctx):
    z = sp.symbols('z')

    mat = rs.invertMatrix(Wm.W_matrix(L,hodge_type)) 

    np_mat = np.array(mat.tolist(), dtype=object)
    np_mat = np.vectorize(lambda f: rs.fmpq_expand(f,z, 0, nmax, ctx))(np_mat)  

    return np_mat


def get_period_matrices(cd):
    period_mat = period_matrix(cd)

    short_pm = np.vectorize(lambda x: x.truncate(sp.ceiling(frac(cd.nmax,cd.p))+1))(period_mat) 

    prod1=np.matmul(short_pm, inverse_W_series_mat(cd.L,cd.hodge_type,sp.ceiling(frac(cd.nmax,cd.p))+1,cd.ctx))
    inverse_short_pm=np.matmul(sigma_mat(cd.hodge_type),prod1).T 

    return period_mat, inverse_short_pm


def U0_matrix(gamma,cd):
    if(cd.hodge_type == [1,1,1]):
        return np.array([[1,0,0],[0,cd.p,0],[0,0,cd.p**2]]) 
    elif(cd.hodge_type == [1,1,1,1]):
        return np.array([[1,0,0,0],[0,cd.p,0,0],[0,0,cd.p**2,0],[gamma*cd.p**3,0,0,cd.p**3]])
    elif(cd.hodge_type == [1,1,1,1,1]):
        return np.array([[1,0,0,0,0],[0,cd.p,0,0,0],[0,0,cd.p**2,0,0],[gamma*cd.p**3,0,0,cd.p**3,0],[0,gamma*cd.p**4,0,0,cd.p**4]])    
    else:    
        raise Exception("The Hodge type currently not supported.")   
    

def U_matrix(gamma,cd):
    period_mat, inv_period_mat = get_period_matrices(cd)

    inv_period_mat_pow = np.vectorize(lambda f: f.monomial_pow(cd.p))(inv_period_mat) 

    prod1=np.matmul(inv_period_mat_pow,U0_matrix(gamma,cd))
    prod2=np.matmul(prod1,period_mat)

    return prod2


def U_denominator(cd):
    z = sp.symbols('z')
    if(cd.pacc > len(cd.hodge_type)-1):
            return (cd.conifold_locus).subs(z, z**cd.p)**(cd.pacc-len(cd.hodge_type))*(cd.other_sing_locus).subs(z, z**cd.p)**(cd.pacc-2)*(cd.apparent_sing_locus).subs(z, z**cd.p) 
    else:
        return 1

def U_numerator(gamma,cd):   
        z = sp.symbols('z')

        #WARNING: I the denominator is computed as series to 10*p terms here. This will couse problems if the real degree is greater. 
        denom_series = rs.fmpq_expand(U_denominator(cd),z, 0, 10*cd.p, cd.ctx)
        
        denom_mat = np.diag([denom_series] * sum(cd.hodge_type)).astype(object)

        numer_mat = np.matmul(denom_mat,U_matrix(gamma,cd))
        numer_padic = np.vectorize(lambda f: f.to_p_adic_series_round(cd.p,cd.pacc))(numer_mat)

        #Test whether the series really terminates. If not, log this. Log the degree of the series, and nmax for every prime.
        log(f"({cd.p},{np.vectorize(lambda f: f.degree())(numer_padic).max()},{cd.nmax})",cd)

        if(cd.nmax-np.vectorize(lambda f: f.degree())(numer_padic).max() < cd.nadd):
           log(f"The numerator did not converge to specified accuracy for p={cd.p}",cd)
           log(f"Nmax: {cd.nmax}, series terminated at {np.vectorize(lambda f: f.degree())(numer_padic).max()}",cd)

        #This takes into account that we will be substituting in Teichmüller representatives for which a=a**p
        numer_teich = np.vectorize(lambda f: f.monomial_power_mod(cd.p-1))(numer_padic)

        return numer_teich


#############
## Evaluate the matrices and compute the polynomials
#############

#Evaluate the U-matrix at given Teichmüller coordinates.
def U_numerators_evaluated(gamma,cd,teich_coords_list):

    mod_ctx = flint.fmpz_mod_poly_ctx(cd.p**cd.pacc)
    U_numerator_poly = np.vectorize(lambda f: f.to_fmpz_mod_poly(mod_ctx))(U_numerator(gamma,cd))
    
    return mpe.multipoint_evaluate_mat(U_numerator_poly,teich_coords_list,mod_ctx)
    
#Compute the coefficients of the polynomial R(T;X) from the traces of the U-matrix products.
def R_poly_coeffs(tr_list,cd):
    a_coeff = -tr_list[0]

    b_coeff = p_utils.rational_to_padic_round((tr_list[0]**2-tr_list[1])/flint.fmpq(2),cd.p,cd.pacc)

    if(not (a_coeff.denom() == 1 and b_coeff.denom() == 1)):
        raise Exception("Something went wrong - the Zeta function coefficient has non-trivial denominator in p.")

    return [p_utils.cmod(a_coeff.numer(),flint.fmpz(cd.p)**(cd.pacc)),p_utils.cmod(flint.fmpq(b_coeff).numer(),flint.fmpz(cd.p)**(cd.pacc))]

def compute_coefficient_list_new(cd):
    pgamma = p_utils.rational_to_padic(-cd.chiratio*p_utils.padic_zeta3(cd.p,cd.pacc+1),cd.p,cd.pacc+1)

    z = sp.symbols('z')

    mod_ctx = flint.fmpz_mod_poly_ctx(cd.p**cd.pacc)

    denom_series = rs.fmpq_expand(U_denominator(cd),z, 0, 10*cd.p, cd.ctx).monomial_power_mod(cd.p-1) 
    denom = denom_series.to_fmpz_mod_poly(mod_ctx)

    teich_coords_list = [p_utils.teich(i,1,cd.p,cd.pacc + 1) for i in range(1, cd.p)]

    
    U_numerator_list = U_numerators_evaluated(pgamma,cd,teich_coords_list)   
    denom_eval_list = mpe.multipoint_evaluate(denom,teich_coords_list,mod_ctx)  



    result_list = []

    for zval in range(1,cd.p):

        #Apparent or other singularities
        if ((cd.apparent_sing_locus * cd.other_sing_locus).subs(z,zval) % cd.p == 0):
            result_list.append([cd.p,zval,0])
        #Conifolds
        elif((cd.conifold_locus).subs(z,zval) % cd.p == 0):
            trU = np.trace(U_numerator_list[zval-1])
            trUU = np.trace(np.matmul(U_numerator_list[zval-1], U_numerator_list[zval-1]))
            
            tr_list = [p_utils.rational_to_padic_round(flint.fmpq(trU,denom_eval_list[zval-1]),cd.p,cd.pacc),p_utils.rational_to_padic_round(flint.fmpq(trUU,denom_eval_list[zval-1]**2),cd.p,cd.pacc)]

            result_list.append([cd.p,zval,R_poly_coeffs(tr_list,cd),"C"]) 
        #Smooth
        else:
            trU = np.trace(U_numerator_list[zval-1])
            trUU = np.trace(np.matmul(U_numerator_list[zval-1], U_numerator_list[zval-1]))
            
            tr_list = [p_utils.rational_to_padic_round(flint.fmpq(trU,denom_eval_list[zval-1]),cd.p,cd.pacc),p_utils.rational_to_padic_round(flint.fmpq(trUU,denom_eval_list[zval-1]**2),cd.p,cd.pacc)]
            
            result_list.append([cd.p,zval,R_poly_coeffs(tr_list,cd)]) 

    return result_list


##################
##Logs
##################

def log(message, cd):
    folder = "logs"
    file_path = os.path.join(folder, "log_"+cd.label)

    os.makedirs(folder, exist_ok=True)

    with open(file_path, "a") as f:
        f.write(message+"\n")

##################
## Main function
##################

def L_functions(L,C,p,conifold_locus,apparent_sing_locus,other_sing_locus,chiratio,label,pacc_init=0,nadd=0):
    # pacc_init is the variable A in the paper.

    CTX = flint.fmpq_mpoly_ctx.get(('z', 1), 'lex')

    theta = sp.symbols('theta')
    b = sp.degree(L,theta)
    hodge_type = [1] * b

    nmax = math.ceil(C * p) + nadd

    k = sp.ceiling(b/2)
    pacc = int(sp.ceiling(sp.log(sp.binomial(b, k), p) + k * (b - 1) / 2))

    if(not(nadd==0) and pacc_init ==0):
        raise Exception(f"Initial p-adic accuracy A must be provided explicitly if additional terms in the period series are included.")  

    #TODO: We should check the denominator of W
    if(pacc_init == 0):
        if chiratio == 0:
            pacc_init = int(pacc + (2*b-1)*math.ceil(C)-b+1)
        else:
            pacc_init = int(pacc + (2*b-1)*math.ceil(C)-b+1-p_utils.ord(chiratio,p))  

    cd = computation_data(L,hodge_type,nmax,p,pacc,pacc_init,conifold_locus,apparent_sing_locus,other_sing_locus,chiratio,CTX,nadd,label)
    return compute_coefficient_list_new(cd)
