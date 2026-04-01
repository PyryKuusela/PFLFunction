import sys
import fcntl
import json
import sympy
import flint
import math
import os

from pflfunction.one_parameter_zeta_functions import L_functions

z,theta = sympy.symbols('z,theta')

def build_PF(PF_arr):
    PF_eqs = []
    for PF_eq in PF_arr:
        PF = 0
        for z_pow in range(len(PF_eq)):
            theta_poly = array_to_poly(PF_eq[z_pow],theta)
            PF += sympy.factor(z**z_pow * theta_poly)
        PF_eqs.append(PF)
    return PF_eqs

def array_to_poly(array,var):
    if array == []: # if nothing... for discriminants
        return 1
    else: # if something
        poly = 0
        for i in range(len(array)):
            poly += array[i]*var**i
        return sympy.factor(poly)

def generate_data(PF_arr,p,settings,geom,disc,nadd,filename):
    
    #
    # Do computations
    #
    
    PF = build_PF(PF_arr)
    p = p
    padic_acc = settings[0]
    scaling = settings[1]
    geom_param = flint.fmpq(-geom[0],geom[1]) # data taken in is chi of mirror; is always 0 for K3s.
    
    if len(disc[0]) > 0:
        coni = math.prod([array_to_poly(coni,z) for coni in disc[0]])
    else:
        coni = array_to_poly(disc[0],z) # will give 1, meaning polynomial 1. No roots, so no coni singularities.
    if len(disc[1]) > 0:
        app = math.prod([array_to_poly(app,z) for app in disc[1]])
    else:
        app = array_to_poly(disc[1],z) # will give 1, meaning polynomial 1. No roots, so no app singularities.
    if len(disc[2]) > 0:
        rest = math.prod([array_to_poly(rest,z) for rest in disc[2]])
    else:
        rest = array_to_poly(disc[2],z) # will give 1, meaning polynomial 1. No roots, so no 'rest' singularities.
    
    results = L_functions(PF[0],scaling,p,coni,app,rest,geom_param,filename,pacc_init=padic_acc,nadd=nadd)

    filename = "outputs/outputs_" + filename
    # print("Current working directory:", os.getcwd())
    # print("Trying to write to:", os.path.abspath(filename))
    
    # only want to do this after computations are done
    # I think this is because flock will wait until its turn to call
    if not filename.endswith(".txt"):
        filename = filename + ".txt"
    
    with open(filename,"a") as f: 
        fcntl.flock(f, fcntl.LOCK_EX)
        
        try:
            for result in results:
                print(result,file=f)
            f.flush()

        finally:
            fcntl.flock(f, fcntl.LOCK_UN)
            
        
if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python bash_test.py <var1> <var2> <var3> <var4> <var5> <var6> <var7>")
        sys.exit(1)

    p = json.loads(sys.argv[1]) # integer prime
    PF = json.loads(sys.argv[2]) # array of [PF equations]
    settings = json.loads(sys.argv[3]) # array of [padic_acc, scaling]
    geom = json.loads(sys.argv[4]) # array of [euler_char, int_num]
    disc = json.loads(sys.argv[5]) # array of [[coni],[apparent],[rest]] excluding inf
    nadd = json.loads(sys.argv[6]) # number of additional terms; 0 by default
    filename = "" + sys.argv[7] # hardcode filename
    
    generate_data(PF,p,settings,geom,disc,nadd,filename)
