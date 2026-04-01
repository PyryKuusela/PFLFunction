import sympy
import math
import flint
import ast
import re
import json
import fractions
import decimal
import sys
import sympy

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

def get_data(operator_name):
    with open("input_data.txt","r") as infile:
        for line in infile:
            if line.startswith(operator_name):
                operator_data_in = line.strip()
                break
        else:
            print(f"The CY operator {operator_name} was not found.")
            operator_data_in = None
    return operator_data_in

def extract_top_level_brackets(s):
    parts = []
    depth = 0
    current = ''
    for char in s:
        if char == '[':
            if depth == 0:
                current = ''  # Start a new part
            depth += 1
        if depth > 0:
            current += char
        if char == ']':
            depth -= 1
            if depth == 0:
                parts.append(current.strip())
    return parts

def create_inputs_range(primes,operator_name,filename,scaling,style="enum"):
    operator_data_init = get_data(operator_name)
    assert operator_data_init is not None, "CY operator not found."

    # Data looks like e.g. ---> 4.1.1 [[[0,0,0,0,1],[-120,-1250,-4375,-6250,-3125]]] [-200,5] [[[1,-3125]],[],[]]
    # splits at first space (discarding 4.1.1) and splits remaining data according to outer brackets
    PF,geom,sings = extract_top_level_brackets(operator_data_init.split(" ",1)[1])
    operator_data = PF + " " + "[acc,scaling]" + " " + geom + " " + sings
    
    assert type(primes) is int or type(primes) is list, "Type of primes not int or list"
    if type(primes) is list:
        assert len(primes) == 2, "List needs to be in form [p_min,p_max]"

    if type(primes) is int and style == "enum":
        p_min = 1
        p_max = p    
        
    elif type(primes) is int and style == "range":
        p_min = 2
        p_max = p
        
    else:
        p_min = primes[0]
        p_max = primes[1]    
    
    with open(filename,"w") as f: 
        if style == "enum":
            for i in range(p_min,p_max+1):
                p = sympy.prime(i)
                print(str(p) + " " + operator_data.replace("scaling",f"{scaling}").replace("acc",f"{padic_acc}"),file=f)
                
        elif style == "range":
            for p in list(sympy.primerange(a=p_min,b=p_max+1)):
                print(str(p) + " " + operator_data.replace("scaling",f"{scaling}").replace("acc",f"{padic_acc}"),file=f)
                
        else:
            raise Exception("Unsupported style")
        
        f.close()
        
if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python bash_test.py <var1> <var2> <var3> <var4> <var5>")
        sys.exit(1)

    primes = json.loads(sys.argv[1]) # integer array
    operator_name = sys.argv[2] # string
    filename = "inputs/inputs_" + sys.argv[3] # string
    scaling = decimal.Decimal(sys.argv[4]) # rational
    padic_acc = int(sys.argv[5])
    
    if not filename.endswith(".txt"):
        filename = filename + ".txt"
    
    create_inputs_range(primes,operator_name,filename,scaling)
