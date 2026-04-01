import sympy as s
import math
import time

import os
import sys

from ..pAdic import p_adic_utilities as p_adic

from fractions import Fraction as frac


z = s.symbols('z')
theta = s.symbols('theta')


def generate_recurrence_from_operator(L):
    """
    Derives the recurrence relation for the power series solution to a one-parameter differential operator L.
    
    Parameters:
    - L: A differential operator in terms of theta=(z*d/dz) and z.
    
    Returns:
    - recurrence_relation: A SymPy expression representing the recurrence relation in terms of c(n+s_i) and n.
    """
    # Define symbols
    n = s.symbols('n')
    c = s.Function('c')
    z = s.symbols('z')
    theta = s.symbols('theta')

    # Expand operator to avoid silly errors
    L = L.expand()

    # Start with the operator L and substitute theta^k z^n = n^k z^n
    recurrence_expr = 0
    z_power_list = []
    for term in L.as_ordered_terms():
        z_power = s.Poly(term,z).as_poly().degree() 
        z_power_list.append(z_power)
        max_z_power = max(z_power_list)

    for term in L.as_ordered_terms():
        coeff = term.as_coeff_Mul()[0]    # Coefficient of the term
        theta_power = s.Poly(term,theta).as_poly().degree()  # Degree of theta
        z_power = s.Poly(term,z).as_poly().degree()
        recursion_depth = max_z_power-z_power 
        # Substitute theta^k * z^n with n^k * z^n
        recurrence_expr += coeff * (n+recursion_depth)**theta_power * c(n+recursion_depth)

    # now recurrence_relation = 0 is the desired equation
    recurrence_relation = recurrence_expr

    return recurrence_relation


def rational_periods(L, k, max_n):
    """
    Finds the power series parts d_a(n) of the period vector up to order k
    for a Fuchsian differential operator L of the form 
    L=sum_k P_k(z)*theta^k

    Parameters:
    - L: A differential operator in terms of theta=(z*d/dz) and z.
    - k: The number of derivatives of the recurrence relation (k => 0)
    - max_n: The maximum value of n to compute.
    
    Returns:
    - dk_dict_list: A list of dictionaries with sympy.Rational values for d_a(n) 
                    for a in {0,...,k} up and n in {0,...,max_n},
                    with all values as rational numbers.
    """
    
    c = s.Function('c')
    d = s.Function('d')
    n = s.Symbol('n')
    w = s.Wild('w')

    recurrence_relation = generate_recurrence_from_operator(L)

    # Get the shifts s_i from the recurrence relation
    # They are the arguments of each c(n + s_i) in the relation
    shifts = sorted({arg.args[0] - n for arg in recurrence_relation.atoms(c)}, key=lambda x: x.evalf())
    max_shift = max(shifts)
    
    initial_conditions = {0:1}
    dk_dict_list = []

    for b in range(k+1):
        if b > 0:    
            initial_conditions = {0:0}
        # Initialize the dictionary to store d_b(n) values with the initial conditions
        dk_values = {key: s.Rational(value) for key, value in initial_conditions.items()}
        for i in range(1,max_shift):
            dk_values[-i] = 0

        # Obtain the recurrence relation for d_b(n) by taking b derivatives of the recurrence relation for L
        log_k_recurrence_relation = s.diff(recurrence_relation,n,b)
        log_k_recurrence_relation = log_k_recurrence_relation.replace(c(w),d(0,w)).doit() 
        for a in range(1,b+1):    
            log_k_recurrence_relation = log_k_recurrence_relation.replace(s.Derivative(d(0,w),w,a),d(a,w)).doit() 


        # Reorganize recurrence relation for solving d_b(n + max_shift)
        # Isolate d_b(n + max_shift) to express it in terms of lower values
        recurrence_solved = s.solve(log_k_recurrence_relation, d(b,n + max_shift))[0]
        
        # Define the recurrence function by lambdifying n and each c(n + s_i),d_1(n + s_i),d_2(n + s_i),...,d_b(n + s_i) term
        shifted_ds = []
        if b > 0 :
            for a in range(b):
                shifted_ds.append([d(a,n + shift) for shift in shifts])
        # The solved recurrence relation for the current derivative b does not depend on d_b(n+max_shift)
        shifted_ds.append([d(b,n + shift) for shift in shifts if shift < max_shift])

        recurrence_func = s.lambdify((n, *shifted_ds), recurrence_solved, 'sympy')

        # Iteratively compute d_b(i) with lambdified function
        for i in range(-max_shift+1,max_n-(max_shift-1)):
            # Gather the values of d_b(n + s_i) from previous results
            previous_ds = []
            if b > 0:
                for a in range(b):
                    previous_ds.append([dk_dict_list[a][i + shift] for shift in shifts])
            previous_ds.append([dk_values[i + shift] for shift in shifts if shift < max_shift])
            # Evaluate d_b(i + 1) using the recurrence function and convert to Rational
            result = recurrence_func(i, *previous_ds)

            # Check if result is finite and valid
            is_valid = (
                isinstance(result, (int, float)) and math.isfinite(result)
                ) or (
                isinstance(result, s.Basic) and result.is_finite
                )
            
            if is_valid:
                dk_values[i + max_shift] = s.Rational(result)    
            else:
                dk_values[i + max_shift] = 0

        dk_dict_list.append(dk_values)
    return dk_dict_list


def padic_periods(L, k, max_n, p, acc):
    """
    Finds the power series parts d_a(n) of the period vector up to order k
    for a Fuchsian differential operator L of the form 
    L=sum_k P_k(z)*theta^k

    Parameters:
    - L: A differential operator in terms of theta=(z*d/dz) and z.
    - k: The number of derivatives of the recurrence relation (k => 0)
    - max_n: The maximum value of n to compute.
    
    Returns:
    - dk_dict_list: A list of dictionaries with frac values for d_a(n) 
                    for a in {0,...,k} up and n in {0,...,max_n},
                    with all values as p-adic approximations.
    """
    
    c = s.Function('c')
    d = s.Function('d')
    n = s.Symbol('n')
    w = s.Wild('w')

    recurrence_relation = generate_recurrence_from_operator(L)

    # Get the shifts s_i from the recurrence relation
    # They are the arguments of each c(n + s_i) in the relation
    shifts = sorted({arg.args[0] - n for arg in recurrence_relation.atoms(c)}, key=lambda x: x.evalf())
    max_shift = max(shifts)
    
    initial_conditions = {0:1}
    dk_dict_list = []

    for b in range(k+1):
        if b > 0:    
            initial_conditions = {0:0}
        # Initialize the dictionary to store d_b(n) values with the initial conditions
        dk_values = {key: s.Rational(value) for key, value in initial_conditions.items()}
        for i in range(1,max_shift):
            dk_values[-i] = 0

        # Obtain the recurrence relation for d_b(n) by taking b derivatives of the recurrence relation for L
        log_k_recurrence_relation = s.diff(recurrence_relation,n,b)
        log_k_recurrence_relation = log_k_recurrence_relation.replace(c(w),d(0,w)).doit() 
        for a in range(1,b+1):    
            log_k_recurrence_relation = log_k_recurrence_relation.replace(s.Derivative(d(0,w),w,a),d(a,w)).doit() 


        # Reorganize recurrence relation for solving d_b(n + max_shift)
        # Isolate d_b(n + max_shift) to express it in terms of lower values
        recurrence_solved = s.solve(log_k_recurrence_relation, d(b,n + max_shift))[0]
        
        # Define the recurrence function by lambdifying n and each c(n + s_i),d_1(n + s_i),d_2(n + s_i),...,d_b(n + s_i) term
        shifted_ds = []
        if b > 0 :
            for a in range(b):
                shifted_ds.append([d(a,n + shift) for shift in shifts])
        # The solved recurrence relation for the current derivative b does not depend on d_b(n+max_shift)
        shifted_ds.append([d(b,n + shift) for shift in shifts if shift < max_shift])


        # num, den = s.fraction(recurrence_solved)
        # num = s.factor(num)
        # #den = s.factor(den)
        # print(s.factor(num))
        #recurrence_solved = num/den

        recurrence_func = s.lambdify((n, *shifted_ds), recurrence_solved, 'sympy')

        # Iteratively compute d_b(i) with lambdified function
        for i in range(-max_shift+1,max_n):

            # Gather the values of d_b(n + s_i) from previous results
            previous_ds = []
            if b > 0:
                for a in range(b):
                    previous_ds.append([dk_dict_list[a][i + shift] for shift in shifts])
            previous_ds.append([dk_values[i + shift] for shift in shifts if shift < max_shift])

            # Evaluate d_b(i + 1) using the recurrence function and convert to Rational
            result = recurrence_func(i, *previous_ds)

            # Check if result is finite
            is_valid = (
                isinstance(result, (int, float, frac)) and math.isfinite(result)
                ) or (
                isinstance(result, s.Basic) and result.is_finite
                )

            #is_valid = True
            
            if is_valid:
                dk_values[i + max_shift] = p_adic.rational_to_padic_round(frac(result),p,acc)   
                #dk_values[i + max_shift] = p_adic.numden_to_padic_round(num_f(i, *previous_ds),den_f(i, *previous_ds),p,acc)   
            else:
                dk_values[i + max_shift] = frac(0)
        dk_dict_list.append(dk_values) 
    return dk_dict_list
