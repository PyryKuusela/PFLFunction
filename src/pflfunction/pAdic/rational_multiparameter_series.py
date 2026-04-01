from fractions import Fraction as frac
import sympy as sp
import numpy as np
import flint 
from functools import reduce
import time
from . import p_adic_utilities
import math
from collections import defaultdict

##############################################################
###Series and polynomial operations
##############################################################

def fmpq_expand(sp_f, x, x0: int, nmax: int, ctx):
    if sp_f == 0:
        return fmpq_mseries({(0,):0},nmax-1,ctx)
    elif isinstance(sp_f, (frac,int)):
        return fmpq_mseries({(0,):p_adic_utilities.frac_to_fmpq(sp_f)},nmax-1,ctx)
    else:
        terms_dict = sp_f.series(x, x0, nmax).removeO().as_coefficients_dict()

        coeff_dict = {}
        for term, coeff in terms_dict.items():
            if term == 1:
                exponents = (0,)
            else:
                exponents = tuple(int(term.as_powers_dict().get(v, 0)) for v in (x,))

            coeff_dict[exponents] = p_adic_utilities.frac_to_fmpq(coeff)

        return fmpq_mseries(coeff_dict,nmax-1,ctx)    

# function for truncating the flint polynomials
def truncatePoly(poly,nmax,ctx):
    polyDict = dict(poly.terms())

    delete = [key for key in polyDict if sum(key) > nmax]

    # Delete the keys
    for key in delete:
        del polyDict[key]

    return ctx.from_dict(polyDict)

#Also deletes any zero values
def truncate_dict(original_dict,nmax):
    dict_to_modify = original_dict.copy()
    delete_too_large = [key for key in dict_to_modify if sum(key) > nmax]
    delete_zero = [key for key in dict_to_modify if dict_to_modify[key] == 0]    

    # Delete the keys
    for key in delete_too_large:
        del dict_to_modify[key]

    for key in delete_zero:
        del dict_to_modify[key]    

    return dict_to_modify 

#Shifts key -> key + shift
def shift_dict(dict,shift):
    return {tuple(a + b for a, b in zip(key, shift)): value for key, value in dict.items()}    

def multiply_dict_vals(dict,mult):
    return {key: mult*value for key, value in dict.items()}    


#Converts a dictionart with rational (frac or fmpq) entries to rational (frac or fmpq) representations of p-adic numbers to the specified accuracy.
def dict_to_padic(dict,p,padicacc):
    return {key: p_adic_utilities.rational_to_padic(value,p,padicacc) for key, value in dict.items()}

def dict_to_padic_round(dict,p,padicacc):
    return {key: p_adic_utilities.rational_to_padic_round(value,p,padicacc) for key, value in dict.items()}

#Converts a polynomial with rational (frac or fmpq) entries to rational (frac or fmpq) representations of p-adic numbers to the specified accuracy.
def polyTopadic(poly,p,padicacc,ctx):
    polyDict = dict(poly.terms())
    newPolyDict = {key: p_adic_utilities.rational_to_padic(value,p,padicacc) for key, value in polyDict.items()}

    return ctx.from_dict(newPolyDict)

#The new series class uses polynomials only for multiplication and addition. The data is stored purely as a dictionary
class fmpq_mseries:

#TODO: We coerce the degrees to int here. This might be dangerous.

    #WARNING: This assumes that shift is compatible with the polynomial. Only use shift=0 if literally converting a polynomial to a series.
    def __init__(self,dict,max_order,ctx):

        #truncate_dict also gets rid of zero values, so a zero series should be represented by an empty dictionary.
        self.terms = truncate_dict(dict,max_order)
        self.ctx=ctx
        self.nmax = int(max_order)

    @classmethod
    def from_poly(cls,fmpq_mpoly,nmax,ctx):

        poly_dict = dict(fmpq_mpoly.terms())

        return cls(poly_dict, nmax, ctx)

    #######
    # Get basic information
    #######

    #Returns the minimal degree with which each variable appears.
    def min_degs(self):
        if(self.terms == {}):
            return tuple(0 for i in self.ctx.gens())
        
        comps = zip(*self.terms.keys())
        return tuple(min(x) for x in comps) 

    #Returns the minimal total degree of a monomal appearing in the series.
    def min_degree(self):
        if(self.terms == {}):
            return 0
        
        return int(min(sum(t) for t in list(self.terms.keys())))
    
    #Returns the maximal total degree of a monomal appearing in the series.
    def degree(self):
        if(self.terms == {}):
            return 0
        
        return int(max(sum(t) for t in list(self.terms.keys())))

    #Returns the dictionary where the exponents have been shifted corresponding to multiplication by a prefactor
    def to_dictionary(self,prefactor_exponents=0):
        if(prefactor_exponents == 0):
            return self.terms
        else: 
            return {tuple(a - b for a, b in zip(key, prefactor_exponents)): value for key, value in self.terms.items()}

    def to_poly(self,prefactor_exponents=0):
        if prefactor_exponents == 0:
            return self.ctx.from_dict(self.to_dictionary(self.min_degs()))
        else:
            return self.ctx.from_dict(self.to_dictionary(prefactor_exponents))    


    #######
    # Basic operations
    #######

    #This does not check whether the ctx is the same for the series to be multiplied. Maybe checked by python-flint though?
    def __rmul__(self, other):
        return self.__mul__(other)

    def __mul__(self, other):
        if isinstance(other, fmpq_mseries): 

            new_nmax = min(self.nmax+other.min_degree(),other.nmax+self.min_degree())

            #If one of the series is zero, return zero
            if(self.terms == {} or other.terms == {}):
                return fmpq_mseries({},new_nmax,self.ctx)
            
            new_prefactor = tuple(a + b for a, b in zip(self.min_degs(),other.min_degs()))
            new_poly = self.to_poly(self.min_degs())*other.to_poly(other.min_degs())
            new_dict = shift_dict(dict(new_poly.terms()),new_prefactor)

            return fmpq_mseries(new_dict,new_nmax,self.ctx)

        if isinstance(other, (flint.fmpq,flint.fmpz)): 
            new_dict = multiply_dict_vals(self.terms,other)
            return fmpq_mseries(new_dict,self.nmax,self.ctx)
        
        if isinstance(other, (frac,int)): 
            return p_adic_utilities.frac_to_fmpq(other)*self
        
        else:
            raise TypeError(f"Cannot multiply series with {type(other)}")    
        
    def __add__(self,other):
        if isinstance(other, fmpq_mseries): 
            new_nmax = min(self.nmax,other.nmax)

            #if self is zero return other and vice versa
            if(self.terms == {}):
                return other            
            if(other.terms == {}):
                return self
            
            #Otherwise, calculate the "LCM" of the two Laurent series
            comps = zip(self.min_degs(),other.min_degs())
            common_prefactor = tuple(min(x) for x in comps)

            #Compute the polynomial addition and turn it into a dictonary
            new_poly = self.to_poly(common_prefactor)+other.to_poly(common_prefactor)
            new_dict = shift_dict(dict(new_poly.terms()),common_prefactor)

            return fmpq_mseries(new_dict,new_nmax,self.ctx)
        
        #Adding a constant does nothing if nmax is less than 0.
        if isinstance(other, (flint.fmpq,flint.fmpz)): 
            if self.nmax < 0:
                return self
            else: 
                new_dict = self.terms
                new_dict[(0, 0)] = new_dict.get((0, 0), 0) + other
            return fmpq_mseries(new_dict,self.nmax,self.ctx)
        
        if isinstance(other, (frac,int)): 
            return self + p_adic_utilities.frac_to_fmpq(other) 
        

    def __repr__(self):
        if all(x == 0 for x in  self.min_degs()):
            prefactor = ""
        else:    
            prefactor = " ".join(tuple(f"{base}^{exp}" for base, exp in zip(self.ctx.gens(), self.min_degs())))

        return f"{prefactor}({self.to_poly()}) + O({self.nmax+1})"    
    
    def truncate(self,truncated_nmax):
        new_nmax = int(min(self.nmax,truncated_nmax))
        new_dict = truncate_dict(self.terms,self.nmax)

        return fmpq_mseries(new_dict,new_nmax,self.ctx)


    #######
    # Additional functionality
    #######
    def to_p_adic_series(self,p,padicacc):
        return fmpq_mseries(dict_to_padic(self.terms,p,padicacc),self.nmax,self.ctx)
    
    def to_p_adic_series_round(self,p,padicacc):
        return fmpq_mseries(dict_to_padic_round(self.terms,p,padicacc),self.nmax,self.ctx)

    #Warning: This should only be used when the series is actually a polynomial. Technically speaking it will otherwise reduce the accuracy to zero
    def monomial_power_mod(self,mods):
        if(isinstance(mods,int)):
            new_dict = defaultdict(int)
            
            #This takes care of summing up the items if multiple keys have the same mod
            for key, value in self.terms.items():
                new_key = tuple(x % mods for x in key)
                new_dict[new_key] += value  # Sum values for duplicate keys
            
            new_dict=dict(new_dict)

            return fmpq_mseries(new_dict,self.nmax,self.ctx)

        elif(len(mods)==self.ctx.nvars()):
            new_dict = defaultdict(int)

            #This takes care of summing up the items if multiple keys have the same mod
            for key, value in self.terms.items():
                new_key = tuple(a % b for a, b in zip(key, mods))
                new_dict[new_key] += value  # Sum values for duplicate keys

            new_dict=dict(new_dict)    

            return fmpq_mseries(new_dict,self.nmax,self.ctx)
        else:
            raise Exception(f"Cannot take mod of powers of the variables {self.ctx.names()} by {mods}")            
    
    #Raises the monomials to the power specified by the (multi-)index powers.
    def monomial_pow(self,powers):
        if(isinstance(powers,int)):
            new_dict = {tuple(powers * x for x in key): value for key, value in self.terms.items()}
            new_nmax = self.nmax * powers

            return fmpq_mseries(new_dict,new_nmax,self.ctx)

        elif(len(powers)==self.ctx.nvars()):
            new_dict = {tuple(a * b for a, b in zip(key, powers)): value for key, value in self.terms.items()}
            new_nmax = self.nmax * min(powers)

            return fmpq_mseries(new_dict,new_nmax,self.ctx)
        else:
            raise Exception(f"Cannot raise the variables {self.ctx.names()} to power {powers}")


    def evaluate_at(self,vals):

        if isinstance(vals,dict):
            raise Exception("Not implemented for dictionaries, please use a tuple") #TODO: Could make this work.
        
        elif isinstance(vals,(tuple,list)):
            fmpq_dict = {str(key): p_adic_utilities.frac_to_fmpq(value) for key, value in zip(self.ctx.gens(),vals)}
            prefact = math.prod(tuple(p_adic_utilities.frac_to_fmpq(base)**exp for base, exp in zip(vals, self.min_degs())))

            evaluated_poly = prefact*self.to_poly().subs(fmpq_dict)
            
            return p_adic_utilities.fmpq_to_frac(dict(evaluated_poly.terms())[(0,)*self.ctx.nvars()])
        else:
            raise TypeError(f"Cannot evaluate series at {type(vals)}") 

#TODO: This could be better and even evaluate the last multiplication mod modulus.       
    def evaluate_at_mod(self,vals,p,padic_acc):

        prefact = math.prod(tuple(flint.fmpq(base,1)**exp for base, exp in zip(vals, self.min_degs())))
        padic_acc_intermediate = padic_acc + p_adic_utilities.ord(prefact,p)

        termsdict = self.to_dictionary()
        polysum = sum(evaluate_term_mod(vals,exps,prefact,p,padic_acc_intermediate) for exps, prefact in termsdict.items())

        return polysum  

#TODO: Warning, this does not currently check that any of this makes sense. It will produce a fmpz polynomial out of any rational multivariate series.
    def to_fmpz_poly(self):
        deg = self.degree()
        termsdict = self.to_dictionary()

        coeff_list = [termsdict.get((i,), flint.fmpq(0)).numer() for i in range(deg + 1)]

        return flint.fmpz_poly(coeff_list)

#TODO: Warning, this does not currently check that any of this makes sense. It will produce a fmpz_md_polynomial out of any rational multivariate series.
    def to_fmpz_mod_poly(self,mod_ctx):
        deg = self.degree()
        termsdict = self.to_dictionary()

        coeff_list = [termsdict.get((i,), flint.fmpq(0)).numer() for i in range(deg + 1)]

        return flint.fmpz_mod_poly(coeff_list,mod_ctx)        
        

    #TODO: This has the problem that it does not affect the zero coefficients, so currently this works only for 'multiplicative' functions. Should be sufficient for out purposes.
    def apply_to_coeffs(self,f):
        polyDict = dict(self.poly.terms())
        newPolyDict = {key: f(value) for key, value in polyDict.items()}

        return fmpq_mseries(self.ctx.from_dict(newPolyDict),self.nmax,self.ctx)


#TODO: Apparently e.g. 0*O(5) + ser + O(10) truncates the series ser to order 5. Should not be a problem, but might not be what we want.


##############################################################
###Matrix operations
##############################################################

#Inverts a polynomial matrix (does not work for series!)
def invertMatrix(mat):
    #TODO: Convert the flint poly to sympy poly if need be

    invmat = mat.inv() 
    return invmat.applyfunc(sp.simplify)


#TODO: Warning: this assumes that vals are integers
#TODO: WARNING: This assumes one-parameter case
#TODO: Rename
def evaluate_term_mod(vals,exps,prefactor,p,padic_acc):
    padic_acc_intermediate = padic_acc - p_adic_utilities.ord(prefactor,p)

    val = prefactor*math.prod(tuple(pow(base,exp,p**padic_acc_intermediate) for base, exp in zip(vals, exps)))  

    return val
