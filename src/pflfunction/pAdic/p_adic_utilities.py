#from pyadic import PAdic
from fractions import Fraction as frac
import sympy as sp
import flint
import time
import shelve

from flint import fmpq, fmpz

import sys
sys.setrecursionlimit(10000)


def prime_exponent(n, p):
    if(n == 0):
        raise ValueError(f"the argument is zero")

    exponent = 0
    while n % p == 0:
        n //= p
        exponent += 1
    return exponent


#The p-adic valuation.
def ord(q: frac, p: int):

    if not sp.isprime(p):
        raise ValueError(f"{p} is not a prime number")

    if isinstance(q, frac):
        num_exp = prime_exponent(q.numerator,p)
        den_exp = prime_exponent(q.denominator,p)

    elif isinstance(q,flint.fmpq): 
        num_exp = prime_exponent(q.numer(),p)
        den_exp = prime_exponent(q.denom(),p) 

    elif isinstance(q,flint.fmpz): 
        num_exp = prime_exponent(q,p)
        den_exp = 0  

    elif isinstance(q, int):
        num_exp = prime_exponent(q.numerator,p)
        den_exp = 0          

    else:
        raise ValueError(f"Unknown type {type(q)}")         

    return num_exp - den_exp

#Extended Euclidean algorithm
def extended_gcd(a, b):
    if b == 0:
        return a, 1, 0
    gcd, x1, y1 = extended_gcd(b, a % b)
    x = y1
    y = x1 - (a // b) * y1
    return gcd, x, y

#Modular inverse
def modular_inverse(a, m):
    gcd, x, y = extended_gcd(a, m)
    if gcd != 1:
        raise ValueError(f"No modular inverse exists for {a} modulo {m}")
    return x % m

#Converting between fmpq and frac
def fmpq_to_frac(q : flint.fmpq):
    return frac(int(q.numerator),int(q.denominator))

def frac_to_fmpq(q : frac):
    return flint.fmpq(q.numerator,q.denominator)

#converts a rational number into form int/p^n. Works for both frac and fmpq
def rational_to_padic(q: frac, p: int, acc: int):

    if q == 0:
        return 0

    if isinstance(q, int):
        ord_numer = ord(q,p)

        return ((q//p**ord_numer) % p**acc) * frac(p)**(ord_numer)

    if isinstance(q, flint.fmpz):
        ord_numer = ord(q,p)

        return ((q//fmpz(p)**ord_numer) % fmpz(p)**acc) * fmpq(p,1)**(ord_numer)     

    if isinstance(q, frac):
        ord_numer = ord(q.numerator,p)
        ord_denom = ord(q.denominator,p)

        return (((q.numerator//p**ord_numer) * modular_inverse(q.denominator//p**ord_denom,p**acc)) % p**acc) * frac(p)**(ord_numer-ord_denom)
    
    elif isinstance(q,flint.fmpq):
        ord_numer = ord(q.numer(),p)
        ord_denom = ord(q.denom(),p)

        return (((q.numer()//fmpz(p)**ord_numer) * modular_inverse(q.denom()//fmpz(p)**ord_denom,fmpz(p)**acc)) % fmpz(p)**acc) * fmpq(p,1)**(ord_numer-ord_denom) 
    
    else:
        raise ValueError(f"Unknown type {type(q)}")    


def rational_to_padic_round(q: frac, p: int, acc: int):

    if q == 0:
        return 0

    elif isinstance(q, int):
        ord_numer = ord(q,p)

        return (q % p**acc)


    elif isinstance(q, frac):
        ord_numer = ord(q.numerator,p)
        ord_denom = ord(q.denominator,p)

        acc = acc + ord_denom

        if(ord_numer-ord_denom>=acc):
            return 0

        else:
            pacc = p**acc
            pord_numer = p**ord_numer

            modinv = pow(q.denominator//p**ord_denom, -1, pacc)

            return (((q.numerator//pord_numer) * modinv * pord_numer) % pacc) * frac(p)**(-ord_denom)

    elif isinstance(q,flint.fmpq):
        acc = acc + ord(q.denominator,p)

        ord_numer = ord(q.numerator, p)
        ord_denom = ord(q.denominator, p)
        p = fmpz(p)
        modulus = p ** acc
        num = fmpz(q.numerator) // (p ** ord_numer)
        denom = fmpz(q.denominator) // (p ** ord_denom)
        inv_denom = pow(denom, -1, modulus)
        result_mod = (num * inv_denom * (p ** ord_numer)) % modulus
        return fmpq(result_mod) * fmpq(p) ** (-ord_denom)  
    
    else:
        raise ValueError(f"Unknown type {type(q)}")    

#Compute the modular representative in the range [-N/2, N/2]
def cmod(x, N):
    mod_val = x % N
    if mod_val > N // 2:
        mod_val -= N
    return mod_val


#Gives the Teichmueller representative (=Hansel lift)
def teich(q: frac, exp: int, p: int, acc: int):
    if not sp.isprime(p):
        raise ValueError(f"{p} is not a prime number")
    
    if exp == 0:
        return 1
    elif exp == 1:
        return pow(q,p**(acc-1),p**acc)
    elif exp > 1:
        return teich((q**(exp % (p-1))) % p,1,p,acc)
    

###############Computing p-adic gamma functions, zeta functions, etc.

def flint_factorial(n):
    return fmpz().fac_ui(n)

def dw(n, p):
    total = fmpq(0)
    p = fmpz(p)
    for i in range(n + 1):
        denom = flint_factorial((n - i) * p) * flint_factorial(i) * (p ** i)
        total += fmpq(1, denom)
    return total

def nmax_pgamma(acc: int, p: int):
    if not sp.isprime(p):
        raise ValueError(f"{p} is not a prime number")   
    
    c = acc
    red = acc/p
    while(sp.floor(red)>0):
        c+=sp.floor(red)
        red/=p

    return c   


def padic_gamma1(p: int,acc: int):
    if not sp.isprime(p):
        raise ValueError(f"{p} is not a prime number")

    return rational_to_padic(sum([dw(i+1,p)*fmpz(p)**i*flint_factorial(i) for i in range(0,nmax_pgamma(acc,p)+1)]),p,acc)


def padic_gamma3(p: int,acc: int):
    if not sp.isprime(p):
        raise ValueError(f"{p} is not a prime number") 

    return rational_to_padic(3*sum([dw(i+3,p) * fmpz(p)**i * flint_factorial(i+2) * frac_to_fmpq(sp.harmonic(i+2)**2 - sp.harmonic(i+2,2)) for i in range(0,nmax_pgamma(acc,p)+1)]),p,acc)

def padic_zeta3(p: int, acc: int):
    if not sp.isprime(p):
        raise ValueError(f"{p} is not a prime number")
    
    result = rational_to_padic((padic_gamma1(p,acc+2)**3-padic_gamma3(p,acc+2))/2,p,acc+2)

    return result.numer() % fmpz(p)**acc
