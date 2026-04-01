from .CYnOperatorPeriods import rational_periods, padic_periods
from ..pAdic.rational_multiparameter_series import fmpq_mseries
import sympy as s
from flint import fmpq, fmpq_mpoly_ctx
from math import comb, factorial
import time
import flint

def period_matrix(periods,ctx):
    """
    Construct the corrected E-tilde matrix from periods 
    using Eq. (2.26) from arXiv:2312.07611.
    Here, tilde{E}_ab(z) = theta^a f_b(z) where
    f_b are the power series parts of the periods.
    
    Parameters:
    - periods: List of dictionaries d_j(n), with j from 0 to k.

    Returns:
    - matrix: k x k matrix of FLINT power series.
    """
    k = len(periods)
    n_max = max(periods[0].keys())
    #ctx = fmpq_mpoly_ctx.get(('z',))  # Univariate power series context

    matrix = []
    for a in range(k):  # rows: theta^a
        row = []
        for b in range(k):  # cols: f_b
            series_dict = {}
            j_min = max(0, b - a)
            j_max = b
            for n in range(n_max): 
                coeff =  (
                    sum([comb(a, b-j) *
                         s.Rational(1,factorial(j)) * 
                         periods[j][n] * 
                         n ** (a+j-b) 
                        for j in range(j_min, j_max+1)
                        ])
                )
                key = (int(n),)
                series_dict[key] = fmpq(int(coeff.numerator), int(coeff.denominator))

            #Modified n_max -> n_max-1
            row.append(fmpq_mseries(series_dict,n_max-1,ctx)) # Custom univariate power series 

        matrix.append(row)
    return matrix
