from flint import fmpz, fmpz_poly, fmpz_mod_poly_ctx, fmpz_mod_poly, fmpz_mod, fmpz_mod_ctx
import random
import time
import numpy as np

def build_balanced_subproduct_tree(points, mod_ctx):
    """
    Builds a balanced subproduct tree for arbitrary points using recursive splitting.
    Returns the root product polynomial and a nested tree structure.
    """
    if len(points) == 1:
        # Base case: single node (x - a)
        poly = fmpz_mod_poly([-points[0] % mod_ctx.modulus(), 1], mod_ctx)
        return poly, [poly]
    
    mid = len(points) // 2
    left_poly, left_tree = build_balanced_subproduct_tree(points[:mid], mod_ctx)
    right_poly, right_tree = build_balanced_subproduct_tree(points[mid:], mod_ctx)

    root = left_poly * right_poly

    return root, [root] + [left_tree, right_tree]

def fast_multipoint_eval(f, tree):
    """
    Recursively evaluate f at the points defined by the subproduct tree.
    Returns a list of evaluations at leaves in the order of the original points.
    
    - f: fmpz_mod_poly to evaluate
    - tree: subproduct tree structure [poly, left_tree, right_tree] or leaf poly
    """
    if len(tree) == 1:
        # Leaf node: remainder modulo (x - a) = evaluation at a
        # f mod (x - a) = constant term of remainder
        remainder = f % tree[0]

        if(remainder.is_zero()): # the constant term is 0
            return [0]
        else:
            val = remainder.coeffs()[0]  # constant term
            return [val]
    
    poly = tree[0]
    left_tree = tree[1]
    right_tree = tree[2]
    
    # Compute remainder of f modulo the node's polynomial
    remainder = f % poly
    
    # Recursively evaluate in left and right subtrees
    left_vals = fast_multipoint_eval(remainder, left_tree)
    right_vals = fast_multipoint_eval(remainder, right_tree)
    
    return left_vals + right_vals

#mod_ctx here is really mod_poly_ctx
def multipoint_evaluate(f, points, mod_ctx):
    tree = build_balanced_subproduct_tree(points, mod_ctx)
    results = fast_multipoint_eval(f, tree[1])
    results2 = [fmpz(str(i)) for i in results]
    return results2


def multipoint_evaluate_mat(mat, points, ctx):
    """
    Evaluate a matrix of fmpz_mod_poly entries at many points.
    
    Parameters:
        M: numpy ndarray of shape (m, k) with fmpz_mod_poly entries
        points: list of fmpz points to evaluate at
        ctx: fmpz_mod_ctx

    Returns:
        List of numpy matrices [M(x1), ..., M(xn)] of shape (m, k)
    """
    m, k = mat.shape
    n_points = len(points)

    polys = mat.flatten()

    tree = build_balanced_subproduct_tree(points, ctx)

    #Evaluate each polynomial at all points
    all_values = []
    for poly in polys:
        evals = fast_multipoint_eval(poly, tree[1])
        evals2 = [fmpz(str(i)) for i in evals]
        all_values.append(evals2)

    #Reshape to list of matrices
    result = []
    for i in range(n_points):
        matrix_i = np.array([all_values[j][i] for j in range(m * k)], dtype=object).reshape((m, k))
        result.append(matrix_i)

    return result
