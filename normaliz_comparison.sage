from PyNormaliz import *
import statistics
import time
attach("main.sage")

def weighted_h_star_order_polytope_normaliz(P, weight):
    """
    Uses PyNormaliz' WeightedEhrhartSeries function to compute the weighted h-star vector of the order polytope O(P) of a poset P with respect to a polynomial weight. 

    Note that for weighted Ehrhart series, the computation depends on the correspondence between the variables of the polynomial weight and the coordinates of the ambient space. This function assumes that the ith variable of the parent polynomial ring R of the weight, which can be accessed as the ith element of the tuple yielded by R.gens(), corresponds to the ith coordinate (1-indexed) of the ambient space.

    Uses the cover relations of P to extract a hyperplane representation of O(P), which is then passed into WeightedEhrhartSeries.

    Requires: 
        1 <= len(R.gens()) <= len(P) where R is the parent polynomial ring of the weight.
        The elements of P are labeled with [n]. Note that the elements of SageMath's Posets(n) generator are labeled with {0, ..., n - 1}, NOT [n] as required by this function.

    Example:
    sage: D = DiGraph({1: [2], 2: [6], 3: [7], 4: [7], 5: [7], 6: [7], 7: []})
    sage: P = Poset(D)
    sage: R.<x_1, x_2, x_3, x_4, x_5, x_6> = PolynomialRing(ZZ)
    sage: f = x_1^2 * x_3 * x_5 * x_6 + x_2 * x_3^3 * x_4
    sage: weighted_h_star_order_polytope_normaliz(P, f)
    1620*z^8 + 31360*z^7 + 146445*z^6 + 226362*z^5 + 124622*z^4 + 22896*z^3 + 1089*z^2 + 6*z
    """
    R = PolynomialRing(ZZ, 'z')
    z = R.gens()[0]
    n = len(P)
    weight_dict = weight.dict()
    # Must convert SageMath polynomial into Normaliz polynomial format first (e.g. x[1]^2 * x[2])
    monomial_terms = [""] * len(weight_dict)
    for index, (exp_vec, coeff) in enumerate(weight_dict.items()):
        variable_terms = [""] * (sum(x != 0 for x in exp_vec) + (0 if coeff == 1 else 1))
        if coeff != 1:
            variable_terms[0] = str(coeff)
        var_terms_idx = 0 if coeff == 1 else 1
        for i in range(len(exp_vec)):
            if exp_vec[i] != 0:
                variable_terms[var_terms_idx] = f"x[{i + 1}]" if exp_vec[i] == 1 else f"x[{i + 1}]^{exp_vec[i]}"
                var_terms_idx += 1
        monomial_terms[index] = " * ".join(variable_terms)
    normaliz_weight = " + ".join(monomial_terms)
    ineqs = []
    # Note that given a hyperplane a_1*x_1 + ... + a_n*x_n + c >= 0, SageMath represents it as [c, a_1, ..., a_n] and Normaliz uses [a_1, ..., a_n, c]
    for p_1, p_2 in P.cover_relations():
        row = [0] * (n + 1)
        row[p_2 - 1] = 1
        row[p_1 - 1] = -1
        ineqs.append(row)
    for min_elem in P.minimal_elements():
        row = [0] * (n + 1)
        row[min_elem - 1] = 1
        ineqs.append(row)
    for max_elem in P.maximal_elements():
        row = [0] * (n + 1)
        row[max_elem - 1] = -1
        row[-1] = 1
        ineqs.append(row)
    cone = Cone(inhom_inequalities = ineqs, polynomial = normaliz_weight)
    # Output format of WeightedEhrhartSeries: list of length 4, entry 0 is list of coefficients of h-star vector starting from lowest nonzero coefficient, entry 2 is degree of lowest nonzero coefficient term
    ehrhart_series = cone.WeightedEhrhartSeries()
    h_star = 0
    for i in range(len(ehrhart_series[0])):
        h_star += ehrhart_series[0][i] * z^(ehrhart_series[2] + i)
    return h_star

def compare_to_normaliz(P, weight):
    """
    Returns a tuple (time_1, time_2) where time_1 (resp. time_2) is the time in seconds it took weighted_h_star_order_polytope (resp. weighted_h_star_order_polytope_normaliz) to compute the weighted h-star vector of the order polytope O(P) of the poset P with respect to a polynomial weight.

    Takes the median of three runs for each algorithm and uses time.perf_counter() to record wall time.

    Requires: 
        1 <= len(R.gens()) <= len(P) where R is the parent polynomial ring of the weight.
        The elements of P are labeled with [n]. Note that the elements of SageMath's Posets(n) generator are labeled with {0, ..., n - 1}, NOT [n] as required by this function.

    """
    time_1 = [0] * 3
    for i in range(3):
        start = time.perf_counter()
        weighted_h_star_order_polytope(P, weight)
        end = time.perf_counter()
        time_1[i] = end - start
    time_2 = [0] * 3
    for i in range(3):
        start = time.perf_counter()
        weighted_h_star_order_polytope_normaliz(P, weight)
        end = time.perf_counter()
        time_2[i] = end - start
    return (statistics.median(time_1), statistics.median(time_2))