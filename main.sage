from itertools import combinations

def weighted_h_star_order_polytope(P, weight):
    """
    Computes the weighted h-star vector of the order polytope O(P) of a poset P with respect to a polynomial weight. 

    Note that for weighted Ehrhart series, the computation depends on the correspondence between the variables of the polynomial weight and the coordinates of the ambient space. This function assumes that the ith variable of the parent polynomial ring R of the weight, which can be accessed as the ith element of the tuple yielded by R.gens(), corresponds to the ith coordinate (1-indexed) of the ambient space.

    Requires: 
        1 <= len(R.gens()) <= len(P) where R is the parent polynomial ring of the weight.
        The elements of P are labeled with [n]. Note that the elements of SageMath's Posets(n) generator are labeled with {0, ..., n - 1}, NOT [n] as required by this function.

    Example:
    sage: D = DiGraph({1: [2], 2: [6], 3: [7], 4: [7], 5: [7], 6: [7], 7: []})
    sage: P = Poset(D)
    sage: R.<x_1, x_2, x_3, x_4, x_5, x_6> = PolynomialRing(ZZ)
    sage: f = x_1^2 * x_3 * x_5 * x_6 + x_2 * x_3^3 * x_4
    sage: weighted_h_star_order_polytope(P, f)
    1620*z^8 + 31360*z^7 + 146445*z^6 + 226362*z^5 + 124622*z^4 + 22896*z^3 + 1089*z^2 + 6*z
    """
    R = PolynomialRing(ZZ, 'z')
    z = R.gens()[0]
    n = len(P)
    # Counting descents only works with order-preserving labeling, which the input P may not have. But we don't want to just relabel P, since that changes the polytope by permuting coordinates and so the original weight indices will no longer be accurate. So we need to also relabel the weight indices later inside the for-loop.
    base_lin_ext = P.linear_extensions()[0]
    natural_mapping = {element: j + 1 for j, element in enumerate(base_lin_ext)}
    P = P.relabel(natural_mapping)
    # Algorithm works on monomials, so need to extract those as a dictionary whose values are coefficients. Then convert the monomials into lists of variable indices and store them.
    weight_dict = weight.dict()
    monomial_weights = [0] * len(weight_dict)
    for index, (exp_vec, coeff) in enumerate(weight_dict.items()):
        monomial_weight = [idx for idx, power in enumerate(exp_vec, start=1) for _ in range(power)]
        # Relabeling weight
        monomial_weight = [natural_mapping[i] for i in monomial_weight]
        monomial_weights[index] = (monomial_weight, coeff)
    h_star = 0
    for ext in P.linear_extensions():
        num_descents_before = [0] * n
        for i in range(1, n):
            num_descents_before[i] = num_descents_before[i - 1] + (1 if ext[i - 1] > ext[i] else 0)
        first_sum_term = z^(num_descents_before[-1])
        pi_inverses = [0] * n
        for i in range(n):
            pi_inverses[ext[i] - 1] = i + 1
        h_star_linext = 0
        for (monomial_weight, coeff) in monomial_weights:
            deg_s = len(monomial_weight)
            # Sort weight variable indices in increasing order of pi inverse
            monomial_weight = sorted(monomial_weight, key=lambda x: pi_inverses[x - 1])
            sum_over_subsets = 0
            for k in range(deg_s + 1):
                for subset in combinations(range(deg_s), k):
                    prod_in_subset = 1
                    prod_out_subset = 1
                    for j in range(deg_s):
                        r_j = 0
                        for l in subset:
                            if l < j:
                                r_j += 1
                        if j in subset:
                            prod_in_subset *= (pi_inverses[monomial_weight[j] - 1] + r_j)
                        else:
                            prod_out_subset *= (num_descents_before[pi_inverses[monomial_weight[j] - 1] - 1] + r_j)
                    second_sum_term = z^k * (1 - z)^(deg_s - k) * prod_in_subset * prod_out_subset
                    sum_over_subsets += second_sum_term
            h_star_linext += coeff * sum_over_subsets
        h_star_linext *= first_sum_term
        h_star += h_star_linext
    return h_star