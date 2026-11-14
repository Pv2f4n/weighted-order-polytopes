from PyNormaliz import *
import matplotlib.pyplot as plt
import statistics
import time
attach("main.sage")

def sage_to_normaliz_polynomial(poly):
    """
    Converts a multivariate polynomial in SageMath format to Normaliz format. 
    
    The ith variable of the parent polynomial ring R of poly, which can be accessed as the ith element of the tuple yielded by R.gens(), will be converted to the variable x[i] in the output.

    Example:
    sage: R.<x, y, z> = PolynomialRing(ZZ)
    sage: f = -3*x*y^2*z - 5*y^3*z + x*y^2 - 2*x^2 + 3*z + 1
    sage: sage_to_normaliz_polynomial(f)
    '-3 * x[1] * x[2]^2 * x[3] - 5 * x[2]^3 * x[3] + x[1] * x[2]^2 - 2 * x[1]^2 + 3 * x[3] + 1'
    """
    if poly.is_zero():
        return "0"
    terms = []
    for idx, (coeff, exp_tuple) in enumerate(zip(poly.coefficients(), poly.exponents())):
        var_factors = []
        for i, exp in enumerate(exp_tuple, start=1):
            if exp == 1:
                var_factors.append(f"x[{i}]")
            elif exp > 1:
                var_factors.append(f"x[{i}]^{exp}")
        is_negative = coeff < 0
        abs_coeff = abs(coeff)
        if not var_factors:
            term_body = str(abs_coeff)
        elif abs_coeff == 1:
            term_body = " * ".join(var_factors)
        else:
            term_body = f"{abs_coeff} * " + " * ".join(var_factors)

        if idx == 0:
            terms.append(f"-{term_body}" if is_negative else term_body)
        else:
            sign = " - " if is_negative else " + "
            terms.append(f"{sign}{term_body}")
    return "".join(terms)

def weighted_h_star_order_polytope_normaliz(P, weight):
    """
    Uses PyNormaliz' WeightedEhrhartSeries function to compute the weighted h-star polynomial of the order polytope O(P) of a poset P with respect to a polynomial weight. 

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
    normaliz_weight = sage_to_normaliz_polynomial(weight)
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
    Returns a tuple (time_1, time_2) where time_1 (resp. time_2) is the time in seconds it took weighted_h_star_order_polytope (resp. weighted_h_star_order_polytope_normaliz) to compute the weighted h-star polynomial of the order polytope O(P) of the poset P with respect to a polynomial weight.

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

def plot_execution_times(data):
    """
    Plots a log-log scatter plot comparing execution times of our formula to that of Normaliz, where data is a list of the form (num_lin_ext, (our_formula_time, normaliz_time)) with times in seconds.
    """
    # Unpack the data into separate lists for plotting
    x_vals = [entry[0] for entry in data]
    our_formula_times = [entry[1][0] for entry in data]
    normaliz_times = [entry[1][1] for entry in data]

    # Create the figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot the scatter points
    ax.scatter(x_vals, our_formula_times, color='blue', label='Our formula', alpha=0.8)
    ax.scatter(x_vals, normaliz_times, color='orange', label='Normaliz', alpha=0.8)

    ax.tick_params(axis='both', labelsize=14)

    # Set both axes to logarithmic scale
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Set axis labels
    ax.set_xlabel('Number of Linear Extensions', fontsize=14)
    ax.set_ylabel('Execution Time (seconds)', fontsize=14)

    # Add a legend (key) and grid for readability
    ax.legend(fontsize=14)
    ax.grid(True, which="both", linestyle='--', linewidth=0.5, alpha=0.7)

    # Adjust layout and display the plot
    plt.tight_layout()
    plt.show()