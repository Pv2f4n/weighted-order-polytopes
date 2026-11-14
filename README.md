# Computing weighted h-star polynomials of order polytopes
A SageMath implementation of the combinatorial formula for the weighted h-star polynomial of an order polytope presented in forthcoming work with Daniel Hwang, David Jin, Juliet Whidden, and Katelyn Zhou.

main.sage contains the implementation, and normaliz_comparison.sage contains the code that was used to benchmark it against Normaliz' WeightedEhrhartSeries function. 

AI disclosure: Gemini 3.1 Pro was used to write the functions `sage_to_normaliz_polynomial()` and `plot_execution_times()`. All other code was written by hand without the assistance of AI tools.