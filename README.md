# EM_alg_for_MSN

EM_dim_1.R and EM_dim_morethan_1.R are two example codes for fitting modified skew-normal (MSN) distribution to simulated observations from the MSN distribution with location parameter xi, scale parameter Psi, and skewness parameter eta, for dimension one and dimension more than on, using the
EM algorithm.

Users can edit the stopping criteria of the EM algorithm in the "while" loop inside the code at their convenience.

The variables fin_xi, fin_Psi, and fin_eta return the final estimates of the parameters xi, Psi, and eta, respectively.

For fitting the MSN model to a dataset, replace the matrix X with the data matrix.
