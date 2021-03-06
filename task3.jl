using Random,Distributions,Plots,QuadGK

# Define functions for calculating ll, prior, posterior distribution
# likelihood(t) = prod(pdf.(Normal(t, sigma_pop), x_sample))
# log_likelihood(t) = loglikelihood(Normal(t, sigma_pop), x_sample)
# prior_pdf(t) = pdf.(Normal(mean_prior, std_prior), t)
# ll_prior(x) = likelihood(x) * prior_pdf(x)
# post_pdf(t) = (1/k_val) .* ll_prior(t)

## Input data, units kN, cm, kPa
# Measured deflection
x_sample = 5 
d_unc = 1

# Prior distribution parameters (guess) in kPa
prior_mean = 3000
prior_std = 6000

# Deflection as a function of E
Q = 100
L = 10*100
I = 1e5 

d(ee) = Q * L^3 / (48 * ee * I)
E_mod(dd) = (Q*L^3) / (48 * I) * dd

# Likelihood function of deflection based on measurement
tt = range(0, stop=25000, length=1000)
ll_d(t) = prod(pdf.(Normal(d(t), d_unc), x_sample))

# MLE estimator
log_ll_d(t) = -loglikelihood(Normal(d(t), d_unc), x_sample)
E_mle = E_mod(minimum(log_ll_d.(tt)))

bbb = ll_d.(tt)
# Plot likelihood function
begin
    gr()
    plot(tt, ll_d.(tt),label="likelihood")
    xaxis!("E_mod d [kPa]")
end;

# Plot log-likelihood function
begin
    gr()
    plot!(tt, log_ll_d.(tt))
end;

# Define prior dist. and integration function
prior_pdf(t) = pdf.(Normal(prior_mean, prior_std), t)
ll_prior(x) = ll_d(x) * prior_pdf(x)

# Integrate normalization coeff.
k_val, _ = quadgk(ll_prior, 0, 100000)

# Calculate posterior distribution
post_pdf(t) = (1/k_val) .* ll_prior(t)

begin
    gr()
    plot(tt, prior_pdf.(tt), label="prior")
    plot!(tt, post_pdf.(tt), label="posterior")
end;

# Calculate posterior mean
post_pdf_t(t) = t .* post_pdf(t)
mean_post, _ = quadgk(post_pdf_t, 0, 1e5)
println("Posterior mean: $(mean_post)")
