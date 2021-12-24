using Random, Distributions,Plots, QuadGK

# Input variables
x_sample = [30.1, 28.7, 21.2, 28.6, 25.8, 25.1, 24.2, 27.3, 20.0, 34.8];
sigma_pop = 5;
mean_prior_mean = 40;
mean_prior_std = 100;

# Define log-likelihood function
log_likelihood(t) = loglikelihood(Normal(t, sigma_pop), x_sample)

# Calcualte the MLE estimator
tt = range(10, stop=50, length=100);
log_ll = -log_likelihood.(tt);
mle = minimum(log_ll)
display("MLE gives $(mle)")

# Display Log Likelihood and MLE value
begin
    gr()
    scatter(tt, log_ll, label="log_ll");
    vline!([mle], label="MLE")        
end

# Likelihood function
likelihood(t) = prod(pdf.(Normal(t, sigma_pop), x_sample))
begin
    gr()
    plot(tt, likelihood.(tt), label="ll")
end

# Define prior and integration function
prior_pdf(t) = pdf.(Normal(mean_prior_mean, mean_prior_std), t)
ll_prior(x) = likelihood(x) * prior_pdf(x)

# Calculate normalization coefficient k_val
k_val, _ = quadgk(ll_prior, -Inf, Inf)

# Calculate posterior definition
post_pdf(t) = (1/k_val) .* ll_prior(t)
begin
    gr()
    plot(tt, prior_pdf.(tt), label="prior")
    plot!(tt, post_pdf.(tt), label="posterior")
end

# Calculate posterior mean 
post_mean(t) = t .* post_pdf(t)
mean_post_mean, _ = quadgk(post_mean, -Inf, Inf)
println("Posterior mean: $(mean_post_mean)")
