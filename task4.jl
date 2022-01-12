using Random
using Distributions
using Plots
using QuadGK

## Input variables
# First sensor, midspan
x_sample_1 = 5

# Second sensor, x = L/4 (from right)
x_sample_2 = 3.5

# Structure Input, units kN, cm, kPa
Q = 100
L = 10 .* 100
I = 1e5

prior_E_mean = 5e3
prior_E_std = 10e3

xL = L / 4

# Deflection functions
d1_f(ee) = Q * L^3 / (48 * ee * I)
d2_f(ee, xL) = Q .* xL / (48 .* ee .* I) .* (3*L^2 - 4*xL^2)

## Multivariate normal distribution of model uncertanties
# Means array

# Covariance array
sig = [1 0.5*1*1; 0.5*1*1 1]
means(t) = [d1_f(t), d2_f(t, xL)]

## Likelihood 
tt = range(0, stop=25000, length=1000)
ll_d(t) = pdf(MvNormal(means(t), sig), [x_sample_1; x_sample_2])

# Plot likelihood function
begin
    gr()
    plot(tt, ll_d.(tt),label="likelihood")
    xaxis!("E_mod d [kPa]")
end

# Prior 
prior_E_pdf(t) = pdf.(Normal(prior_E_mean, prior_E_std), t)
ll_prior(t) = ll_d(t) * prior_E_pdf(t)

# Normalizaiton coeff.
k_val,_ = quadgk(ll_prior, 0, 10e3)

# E posterior dist. function
post_pdf(t) = (1/k_val) .* ll_prior(t)

begin
    gr()
    plot(tt, prior_E_pdf.(tt), label="prior")
    plot!(tt, post_pdf.(tt), label="posterior")
end


# Exploring different correlations
rho_corr(dx, c) = exp(-dx/c)
xi = range(-3, stop=1, length=1000)
xi_i = xi' .* ones(length(xi))
c = -1:0.2:0.8
c_i = ones(length(c))' .* c
begin
    gr()
    for c_val in [0.5]
        plot(xi, rho_corr.(xi, c_val), label="c= $(c_val)")
    end
    ylabel!("rho")
    xlabel!("sensor distance")
end