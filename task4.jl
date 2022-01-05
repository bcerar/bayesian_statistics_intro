using Random, Distributions, Plots, QuadGK

## Input variables
# First sensor, midspan
x_sample_1 = 5

# Second sensor, x = L/4 (from right)
x_sample_2 = 3.5

# Structure Input, units kN, cm, kPa
Q = 100
L = 10 .* 100
I = 1e5

xL = L / 4

# Deflection functions
d1_f(ee) = Q * L^3 / (48 * ee * I)
d2_f(ee, xL) = Q .* xL / (48 .* ee .* I) .* (3*L^2 - 4*xL^2)

## Multivariate normal distribution of model uncertanties
# Means array
mu(t) = [d1_f(t), d2_f(t, xL)]

# Covariance array
sig = [1 0.5*1*1; 0.5*1*1 1]

## Likelihood 
tt = range(0, stop=25000, length=1000)
ll_d(t) = prod(pdf.(MultivariateNormal(mu(t), sig), [x_sample_1; x_sample_2]))

# Plot likelihood function
begin
    gr()
    plot(tt, ll_d.(tt),label="likelihood")
    xaxis!("E_mod d [kPa]")
end;
