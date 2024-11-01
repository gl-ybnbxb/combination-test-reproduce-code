source("simulation_helper.R")


# Reproduce all tables
# Table 2 & S3
set.seed(123)
validity_different_methods_simulator(p=5,rho = 0, alpha = 0.05, reptimes = 1e6, copula='t', dof=2, one_sided = T)
validity_different_methods_simulator(p=5,rho = 0.5, alpha = 0.05, reptimes = 1e6, copula='t', dof=2, one_sided = T)
validity_different_methods_simulator(p=5,rho = 0.9, alpha = 0.05, reptimes = 1e6, copula='t', dof=2, one_sided = T)
validity_different_methods_simulator(p=5,rho = 0.99, alpha = 0.05, reptimes = 1e6, copula='t', dof=2, one_sided = T)

validity_different_methods_simulator(p=5,rho = 0, alpha = 5e-4, reptimes = 1e6, copula='t', dof=2, one_sided = T)
validity_different_methods_simulator(p=5,rho = 0.5, alpha = 5e-4, reptimes = 1e6, copula='t', dof=2, one_sided = T)
validity_different_methods_simulator(p=5,rho = 0.9, alpha = 5e-4, reptimes = 1e6, copula='t', dof=2, one_sided = T)
validity_different_methods_simulator(p=5,rho = 0.99, alpha = 5e-4, reptimes = 1e6, copula='t', dof=2, one_sided = T)

validity_different_methods_simulator(p=100,rho = 0, alpha = 0.05, reptimes = 1e6, copula='t', dof=2, one_sided = T)
validity_different_methods_simulator(p=100,rho = 0.5, alpha = 0.05, reptimes = 1e6, copula='t', dof=2, one_sided = T)
validity_different_methods_simulator(p=100,rho = 0.9, alpha = 0.05, reptimes = 1e6, copula='t', dof=2, one_sided = T)
validity_different_methods_simulator(p=100,rho = 0.99, alpha = 0.05, reptimes = 1e6, copula='t', dof=2, one_sided = T)

validity_different_methods_simulator(p=100,rho = 0, alpha = 5e-4, reptimes = 1e6, copula='t', dof=2, one_sided = T)
validity_different_methods_simulator(p=100,rho = 0.5, alpha = 5e-4, reptimes = 1e6, copula='t', dof=2, one_sided = T)
validity_different_methods_simulator(p=100,rho = 0.9, alpha = 5e-4, reptimes = 1e6, copula='t', dof=2, one_sided = T)
validity_different_methods_simulator(p=100,rho = 0.99, alpha = 5e-4, reptimes = 1e6, copula='t', dof=2, one_sided = T)


# Table S1
set.seed(1)
validity_different_methods_simulator(p=2,rho = -0.5, alpha = 0.05, reptimes = 5e4, one_sided = T)
validity_different_methods_simulator(p=2,rho = -0.9, alpha = 0.05, reptimes = 5e4, one_sided = T)
validity_different_methods_simulator(p=2,rho = -0.99, alpha = 0.05, reptimes = 5e4, one_sided = T)


# Table S2
# alpha decreases from 5e-2 to 5e-8
a = 0.05
N=1.3e10
za = qnorm(1-a/2)
error_tab = validity_different_methods_simulator(p=2,rho = -0.9, alpha = a, reptimes = N, one_sided = T)

# Wilson binomial confidence interval
p_hat = error_tab[1,]
first_term = round((p_hat+za^2/2/N)/(1+za^2/N),digits = 3)
second_term = round(za*sqrt(p_hat*(1-p_hat)/N+za^2/4/N^2)/(1+za^2/N),digits = 3)
c(first_term-second_term, first_term+second_term)

