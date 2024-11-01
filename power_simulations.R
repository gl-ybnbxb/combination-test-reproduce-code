source("simulation_helper.R")

#' This is a function to generate simulated p-values and compute power
#'
#' @param signal_seq a sequence of number describing the signal level.
#' @param prop the number of 0s and nonzero mus in the mean vector.
#' @param rho the correlation in the correlation matrix
#' The default value is rho=0.
#' @param alpha the significance level for the global test.
#' @param copula dependence structure. 'gaussian' or 't'.
#' @param dof degree of freedom of multivariate t.
#' @param one_sided one-sided p-values or not.
#' @param reptimes the repeated times for simulations.
#' @return a vector contains empirical power of each method.
#'
power_simulator = function(signal_seq = seq(0, 6, by = 0.5),
                           prop = c(3,2),rho = 0, alpha = 5e-4,
                           copula = 't', dof = 2,
                           one_sided = T, reptimes = 1e5,){

  require(reshape2)
  dim = sum(prop)

  # compute power
  power.mat = t(sapply(signal_seq, function(u){
                    mu_vec = c(rep(0,prop[1]),rep(u,prop[2]));
                    power.vec=power_different_methods_simulator(mu_vec, rho = rho,
                                                                alpha = alpha,
                                                                copula = copula, dof = dof,
                                                                one_sided = one_sided,
                                                                reptimes = reptimes);
                    return(power.vec)}))
  row.names(power.mat) = signal_seq
  print('Power table done!')

  # melt data
  num_method = ncol(power.mat)
  power.dat = melt(power.mat)
  colnames(power.dat) = c("Mu","Method","Power")

  # save the power data
  if(prop[1]==0){
    power.filename = paste0("power_alpha_",alpha,"_cor_",rho,"_mean_(u_", prop[2],")_copula_",copula,".csv")
  }else{
    power.filename = paste0("power_alpha_",alpha,"_cor_",rho,"_mean_(0_",prop[1],",u_", prop[2],")_copula_",copula,".csv")
  }
  write.csv(power.dat, power.filename)
}

# d1: 0  1  0   18   0  95
# d2: 5  4  20  2    100   5
# alpha: 5e-2, 5e-4
# rho: 0, 0.5, 0.9, 0.99


# example code for Figure 3
power_simulator(prop = c(0,5),rho = 0.5, alpha = 0.05, copula = 'gaussian', reptimes = 1e6)

# example code for Figure 4
power_simulator(prop = c(0,5),rho = 0.5, alpha = 0.05, copula = 't', dof = 2, reptimes = 1e6)

# example code for Figure S2
power_simulator(signal_seq = seq(0, 0.5, by = 0.01), prop = c(4,1),
                rho = -0.2, alpha = 0.05, copula = 'gaussian', reptimes = 1e6)


############################################


