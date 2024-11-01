source("../utils.R")
source("../simulation_function.R")

#' Simulation function for investgating connection between tail index and validity.
#'
#' @param p number of base hypotheses
#' @param rho the correlation in the correlation matrix
#' The default value is rho=0.
#' @param alpha the significance level for the global test.
#' @param gamma_list the tail index sequence.
#' @param reptimes the repeated times for simulations.
#' @return a data frame whose each row has type-I error and method name of each method.
#'
validity_vs_tail_index_simulator = function(p, rho = 0, alpha = 5e-4,
                                            reptimes = 1e5, gamma_list = seq(1,2,by=0.01)){
  require(mvtnorm)
  require(rmutil)
  require(invgamma)

  mu = rep(0,p)
  test.dat = data_Gen(n = reptimes, p = p, mu = mu, rho = rho, copula = 'gaussian')
  p.mat = test.dat$p.mat

  # frechet with different alpha
  p_frechet = sapply(gamma_list, function(y){
    S.frechet = apply(p.mat, 1, function(x) sum(qfrechet(1-x,alpha=y)))
    p.frechet = p*pfrechet(S.frechet, lower.tail = F,alpha=y)
    p.frechet = sapply(p.frechet, function(x) min(x,1))
  })

  # t with different dof
  p_t = sapply(gamma_list, function(y){
    S.t = apply(p.mat, 1, function(x) sum(qt(1-x,df=y)))
    p.t = p*pt(S.t,df = y,lower.tail = F)
    p.t = sapply(p.t, function(x) min(x,1))
  })

  # Pareto with different gamma
  p_pareto = sapply(gamma_list, function(y){
    S.pareto = apply(p.mat, 1, function(x) sum(qpareto(1-x,alpha=y)))
    p.pareto = p*ppareto(S.pareto, lower.tail = F,alpha=y)
    p.pareto = sapply(p.pareto, function(x) min(x,1))
  })

  # inverse gamma with different dof
  p_invgamma = sapply(gamma_list, function(y){
    S.invgamma = apply(p.mat, 1, function(x) sum(qinvgamma(1-x,shape = y)))
    p.invgamma = p*pinvgamma(S.invgamma,shape = y,lower.tail = F)
    p.invgamma = sapply(p.invgamma, function(x) min(x,1))
  })

  # Cauchy
  T.cauchy = apply(p.mat, 1, function(x) sum(tan((0.5-x)*pi)/p))
  p.cauchy = pcauchy(T.cauchy,lower.tail = F)
  p.cauchy = sapply(p.cauchy, function(x) min(x,1))

  # Truncated t with the threshold 0.5
  p_tt05 = sapply(gamma_list, function(y){
    S.tt05 = apply(p.mat, 1, function(x){sum(qt(1-0.5*x,df=y))})
    p.tt05 = p*pt(S.tt05,df = y,lower.tail = F)/0.5
    p.tt05 = sapply(p.tt05, function(x) min(x,1))
  })

  # Truncated t with the threshold 0.7
  p_tt07 = sapply(gamma_list, function(y){
    S.tt07 = apply(p.mat, 1, function(x){sum(qt(1-0.7*x,df=y))})
    p.tt07 = p*pt(S.tt07,df = y,lower.tail = F)/0.7
    p.tt07 = sapply(p.tt07, function(x) min(x,1))
  })

  # Truncated t with the threshold 0.9
  p_tt09 = sapply(gamma_list, function(y){
    S.tt09 = apply(p.mat, 1, function(x){sum(qt(1-0.9*x,df=y))})
    p.tt09 = p*pt(S.tt09,df = y,lower.tail = F)/0.9
    p.tt09 = sapply(p.tt09, function(x) min(x,1))
  })

  # p value matrix of aggregated p-values
  p.summary = cbind(p_frechet,
                    p_t,
                    p_pareto,
                    p_invgamma,
                    p_tt05,
                    p_tt07,
                    p_tt09,
                    p.cauchy)

  # type I error
  pnt.error = apply(p.summary, 2, function(x) mean(x<=alpha))
  sd.error = apply(p.summary, 2, function(x) sd(x<=alpha)/sqrt(reptimes))

  # table of type I error
  error_df = data.frame(Error=pnt.error, SD=sd.error)
  error_df$gamma = c(rep(gamma_list, 7),1)
  name_list = c('Frechet', 't', 'Pareto', 'Inverse_gamma',
                'Truncated t 0.5', 'Truncated t 0.7', 'Truncated t 0.9')
  error_df$Method = c(rep(name_list, each = length(gamma_list)), 'Cauchy')

  return(error_df)

}


#' The function for testing validity of global tests
#'
#' @param p the number of p-values in one p-value vector. \code{p} is the
#' number of hypotheses.
#' @param rho the correlation in the correlation matrix
#' The default value is rho=0.
#' @param alpha the significance level for the global test.
#' @param gamma_list the tail index sequence.
#' @param reptimes the repeated times for simulations.
#' @return a matrix containing empirical type-I errors (with standard error) of different methods.
#'
validity_different_methods_simulator = function(p, rho = 0, alpha = 5e-4,
                                                copula = 'gaussian', dof = 5,
                                                one_sided = F, reptimes = 1e5){
  require(mvtnorm)
  require(rmutil)
  require(invgamma)

  mu = rep(0,p)
  test.dat = data_Gen(n = reptimes, p = p, mu = mu, rho = rho,
                      copula = copula, dof = dof,
                      one_sided = one_sided)
  p.mat = test.dat$p.mat

  # Cauchy
  T.cauchy = apply(p.mat, 1, function(x) sum(tan((0.5-x)*pi)/p))
  p.cauchy = pcauchy(T.cauchy,lower.tail = F)
  p.cauchy = sapply(p.cauchy, function(x) min(x,1))

  # Pareto: alpha=1
  S.pareto = apply(p.mat, 1, function(x) sum(qpareto(1-x)))
  p.pareto = p*ppareto(S.pareto,lower.tail = F)
  p.pareto = sapply(p.pareto, function(x) min(x,1))

  # Truncated Cauchy:
  T.tcauchy = apply(p.mat, 1, function(x) sum(tan((0.5-0.9*x)*pi)/p))
  p.tcauchy = pcauchy(T.tcauchy,lower.tail = F)/0.9
  p.tcauchy = sapply(p.tcauchy, function(x) min(x,1))

  # Frechet:
  S.frechet = apply(p.mat, 1, function(x) sum(qfrechet(1-x)))
  p.frechet = p*pfrechet(S.frechet, lower.tail = F)
  p.frechet = sapply(p.frechet, function(x) min(x,1))

  # Levy: m=0, s=1
  S.levy = apply(p.mat, 1, function(x) sum(qlevy(1-x)))
  p.levy = p*(1-plevy(S.levy))
  p.levy = sapply(p.levy, function(x) min(x,1))

  # Bonferroni
  p.bon = apply(p.mat, 1, function(x) min(p*min(x),1))

  # Fisher
  T.fisher = apply(p.mat, 1, function(x) -2*sum(log(x)))
  p.fisher = pchisq(T.fisher,df = 2*p,lower.tail = F)

  # p value matrix of aggregated p-values
  p.summary = cbind(Cauchy = p.cauchy,
                    Pareto = p.pareto,
                    Truncated_Cauchy = p.tcauchy,
                    Frechet  = p.frechet,
                    Levy = p.levy,
                    Bonferroni = p.bon,
                    Fisher = p.fisher)

  # type I error
  pnt.error = apply(p.summary, 2, function(x) mean(x<=alpha))
  sd.error = apply(p.summary, 2, function(x) sd(x<=alpha)/sqrt(reptimes))
  error_tab = rbind(Error=pnt.error, SD=sd.error)

  return(error_tab)

}


#' The function for investigating power of global tests
#'
#' @param mu the mean vector for the multivariate gaussian/t
#' @param rho the correlation in the correlation matrix
#' The default value is rho=0.
#' @param alpha the significance level for the global test.
#' @param copula dependence structure. 'gaussian' or 't'.
#' @param dof degree of freedom of multivariate t.
#' @param one_sided one-sided p-values or not. The default value is False.
#' @param reptimes the repeated times for simulations.
#' @return a data frame whose each row has type-I error and method name of each method.
#'
power_different_methods_simulator = function(mu, rho = 0, alpha = 5e-4,
                                            copula = 'gaussian', dof = 5,
                                            one_sided = F, reptimes = 1e5){
  require(mvtnorm)
  require(rmutil)
  require(invgamma)

  p = length(mu)
  test.dat = data_Gen(n = reptimes, p = p, mu = mu, rho = rho,
                      copula = copula, dof = dof,
                      one_sided = one_sided)
  p.mat = test.dat$p.mat

  # Cauchy
  T.cauchy = apply(p.mat, 1, function(x) sum(tan((0.5-x)*pi)/p))
  p.cauchy = pcauchy(T.cauchy,lower.tail = F)
  p.cauchy = sapply(p.cauchy, function(x) min(x,1))

  # Pareto: alpha=1
  S.pareto = apply(p.mat, 1, function(x) sum(qpareto(1-x)))
  p.pareto = p*ppareto(S.pareto,lower.tail = F)
  p.pareto = sapply(p.pareto, function(x) min(x,1))

  # Truncated Cauchy:
  T.tcauchy = apply(p.mat, 1, function(x) sum(tan((0.5-0.9*x)*pi)/p))
  p.tcauchy = pcauchy(T.tcauchy,lower.tail = F)/0.9
  p.tcauchy = sapply(p.tcauchy, function(x) min(x,1))

  # Frechet:
  S.frechet = apply(p.mat, 1, function(x) sum(qfrechet(1-x)))
  p.frechet = p*pfrechet(S.frechet, lower.tail = F)
  p.frechet = sapply(p.frechet, function(x) min(x,1))

  # Levy: m=0, s=1
  if(copula=='gaussian'){
    S.levy = apply(p.mat, 1, function(x) sum(qlevy(1-x)))
    p.levy = p*(1-plevy(S.levy))
    p.levy = sapply(p.levy, function(x) min(x,1))
  }

  # Bonferroni
  p.bon = apply(p.mat, 1, function(x) min(p*min(x),1))

  # Fisher
  T.fisher = apply(p.mat, 1, function(x) -2*sum(log(x)))
  p.fisher = pchisq(T.fisher,df = 2*p,lower.tail = F)

  # p value matrix of aggregated p-values
  if(copula=='gaussian'){
    p.summary = cbind(Cauchy = p.cauchy,
                      Pareto = p.pareto,
                      Truncated_Cauchy = p.tcauchy,
                      Frechet  = p.frechet,
                      Levy = p.levy,
                      Bonferroni = p.bon,
                      Fisher = p.fisher)
  }else{
    p.summary = cbind(Cauchy = p.cauchy,
                      Pareto = p.pareto,
                      Truncated_Cauchy = p.tcauchy,
                      Frechet  = p.frechet,
                      Bonferroni = p.bon,
                      Fisher = p.fisher)
  }

  # compute power
  pnt.power = apply(p.summary, 2, function(x) mean(x<=alpha))

  return(pnt.power)

}

