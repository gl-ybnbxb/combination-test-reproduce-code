#' The main function for the tranformation-based global testing method
#'
#' @param p.vector a vector of p-values with length \code{d} of p-values. \code{d} is the
#' number of hypotheses.
#' @param method one of "Cauchy", "Pareto", "log_Cauchy", "Frechet", "Levy", with the default method "Cauchy". 
#'
#' @param gamma when the method is 'Pareto' or 'Frechet', one paramater of the distribution is needed. 
#' The default value is gamma=1.
#' @return an aggregated p-value of the global test.
#'


transformation_global_test = function(p.vector, method='Cauchy',...){
  require(rmutil)
  
  params <- list(...)
  
  d <- length(p.vector)
  
  if(1 %in% p.vector){
    # cat('There are p-values that are exactly 1!')
    p.vector[p.vector==1]=max(0.999, 1-1/length(p.vector))
  }
  
  if(method =='Fisher'){
    chisq <- sum(-2*log(p.vector))
    p.global <- 1-pchisq(chisq, df = 2*d)
  }
  
  if(method == 'Bonferroni'){
    p.global <- min(1,min(p.vector)*d)
  }
  
  if(grepl('truncated_',method)){
    truncate.threshold=as.numeric(unlist(strsplit(method,'_'))[2])
    S = sum(1/length(p.vector)*tan((0.5-truncate.threshold*p.vector)*pi))
    p.global = min(pcauchy(S,lower.tail = F)/truncate.threshold,1)
  }
  
  if(method == 'Cauchy'){
    Sd <- sum(tan((0.5-p.vector)*pi)/d)
    p.global <- pcauchy(Sd,lower.tail = F)
    p.global <- min(p.global,1)
  }
  
  if(method == 'Pareto'){
    missing_gamma <- (length(params)==0)
    if(missing_gamma){
      gamma=1
    }else{
      gamma = params$gamma
    }
    
    # cdf of pareto
    ppareto <- function(x, x_m = 1, alpha = 1, lower.tail = T){
      if (lower.tail){
        y = 1-(x_m/x)^alpha
      }
      else{
        y = (x_m/x)^alpha
      }
      return(y)
    }
    
    # quantile function of pareto
    qpareto <- function(x, x_m = 1, alpha = 1){
      y = x_m/(1-x)^(1/alpha)
      return(y)
    }
    
    Sd <- sum(qpareto(1-p.vector, alpha=gamma))
    p.global <- d*ppareto(Sd,lower.tail = F, alpha=gamma)
    p.global <- min(p.global,1)
  }
  
  if(method == 'log Cauchy'){
    lcauchy_mid <- function(x){
      tan(pi*(x-0.5))
    }
    plcauchy_mid <- function(x){
      0.5-atan(x)/pi
    }
    lcauchy_trans <- function(x){
      p = length(x)
      trans_x = sapply(x,function(s) lcauchy_mid(1-s))
      m = max(trans_x)
      mid = m+log(sum(exp(trans_x-m)))
      p.lcauchy  = min(1,p*plcauchy_mid(mid))
      
      return(p.lcauchy)
    }
    
    p.global  = lcauchy_trans(p.vector)
  }
  
  if(method == 'Frechet'){
    missing_gamma <- (length(params)==0)
    if(missing_gamma){
      gamma=1
    }else{
      gamma = params$gamma
    }
    
    # cdf of frechet
    pfrechet <- function(x, m=0, s=1, alpha=1, lower.tail = T){
      if (lower.tail){
        y = exp(-((x-m)/s)^(-alpha))
      }
      else{
        y = 1-exp(-((x-m)/s)^(-alpha))
      }
      return(y)
    }
    # quantile function of frechet
    qfrechet <- function(x, m=0, s=1, alpha=1){
      m+s/(-log(x))^(1/alpha)
    }
    
    Sd = sum(qfrechet(1-p.vector,alpha=gamma))
    p.global = d*pfrechet(Sd, lower.tail = F,alpha=gamma)
    p.global = min(p.global,1)
  }
  
  if(method == 'Levy'){
    # Levy: m=0, s=1
    Sd = sum(qlevy(1-p.vector))
    p.global = d*(1-plevy(Sd))
    p.global = min(p.global,1)
  }
  
  return(p.global)
}
