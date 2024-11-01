#' This is a function to compute maximum power gain
#'
#' @param prop the number of 0s and nonzero mus in the mean vector.
#' @param rho the correlation in the correlation matrix
#' @param alpha the significance level for the global test.
#' @param copula dependence structure. 'gaussian' or 't'.
#' @return Figure S2
#'
plot_figure_s2 = function(prop = c(4,1), rho = -0.2, alpha = 5e-4, copula='gaussian'){

  require(ggplot2)
  require(reshape2)
  dim = sum(prop)

  # read and organize the data
  if(prop[1]==0){
    power.filename = paste0("power_alpha_",alpha,"_cor_",rho,"_mean_(u_", prop[2],")_copula_",copula,".csv")
  }else{
    power.filename = paste0("power_alpha_",alpha,"_cor_",rho,"_mean_(0_",prop[1],",u_", prop[2],")_copula_",copula,".csv")
  }

  power.mat = read.csv(power.filename)
  power.mat = power.mat[,-1] # remove the first column
  power.mat = power.mat[power.mat$Method!='Frechet',]
  power.mat = power.mat[power.mat$Method!='Levy',]
  power.mat = power.mat[power.mat$Method!='Fisher',]
  power.mat$Method = as.factor(power.mat$Method)
  levels(power.mat$Method) = c('Bonferroni','Cauchy','Pareto 1 (Frechet 1)','Truncated t1 0.9')
  power.mat$Method = as.character(power.mat$Method)

  # plot the power
  p = ggplot(power.mat, aes(Mu, Power, group = Method, color = Method)) +
    geom_line(aes(linetype=Method),lwd = 1.5) +
    theme_bw() +
    scale_linetype_manual(values=c(rep(1,3), 6))+
    scale_color_manual(values=c('black','#AD0826',"#C77CFF","#EA5839"))+
    theme(text = element_text(size = 25),
          axis.text = element_text(size = 25),
          legend.direction = 'horizontal',
          #legend.position = 'none',
          legend.position = "bottom")+
    labs(title = paste0("Correlation = ",rho), x=expression(paste(mu)))+
    geom_hline(yintercept = alpha, linetype = 2)

  return(p)
}

# Example code to generate Figure S2
power_plot(prop = c(4,1), rho = -0.2, alpha = 5e-2)
power_plot(prop = c(4,1), rho = -0.2, alpha = 5e-4)
power_plot(prop = c(0,5), rho = -0.2, alpha = 5e-2)
power_plot(prop = c(0,5), rho = -0.2, alpha = 5e-4)


