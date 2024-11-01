library(ggplot2)

#######################################
###      Aggregate data for plot    ###
#######################################

#' This is a function to compute maximum power gain
#'
#' @param alpha the significance level for the global test.
#' @param d1 the number of 0s in the mean vector.
#' @param d2 the number of nonzero mus in the mean vector.
#' @param rho the correlation in the correlation matrix
#' @param copula dependence structure. 'gaussian' or 't'.
#' @return a data frame contains maximum power gain of each method over Bonferroni's test.
#'
summary_data = function(alpha,d1,d2,rho,copula='gaussian'){
  dim=d1+d2
  if(d1==0){
    df.name = paste0("power_alpha_",alpha,"_cor_",rho,"_mean_(u_", d2,")_copula_",copula,".csv")
    type='dense'
  }else{
    df.name = paste0("power_alpha_",alpha,"_cor_",rho,"_mean_(0_",d1,",u_", d2,")_copula_",copula,".csv")
    type='sparse'
  }

  df = read.csv(df.name)
  df = df[df$Mu!=0,]
  power_diff_cauchy = max(df$Power[df$Method=='Cauchy'] - df$Power[df$Method=='Bonferroni'])
  power_diff_pareto = max(df$Power[df$Method=='Pareto'] - df$Power[df$Method=='Bonferroni'])
  power_diff_frechet = max(df$Power[df$Method=='Frechet'] - df$Power[df$Method=='Bonferroni'])
  power_diff_levy = max(df$Power[df$Method=='Levy'] - df$Power[df$Method=='Bonferroni'])
  power_diff_tcauchy = max(df$Power[df$Method=='Truncated_Cauchy'] - df$Power[df$Method=='Bonferroni'])

  if(copula=='gaussian'){
    summary_df = data.frame(n = rep(dim, 5),
                            max_power_diff = c(power_diff_cauchy, power_diff_pareto,
                                               power_diff_frechet, power_diff_levy,
                                               power_diff_tcauchy),
                            rho = rep(rho, 5),
                            Method = c('Cauchy','Pareto','Frechet','Levy','Truncated Cauchy 0.9'),
                            alpha = rep(alpha, 5),
                            signal = rep(type, 5))
  }else{
    summary_df = data.frame(n = rep(dim, 4),
                            max_power_diff = c(power_diff_cauchy, power_diff_pareto,
                                               power_diff_frechet, power_diff_tcauchy),
                            rho = rep(rho, 4),
                            Method = c('Cauchy','Pareto','Frechet','Truncated Cauchy 0.9'),
                            alpha = rep(alpha, 4),
                            signal = rep(type, 4))
  }


  return(summary_df)

}

# Example code to aggregate data from different settings
dim_matrix = matrix(c(0,5,0,20,0,100,4,1,19,1,95,5),ncol=2,byrow = T)

summary_table = data.frame()
for (a in c(5e-2,5e-4)) {
  for (i in 1:nrow(dim_matrix)) {
    d1 = dim_matrix[i,1]
    d2 = dim_matrix[i,2]
    if(d1+d2==5){
      rho_seq = c(-0.2,0,0.5,0.9,0.99)
    }else{
      rho_seq = c(0,0.5,0.9,0.99)
    }
    for (r in rho_seq) {
      summary_table = rbind(summary_table, summary_data(alpha=a,d1=d1,d2=d2,rho=r,copula='t'))
    }
  }
}
View(summary_table)
colnames(summary_table) = c('n','Maximum of Power Difference',
                            'Correlation', 'Method', 'Alpha', 'Signal Type')
summary_table$Method[summary_table$Method=='Frechet'] = 'Frechet 1'
summary_table$Method[summary_table$Method=='Pareto'] = 'Pareto 1'
summary_table$Method[summary_table$Method=='Truncated Cauchy 0.9'] = 'Truncated t1 0.9'


#######################################
###           Figure 3 & 4          ###
#######################################

#' This is a function to generate Figure 3 and 4
#'
#' @param alpha the significance level for the global test.
#' @param type the signal type. 'dense' or 'sparse'.
#' @param data the summary data of maximum power gain in the previous step.
#' @param copula dependence structure. 'gaussian' or 't'.
#' @return Figure 3 or 4.
#'
summary_power_plot = function(alpha, type, data = summary_table, copula = 'gaussian'){
  subdf = data[data$Alpha==alpha & data$`Signal Type`==type,]
  if(copula=='gaussian'){
    y_lim_s = 0
    if(alpha == 0.05){
      y_scale_int = seq(0,0.6,by = 0.1)
      y_lim_l = 0.6
    }else{
      y_scale_int = seq(0,0.4,by = 0.1)
      y_lim_l = 0.4
    }
  }else{
    y_scale_int = seq(0,1,by = 0.2)
    y_lim_l = 1
    y_lim_s = 0
  }

  if(copula == 'gaussian'){
    pd = position_dodge(0.4)
    p = ggplot(subdf, aes(x=as.factor(Correlation), y=`Maximum of Power Difference`, group=Method, color=Method, shape=Method, fill=Method)) +
      geom_point(size = 8,alpha=0.5) +
      geom_line(lwd = 3,alpha=0.5)+
      facet_wrap(. ~ n, scales = "free_x")+
      xlab('Correlation')+
      ylab('Max. power gain\n compared to Bonferroni')+
      theme_minimal()+
      theme(axis.line = element_line(color="black"),
            text = element_text(size = 35),
            axis.text = element_text(size = 30),
            axis.title.y = element_text(size=40),
            strip.text = element_text(size = 40),
            #legend.direction = "horizontal",
            legend.position = "none",
            panel.spacing = unit(5, "lines"))+
      scale_y_continuous(breaks = y_scale_int,limits = c(y_lim_s,y_lim_l))+
      scale_shape_manual(values = 21:25)+
      scale_color_manual(values=c("#AD0826", "#7CAE00", "#00BFC4","#C77CFF","#EA5839"))
  }else{
    pd = position_dodge(0.4)
    p = ggplot(subdf, aes(x=as.factor(Correlation), y=`Maximum of Power Difference`, group=Method, color=Method, shape=Method, fill=Method)) +
      geom_point(size = 8,alpha=0.5) +
      geom_line(lwd = 3,alpha=0.5)+
      facet_wrap(. ~ n, scales = "free_x")+
      xlab('Correlation')+
      ylab('Max. power gain\n compared to Bonferroni')+
      theme_minimal()+
      theme(axis.line = element_line(color="black"),
            text = element_text(size = 35),
            axis.text = element_text(size = 30),
            axis.title.y = element_text(size=40),
            strip.text = element_text(size = 40),
            legend.direction = "horizontal",
            legend.position = "bottom",
            panel.spacing = unit(5, "lines"))+
      scale_y_continuous(breaks = y_scale_int,limits = c(y_lim_s,y_lim_l))+
      scale_shape_manual(values = c(21,22,24,25))+
      scale_color_manual(values=c("#AD0826", "#7CAE00","#C77CFF","#EA5839"))
  }

  return(p)
}

# Example code for Gaussian
summary_power_plot(alpha=0.05, type='dense', data = summary_table)
summary_power_plot(alpha=0.05, type='sparse', data = summary_table)

# Example code for t
summary_power_plot(alpha=5e-4, type='dense', data = summary_table, copula = 't')
summary_power_plot(alpha=5e-4, type='sparse', data = summary_table, copula = 't')

