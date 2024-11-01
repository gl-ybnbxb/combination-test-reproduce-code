# plot function for Figure 2 and S1
plot_figure2_s1 = function(data_plot, rho, alpha, dim){
  require(ggplot2)

  n = nrow(data_plot)
  pnt_data = data_plot[n,]
  data_plot = data_plot[-n,]

  if(dim == 5){
    if(alpha==0.05){
      y_lim_s = -1.6
      y_lim_l = -0.9
      y_scale_int = seq(-1.6,-0.9,by = 0.3)
    }else{
      y_lim_s = -3.6
      y_lim_l = -2.9
      y_scale_int = seq(-3.6,-1.3,by=0.3)
    }
  }else if(dim == 100){
    if(alpha==0.05){
      y_lim_s = -1.9
      y_lim_l = 0
      y_scale_int = seq(-1.9,0,by = 0.6)
    }else{
      y_lim_s = -3.9
      y_lim_l = -2.25
      y_scale_int = seq(-3.9,-2.3, by = 0.4)
    }
  }

  # plot the error
  p = ggplot(data_plot, aes(gamma, log10(Error), group = Method, color = Method)) +
    geom_line(lwd = 2) +
    geom_point(data = pnt_data, mapping=aes(x=gamma,y=log10(Error),group=Method, color=Method), size=30, shape='*')+
    theme_bw() +
    scale_x_continuous(breaks = seq(0.7, 1.5, by = 0.2))+
    scale_y_continuous(breaks = y_scale_int,
                       labels=format(10^y_scale_int,scientific = T,digits=2),
                       limits = c(y_lim_s,y_lim_l))+
    labs(title = paste0("Correlation = ",rho), x=expression(paste(gamma)), y='Error')+
    geom_hline(yintercept = log10(alpha), linetype = 2)+
    theme(text = element_text(size = 40),
          axis.text = element_text(size = 35),
          axis.title.x = element_text(size = 40),
          legend.direction = 'horizontal',
          legend.position = 'bottom')+
    scale_color_manual(values=c('#AD0826', "#7CAE00", "#8494FF", "#C77CFF",
                                '#AD0826', "#FDAE61", "#F7824D", "#EA5839"),
                       guide = guide_legend(override.aes = list(linetype = c('blank',rep('solid',7)),
                                                                shape = c('*',rep('',7)))))

  # save the plot
  filename = paste0("validity_n_",dim,"_error_",alpha,"_correlation_",rho,'.pdf')
  pathname = getwd()
  ggsave(filename, p, width = 20, height = 10, path = pathname)
}




