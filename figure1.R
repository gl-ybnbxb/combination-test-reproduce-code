source('../utils.R')
library(Hmisc)
library(colorspace)
library(scales)
cols = hue_pal()(100)

# This is the function which generates Figure 1
plot_rejection = function(alpha,
                          x=seq(-5,5,length=1000), y=seq(-5,5,length=1000),
                          method_list = c('Cauchy', 'Frechet'),
                          pareto_gamma_list = c(1/2,1),
                          frechet_gamma_list = seq(1,2,by=0.1)){
  require(rmutil)

  p1 = function(x) 2*pnorm(abs(x),lower.tail = F)
  p2 = function(y) 2*pnorm(abs(y),lower.tail = F)
  p_inv = function(x) 1-2*pnorm(abs(x),lower.tail = F)

  # Bonferroni
  z.bon<-outer(x,y,function(x,y) 4*pnorm(pmax(abs(x),abs(y)),lower.tail = F)-alpha)
  contour(x,y,z.bon,levels=0, col = 'black',drawlabels = F, lwd=4,cex.axis=2.5)

  # Fisher
  z.fisher<-outer(x,y,function(x,y) -2*log(p1(x))-2*log(p2(y))-qchisq(1-alpha,df = 4))
  contour(x,y,z.fisher,levels=0, col = '#619CFF',add = T,drawlabels = F, lwd=4)

  # Cauchy
  z.cauchy<-outer(x,y,function(x,y) {tan((2*pnorm(abs(x))-1.5)*pi) + tan((2*pnorm(abs(y))-1.5)*pi)-2*qcauchy(1-alpha)})
  contour(x,y,z.cauchy,levels=0,col = '#AD0826', add = T, drawlabels = F, lwd=4)

  # Frechet, gamma=1
  z.frechet<-outer(x,y,function(x,y) {qfrechet(p_inv(x))+qfrechet(p_inv(y))-qfrechet(1-alpha/2)})
  contour(x,y,z.frechet,levels=0, col = '#7CAE00',add = T,drawlabels = F, lwd=4)

  # title of the plot
  title(main = bquote(alpha==.(alpha)),cex.main=3)
}

# Figure 1
plot_rejection(5e-2)
plot_rejection(5e-3)
plot_rejection(5e-4)
