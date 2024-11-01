source("simulation_helper.R")
source("plot_helper.R")

args <- commandArgs(trailingOnly = TRUE)
d = as.integer(args[1])
r = as.numeric(args[2])
a = as.numeric(args[3])

cat('Number of hypotheses:',d,'\n')
cat('Correlation:',r,'\n')
cat('Significance level:',a,'\n\n')

# run simulations
set.seed(12)
error_df = validity_vs_tail_index_simulator(p=d, rho = r, alpha = a,
                                             reptimes = 1e6, gamma_list = seq(0.7,1.5,by=0.01))

# generate Figure 2 & S1
plot_figure2_s1(error_df, rho=r, alpha=a, dim=d)
