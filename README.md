# Code for Reproducibility

This is code for reproducing all figures and tables in the paper *Aggregating Dependent Signals with Heavy-Tailed Combination Tests*.


## Two helper files
- `simulation_helper.R` contains basic simulation functions to generate p-values and compute type one error or power of different methods.
- `plot_helper.R` contains the plotting function for all figures.

## Figures 
`figure1.R` is used to generates Figure 1.

To reproduce Figure 2 and S1. Remerber first change the parameters in `run_figures_s1.sh`. Then run
```
cd ./R/figure/error/
chmod +x run_figures_s1.sh
./run_figures_s1.sh
```
Then the plot will be automatically saved in this `error` folder.

To reproduce Figure 3 and 4, first, use `power_simulations.R` to generate data frames containing power of different methods varying signal levels. The data frames will be saved in the `figure` folder as csv files. Then, use `figure_3_4.R`:
- `summary_data` function computes the maximum power gain of each method over the Bonferroni's test. Example code below this function can aggregate power gain results under different settings.
- `summary_power_plot` function generates Figure 3 and 4. Example code below this function can generate Figure 3 and 4 for Gaussian and t copla separately.

To reproduce Figure S2, first use example code in `power_simulations.R` to generate data frames containing power of different methods varying signal levels. The data frames will be saved in the `figure` folder as csv files. Then, use example code in `figure_s2.R` to generate each sub figure in Figure S2. 

## Tables

`tables.R` can reproduce Table 2 and S1-S3 in the paper.
