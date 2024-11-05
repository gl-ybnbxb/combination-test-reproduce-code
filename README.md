# Code for Reproducibility

This is code for reproducing all figures and tables in the paper *Aggregating Dependent Signals with Heavy-Tailed Combination Tests*.

First, use the following code to install the package `heavytailcombtest`
```
library(devtools)
install_github("gl-ybnbxb/heavytailcombtest")
```

Note that the default directory of all commands below is the working directory.

## Two helper files
- `simulation_helper.R` contains basic simulation functions to generate p-values and compute type one error or power of different methods.
- `plot_helper.R` contains the plotting function for all figures.

## Figures 
`figure1.R` is used to generates Figure 1.

To reproduce Figure 2 and S1. Remember first change the parameters in `run_figures_s1.sh`. Then run
```
mkdir figure
chmod +x run_figures_s1.sh
./run_figures_s1.sh
```
Then the plot will be automatically saved in this `figure` folder.

To reproduce Figure 3 and 4, first make the folder `data` and then use `power_simulations.R` to generate data frames containing power of different methods varying signal levels. The data frames will be saved in the `data` folder as csv files. Then, use `figure_3_4.R`:
- `summary_data` function computes the maximum power gain of each method over the Bonferroni's test. Example code below this function can aggregate power gain results under different settings.
- `summary_power_plot` function generates Figure 3 and 4. Example code below this function can generate Figure 3 and 4 for Gaussian and t copla separately.

To reproduce Figure S2, first make the folder `data` and then use example code in `power_simulations.R` to generate data frames containing power of different methods varying signal levels. The data frames will be saved in the `data` folder as csv files. Then, use example code in `figure_s2.R` to generate each sub figure in Figure S2. 

To reproduce Figures 5 and 6, first make the folder `figure` and then use scripts `Figure_5_circadian_v1.R` and `Figure_6_GWAS_v3.R`. Data and additional utility and helper functions are in folder data_scripts.

## Tables

`tables.R` can reproduce Table 2 and S1-S3 in the paper.
