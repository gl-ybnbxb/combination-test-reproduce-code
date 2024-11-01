#!/bin/bash


# change this part
DIM=5
RHO=0.5
ALPHA=0.05

# run the r script
Rscript --vanilla figure2_s1.R $DIM $RHO $ALPHA
