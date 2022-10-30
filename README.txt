Copyright (C) 2020-2022 Justin A Fishbone

This package contains scripts to reproduce the figures and tables in the associated paper:
Fishbone, J.A., and Mili L. (2023?) New Highly Efficient High-Breakdown Estimator of Multivariate Scatter and Location for Elliptical Distributions. Canadian Journal of Statistics.


INSTALLATION
1) Install SEstimator estimator package from:
    https://github.com/JAFishbone/SEstimator
2a) Download KSD estimator from:
    https://doi.org/10.1016/j.csda.2016.11.006
        Supplementary material: 
            MMC S2. “Codes”: Contains Matlab code for the proposed estimators.
2b) Add the directory containing these files to the path: addpath(...);
3) Run add_paths from this directory


EXECUTION
0a) To generate Figure 1, execute: plot_rhos_and_ws
0b) To generate Figure 2, execute: plot_ws
This takes less than 1 minute to run
This script outputs plot files to the current directory

1) To generate theoretical asymptotic efficiency plots, execute:
    run_plot_asymptotic_efficiencies;
This takes approximately 10 minutes to complete on an Intel Core i7-8850h
This script outputs plot files to the current directory

2a) To generate low-resolution (lower sample spacings and number of trials) versions of Figures 5b and 8, execute:
    quick_run_simulations;
This takes less than an hour to complete on an Intel Core i7-8850h
This script outputs plot files to the current directory

2b) To generate full-resolution versions of Figures 5b and 8, execute:
    run_simulations;
This takes approximately 2 days to complete on an Intel Core i7-8850h
This script outputs plot files to the current directory

3) To generate influence function figure, execute:
    run_plot_influence_functions;
This takes less than 1 minute to run
This script outputs plot files to the current directory

4) To generate Table 4, execute:
    run_computational_sims;
This takes approximately 2.5 hours to complete on an Intel Core i7-8850h
This script outputs results in the file ./computational_simulation_results.txt

5) To generate Tables 5 and 6, execute:
    run_portfolio_optimization;
This takes less than 1 minute to run
Results are printed to the screen


NOTES ON FINANCIAL DATA
The stock return data in linear_returns.csv are the gross (i.e. linear) daily returns of the component stocks of the DJIA.
These were calculated using the adjusted close price, return(i) = adjClosePrice(i)/adjClosePrice(i-1).
Original price data was downloaded from Yahoo Finance ( https://finance.yahoo.com )
Only the 25 component stocks that were in the index for the entirety of 2016-2020 were included in this experiment.