# Bayesian-MIDAS

Matlab code for the BMIDAS models proposed in Kohns & Potagailo (2023) "Flexible Bayesian MIDAS: time‑variation, group‑shrinkage and sparsity", a working paper version can be found [here](https://www.bankofengland.co.uk/-/media/boe/files/working-paper/2023/flexible-bayesian-midas-time-variation-group-shrinkage-and-sparsity.pdf). The code runs the nowcasting exercise for UK GDP growth using the T-SVt-BMIDAS model with GIGG prior and ex-post sparsification. It also incorporates alternative versions of the model without time-varying components or without sparsification step.

If you use the code in your work please cite as: 

Kohns, D., & Potjagailo, G. (2023). Flexible Bayesian MIDAS: time‑variation, group‑shrinkage and sparsity. Bank of England Working Paper.

-----

This code is free to use for academic purposes only, provided that the paper is cited appropriately. This code comes without technical support of any kind. It is expected to reproduce the results reported in the paper. Under no circumstances will the authors be held responsible for any use (or misuse) of this code in any way.

-----

# Run the Code

In order to run the code, makes sure the following is accessible in your directory according to the file structure on this git:
* Estimation file: `Main_estimation_bmidas.m`
* Excel files: 1) `UK_data_bmidas.xlsx` (contains the UK macro data with the transformation codes), 2) `publication_calendar.xlsx` (contains pseudo publication calendar for nowcasting application)
* Functions folder which include sampling functions and scripts for data handling.

## Example

To run 'Main_estimation_bmidas.m", one can set various modelling and specification choices at the beginning of the code, which will then be executed in the automatic part of the code.
- define correct directory
- estimation sample period and out-of-sample evaluation period
- groups of indicators can be included or dropped
- the nowcast evaluation can be shut off and only the latest quarter and the latest nowcast period can be evaluation
- model components can be shut off: time-varyind trend, stochastic volatility, t-distributed errors
- the ex-post sparsification step can be shut off
- the MIDAS frequency mismatch and number of high-frequency lags can be defined, and Almon lag polynomial restrictions can be imposed
  


