# Bayesian-MIDAS

Matlab code for the BMIDAS models proposed in Kohns & Potagailo (2023) "Flexible Bayesian MIDAS: time‑variation, group‑shrinkage and sparsity", a working paper version can be found [here](https://www.bankofengland.co.uk/-/media/boe/files/working-paper/2023/flexible-bayesian-midas-time-variation-group-shrinkage-and-sparsity.pdf). If you use the code in your work please cite as: 

Kohns, D., & Potjagailo, G. (2023). Flexible Bayesian MIDAS: time‑variation, group‑shrinkage and sparsity. Bank of England Working Paper.

-----

This code is free to use for academic purposes only, provided that the paper is cited appropriately. This code comes without technical support of any kind. It is expected to reproduce the results reported in the paper. Under no circumstances will the authors be held responsible for any use (or misuse) of this code in any way.

-----

# Run the Code

In order to run the code, makes sure the following is accessible in your directory according to the file structure on this git:
* Estimation file: `Main_estimation_bmidas.m`
* Excel files: 1) `MF_DFM_FAMEdata.xlsx` (contains the macro data downloaded from fame with the transformation codes), 2) `publication_calendar.xlsx` (contains pseudo publication calendar for nowcasting application)
* Functions folder which include sampling functions and scripts for data handling.

## Example

