# Bayesian-MIDAS

MATLAB code for Bayesian MIDAS models from Kohns and Potjagailo (2025), Flexible Bayesian MIDAS: time-variation, group-shrinkage and sparsity.

Running the main script reproduces the RMSFE patterns, the main figures from the paper, and a summary table with RMSFE and CRPS.

Figures are generated for each selected model specification separately, rather than as combined model-comparison figures.

## Quick Start

1. Open Main.m in MATLAB.
2. Set your model options and MCMC settings.
3. Run Main.m.
4. Results and figures are saved under output/ in a folder named by the selected specification.

## Notes

- The code supports alternative prior choices and state specifications through options in Main.m.
- You can select variable groups to include in estimation (survey, activity/trade, labour, mortgages). This allows quick inclusion or exclusion of groups, for example dropping survey data.
- Selected groups are reflected in the output folder name via a data tag (for example, data_suractlabmort or data_actlabmort).
- You can also change which series are included inside each group. This within-group composition is currently not reflected in the output folder title.
- Absolute evaluation tables (RMSFE and CRPS) are always exported alongside model outputs. If a compatible benchmark output is available (currently the Mogliani and Simoni BMIDAS), additional benchmark-relative and DM-based comparison tables are exported.
- The provided data are a public subset of the original dataset and exclude survey data, which cannot be publicly shared. As a result, reproduced results will not exactly match the original paper.

## Citation

If you use these codes, please cite:

Kohns, D., & Potjagailo, G. (2025). Flexible Bayesian MIDAS: time-variation, group-shrinkage and sparsity. *Journal of Business & Economic Statistics*, *43*(4), 1034-1050.




