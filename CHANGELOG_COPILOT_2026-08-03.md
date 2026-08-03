# Copilot Change Log (2026-08-03)

## Scope
This note documents all Copilot code edits from this round, not only the horseshoe fix.

## Files Updated
- Main_new.m
- functions/dmtest_modified.m
- functions/bmidas_horseshoe.m

## Documentation / Meta
- Created this log file: CHANGELOG_COPILOT_2026-08-03.md
- No additional source-code files were edited by Copilot beyond the three listed above.

## 1) Horseshoe plus PC slice-sampling fix

### File
- functions/bmidas_horseshoe.m

### Problem
- Runtime error path ended at uni_slice with Invalid slice sampling argument.
- In PC-SV updates, the sampler support for V is (0.0001, 0.5), while the incoming default level for V_h0 can be outside this range unless scaled before sampling.

### Change
- Added state-variance scaling in horseshoe initialization to match the existing gigg path:
  - V_g0 = V_g0/10000
  - V_tau0 = V_tau0/10000
  - if sv_obs_type is PC: V_h0 = V_h0/100, else V_h0 = V_h0/10
  - V_omegah = V_omegah*100

### Effect
- Horseshoe with PC-SV now starts from admissible values for the slice update and avoids immediate argument-validation failure.

## 2) Pre-defining and parfor-safe output storage

### File
- Main_new.m

### Change set
- Pre-allocated date storage for parfor-safe indexed writes:
  - dq_nfor changed from dynamic append to NaT(1, nfor)
  - Write changed from concatenation to indexed assignment dq_nfor(tperiod) = ...
- Pre-allocated additional output arrays before loop:
  - cyc_pred_all
  - tau_all and sv_trend_all (only if trend component is active)
  - sv_all (only if observation SV is active)
- Added guarded output writes for MAL prior:
  - Only writes cycle/trend/SV arrays when prior.midas is not MAL
  - Only attaches output.cyc_pred_all when prior.midas is not MAL

### Effect
- Removes dynamic growth patterns inside parfor writes and avoids shape/runtime issues from storing optional arrays for model variants that do not produce them.

## 3) DM statistic and evaluation-table updates

### Files
- functions/dmtest_modified.m
- Main_new.m

### Change in dmtest_modified
- Extended function signature to accept an optional loss_type:
  - squared (default): d = e1.^2 - e2.^2
  - raw: d = e1 - e2
- Kept HLN-corrected DM statistic and Newey-West HAC variance framework.

### Change in Main_new
- Added a post-estimation block that:
  - Computes WQS for the current model
  - Loads MAL benchmark results (if available)
  - Builds subsample ranges (Full, GFC, Tranquil, Pandemic, Pre)
  - Computes RMSFE, CRPS, WQS summary metrics
  - Runs DM tests for relative comparisons to benchmark
  - Formats ratios with significance stars
  - Exports a summary table to Excel in the model output folder

### Effect
- Provides an automated, reproducible evaluation report versus the benchmark model with significance annotation.

## 4) Additional robustness update in Main_new

### File
- Main_new.m

### Change
- Made evaluation start-date handling robust:
  - Parses beg_eval_per as datetime
  - Falls back to first valid date at or after requested start
  - Falls back to final available date if needed
- Updated Tq assignment to explicitly use size(y,1).

### Effect
- Avoids failures when the requested evaluation start label does not exactly match d_q formatting/content.

## Cleans / Non-Changes
- No destructive git operations were run.
- No bulk formatting/refactor sweep was applied.
- No output artifacts were deleted by Copilot in this edit pass.
- Copilot did not execute MATLAB runs; output folder changes in the working tree are runtime artifacts, not code edits applied in this patch pass.
- Existing static-analysis warnings unrelated to these changes were left unchanged.

## Verification Performed
- Reviewed stack-trace files and compared related code paths.
- Checked diffs for Main_new, dmtest_modified, and bmidas_horseshoe.
- Ran per-file diagnostics after patching bmidas_horseshoe; no syntax errors introduced by the new lines.

## Notes
- Full MATLAB runtime execution was not run in this environment.
- Recommended follow-up: run the intended model configurations once and verify:
  - horseshoe plus PC path runs past prior failure point
  - parfor outputs are populated with expected shapes
  - exported DM Excel table matches expected benchmark comparisons
