# FLOCK Formation × Adaptive Metropolis: controlled benchmark

This experiment compares whether organizing MCMC exploration into stages helps **even when the local MCMC transition algorithm is kept the same**. It also tests the separate hypothesis that FLOCK can use simple random-walk MH local updates without requiring Adaptive Metropolis.

## Downloaded files

- `flock_am_formation_comparison.R`: **standalone base R script**; no packages, external files, or other FLOCK scripts required.
- This README.

The script reproduces **the same target likelihoods and priors** as the earlier six-scenario `DiagnosticLCA` (5-dimensional, two-test imperfect diagnostics with disease-label ambiguity) and `HighDimCorrelated` (40-dimensional, high-collinearity Gaussian regression with a proper sparse Gaussian-mixture prior). Model data are held fixed across repeated sampler trials; data-generating truths are not substituted for Bayesian reference posterior quantities.

## Methods

| Method | Exploration | Posterior estimator |
|---|---|---|
| `IndependentAM` | Eight independently initialized AM chains | Pool post-burn-in draws **with diagnostics**; pooled estimates can be wrong if chains remain trapped. |
| `StagedAM` | One AM scout + six exploratory AM chains with *undirected* starting points | One corrected final sampler, using 20% final warmup for local covariance adaptation then frozen local/global proposals |
| `FLOCK_AM` | One AM scout + three novelty-oriented followers + three softly orthogonal pioneers | Exactly the same final sampler as `StagedAM` |
| `FLOCK_RWMH` | Same FLOCK allocation and starts, **fixed** random-walk local MH updates | One final sampler using fixed local random-walk MH **plus occasional corrected independence-MH proposals** |

All methods use the **same target density** and a fixed total number of **log-posterior evaluations**, including initialization, scouting, candidate evaluations, all middle-stage transitions, and final sampling. They do not necessarily consume equal wall-clock time or have equal retained sample counts. Every method has a 2% prior-independence proposal component, handled with its correct MH acceptance ratio.

For the three staged methods, evaluation allocation is **15% scout / 30% six explorers / 55% final sampling**, with final warmup = 20% of final transitions. `IndependentAM` uses eight chains, each dropping the first 35% as warmup. Exploratory chains are **not** counted as posterior draws for staged methods.

`StagedAM` and `FLOCK_AM` use identical basic kernels, proposal fitting, and final inference; the intended difference is how their six middle-stage starting locations are chosen. This is the principal formation ablation. `FLOCK_RWMH` is a separate kernel ablation and should not be described as a pure local RWMH implementation because it also uses global independence proposals in the final stage.

## Run

Save the `.R` file to your working directory, then in a fresh R session:

```r
Sys.unsetenv(c("REPS","BUDGET","OUTDIR","METHODS","SCENARIOS"))
Sys.setenv(RUN_MODE="smoke",SCENARIOS="DiagnosticLCA")
source("flock_am_formation_comparison.R")
```

This is **only** a code smoke test (2 repetitions, 4,000 target calls per method). Do not draw performance conclusions from it.

For a substantive exploratory pilot using both models:

```r
Sys.setenv(RUN_MODE="pilot", SCENARIOS="all",
           REPS="12", BUDGET="30000",
           OUTDIR="flock_am_pilot_30k")
source("flock_am_formation_comparison.R")
```

For a larger comparison:

```r
Sys.setenv(RUN_MODE="full", SCENARIOS="all",
           REPS="60", BUDGET="60000",
           OUTDIR="flock_am_full_60k")
source("flock_am_formation_comparison.R")
```

You can omit the fourth method with `Sys.setenv(METHODS="IndependentAM,StagedAM,FLOCK_AM")` **before** sourcing. Choose a fresh `OUTDIR` whenever you change model, method, seed, or budget settings. The script maintains a manifest to prevent incompatible checkpoints from being mixed.

## Outputs

- `method_summary.csv`: aggregated RMSE, region-probability errors, elapsed time, transitions.
- `summary.csv`, `metrics.csv`, `groups.csv`, `chains.csv`, `parameters.csv`: detailed per-replicate metrics, posterior parameter mean/SD/95% intervals, group occupancy, per-chain scalar ESS and transitions.
- `paired_flock_am_vs_staged_am.csv`: **key formation comparison** (negative mean RMSE difference favors FLOCK-AM); descriptive paired bootstrap intervals.
- `paired_vs_independent_am.csv`: each staged method versus conventional AM; descriptive paired bootstrap intervals.
- `stage_coverage.csv`: coverage during the scout, middle, and final phases; records whether exploratory stages found a group that the final chain later missed.
- `reference_metrics.csv`: analytic/exact Bayesian posterior reference quantities.
- `posterior_diagnostics.pdf`: reference density overlays, retained trace, parameter relationship, and accuracy bar plot.
- `DiagnosticLCA_FLOCK_AM_FINAL_posterior.csv`, etc.: first-replicate **final retained posterior** draws only.
- `DiagnosticLCA_IndependentAM_eight_retained_chains.csv`, etc.: first-replicate eight-chain post-warmup draws with chain IDs.
- `*_scout_eigenvalues.csv`: scout covariance eigenvalues and orthogonal direction probabilities.
- `checkpoints/`: summary checkpoint per completed method/replicate, supporting resume with identical settings.

## How to interpret results

**First selling point:** compare `FLOCK_AM` with `StagedAM`. Are consequential-region occupancy and posterior functionals more accurate because the middle chains start strategically?

**Second selling point:** compare `FLOCK_RWMH` with `IndependentAM`. Can the organized formation with simple fixed-covariance local MH be competitive with independently run Adaptive Metropolis? Remember both staged variants also have a corrected global independence-MH proposal.

Examine **acceptance rates, distinct final states, effective sample size, and region occupancy**, not just posterior mean error. A one-chain final sampler is easier to report but still cannot prove convergence or guarantee that other modes do not exist.

**Important limitations:**

- The LCA reference is a disease-label-symmetric, partially identified posterior. Equal label orientations do not automatically represent different scientifically meaningful conclusions.
- The high-dimensional benchmark has a known analytically expressible posterior; its group indicator is a coefficient contrast, not necessarily a distinct density mode.
- The final sampler's global Gaussian-mixture proposal is still a learned proposal. A poorly fitted high-dimensional mixture can have very low global acceptance. The script reports this instead of hiding it.
- The numerical experiments test methods, not general convergence guarantees or publication-ready superiority. The proposal is frozen, and adaptive local transitions freeze after final warmup.
- Every target call is counted, but proposal fitting, matrix algebra, and I/O have additional wall-clock costs.
- This file has been **statically checked**, but **not executed end to end in R** in the present environment. Run the smoke test before a long pilot.