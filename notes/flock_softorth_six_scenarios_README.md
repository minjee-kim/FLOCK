# FLOCK-Orthogonal-Soft — six-scenario R benchmark

A standalone, **base R** research prototype comparing:

- `AdaptiveMetropolis`: eight independent local Adaptive Metropolis chains, adaptation during the first 35% of each chain, then frozen proposals.
- `FLOCK_Novelty`: two exploratory scouts, three followers, three novelty pioneers; one final frozen global mixture-proposal MH chain.
- `FLOCK_OrthogonalSoft`: one adaptive scout, three followers beginning in low-coverage plausible locations, and three pioneers encouraged (not forced) into low-variance eigendirections of the scout's local covariance; one final corrected MH chain.

**Only the final post-burn-in chain contributes posterior draws for each FLOCK method.** Adaptive Metropolis retains post-warmup draws from all eight chains, but pooling them is reliable only if they have adequately explored a common posterior.

## Scenarios

| Name | Parameters | Main challenge | Reference |
|---|---:|---|---|
| `ProductRidge` | 2 | Product nonidentifiability, sign symmetry | 1D numerical integration with Gaussian conditionals |
| `WeakSignalFunnel` | 7 | Hierarchical funnel under weak data | 1D numerical integration with Gaussian conditionals |
| `ThreeComponentMixture` | 5 | Six label-switched configurations | Exact enumeration of latent allocations (3^9) |
| `CorrelatedPredictors` | 2 | Collinear predictors, spike/slab ambiguity | Exact 4-component Gaussian-mixture posterior |
| `DiagnosticLCA` | 5 | Two tests without gold standard; latent-disease label symmetry | Exact enumeration of latent counts, conjugate Beta integration |
| **`HighDimCorrelated`** | **40** | 40 regression coefficients, nearly indistinguishable leading predictors, sparse competing explanations | Exact 4-component, 40-dimensional Gaussian-mixture posterior |

The high-dimensional case uses a **correctly specified Gaussian regression likelihood with mixture priors**, not a manufactured target density. Other coefficients include modestly correlated nuisance predictors. Posterior expectations are compared against an analytically calculated reference rather than parameter-generating truth.

## Commands

Save the script in your working directory and run the high-dimensional scenario **by itself** first:

```r
Sys.setenv(
  SCENARIOS = "HighDimCorrelated",
  RUN_MODE = "smoke",
  REPS = "2",
  BUDGET = "4000",
  OUTDIR = "flock_highdim_smoke"
)
source("flock_softorth_six_scenarios.R")
```

Once it succeeds, use a meaningful pilot:

```r
Sys.setenv(
  SCENARIOS = "HighDimCorrelated",
  RUN_MODE = "pilot",
  REPS = "12",
  BUDGET = "30000",
  OUTDIR = "flock_highdim_pilot"
)
source("flock_softorth_six_scenarios.R")
```

For all six scenarios use `SCENARIOS = "all"` and a **fresh** `OUTDIR`. You may select two or more with `SCENARIOS = "ProductRidge,DiagnosticLCA,HighDimCorrelated"`.

Built-in defaults: smoke 2 replicates × 4,000 calls; pilot 12 × 30,000; full 60 × 60,000. These are **per method per scenario**. In high dimension, map construction and proposal-density evaluations can be expensive despite equal target-call budgets. Compare wall-clock time as well.

## Outputs

- `summary.csv`: errors, retained draw counts, wall-clock seconds, transition/coverage diagnostics by run.
- `method_summary.csv`: aggregate RMSE, occupancy error and time.
- `metrics.csv`: each decision-relevant posterior functional and its known reference.
- `parameters.csv`: posterior means, standard deviations and 95% intervals for all parameters.
- `groups.csv`: exploratory region-occupancy summaries (not proof of unique density modes).
- `chains.csv`: per-chain trace diagnostics; one final chain for FLOCK, eight for AM.
- `paired_comparisons.csv`: paired bootstrap comparisons if at least eight replicates are available.
- `reference_metrics.csv`, `reference_marginal_intervals.csv`: analytic/integrated truth for scoring.
- `<Scenario>_FLOCK_Novelty_FINAL_posterior.csv`: **only retained final chain**, replicate 1.
- `<Scenario>_FLOCK_OrthogonalSoft_FINAL_posterior.csv`: likewise.
- `<Scenario>_AdaptiveMetropolis_eight_chains.csv`: retained draws with chain ID, replicate 1.
- `<Scenario>_FLOCK_OrthogonalSoft_scout_geometry.csv`: covariance eigensystem and direction preference, replicate 1.
- `posterior_illustrations.pdf`: reference density, first-replicate estimated density and traces.
- `flock_stage_exploration.pdf`: scout, follower, pioneer and final trajectories projected to two dimensions.
- `checkpoints/`: automatically saved after **each scenario × replicate × method**. Rerunning with the same settings resumes from saved checkpoints. The manifest disallows incompatible reuse of an output directory.

## Scientific cautions

1. An accurate **single** FLOCK final chain does not establish convergence. Use independent runs and mode-specific accuracy diagnostics. If the earlier stages miss a consequential region, the final chain may still miss it.
2. The 40-dimensional benchmark challenges AM's covariance-learning efficiency, **but FLOCK itself uses full-covariance eigenanalysis and Gaussian mixture density evaluations**. Its complexity is not inherently lower. A low-rank variant would be a different algorithm and deserves a separate ablation.
3. No algorithms are given analytical posterior quantities during sampling. These are used only for evaluation.
4. The comparisons have identical target evaluation counts but not necessarily identical proposal-tuning overhead. They are research prototypes, not state-of-the-art or automatically well-tuned implementations.
5. The script has been statically checked but **not executed in R in this environment**. Start with smoke mode and inspect diagnostics before a long run.
6. The very highly correlated predictor pair is nearly observationally equivalent: different coefficients can produce almost identical in-sample fits but different extrapolation under covariate shifts. The contrast is a **consequential** quantity to estimate.