# FLOCK-AM + Dog Walker: two-model benchmark

This benchmark is a standalone **base R** extension of `flock_am_formation_comparison.R`, leaving the original unchanged. It tests the same **DiagnosticLCA** (5 parameters, two imperfect diagnostic tests, no gold standard, disease-label ambiguity) and **HighDimCorrelated** (40 regression coefficients, near-collinear predictors, spike/slab Gaussian-mixture prior) target posteriors and their reference posterior functionals.

## Algorithms

| Method | Exploration | Retained posterior |
|---|---|---|
| `IndependentAM` | Eight separately adapted AM chains | Pool separately burned-in draws; look at each chain too |
| `StagedAM` | One AM scout + six undirected AM middle chains | One final AM/local-global MH chain after final warmup |
| `FLOCK_AM` | One AM scout + three novelty starts + three soft-orthogonal starts | One final AM/local-global MH chain after final warmup |
| `FLOCK_RWMH` | Same formation with fixed random-walk local updates | One final local/global MH chain after burn-in |
| `DogWalker65` | Seven mutually repelling walkers; 65% exploration | One final AM/local-global MH chain after final warmup |
| `DogWalker65_NoRepulsion` | Same as `DogWalker65` but repulsion strength zero | Same final sampler |
| `DogWalker90` | Seven mutually repelling walkers; 90% exploration | One final AM/local-global MH chain after final warmup |
| `DogWalker90_NoRepulsion` | Same as `DogWalker90` but repulsion strength zero | Same final sampler |

**Dog Walker exploration differs from the earlier one-dimensional illustration**: each walker now learns a regularized local AM covariance while exploring; it uses occasional larger symmetric random-walk steps and 2% prior-independence proposals, with the correct proposal correction. This ensures the comparison does not disadvantage Dog Walker by using completely fixed local proposals, but results should **not** be treated as identical to those from the original `flock_dogwalker_theta2.R` implementation.

### Repulsion

For walker states `x_1,...,x_7`, the exploration-only joint target is proportional to the product of original posterior densities multiplied by a repulsion factor over walker pairs. Distances are standardized using fixed coordinate scales estimated from **prior draws**, with dimension adjustment. Conditional walker updates use the full repulsion log-ratio; prior-independence jumps also use the correct prior proposal ratio. Repelled walker histories are **not posterior draws**. The final chain uses the original posterior and a frozen global mixture proposal plus a covariance-tuned local MH kernel, with the appropriate correction.

## Fairness and limitations

All methods use exactly the same **total log-posterior evaluation budget** per method and replicate. Matrix operations, proposal density evaluations, and wall-clock time are recorded separately via total runtime. Total budget is *not* the same as number of retained posterior draws: the 90% exploration version will have far fewer final samples. `StagedAM` and `FLOCK_AM` share identical local kernels and final MH kernels but differ in exploratory starting points. Dog Walker adds a distinct repulsive exploratory target, so its comparison with FLOCK is a **whole-algorithm** comparison. The no-repulsion controls isolate the added repulsion at the same exploration budget.

**One final posterior chain cannot establish convergence by itself.** `stage_coverage.csv`, `chains.csv`, `summary.csv` and retained posterior outputs are provided to assess missed regions, transitions, final global acceptance, and effective sample sizes. For the 40D regression, the coefficient-sign groups are not guaranteed separate density modes, and the LCA example's label-switching symmetry should not be mistaken for a separate scientific explanation.

## Running

Save the `.R` file in your R working directory.

Small execution check (NOT a meaningful performance comparison):

```r
Sys.unsetenv(c("REPS","BUDGET","OUTDIR","METHODS","SCENARIOS"))
Sys.setenv(RUN_MODE="smoke", SCENARIOS="all")
source("flock_am_dogwalker_comparison.R")
```

Run a substantive comparison:

```r
Sys.setenv(
  RUN_MODE="pilot",
  SCENARIOS="all",
  REPS="12",
  BUDGET="30000",
  OUTDIR="flock_am_dogwalker_pilot_30k"
)
source("flock_am_dogwalker_comparison.R")
```

If you only want the primary four methods for speed:

```r
Sys.setenv(METHODS="IndependentAM,FLOCK_AM,DogWalker65,DogWalker65_NoRepulsion")
Sys.setenv(OUTDIR="flock_am_dogwalker_primary_30k")
source("flock_am_dogwalker_comparison.R")
```

Choose a **fresh `OUTDIR` whenever method selections or settings change**. A manifest prevents incompatible checkpoints from being mixed. The script creates one checkpoint per completed scenario–replicate–method combination and can resume with identical settings.

## Main outputs

- `method_summary.csv`: average normalized RMSE, group L1 error, seconds, transitions, and, when applicable, walker/global acceptance.
- `summary.csv`, `metrics.csv`, `groups.csv`, `chains.csv`, `parameters.csv`: per-run posterior inference and diagnostics.
- `stage_coverage.csv`: scout/middle or walker stage region coverage vs final retained posterior coverage.
- `paired_vs_independent_am.csv`: paired errors vs independent AM.
- `paired_flock_am_vs_staged_am.csv`: isolated formation ablation.
- `paired_dogwalker_repulsion.csv`: *repelled minus no-repulsion* errors (negative favors repulsion), separately for 65% and 90% exploratory allocations.
- `posterior_diagnostics.pdf`: reference densities, retained traces, parameter relationships, and functional errors.
- `*_FINAL_posterior.csv`: retained final-chain draws from first replicate (FLOCK and Dog Walker); `*_eight_retained_chains.csv` for Independent AM.
- `*_walkers_final_positions.csv`: final walker states for first replicate; dog walker exploration **not** used for posterior inference.
- `reference_metrics.csv` and `checkpoints/`.

### Important runtime caveat

The file was statically checked in this environment, but **R is not installed here**. Run the smoke test and check its outputs before committing to extensive replications. The earlier 40D benchmark showed practically zero accepted global proposals: the new design reports `global_accept` and `final_distinct` rather than assuming its final sampler has mixed.