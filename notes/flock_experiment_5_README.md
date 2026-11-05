# FLOCK Experiment 5

## Run

This is a standalone base-R script. Do not source earlier FLOCK scripts first.
Use a fresh R session and place `flock_experiment_5.R` in your working directory.

```r
# Clear overrides left behind by earlier experiments.
Sys.unsetenv(c("REPS", "BUDGET", "OUTDIR", "METHODS", "FLOCK_AUTORUN"))
Sys.setenv(RUN_MODE = "smoke")
source("flock_experiment_5.R")
```

| Mode | Replicates | Target calls per method and replicate |
|---|---:|---:|
| `smoke` | 3 | 5,000 |
| `pilot` | 10 | 20,000 |
| `full` | 100 | 30,000 |

After the smoke test succeeds, run the pilot. Inspect the map and merge audits
before committing to the full study.

```r
Sys.setenv(RUN_MODE = "pilot")
source("flock_experiment_5.R")

# After reviewing the pilot:
Sys.setenv(RUN_MODE = "full")
source("flock_experiment_5.R")
```

`REPS` and `BUDGET`, when set, override the mode defaults. The script prints the
actual settings. Each mode/budget combination has its own default output folder.

For an individual-method debug run in a fresh folder:

```r
Sys.setenv(RUN_MODE = "smoke", REPS = "1", BUDGET = "5000",
           METHODS = "FLOCK_Map,FLOCK_MapMerge", OUTDIR = "flock_map_debug")
source("flock_experiment_5.R")
```

Clear those overrides before returning to the full eight-method comparison.

## Methods

1. **FLOCK_Novelty:** earlier density-penalized distance rule for pioneer starts.
2. **FLOCK_Unexplored:** random choice among the more distant half of a candidate
   pool relative to scout coverage. This is a relative coverage heuristic, not
   a guarantee that a region is truly unexplored.
3. **FLOCK_Map:** provisional scout clusters guide follower starts and covariance;
   pioneers use the Unexplored rule; the cumulative map supplies the final proposal.
4. **FLOCK_MapMerge:** the same map design, with conservative merging before
   follower assignment and after middle-stage exploration.
5. **AdaptiveMixtureMH:** custom eight-chain, full-covariance adaptive exploratory
   baseline followed by the same frozen mixture proposal family and final MH kernel.
6. **ParallelTempering:** eight fixed temperatures, warmup scale adaptation,
   adjacent swaps, and retained cold-chain draws only.
7. **SAChES_inspired:** custom differential-evolution exploration followed by
   separately adaptive local chains, frozen for retained sampling.
8. **RWMH_4chains:** four independent random-walk MH chains.

The custom adaptive-mixture comparator is **not AIMM**. `SAChES_inspired` is **not a
reproduction of the published SAChES algorithm**. Performance against either
custom implementation is not evidence of superiority over the published method.

## Target and reference diagnostics

The target is a correlated 10-dimensional, six-component Gaussian mixture with
weights `(0.30, 0.25, 0.20, 0.15, 0.08, 0.02)`. A/B and D/E are overlapping pairs.
The sampler sees only a log-density function and generic dispersed starts.

Analytical reference means, covariance matrices and marginal quantiles are used
only after sampling. Component weights are estimated by average true-component
responsibilities; these estimates have the correct reference expectation even
when components overlap.

**Components, map clusters, and density modes are not interchangeable.** The
script reports visits to 95% Mahalanobis component cores and to prespecified
component groups AB, C, DE, F. It does not claim to know the exact attraction
basins of the mixture density. Core visits can overlap and do not establish
accurate posterior mass estimation.

## The map comparison

FLOCK_Map and FLOCK_MapMerge share scouts, stage-specific seeds, middle-chain
budgets, pioneer policy, covariance regularization, and final MH machinery.
Their deliberate difference is merging and its downstream consequences.

Both reserve 2% of the total budget **inside** the 45% middle/preparation stage.
Density-bridge checks spend that reserve. Any unused reserve becomes discarded
MH settling transitions using the frozen global proposal. Thus, merge checks are
not free, and both variants still have the same 35% final-transition budget.

The map is refitted to cumulative histories at stage boundaries. Middle chains
receive a frozen scout map and do not communicate while running. Execution is
serial in this script; this is not multicore benchmarking.

Merging requires sufficient unique exploratory observations, Mahalanobis center
distance at most 1.5, symmetric Gaussian KL at most 1, and (by default) a
nine-point straight-line density check with log-density drop at most 2 relative
to the lower endpoint. These are fixed pilot heuristics, not calibrated tests or
proofs of connectivity. `BRIDGE_CHECK=FALSE` disables the density check and is
recorded in the run manifest. Use a new output folder when changing it.

The audit records both rejected and accepted pairs and the cost of each check.
Retrospective flags for merging entries with different dominant reference groups
are **not** reported as merge error rates or precision/recall.

## Evaluation budgets and time

Every method must use exactly the requested number of target-density calls.
Initialization, candidate scoring and density-bridge checks are counted.
Diagnostics are oracle-assisted postprocessing, excluded from the algorithm's
sampling budget and elapsed time; diagnostic time is reported separately.

Equal target-call budgets do not imply equal step sizes, retained sample counts,
tuning strength, or wall-clock cost. Competitor-specific settings are fixed in
the script and saved in the output folder. Warmup and adaptation are discarded.
A correct fixed MH kernel does not guarantee adequate finite-budget mixing.

## Files and restart behavior

Default output directory: `flock_exp5_<mode>_B<budget>/`.

Key outputs:

- `flock_experiment_5_results.csv`: replicate-level errors, estimates and timing.
- `flock_experiment_5_summary.csv`: means, medians and between-run variability.
- `flock_experiment_5_paired.csv`: paired differences and exploratory bootstrap intervals.
- `flock_experiment_5_component_discovery.csv`: visits by method, stage and source.
- `flock_experiment_5_discovery_summary.csv`: stage-specific visit rates.
- `flock_experiment_5_first_evaluated.csv`: first evaluated point in each component core.
- `flock_experiment_5_budget_audit.csv`: exact target calls by stage.
- `flock_experiment_5_map_summary.csv`: map size and regional summaries.
- `flock_experiment_5_merge_audit.csv`: accepted and rejected merge decisions.
- `flock_experiment_5_merge_reference_diagnostics.csv`: retrospective group-overlap flags.
- `flock_experiment_5_follower_assignments.csv`: middle-chain assignments.
- `flock_experiment_5_diagnostics.pdf`: multipage plots.
- `settings.txt`, `session_info.txt`, `script_used.R`, `INTERPRETATION.txt`: provenance.
- `checkpoints/`: a checkpoint after each completed replicate-method pair.

Rerun the same settings and script with `RESUME=TRUE` (default) to skip completed
pairs. The manifest rejects incompatible script/settings/R-version changes.
Increasing `REPS` while keeping other settings fixed is allowed. Changing the
algorithm requires a new output directory; do not overwrite historical runs.

First-replicate retained draws are saved in checkpoints. Set `SAVE_DRAWS=TRUE`
before a run to retain every replicate's draws. Reusing existing checkpoints
will not retroactively create missing draws.

## Validation status and an earlier correction

R is not installed in the authoring environment. The script has been checked for
balanced delimiters/strings and reviewed, but **has not been executed end-to-end
in R here**. The first local smoke test remains essential.

Built-in startup checks test Gaussian covariance-factor conventions, mixture
sampling moments, a pointwise MH flux identity, basic merge geometry, singleton
candidate selection, and integer budget splitting. Every algorithm also has
runtime budget and finite-draw assertions. Passing these checks is not a proof
that the whole method produces accurate finite-run posterior inference.

An earlier Experiment 2/3 helper returned a generic covariance factor `L` with
`L %*% t(L) = S`, but generated proposals using `t(L) %*% z` while evaluating their
density with `solve(L, ...)`. In general those describe different covariances.
This script instead uses the upper Cholesky factor `R = chol(S)` consistently:
`mu + t(R) %*% z` for sampling and `backsolve(R, x-mu, transpose=TRUE)` for density
evaluation. Do not use this version to append to earlier result tables, and do
not treat earlier results affected by that mismatch as validated MH benchmarks.

## Reading the outcome

The main ablation is MapMerge versus Map. Also examine Map versus Unexplored,
posterior covariance/quantiles, rare-component mass, and wall-clock overhead.
A better mean alone does not establish reliable posterior recovery. If almost
no merge occurs, this configuration has not meaningfully tested merging.
Bootstrap intervals quantify variation across seeds for this one fixed target;
they are exploratory and not adjusted for multiple comparisons. Avoid repeatedly
retuning the same benchmark until FLOCK wins.

## Methodological sources

- Sandia, *SAChES: Scalable Adaptive Chain-Ensemble Sampling* (2017):
  https://www.sandia.gov/research/publications/details/saches-scalable-adaptive-chain-ensemble-sampling-2017-08-01/
- Maire et al., *Adaptive Incremental Mixture Markov Chain Monte Carlo*:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7224357/

These references motivate comparison choices; this script does not reproduce
either published method.