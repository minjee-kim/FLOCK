# FLOCK multimodal benchmark — how to run

This script (`flock_hard_multimodal_benchmark.R`) is a standalone R benchmark that preserves the **three targets** from the earlier FLOCK-Orthogonal experiment and provides harder versions of each. It also offers a fourth, optional posterior where different modes have different predictive consequences.

**Important:** This is a benchmark of **difficult but correctly specified Bayesian posteriors**, not a demonstration that MCMC repairs model misspecification. Nonidentifiability, label switching, and a weak hierarchical scale are different statistical phenomena.

## Install and run

Only base R is needed for the six always-available methods. To include published DREAMzs and DEzs implementations, install the maintained CRAN package once:

```r
install.packages("BayesianTools")
```

The script detects it automatically (when `INCLUDE_DREAM="auto"`). To force it to run those competitors, set `INCLUDE_DREAM="TRUE"`; to skip them, `INCLUDE_DREAM="FALSE"`.

Run the code-check version first:

```r
Sys.setenv(RUN_MODE="smoke", HARD="TRUE", INCLUDE_DREAM="FALSE")
Sys.unsetenv(c("REPS","BUDGET","OUTDIR","SCENARIOS"))
source("flock_hard_multimodal_benchmark.R")
```

Then a substantive pilot with **10 independent replicates and 20,000 evaluations per method**, if the smoke test succeeds:

```r
Sys.setenv(RUN_MODE="pilot", HARD="TRUE", INCLUDE_DREAM="auto")
source("flock_hard_multimodal_benchmark.R")
```

Or a larger study:

```r
Sys.setenv(RUN_MODE="full", HARD="TRUE", INCLUDE_DREAM="auto",
           REPS="40", BUDGET="50000")
source("flock_hard_multimodal_benchmark.R")
```

The default scenarios are ProductRidge, WeakSignalFunnel and ThreeComponentMixture. Include the predictive-decision example by running:

```r
Sys.setenv(SCENARIOS="all", OUTDIR="flock_hard_all_four")
source("flock_hard_multimodal_benchmark.R")
```

You can select scenarios by name, e.g. `SCENARIOS="ProductRidge,ConsequentialModes"`. To use the original easier parameter settings, set `HARD="FALSE"` and a new output directory. When overriding options, use `Sys.unsetenv()` to clear previous settings.

## Methods

- `AdaptiveMH`: eight independent Adaptive Metropolis chains; 35% warmup is discarded in each.
- `FLOCK_Novelty`: two aggressive scouts, three local followers, three novelty-seeking pioneers, one final MH sampler.
- `FLOCK_OrthogonalHard`: one adaptive scout, three underexplored-region followers, three pioneers initialized along low-variance eigendirections, then one final sampler.
- `FLOCK_OrthogonalSoft`: like Hard but with probabilistically weighted eigendirections and perturbed initialization; other sampling machinery is unchanged.
- `UndirectedMixtureMH`: an undirected staged exploration control using the same frozen global-mixture final kernel as FLOCK, without novelty- or orthogonal-directed starts.
- `ParallelTempering`: base-R likelihood-tempered parallel tempering with eight temperatures, adjacent swaps, and post-warmup **cold-chain only** estimates.
- `DREAMzs_Published` and `DEzs_Published` (optional): authentic algorithms via `BayesianTools::runMCMC()` rather than homemade approximations. Actual log-likelihood evaluations and wall-clock cost are measured; internal initialization may yield budgets different from the nominal value.

For each FLOCK variant, **only final-sampler post-burn-in draws** enter inference. Exploratory points can be intentionally biased toward novel directions; they are not pooled as posterior samples. All FLOCK variants have a common final sampling kernel: 80% globally learned independence proposals with a correct MH ratio and 20% local symmetric RWMH steps. Final proposals are frozen before sampling.

## Scenarios and targets

1. `ProductRidge`: posterior with sign-separated branches and a curved, weakly identified parameter ridge; reference obtained by integrating out one normal parameter.
2. `WeakSignalFunnel`: weakly observed hierarchical effects and a strongly changing conditional scale; reference moments integrated over the one-dimensional log-scale.
3. `ThreeComponentMixture`: Gaussian mixture with three unknown means and weights, symmetric label permutations, and an **exact** reference from summing all `3^9` latent allocations. Note that label permutations are observationally equivalent.
4. `ConsequentialModes` (optional): nonidentifiable sign of a mean under training data, asymmetric prior, and different future predictions for the two modes. A one-dimensional quadrature gives the reference; this is the strongest test of *inferential consequences*.

Set `HARD="TRUE"` for higher mode separation, a sharper product ridge, and a weaker-data funnel. This increases difficulty but does not guarantee that the methods rank differently.

## Outputs

Each run creates an output folder `flock_hard_<mode>_B<budget>_<hard/original>` with:

- `reference_values.csv`: reference values for each scientific or probability summary.
- `replicate_results.csv`: errors, posterior mode/event probabilities, transitions, retained draws, elapsed time, and actual target-density calls for each completed method/replicate.
- `method_summary.csv`: means and medians of accuracy metrics and failure frequencies.
- `paired_comparisons.csv`: descriptive paired error differences and bootstrap intervals (interpret only with adequate replicates).
- `posterior_reference_overlays.pdf`: retained density curves against the analytic/numerical reference marginal for each method on replicate 1.
- `<scenario>__<method>__retained_posterior.csv`: retained posterior parameter draws for replicate 1; useful for trace plots or additional diagnostics.
- `checkpoints/*.rds`: saves **every completed method-by-replicate**, so rerunning exactly the same configuration resumes completed work.

The manifest prevents accidental reuse of incompatible output folders. When using `BayesianTools`, an external sampler may consume a slightly different number of posterior evaluations because of initialization; `comparable_budget` flags results more than 5% away from the nominal budget. Do not make a strict fixed-budget claim from unmatched runs.

## What constitutes success?

The primary test is **posterior estimation error** on reference-known functionals and probabilities of substantively different regions, not merely mode discoveries. Compare paired differences versus the undirected mixture control and the **published** DREAMzs/DEzs methods when available. Plot errors against actual evaluations and wall-clock times. A single FLOCK final chain can still miss modes even if the exploratory chains found them.

The pilot is a first implementation, not a validation of FLOCK correctness or a reproduction of every state-of-the-art sampler. The code has been statically checked but **not executed in R in this environment**; do a smoke run before the expensive study.

## Published comparator documentation

- BayesianTools R package: https://florianhartig.github.io/BayesianTools/reference/runMCMC.html
- DREAMzs: https://florianhartig.github.io/BayesianTools/reference/DREAMzs.html
- The adaptive-incremental-mixture MCMC literature (AIMM) is relevant, but AIMM is **not** reproduced by our undirected-mixture control: https://pmc.ncbi.nlm.nih.gov/articles/PMC7224357/