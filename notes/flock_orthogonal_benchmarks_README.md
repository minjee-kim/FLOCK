# FLOCK orthogonal benchmark (prototype)

Standalone base-R study comparing **FLOCK-Novelty**, **FLOCK-Orthogonal**, and **Adaptive Metropolis** on three Bayesian targets. This is a controlled **total log-posterior evaluation budget** comparison, not a published-method reproduction.

## Run

Put `flock_orthogonal_benchmarks.R` in your R working directory and run:

```r
# Syntax/runtime debugging only. DO NOT interpret a smoke-test ranking.
Sys.setenv(RUN_MODE="smoke")
source("flock_orthogonal_benchmarks.R")

# Scientific pilot
Sys.setenv(RUN_MODE="pilot", REPS="10", BUDGET="30000")
source("flock_orthogonal_benchmarks.R")

# Larger experiment, preferably after inspecting pilot behavior
Sys.setenv(RUN_MODE="full", REPS="100", BUDGET="100000")
source("flock_orthogonal_benchmarks.R")
```

When changing `REPS` or `BUDGET`, use a new `OUTDIR` or the default folder name, which includes the run mode and budget. An incompatible run manifest causes an error rather than silently mixing checkpoints. If resuming, preserve the same script version and manifest settings.

Default budgets: smoke 2 replicates × 3,000 evaluations; pilot 10 × 20,000; full 100 × 60,000. User environment overrides are authoritative. **The smoke run is for finding implementation errors, not measuring convergence.**

## Targets and reference values

1. **ProductRidge (2D).** Six observations from `N(theta1*theta2, 0.7^2)` with Gaussian priors. The nonlinear product creates curved, sign-symmetric posterior branches. Reference expectations use a one-dimensional numerical integral after analytically integrating out theta2.
2. **WeakSignalFunnel (7D).** Six weakly observed random effects with a shared log-scale. The reference posterior for the log-scale is computed via a one-dimensional integral of the exactly marginalized normal model; latent-effect moments use Gaussian conditional formulas. The `log(tau)>0` event is a *tail event*, not a mode.
3. **ThreeComponentMixture (5D).** Nine observations fit with a 3-component normal mixture, unit observation variance, normal priors on means, and a Dirichlet(1,1,1) prior on weights. The unconstrained parameterization uses log(w1/w3), log(w2/w3), **including the correct Jacobian in the prior**. Reference functionals use exact enumeration of `3^9 = 19,683` allocations and normal/Dirichlet conjugacy. The six label permutations have equal posterior probability by exchangeability; labeled posterior mean comparisons alone are misleading.

References are used ONLY for retrospective scoring—not initialization, proposal selection, or adaptation.

## Algorithms and burn-in

**AdaptiveMetropolis:** Eight independent random-walk chains. Each receives approximately B/8 target calls, adapts empirical covariance during the first 35% of its run, freezes thereafter, and retains only the post-warmup samples. Their pooled estimate is reported for illustration, but if they remain trapped in different modes, pooling cannot repair nonconvergence. The benchmark records mode switching separately in each chain.

**FLOCK_Novelty:** Two exploratory scouts with large local proposals and occasional prior-independent MH jumps, then three follower chains starting from scout discoveries and three exploratory pioneers initialized at prior-drawn candidates ranked using novelty and density. Only the 40% final-stage independence/global mixture MH chain (minus 10% of its final transitions as burn-in) contributes posterior samples.

**FLOCK_Orthogonal:** One covariance-adaptive scout, then three *followers* started from plausible, less-visited prior candidate regions, plus three *pioneers* started around high-density scout discoveries and displaced along low-variance scout eigenvectors. Pioneer proposals use a regularized/capped inverse-eigenvalue geometry; followers use the learned covariance direction. The final sampler and posterior-draw handling are the same as Novelty.

FLOCK exploratory chains need not themselves approximate the target: their only job is to learn the proposal. Both FLOCK final samplers freeze their learned global proposal. Local RW MH and global independence MH transitions use **separate correct MH ratios** and the same mixing probabilities. The global proposal includes a prior-based component with full support.

## Budget and diagnostics

All methods use exactly B log-posterior calls per **replicate**; the initialization and pioneer-start candidate scoring are counted. FLOCK allocates roughly 20% scouting, 40% middle, and 40% final. AM uses all its budget across eight local chains. All draws used for FLOCK inference come from the final chain *only*.

The output folder contains:

- `reference_values.csv`: numerical/analytic reference expectations.
- `replicate_results.csv`: functional accuracy, occupancy, transitions, accept rate, cost.
- `metric_estimates.csv`: **every metric estimate and reference** for all methods and replicates.
- `group_occupancy.csv`: positive/negative branches, scale tail event, or six mixture label permutations.
- `method_summary.csv`: summary over replicates (not one-run rankings).
- `paired_comparisons.csv`: paired mean standardized RMSE differences and descriptive bootstrap intervals.
- `posterior_diagnostics.pdf`: reference marginals versus chain KDEs, example traces, and error distributions.
- `checkpoints/`: atomic RDS for each completed scenario–replicate–method.

**Important:** the `standardized_rmse` is a descriptive aggregation of preselected moment/tail-probability errors using fixed problem-specific scales. Inspect individual metrics and group frequencies before ranking. High acceptance or frequent transitions are not proof of posterior accuracy.

## Interpret with care

- The three scenarios test **different** difficulties. A method that wins all three would be interesting, but unlikely to establish broad dominance.
- The two FLOCK variants differ in more than the orthogonal start rule: one vs. two scouts, scout adaptation, and middle-stage proposal geometry. This is a whole-method comparison, not a clean causal ablation. If Orthogonal wins, isolate its eigenvector strategy with a matched-scout control.
- AM is a **local** covariance-adaptive baseline, not a state-of-the-art global adaptive mixture MH implementation. If FLOCK wins, follow with adaptive mixture MH, HMC/NUTS where applicable, and appropriately tuned parallel tempering.
- Weak-data hierarchical posteriors can often be handled more simply by noncentered reparameterization. FLOCK should be compared against that solution before claiming its exploratory geometry is necessary.
- Mixture label switching reflects model symmetry. Judge scientific inference using invariant quantities and predictive distributions, not only whether all six labelings were visited.
- Budget matching by log-target evaluations does not match wall time; proposal construction and density evaluation add overhead and wall time is also reported.
- **The script could not be executed in R in this environment.** It includes numerical self-tests for Gaussian covariance and target/reference finiteness. Run smoke before investing in larger experiments and report any errors.

`reference_marginal_ci.csv` provides numerically integrated 95% credible-interval
bounds for the plotted coordinate. `replicate_results.csv` includes the
corresponding sampler-based posterior mean and 95% credible interval for each
method/replicate. The order-label-dependent marginal for the mixture is
illustrative; focus on label-invariant functional metrics for inference.