# FLOCK: final decision benchmark

**Objective:** Is the FLOCK scout/follower/pioneer *formation* worth further work when all algorithms use the same underlying Adaptive Metropolis (AM) local sampler and the staged algorithms use the same final sampler?

This experiment is a **go/no-go pilot**, not a proof of convergence, superiority, or novelty. It is **not** designed to favor FLOCK. The earlier ARPM results with a mismatched Gaussian proposal factor must not be used as validated evidence.

## Files

- `flock_final_decision.R` — standalone, base-R-only implementation.
- This README — mathematical assumptions, fairness controls, output descriptions, and interpretation plan.

Run using `source('flock_final_decision.R')` in R after placing the file in your working directory. R >= 4.0 is recommended. No extra packages required.

## Five methods and why they're here

| Method | Exploration and final estimation |
|---|---|
| **IndependentAM** | Eight independently initialized chains; each adapts only during its first 35% of transitions, freezes the covariance, and contributes retained posterior draws. No jumps. |
| **IndependentAM_Shifts** | Eight AM chains each spend 45% of their individual evaluation budget on warmup, *share those histories* to learn fixed symmetric region-shift proposals, then continue to collect posterior draws using the same fixed local+jump kernel used in the staged methods. This is the **stronger information-sharing baseline**. |
| **StagedAM** | One AM scout -> six AM explorers with undirected choices from identical candidate banks -> one final chain. The final chain uses local AM (during warmup, then frozen) plus 12% fixed symmetric region-shift proposals. |
| **FLOCK_AM** | **Same AM kernels, same budget, same scout seed, same candidate locations, same final local+jump kernel as StagedAM.** Three middle runs start near high-density scout discoveries; three select novelty-favored plausible starts. Only the starting-point assignment policy changes. |
| **FLOCK_AM_LocalOnly** | Identical FLOCK formation and same posterior final-chain initialization/covariance rule, but disables the informed region-shift moves in the final chain. This shows whether the final sampler actually exploits the exploration. |

The StagedAM/FLOCK-AM final stage is *not* a custom independence Gaussian mixture proposal. It uses **mostly Adaptive Metropolis local updates** and a small fixed proportion of **paired symmetric shifts between discovered centers**. A shift draws `x' = x ± delta + small Gaussian jitter`. For every `delta`, the reverse `-delta` has the same probability. Consequently, the jump proposal is symmetric, and the standard log-MH acceptance ratio `logpi(x') - logpi(x)` applies. The shifts are fixed before the retained final stage. The local covariance adapts **during warmup only** and then freezes.

**This is not a claim of a fundamentally new MCMC transition.** It deliberately isolates whether *strategic allocation of exploratory computation* adds value when using established kernels.

## Three target distributions (identical across all methods)

1. **KnownMixture5D:** A 5-dimensional mixture of four Gaussian components with unequal weights `(0.45, 0.30, 0.20, 0.05)` and different centers and variances. The target is exactly known. The key metric is the estimated rare-component posterior mass (5%); component responsibilities are used **only in evaluation**.
2. **DiagnosticLCA:** A 5-parameter, two-test/one-population latent-class model with no gold standard, Beta(2,2) priors on prevalence, sensitivities and specificities, and counts for `(11,10,01,00)` of `(7,2,2,19)`. Exact posterior reference expectations are computed by enumerating all possible latent disease counts and integrating the Beta priors. The key metric is `P(prevalence > 0.65)`; the disease-label orientation is symmetric and scientifically not identifiable from these observations.
3. **HighDimCorrelated:** A 40-parameter linear regression with highly correlated first two predictors and spike/slab *mixture* priors for those two coefficients. The exact posterior is a mixture of four conjugate multivariate normals, so reference moments and contrast probabilities are analytic. The key metric is `P(beta1 > beta2)`; the contrast sign is **not itself proof of two density modes**.

**Important:** The Gaussian-mixture target is a *known posterior target*, while the latter two are posteriors from actual likelihoods and priors. Treat them differently in interpreting findings. The true parameters and component memberships are not provided to any sampler.

## Runs

Smoke test (tests execution only, NOT scientific evidence):

```r
Sys.unsetenv(c('REPS','BUDGET','METHODS','SCENARIOS','OUTDIR'))
Sys.setenv(RUN_MODE='smoke')
source('flock_final_decision.R')
```

Pilot (recommended first substantive test):

```r
Sys.setenv(RUN_MODE='pilot',REPS='10',BUDGET='30000',SCENARIOS='all')
source('flock_final_decision.R')
```

Predeclared full experiment (requires substantial time):

```r
Sys.setenv(RUN_MODE='full',REPS='60',BUDGET='60000',SCENARIOS='all')
source('flock_final_decision.R')
```

The full experiment uses **54 million** target-density evaluations (60 x 60,000 x 3 scenarios x 5 methods), plus proposal-fitting overhead. Do not assume this will be quick. If computational limits are restrictive, a useful first targeted full comparison is:

```r
Sys.setenv(RUN_MODE='full',REPS='40',BUDGET='60000',
           SCENARIOS='KnownMixture5D,DiagnosticLCA',
           METHODS='IndependentAM,IndependentAM_Shifts,StagedAM,FLOCK_AM',
           OUTDIR='flock_final_targeted')
source('flock_final_decision.R')
```

Set `OUTDIR` to a **new empty directory** when changing settings or script version. The manifest rejects incompatible reruns in the same directory.

## Budget and burn-in

- Each method uses **exactly B evaluations of the identical unnormalized log posterior**, including all initial state and candidate-screening evaluations. Reference computations for evaluating ground truth are separate, not charged to any method.
- Staged/FLOCK methods allocate **15% scout, 30% middle runs, 55% final**; all exploratory histories are discarded from the posterior estimator. The first **18% of final transitions** are discarded for final AM adaptation and warmup; only the **single remaining final chain** is used for posterior summaries.
- IndependentAM uses **eight chains** with 35% of their own transitions discarded; their retained draws are pooled **for estimation only**. Mode transitions count within each chain and never count artificial transitions at concatenation boundaries.
- IndependentAM_Shifts uses **45% per chain for exploration/warmup**, followed by 55% retained draws with covariance frozen, and shared symmetric-shift proposals.
- Staged AM and FLOCK AM receive the **same scout seed** and the **same candidate sets** at every middle exploration slot; only the decision about which candidate to initialize from is different. Their learned proposals may differ, which is the intended consequence of the formation.
- All algorithms may have different numbers of retained draws after accounting for exploration; this is part of the fixed-budget comparison. Report MSE per total evaluation, plus wall time.

## Outputs

The output directory defaults to `flock_final_<mode>_B<budget>` and contains:

- `reference_metrics.csv`: ground-truth posterior quantities for every scenario.
- `replicate_results.csv`: each replicate/method's normalized RMSE, consequential-event error, threshold failure, CI endpoint error, time, transitions, warmup updates, jump acceptances, etc.
- `metric_estimates.csv`: every posterior quantity estimated and its exact reference.
- `method_summary.csv`: averages, medians and failure rates.
- `paired_bootstrap.csv`: paired differences with bootstrap 95% intervals. **Negative** `mean_diff` favors `method_A` for an error metric.
- `diagnostics.pdf`: error distributions, wall times and within-chain region transitions.
- `checkpoints/*.rds`: restartable results saved after each successful *method-replicate* combination. For replicate 1 only, includes retained posterior draws; explorer histories are never posterior draws.
- `manifest.rds`: prevents accidentally mixing incompatible settings in one output directory.

Outputs update at the end of every completed replicate. If R stops, source the same script with identical settings to resume from checkpoints. The checkpoint folder alone also preserves completed methods if a replicate is interrupted.

## Predeclare decision criteria **before** reviewing the full run

The **primary metric for each scenario** is the absolute error in:

- KnownMixture5D: 5% rare component mass;
- DiagnosticLCA: `P(prevalence > 0.65)`;
- HighDimCorrelated: `P(beta1 > beta2)`.

Secondary metrics: normalized posterior-functional RMSE, 95% interval endpoint error, significant failure rates (`important_failure`), *within-chain* region-switch counts, number of accepted informed jumps, and wall-clock seconds (including adaptation and proposal construction).

Pragmatic, **not formal statistical test**, decision rule:

**Continue developing FLOCK as a methods project only if**: (a) FLOCK-AM beats StagedAM on mean primary error in at least **2 of 3** scenarios, (b) it is also at least competitive with IndependentAM_Shifts, preferably outperforming it in one substantive case, (c) its improvement does not come with increased serious failure rates, and (d) sampling time is not prohibitively higher. Favor results with paired confidence intervals excluding zero and reproducible improvements across budgets.

**Deprioritize the general FLOCK claim** if FLOCK-AM repeatedly fails to beat StagedAM or IndependentAM_Shifts, or if final sampling has negligible informed-jump acceptance and cannot represent important regions despite exploration, particularly at the full budget. A 40D loss alone need not refute potential benefit on low-D multimodal or weakly identified posteriors; scope claims accordingly.

The retrospective `important_failure` threshold is `0.035` for 5D rare weight error, `0.10` for LCA prevalence-tail probability error, and `0.20` for high-D coefficient-contrast probability error. These are prechosen practical tolerances, not universally meaningful thresholds.

## Limitations

- The algorithm is a **prototype** and the benchmark is not a direct reproduction of published algorithms (AIMM, SAChES, etc.). Success here would justify a stronger literature-method comparison, not replace it.
- The final sampler returns a single chain, but this does not eliminate the scientific need for independent convergence checks. `final_covered_both`, acceptance diagnostics and exact-reference metrics are **test-only** checks, not general-purpose proof of full posterior exploration.
- Large exploratory distances and covariance eigenvectors cannot certify the existence of an undiscovered mode.
- The learned symmetric jumps can have low acceptance in high dimension; their real acceptance rates are explicitly recorded. Avoid retuning them based on the full benchmark results and presenting that as an untouched test.
- For known-normal-mixture overlaps, individual generating components are not necessarily distinct density modes. For LCA the orientation is label-symmetric; scientific conclusions may require identification constraints. For high-D regression the contrast sign is a region, not necessarily a mode.
- This R file was structurally inspected but **not run end-to-end in this environment**, where R is unavailable. The built-in numerical covariance/proposal symmetry checks help catch some errors but do not replace running and reviewing the smoke test.