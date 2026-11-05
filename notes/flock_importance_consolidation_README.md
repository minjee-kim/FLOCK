# FLOCK: aggressive exploration + importance consolidation

This is a **research prototype**, not a validated sampler or a claim of a new method. The goal is to isolate whether the **scout / follower / pioneer roles** make a better final importance-sampling proposal than matched undirected exploration.

## Mathematical objective

Given a target density proportional to `pi_star(theta)`, the exploration stage constructs a *frozen*, evaluable Gaussian-mixture proposal `q(theta)` that includes a broad defensive component. The final stage generates fresh independent `theta_i ~ q` and computes normalized weights proportional to `pi_star(theta_i)/q(theta_i)`. The output is **one weighted posterior approximation**, **not** a conventional MCMC chain. The exploratory trajectories are **not** pooled into the posterior estimator, and they do **not** need to converge to the posterior.

The target is a fully known mixture of four multivariate normal densities (the posterior analogue of a difficult multimodal target). Available dimensions: `2d`, `5d`, `10d`. Component weights are 0.44, 0.30, 0.21 and 0.05; some modes have correlated coordinates and different widths. Ground-truth means, covariances, marginal intervals, and component masses are computed from the target **only for post-hoc evaluation**. The algorithms use only the counted `logp` oracle and shared broad initial-reference sampler.

## Methods

* **FLOCK_Roles** — 1 scout (20% budget, 40% multiscale MH / 40% tempered MH / 20% random restarts), 3 followers starting near distinct high-log-density scout discoveries, 3 pioneers chosen among candidate locations using novelty plus density plausibility. Exploratory budget overall is 65%; remaining 35% is fresh importance-sampling proposals. Pioneers can learn coverage information from earlier followers and pioneers.
* **Undirected_Matched** — EXACT SAME scout schedule, three local exploration chains, three pioneer exploration chains, **and identical screening-evaluation budget**, but follower starting points are independent broad draws and pioneer candidate selection is random. This is the *primary ablation* of role-directed information transfer.
* **Undirected_RWMH** — Seven independently initialized conventional random-walk MH exploratory chains sharing the same 65% budget; 35% importance consolidation. This is a simple secondary baseline, not as tightly controlled as Undirected_Matched.

All three use the **identical** final proposal-construction algorithm (k-means local Gaussians + 15% broad Gaussian) and importance-sampling weighting. Selection and clustering do not use known mixture component locations. Exploratory `logp` calls, candidate screening and final samples all count toward the same target-evaluation budget.

### Important limitations

* The scout's changing temperature and its uncorrected random restarts are **not** a stationary MCMC algorithm; that is intentional because it is only proposal-learning exploration.
* Both baselines also receive the same final importance-sampling machinery, so results measure the *incremental utility* of FLOCK's roles, not just whether importance sampling beats MCMC.
* Importance weighting is valid only when `q` has support wherever the target does; the broad defensive Gaussian ensures support here, but **not necessarily practical finite-sample accuracy**. In 10 dimensions, importance weights may become highly concentrated.
* Importance-sampling ESS does **not** certify discovery of unseen components. Check component-weight errors, missed-consequential-component diagnostics, and credible-interval errors as well.
* Component responsibilities are **retrospective scoring tools**; overlapping mixture components do not necessarily correspond to distinct density modes.
* No R interpreter was available to execute the script here. The script includes basic algebraic numerical checks; run `smoke` before any substantive run.
* Published adaptive multiple-importance-sampling, SMC, and regional methods are **not** implemented in this limited ablation. Positive results would motivate adding those competitors, not establish state-of-the-art performance.

## Run

Save the `.R` file somewhere accessible in R. To avoid inheriting settings from earlier FLOCK experiments:

```r
Sys.unsetenv(c('REPS','BUDGET','OUTDIR','TARGET'))
Sys.setenv(RUN_MODE='smoke', TARGET='2d')
source('flock_importance_consolidation.R')
```

The smoke test defaults to 2 replicates x 3,000 target calls per method and is **only a code check**. Next:

```r
Sys.setenv(RUN_MODE='pilot', TARGET='5d', REPS='10', BUDGET='12000')
source('flock_importance_consolidation.R')
```

Larger replication:

```r
Sys.setenv(RUN_MODE='full', TARGET='5d', REPS='50', BUDGET='30000')
source('flock_importance_consolidation.R')
```

Optional 10D stress test: set `TARGET='10d'`. This is not guaranteed to perform well with Gaussian importance sampling. Leave `OUTDIR` unset for automatic budget- and case-specific output folders. If you edit the algorithms, use a *new output directory* to avoid reusing checkpoints from an incompatible script.

## Output files

* `replicate_results.csv`: per-replicate posterior-mean RMSE, component-mass error, covariance error, 95% marginal interval endpoint error, rare-component mass error, IS ESS, max importance weight, target-call budget and wall clock time.
* `method_summary.csv`: mean/median accuracy, ESS ratio, failure rate and seconds.
* `paired_comparisons.csv`: paired FLOCK-minus-control differences plus descriptive bootstrap intervals; negative favors FLOCK for error metrics, positive favors FLOCK for ESS.
* `exploratory_component_coverage.csv`: how often each original exploratory chain entered a known component core (diagnostic only).
* `<METHOD>_weighted_posterior.csv`: **first replicate** fresh final draws and normalized importance weights for each method; no exploratory draws included.
* `flock_importance_diagnostics.pdf`: error/ESS/time boxplots and 2D weighted-posterior plots for the first replicate.
* `checkpoints/*.rds`: completed method / replicate results; restarting with unchanged settings skips completed work.

## Go / no-go interpretation

The key comparison is **FLOCK_Roles vs Undirected_Matched**. Evaluate not just mean RMSE but *probability of missing a component with consequential posterior mass*, relative error in posterior mean and covariance, marginal credible-interval endpoints, IS weight degeneracy, and real elapsed time. Use paired estimates and uncertainty across *many* replicates. If role-specialized exploration doesn't improve the **matched** control, the mechanism isn't yet justified even when it beats undirected local MH.