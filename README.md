# FLOCK — role-specialized posterior exploration

**Archived exploratory research, October–November 2026. Not a validated MCMC algorithm.**

This repository documents an experimental idea: divide a fixed computation budget among **scouts**, **followers**, and **pioneers**. Learn from their discoveries to construct a proposal, then perform separately corrected posterior inference.

![Three-stage FLOCK architecture](figures/flock_structure.svg)



**The final stage was the central unresolved challenge.** Exploratory states are not automatically posterior draws. A properly constructed MH transition can correct a frozen proposal in principle, but the quality of the learned proposal and resulting estimation accuracy still have to be established. Importance-weighted consolidation was examined separately and is not equivalent to a corrected final MH chain.

## Promising results

The exploratory benchmarks identified several targets on which one FLOCK variant achieved lower estimation error than its specified comparator. The models are stated below so the simulations can be interpreted and reproduced; **these are not established general performance improvements**.

| Model | FLOCK configuration | Comparator | Reported mean error | Evidence |
|---|---|---|---|---|
| Hierarchical normal random effects (7D) | FLOCK-Orthogonal | Adaptive Metropolis | **0.062** vs 0.246 | 2 runs per method |
| Nonidentified Gaussian product model (2D) | FLOCK-Novelty | Adaptive Metropolis | **0.161** vs 0.283 | 2 runs per method |
| Sign-ambiguous quadratic Gaussian model (1D, bimodal posterior) | FLOCK-Novelty | Custom AdaptiveMixtureMH | **0.0113** vs 0.0500 mean absolute error in posterior mean | One archived benchmark, 96,000 log-target evaluations per method |
| Five-dimensional multimodal target (Experiment 3) | FLOCK-Unexplored | Custom AdaptiveMixtureMH | **1.901** vs 2.059 mean RMSE | 100 runs; paired interval includes zero |
| Importance-sampling consolidation pilot | FLOCK-Roles | Matched undirected proposal | **0.204** vs 0.325 reported mean RMSE | 10-run previously reported summary, **not independently verified** |

**Metrics differ:** the first two rows summarize standardized errors across selected posterior functionals, while the bimodal row is absolute error in the posterior mean. They are not directly comparable. The small 2-run pilots are descriptive, not robust demonstrations of superiority.

### Model 1: Hierarchical normal random effects (7D)

Six observed values were fixed to
$$
y=(-0.80,\;0.25,\;1.10,\;-0.55,\;0.65,\;-1.05).
$$
The fitted Bayesian model was
$$
t=\log \tau\sim N(-0.6,1^2), \qquad
\eta_i\mid t\stackrel{\mathrm{ind}}{\sim}N(0,e^{2t}),\qquad
y_i\mid\eta_i\stackrel{\mathrm{ind}}{\sim}N(\eta_i,2.3^2),
\quad i=1,\ldots,6.
$$
The unknowns are \(t,\eta_1,\ldots,\eta_6\). When \(t\) is small, the conditional random-effect distribution contracts sharply, creating a funnel-shaped posterior. The reference for \(t\) uses the analytically integrated model \(y_i\mid t\sim N(0,e^{2t}+2.3^2)\), followed by one-dimensional numerical integration. The benchmark scored the posterior mean of \(t\), \(P(t>0\mid y)\), \(E(\eta_1\mid y)\), and \(E(\eta_1^2\mid y)\), each with a fixed error scale. **The event \(t>0\) is a tail event, not a second posterior mode.**

### Model 2: Nonidentified Gaussian product model (2D)

With independent priors \(\theta_1,\theta_2\sim N(0,3^2)\), the likelihood was
$$
y_i\mid\theta_1,\theta_2\stackrel{\mathrm{ind}}{\sim}
N(\theta_1\theta_2,0.7^2),
\quad i=1,\ldots,6.
$$
The data were generated with R seed **101** by \(\texttt{rnorm(6, mean=3.7, sd=0.7)}\). Only the product \(\theta_1\theta_2\) is directly informed by the likelihood, yielding sign-symmetric, curved posterior branches under the prior. The accuracy score aggregates scaled errors for \(E(\theta_1^2\mid y)\), \(E(\theta_2^2\mid y)\), \(E(\theta_1\theta_2\mid y)\), and \(P(\theta_1>0\mid y)\). Reference expectations were obtained by analytically integrating out \(\theta_2\) conditional on \(\theta_1\), then numerically integrating over \(\theta_1\).

### Model 3: Sign-ambiguous quadratic Gaussian model (bimodal posterior, 1D)

This is the **exact Bayesian data-generating and fitted model** for the bimodal row, not a generic two-Gaussian mixture:
$$
\theta_{\mathrm{true}}=3,\qquad
y_i\mid\theta\stackrel{\mathrm{ind}}{\sim}N(\theta^2,1),
\quad i=1,\ldots,12,\qquad
\theta\sim N(2,3^2).
$$
The experiment simulated \(y_i\stackrel{\mathrm{ind}}{\sim}N(9,1)\) using R seed **20261025**. Because \(\theta^2=(-\theta)^2\), the likelihood cannot identify the sign of \(\theta\). The prior favors positive values, so the two posterior modes near \(-3\) and \(+3\) generally have **unequal mass**. Up to a normalizing constant, the posterior is
$$
p(\theta\mid y)\propto
\exp\left[-\frac{12}{2}(\theta^2-\bar y)^2\right]
\exp\left[-\frac{(\theta-2)^2}{2\cdot3^2}\right].
$$
An accurate numerical reference posterior was computed by trapezoidal integration over 12,001 equally spaced grid points on \([-9,9]\). The archived comparison used eight independent final runs per method, each with a **12,000-target-evaluation budget** (96,000 across runs). The table's 0.0113 versus 0.0500 compares **mean absolute error of the estimated posterior mean** versus that reference; it is **not** the RMSE of entire posterior distributions. The corresponding mean absolute errors in \(P(\theta>0\mid y)\) were **0.00188** (FLOCK-Novelty) and **0.00844** (custom AdaptiveMixtureMH). A single seeded simulated dataset and this benchmark do not establish a general win.

The simulator and frozen-final-proposal implementation are preserved in the historical source `scripts/flock_dogwalker_theta2.R` (see [archive inventory](ARCHIVE_INVENTORY.md)). 

### Important limitations

For Models 1 and 2, the archived comparison used **only two replicates per algorithm**, so the apparent improvements are highly uncertain. The Adaptive Metropolis comparator was a local covariance-adaptive implementation, not a fully tuned global mixture sampler. The exploratory importance-sampling pilot used a different final estimator, and its reported summary was not independently verified. See [RESULTS_PROVENANCE.md](RESULTS_PROVENANCE.md).

## Findings

- Direction-aware, role-specialized exploration did **not consistently improve** performance relative to simpler benchmarks.
- The scout/follower/pioneer distinction may improve discovery without producing correctly weighted posterior estimates.
- Importance-based consolidation avoids requiring one final chain to move among all modes, but introduces potentially severe weight degeneracy.
- Mode discovery, mode-weight recovery, estimator error, and computational cost must be assessed separately.

## Source material

See [ARCHIVE_INVENTORY.md](ARCHIVE_INVENTORY.md) for 14 historical R scripts and recovered outputs. Some archival materials have not yet been committed as individual files. The code is historical, not audited or rerun for this publication.

**Manuscript withheld for review:** `main.tex` and `main.pdf` are intentionally absent.

No license has been assigned.
