# FLOCK — role-specialized posterior exploration

**Archived exploratory research, October–November 2026. Not a validated MCMC algorithm.**

This repository documents an experimental idea: divide a fixed computation budget among **scouts**, **followers**, and **pioneers**. Learn from their discoveries to construct a proposal, then perform **separately corrected posterior inference**.

## Conceptual structure

![Three-stage FLOCK architecture](figures/flock_structure.svg)



**The final stage was the central unresolved challenge.** Exploratory states are not automatically posterior draws. A properly constructed MH transition can correct a frozen proposal in principle, but the quality of the learned proposal and resulting estimation accuracy still have to be established. Importance-weighted consolidation was examined separately and is not equivalent to a corrected final MH chain.

## Historical experiment highlights

| Study | Comparison | Interpretation |
|---|---|---|
| Experiment 2 (100 runs) | Targeted mean RMSE 2.192 vs random 2.032 | Directional exploration did not outperform random |
| Experiment 3 (100 runs) | Unexplored 1.901 vs adaptive mixture MH 2.059 | Paired difference interval included zero |
| Experiment 4 (100 runs) | Novelty functional RMSE 0.01477 vs adaptive mixture MH 0.01525 | Small difference; no general superiority |
| Experiment 5 (25 complete visible runs) | Novelty 0.391 vs AM 0.491 | Only 11 of 25 pairs favored Novelty; full results missing |
| Importance pilot (previously reported 10 runs) | Roles 0.204 vs undirected 0.325 | **Unverified** prior report; importance sampling, not MH |

**Warning:** different experiments use different targets and RMSE definitions; do not compare magnitudes across rows. Evidence grades and limitations are in [RESULTS_PROVENANCE.md](RESULTS_PROVENANCE.md).

## Promising results in specific scenarios

Although FLOCK did not consistently beat simpler samplers, **several controlled scenarios showed potentially useful gains**. These are exploratory observations, not evidence of general superiority.

| Scenario | FLOCK configuration | Comparator | Observed result | Evidence |
|---|---|---|---|---|
| **Weak-signal funnel** | Orthogonal exploration | Adaptive Metropolis | Mean error **0.062** vs **0.246** (about 75% lower) | Archived method summary, **2 replicates only** |
| **Product-ridge target** | Novelty exploration | Adaptive Metropolis | Mean error **0.161** vs **0.283** (about 43% lower) | Archived method summary, **2 replicates only** |
| **Bimodal posterior, estimating a target functional** | FLOCK-Novelty | Custom AdaptiveMixtureMH | Mean absolute error **0.0113** vs **0.0500** (about 77% lower) | One archived benchmark configuration, **96,000 target evaluations per method** |
| **Multimodal 5D target (Experiment 3)** | FLOCK-Unexplored | Custom AdaptiveMixtureMH | Mean RMSE **1.901** vs **2.059** (about 8% lower) | **100 replicates**; paired uncertainty interval crossed zero |
| **Importance-sampling consolidation pilot** | FLOCK-Roles | Matched undirected proposal | Reported mean RMSE **0.204** vs **0.325** (about 37% lower) | **10 previously reported runs, not independently verified** |

The funnel and product-ridge comparisons are especially fragile because **n=2**; they are useful leads, not reliable performance claims. The bimodal result comes from one target and does not establish replicated superiority. The importance-sampling pilot evaluates **self-normalized importance estimates**, not a converged MH sampler.

**Interpretation:** FLOCK may help on selected geometries or specific estimands, but observed improvements did not reliably translate across different targets, metrics, and configurations. Larger independent replicated tests and fair target-evaluation budgets would be required to establish a repeatable advantage.

The archived `data/method_summary.csv` and `data/method_posterior_summaries.csv` provide the underlying point summaries for the first three scenarios; see [RESULTS_PROVENANCE.md](RESULTS_PROVENANCE.md) for qualifications.

## Findings

- Direction-aware, role-specialized exploration did **not consistently improve** performance relative to simpler benchmarks.
- The scout/follower/pioneer distinction may improve discovery without producing correctly weighted posterior estimates.
- Importance-based consolidation avoids requiring one final chain to move among all modes, but introduces potentially severe weight degeneracy.
- Mode discovery, mode-weight recovery, estimator error, and computational cost must be assessed separately.

## Source material

See [ARCHIVE_INVENTORY.md](ARCHIVE_INVENTORY.md) for 14 historical R scripts and recovered outputs. Some archival materials have not yet been committed as individual files. The code is historical, not audited or rerun for this publication.

**Manuscript withheld for review:** `main.tex` and `main.pdf` are intentionally absent.

No license has been assigned.