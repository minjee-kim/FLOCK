# FLOCK

**Role-Specialized Exploration and Posterior Consolidation**

*Working research notebook · October 2026 · Exploratory methodological study*

**Manuscript:** [PDF](paper/main.pdf) · [LaTeX](paper/main.tex) · [Bibliography](paper/references.bib) · [Manuscript information](paper/README.md)

FLOCK investigates whether a fixed budget of posterior-density evaluations can be allocated more effectively through **scouts, followers, and pioneers** than through independent, interchangeable chains. A scout maps promising regions; followers refine them; pioneers investigate plausible directions or regions with inadequate coverage. Discoveries inform a final proposal used for **separately corrected posterior inference**.

The project is an exploratory investigation. The benchmark results are mixed, and **FLOCK has not been established as a new, generally superior or validated MCMC algorithm**.

## Method overview

![FLOCK algorithm architecture](figures/flock_structure.svg)

The central distinction is between *discovery* and *inference*. Exploratory chains can deliberately favor novelty, orthogonal directions, or repulsion without representing draws from the target posterior. In the final stage, a proposal frozen after exploration may be used with a valid Metropolis–Hastings correction, or the proposal may instead generate fresh importance-weighted draws. These are two different inference designs.

The investigated mechanisms include novelty-based pioneer placement, covariance/eigenvector-guided searches, repulsive **Dog Walker** exploration, regional proposal summaries, and formation ablations. Mathematical specifications, figures, and pseudocode appear in the [algorithm and variants documentation](docs/algorithm.md).

## Numerical investigations

Experiments included multimodal Gaussian mixtures, a nonidentified Gaussian product likelihood, hierarchical normal random effects, and a sign-ambiguous quadratic likelihood with two posterior modes. We tracked posterior-functional estimation error, mode coverage, and performance under fixed target-evaluation budgets.

Several *specific* comparisons were promising: in a one-dataset sign-ambiguous quadratic model, FLOCK-Novelty had mean absolute posterior-mean error **0.0113** versus **0.0500** for a custom adaptive-mixture MH comparator; in a hierarchical normal pilot, Orthogonal FLOCK had error **0.062** versus **0.246** for Adaptive Metropolis. **The latter comparison had only two runs per method**, and neither observation supports a general performance claim. Other ablations favored undirected exploration or showed inconclusive differences.

See [simulation models, equations, and results](docs/experiments.md) for the exact target distributions and comparison conditions, and [result provenance](RESULTS_PROVENANCE.md) for evidence grades and uncertainty.

## Repository contents

| Path | Contents |
|---|---|
| [`paper/`](paper/) | Working manuscript (PDF, LaTeX, and bibliography) |
| [`docs/algorithm.md`](docs/algorithm.md) | Scout/follower/pioneer designs, geometric and repulsive exploration, MH and importance consolidation |
| [`docs/experiments.md`](docs/experiments.md) | Simulated likelihoods, priors, comparisons, and findings |
| [`docs/reproducibility.md`](docs/reproducibility.md) | Current reproduction status, requirements, and limitations |
| [`data/`](data/) | Historical replicate tables, aggregate summaries, partial console logs, and explicitly unverified pilot data |
| [`figures/`](figures/) | Algorithm and exploration schematics |
| [`scripts/`](scripts/) | 14 original standalone R simulation and comparison prototypes |
| [`analysis/`](analysis/) | Python verification and figure regeneration scripts |
| [`notes/`](notes/) | Experiment-specific historical R run instructions |
| [`ARCHIVE_INVENTORY.md`](ARCHIVE_INVENTORY.md) | Original experiment-script and research-archive inventory |

## Reproducibility and limitations

The repository now includes **all 14 archived R prototypes, experiment-specific notes, recovered replicate CSVs, and Python checking scripts**. These are historical research snapshots: they have **not been independently audited or rerun end to end**. Some raw outputs remain incomplete, including the Experiment 5 console transcript and the independently unverified importance pilot. See [reproducibility notes](docs/reproducibility.md) and [result provenance](RESULTS_PROVENANCE.md).

In particular, a final corrected MH chain need not mix adequately across discovered modes; proposal discovery is not the same as estimating posterior mass. Importance-weighted consolidation may suffer from unstable weights. Custom baseline implementations should not be equated with fully tuned published algorithms.

**Research status:** working research notebook; not peer-reviewed, not a published MCMC package, and not a claim of established algorithmic novelty or convergence. No open-source license has been assigned.
