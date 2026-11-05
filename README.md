# FLOCK

**Role-Specialized Exploration and Posterior Consolidation**

*Working research notebook · October 2026 · Exploratory methodological study*

**Manuscript:** [PDF](paper/main.pdf) · [LaTeX](paper/main.tex) · [Bibliography](paper/references.bib) · [Manuscript information](paper/README.md)

FLOCK investigates whether a fixed budget of posterior-density evaluations can be allocated more effectively through **scouts, followers, and pioneers** than through independent, interchangeable chains. A scout maps promising regions; followers refine them; pioneers investigate plausible directions or regions with inadequate coverage. Discoveries inform a final proposal used for **separately corrected posterior inference**.

The project is an exploratory investigation. The benchmark results are mixed, and **FLOCK has not been established as a new, generally superior or validated MCMC algorithm**.

## Method overview

![FLOCK algorithm architecture](figures/flock_structure.svg)

The central distinction is between *discovery* and *inference*. Exploratory chains can deliberately favor novelty, orthogonal directions, or repulsion without representing draws from the target posterior. In the final stage, a proposal frozen after exploration may be used with a valid Metropolis–Hastings correction, or the proposal may instead generate fresh importance-weighted draws. These are two different inference designs.

## Methodological variants

The following approaches were explored as **proposal-learning or inference designs**, not as independently validated algorithms. Their potential uses describe the problems each was intended to address, rather than proven performance advantages.

### 1. Role-specialized FLOCK

**Introduction.** Scouts conduct broad initial exploration, followers refine plausible posterior regions, and pioneers investigate locations inadequately covered by earlier trajectories. The resulting histories inform a frozen global proposal.

**Potential use.** Exploring a multimodal posterior where independent local chains may repeatedly spend computation in the same basin. Its practical value depends on whether the discovered regions are ultimately assigned accurate posterior mass.

<p align="center">
  <img src="figures/flock_structure.svg" width="520" alt="Three-stage FLOCK architecture: scouts, main-body exploration, and final inference" />
</p>

### 2. Novelty-guided exploration

**Introduction.** Candidate pioneer starts are scored using their distance or distinctiveness relative to previous exploratory states, sometimes combined with density plausibility. The novelty, unexplored, density, and diversity variants tested different selection rules.

**Potential use.** Searching for less-visited parts of a complicated posterior support, particularly when conventional random-walk proposals tend to revisit familiar neighborhoods.

<p align="center">
  <img src="figures/flock_search_mechanisms.svg" width="520" alt="Exploration mechanisms; the left panel illustrates novelty-guided search" />
</p>

### 3. Orthogonal and geometric exploration

**Introduction.** A scout estimates local posterior covariance and its eigendecomposition, $\widehat{\Sigma}=V\Lambda V^\top$. Pioneers then probe selected eigendirections, with regularized step sizes; SoftOrth uses softer directional encouragement.

**Potential use.** Investigating narrow or strongly correlated parameter spaces, where local sampling geometry suggests directions that ordinary covariance-adaptive proposals may rarely explore. A small eigenvalue does **not** imply another mode exists.

<p align="center">
  <img src="figures/geometric_search.svg" width="500" alt="Geometric pioneer search: covariance contours and exploratory directions" />
</p>

### 4. Regional maps and consolidation

**Introduction.** Exploratory histories are summarized as candidate local regions, and the Map / MapMerge variants assess geometric similarities between these summaries before forming a global proposal.

**Potential use.** Avoiding redundant proposal components when several exploratory chains discover overlapping basins. Merging close components must not erase distinct regions or distort inferred posterior probabilities.

<p align="center">
  <img src="figures/flock_method_branches.svg" width="520" alt="Method branches from exploration design to consolidated inference" />
</p>

### 5. Dog Walker: repulsive exploration

**Introduction.** Multiple exploratory walkers interact through an auxiliary repulsion term, encouraging spatial separation. The studies compared seven walkers, different exploration-budget allocations, and matched no-repulsion controls. The exploratory histories are **not posterior draws**.

**Potential use.** Encouraging complementary spatial coverage when multiple trajectories otherwise cluster in one region. Repulsion can also push walkers toward negligible posterior-density areas, so any benefit must be measured against that cost.

<p align="center">
  <img src="figures/dog_walker.svg" width="500" alt="Dog Walker analogy: coordinated repulsive exploration with multiple walkers" />
</p>

### 6. Final inference: corrected MH versus importance weighting

**Introduction.** The original design freezes the proposal learned during exploration and applies a Metropolis–Hastings correction. The later importance-consolidation variant instead generates fresh independent proposal draws and assigns weights proportional to $\pi^*(x)/q(x)$.

**Potential use.** Global-proposal MH aims to obtain posterior draws while preserving the correct stationary law; importance consolidation may avoid requiring one terminal chain to transition among all relevant modes. Neither guarantees adequate mode coverage, and importance weights can degenerate.

<p align="center">
  <img src="figures/flock_method_branches.svg" width="520" alt="Two final inference designs: corrected Metropolis-Hastings and importance-weighted consolidation" />
</p>

**Further details:** [Mathematical formulations, diagrams, pseudocode, and method ablations](docs/algorithm.md). The specific implementations and evaluation results are documented separately in [simulation investigations](docs/experiments.md).

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
