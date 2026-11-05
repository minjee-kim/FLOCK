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
