# FLOCK vs Dog Walker: theta-squared posterior illustration

## Research question

Can seven interacting, mutually repelling exploratory chains teach a frozen global proposal about both posterior modes more effectively than FLOCK-Novelty, adaptive mixture MH, or classical samplers under the same total number of target-density evaluations?

This is a **teaching and exploratory benchmark**, not a validated new sampler or an exact implementation of published SAChES/AIMM.

## Target

Simulated data `y_i ~ N(theta_true^2,1)` with `theta_true = 3` and sample size 12. Fitted likelihood `y_i | theta ~ N(theta^2,1)` and prior `theta ~ N(2, 3^2)`. This is correctly specified but sign-nonidentifiable. The posterior has two separated modes near -3 and +3 with different masses. An accurate reference posterior is calculated by one-dimensional grid integration.

## Run

Save the script as `flock_dogwalker_theta2.R`. In R:

```r
Sys.setenv(BUDGET="12000", OUTDIR="flock_theta2_outputs")
source("flock_dogwalker_theta2.R")
```

For a quick smoke test, set `BUDGET="2500"`; it will still run all methods and produce eight independent runs per method. Re-run with at least 12000 to assess the behavior more meaningfully. `SEED` can also be overridden.

## Methods

- RandomWalkMH: local MH proposal, 25% burn-in.
- AdaptiveMetropolis: learn proposal scale during first 25%, freeze afterwards.
- IndependenceMH: proposals from the prior, correct independent MH ratio, 25% burn-in.
- ParallelTempering: four temperatures, retain only cold chain and discard its first 25%.
- AdaptiveMixtureMH: eight exploratory trajectories (50% budget), global proposal, final MH (10% final burn-in).
- SAChES_inspired: differential-evolution exploration (50%), local final sampler (25% final burn-in). This is **not** a faithful reproduction of SAChES.
- FLOCK_Novelty: two scouts (18%), six middle trajectories (32%), final MH (50%, excluding first 10% of final chain).
- FLOCK_Unexplored: same as above with a random above-median-novelty pioneer selection.
- DogWalker90: seven interacting walkers governed by repulsive joint exploration target (90%), frozen global-proposal MH final (10%, excluding first 10% of final chain).
- DogWalker90_NoRepulsion: identical 90/10 budget but removes repulsion. This **tests whether repulsion itself helps**.
- DogWalker65: same repulsion as DogWalker90, but 65/35 exploration/final allocation. This **tests whether spending 90% exploring starves posterior inference**.

Every independent run within a method uses the same total target evaluation budget. Eight independent runs are performed for every method. Dog Walker's seven interacting trajectories are **within** each one of its eight independent runs.

## Output

The output folder contains `generated_data.csv`, `exact_posterior.csv`, `reference_summary.csv`, `chain_posterior_summaries.csv`, `method_posterior_summaries.csv`, and `stage_mode_coverage.csv` (including elapsed times).

- `eight_chains_vs_exact.pdf`: one page per method, with 8 chain-level posterior KDEs (blue) against the exact posterior (black).
- `eight_chain_traces.pdf`: 8 retained chain traces for each method.
- `pooled_vs_exact_CAUTION.pdf`: pooled densities, **not valid** unless all chains actually mix to the same posterior.
- `flock_exploration_vs_final.pdf`: exploratory trajectories (not posterior draws) versus the final FLOCK chain.
- `dog_walker_exploration_vs_final.pdf`: seven repulsive exploration trajectories versus final Dog Walker chain, plus its controls.
- Per-method `_runs.rds` and `all_runs.rds`: raw draws, exploration histories, and diagnostics.

## Interpreting results

Check `prob_positive` against the exact `P(theta>0|y)` and count negative-to-positive transitions, rather than judging accuracy solely by a smooth density. The true data-generating *parameter* is +3, but the Bayesian posterior is not concentrated at +3 because the sign is unidentified by this likelihood. The correct reference for posterior accuracy is the exact posterior curve, not the vertical parameter-truth line.

Dog Walker exploratory samples target a **repulsive joint exploration distribution**, not the desired posterior; discard them for inference. Its final stage uses a frozen mixture proposal with a valid MH acceptance ratio. Final chain length under 90/10 is short, and a theoretical invariant distribution does not guarantee finite-sample convergence.

## Research note

The idea of interacting repulsive MCMC chains has precedent: Martino et al., *Smelly Parallel MCMC Chains* (ICASSP 2015). Dog Walker should be evaluated as a distinct exploration-and-consolidation **architecture**, not claimed as the invention of repulsive MCMC.