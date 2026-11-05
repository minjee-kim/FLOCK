# Reproducibility and research status

This is an **exploratory research archive**, not a release of a validated MCMC implementation. The manuscript and original experiment summaries record experiments carried out in October 2026.

## What is currently in this repository

- [`paper/main.pdf`](../paper/main.pdf) and [`paper/main.tex`](../paper/main.tex): research notebook and LaTeX source
- [`data/`](../data/): selected summaries and explicitly marked unverified pilot results
- [`figures/`](../figures/): conceptual illustrations
- [`scripts/flock_map_demo.R`](../scripts/flock_map_demo.R): a historical region-map demonstration
- [`RESULTS_PROVENANCE.md`](../RESULTS_PROVENANCE.md): which results come from archived CSV, partial logs, or secondary reports
- [`ARCHIVE_INVENTORY.md`](../ARCHIVE_INVENTORY.md): inventory of the 14 original R prototypes and the research archive

**Incomplete:** Most original experiment-driving R scripts and full run outputs have not yet been committed individually to this GitHub repository. Consequently, the commands below **cannot** reproduce all published-in-repository results from a fresh clone. The archived scripts were not audited or rerun for this publication.

## Planned or recovered benchmark commands

The original archive records the following *historical* workflow, when its corresponding scripts are available:

```r
# Run from the repository root in a fresh R session
setwd("scripts")
Sys.unsetenv(c("REPS", "BUDGET", "OUTDIR", "METHODS", "SCENARIOS"))
Sys.setenv(REPS="5", BUDGET="20000")
source("flock_experiment_3.R")
```

`flock_experiment_3.R` is **not presently committed**; this is provided as a record of the intended invocation, not a currently working reproduction command. Do not run multiple archived benchmark scripts in one R session since they define overlapping global functions.

## Numerical interpretation

Metrics differ by target, functional, and computational budget. The tables in [Simulation targets and results](experiments.md) should not be pooled across rows. Posterior sample accuracy must be distinguished from exploratory discovery and importance-weighted estimation; a correct MH acceptance formula alone does not prove mixing or mode coverage.

## What a defensible release still needs

- Exact source code and environment for each benchmark
- Seeds, input data, reference-integral implementation, and exact target-evaluation budgets
- Independent replication of pilot results
- Comparisons with strong, fairly tuned baselines
- Diagnostic evidence for mode transitions, posterior functionals, and effective sample size

No convergence theorem or general superiority claim has been established.
