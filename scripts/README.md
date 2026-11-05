# Original research R prototypes

This directory contains **14 historical standalone R scripts** recovered from the October 2026 research archive. They include the FLOCK scout/follower/pioneer variants, adaptive Metropolis controls, orthogonal and SoftOrth pioneers, Dog Walker repulsion, regional maps, and importance-sampling consolidation.

These are *research snapshots*, not an installed or audited R package. Several scripts define overlapping global functions; run each **in a fresh R session**.

## Selected historical commands

From the repository root:

```r
setwd("scripts")
Sys.unsetenv(c("REPS", "BUDGET", "OUTDIR", "METHODS", "SCENARIOS"))
Sys.setenv(REPS="5", BUDGET="20000")
source("flock_experiment_3.R")
```

A small smoke run is appropriate for checking execution, but **does not establish sampler accuracy**. See the corresponding [per-script notes](../notes/) for available run modes and settings. Use an output directory that does not overwrite existing results.

The files have been restored from the archive, **not rerun or comprehensively reviewed** as part of this GitHub reorganization. See [reproducibility notes](../docs/reproducibility.md) and [numerical provenance](../RESULTS_PROVENANCE.md).
