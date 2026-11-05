# Reproducibility

The repository now includes **the 14 archived standalone R prototypes**, matching per-script notes, extracted CSV/log outputs, and Python verification and plotting utilities. This is a reconstructed experimental archive, **not a formally validated MCMC software release**.

## Check the recorded summaries

```bash
python -m pip install -r requirements.txt
python analysis/verify_results.py
python analysis/replot_figures.py
```

These commands validate and visualize **saved numerical outputs**; they do not rerun posterior samplers. The original verification script checks Experiment 2–4 summary counts and the 25 complete Experiment 5 replicate blocks. It was run successfully against the recovered archive during repository organization.

## R simulations

The scripts are under [`scripts/`](../scripts/) and their experiment-specific guidance is under [`notes/`](../notes/). Use a new R session for each script because functions and global objects are reused across prototypes.

Example from the repository root:

```r
setwd("scripts")
Sys.unsetenv(c("REPS", "BUDGET", "OUTDIR", "METHODS", "SCENARIOS"))
Sys.setenv(REPS="5", BUDGET="20000")
source("flock_experiment_3.R")
```

This is a short exploratory rerun; use the exact matching seed, target, code revision, evaluation budget, and number of replicates before comparing against `data/`. A smoke check is not a performance result.

**Not independently completed:** a fresh end-to-end rerun and mathematical audit of all historical R experiments. Some benchmark outputs remain incomplete or available only as historical summaries. In particular, the ten-run importance pilot is marked **unverified** and Experiment 5 is based on a partial console transcript.

## Fairness and correctness checks

- Match the total number of **target log-density evaluations**, not just sample count; record wall-clock costs separately.
- Distinguish exploratory states from final corrected MH samples and fresh importance-weighted samples.
- Inspect mode transitions, posterior-functional accuracy, final-chain ESS, and importance-weight ESS as relevant.
- Preserve random seeds, exact priors/likelihoods, warmup rules, baseline specifications, and run environment.
- Treat custom `AdaptiveMixtureMH` and `SAChES_inspired` as **prototype controls**, not faithful published implementations.

The source files and recorded outputs are archived for inspection and replication, **not** evidence of general performance superiority or established convergence.
