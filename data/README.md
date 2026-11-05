# Archived experimental outputs

These tables and logs are restored from the historical research archive. They contain **different targets, estimands, evaluation budgets, and numerical metrics**; do not combine the error values across unrelated experiments.

## Primary experiment files

- `flock_experiment_2_results.csv`, `flock_experiment_2_paired.csv`, `flock_experiment_2_summary.csv`: 100-run formation comparisons.
- `flock_experiment_3_results.csv`, `flock_experiment_3_paired.csv`, `flock_experiment_3_discovery_summary.csv`, `flock_experiment_3_summary.csv`: 100-run policy comparisons.
- `flock_experiment_4_results.csv`, `flock_experiment_4_paired.csv`, `flock_experiment_4_discovery_summary.csv`, `flock_experiment_4_summary.csv`: 100-run posterior-functional comparisons.
- `experiment5_visible_log_records.csv`, `experiment5_console_log_excerpt.txt`, `experiment5_complete_visible_reps_summary.csv`, `experiment5_smoke_summary_from_log.csv`: **incomplete** recovered Experiment 5 transcript; 25 complete replicate blocks.
- `method_summary.csv`, `method_posterior_summaries.csv`, `replicate_results.csv`, `paired_comparisons.csv`, `metric_estimates.csv`, and `reference_summary.csv`: additional prototype comparisons and exact/reference functionals.
- `importance_pilot_reported_prior_analysis_UNVERIFIED.csv`: historical **secondary report**, not raw-archive verified.

## Verification

```bash
python -m pip install -r requirements.txt
python analysis/verify_results.py
```

The checker validates consistency of archived summaries and verifies 25 complete Experiment 5 replicate blocks. It **does not rerun the R algorithms** or prove correctness, mixing, or convergence. See [result provenance](../RESULTS_PROVENANCE.md).
