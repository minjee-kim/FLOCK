# Historical archive inventory

An archived research package was recovered and inspected. It includes **14 original R scripts**, notes, CSV results, plots, and analytical support files. The inventory below records the exact historical filenames. **All 14 R scripts listed here are now committed under `scripts/`.** Their recovery does not imply any fresh correctness audit.

## R source files

- `scripts/flock_am_dogwalker_comparison.R`
- `scripts/flock_am_formation_comparison.R`
- `scripts/flock_dogwalker_theta2.R`
- `scripts/flock_experiment_2.R`
- `scripts/flock_experiment_3.R`
- `scripts/flock_experiment_4.R`
- `scripts/flock_experiment_5.R`
- `scripts/flock_final_decision.R`
- `scripts/flock_hard_multimodal_benchmark.R`
- `scripts/flock_importance_consolidation.R`
- `scripts/flock_map.R`
- `scripts/flock_map_demo.R`
- `scripts/flock_orthogonal_benchmarks.R`
- `scripts/flock_softorth_six_scenarios.R`

## Additional archived materials

- Nine per-script research notes under [`notes/`](notes/)
- Results provenance and original source checksums in `CODE_MANIFEST.csv`
- Recovered experiments 2–4 per-run tabular outputs and the incomplete Experiment 5 console transcript under [`data/`](data/)
- Orthogonal and Dog Walker summary tables and unverified importance-sampling pilot summary
- The four original plotted PNG figures remain in the separate historical archive; the two Python support scripts are committed under [`analysis/`](analysis/)
- `requirements.txt` and `.gitignore`

The full archive is available separately as `FLOCK_GitHub_Research_Package.zip`.

## Manuscript publication status

The working report has since been added to GitHub under [`paper/main.tex`](paper/main.tex) and [`paper/main.pdf`](paper/main.pdf), together with [`paper/references.bib`](paper/references.bib). The historical Overleaf README and the complete original PNG figure set are **not** committed. The current LaTeX source contains figure placeholders where needed, while the archived PDF preserves the original compiled figure content.

## Caution

Original R sources are historical prototypes and were not rerun when the archive was assembled. They may define conflicting global R functions, so run them in separate sessions.