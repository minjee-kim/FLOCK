# Historical archive inventory

An archived research package was recovered and inspected. It includes **14 original R scripts**, notes, CSV results, plots, and analytical support files. The inventory below records the exact historical filenames; it does **not** imply each file is already committed to this repository.

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

- Nine per-script research notes under `notes/`
- Results provenance and original source checksums in `CODE_MANIFEST.csv`
- Recovered experiments 2–4 tabular outputs and incomplete experiment 5 console output
- Orthogonal and Dog Walker summary tables and unverified importance-sampling pilot summary
- Four plotted benchmark figures and two Python support scripts
- `requirements.txt` and `.gitignore`

The full archive is available separately as `FLOCK_GitHub_Research_Package.zip`.

## Manuscript publication status

The working report has since been added to GitHub under [`paper/main.tex`](paper/main.tex) and [`paper/main.pdf`](paper/main.pdf), together with [`paper/references.bib`](paper/references.bib). The historical Overleaf README and the complete original PNG figure set are **not** committed. The current LaTeX source contains figure placeholders where needed, while the archived PDF preserves the original compiled figure content.

## Caution

Original R sources are historical prototypes and were not rerun when the archive was assembled. They may define conflicting global R functions, so run them in separate sessions.