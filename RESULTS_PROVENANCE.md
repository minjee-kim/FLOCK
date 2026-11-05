# Results provenance and interpretation limits

The numbers below are from recovered earlier research records, **not newly rerun experiments**.

## Grades

- **A — archived CSV:** numerical summary recovered from a saved experiment output.
- **B — console excerpt:** numerical values parsed from an incomplete terminal log.
- **C — previously reported:** secondary historical summary not independently verified against raw output.
- **D — code/archive only:** source exists but numerical output not checked.

| Experiment/output | Grade | What can be claimed |
|---|---|---|
| Experiment 2 summary, paired, and replicate CSV | A | 100-run targeted/follower/random comparison |
| Experiment 3 summary, paired, replicate, and discovery CSV | A | 100-run exploration-policy/adaptive-mixture comparison |
| Experiment 4 summary, paired, replicate, and discovery CSV | A | 100-run functional error/discovery comparison |
| Orthogonal benchmark CSV summaries | A | Only two replicates per method/scenario in recovered method summary |
| Dog Walker theta-squared posterior summaries | A | One benchmark configuration with 96,000 target evaluations per method |
| Experiment 5 console log | B | Only 25 complete visible replicate blocks out of requested 100, plus a partial 26th |
| Experiment 5 smoke summary | B | Three replicates; smoke check, not performance evidence |
| Importance pilot comparison, B=12,000 and n=10 | C | Earlier reported, raw ZIP not verified in original preparation |
| Hard benchmark and later-stage pilot ZIPs | D | No defensible aggregate ranking without checking their raw contents |

## Essential distinctions

1. Early 5D target mixture weights (0.45, 0.30, 0.20, 0.05) differ from later importance target weights (0.44, 0.30, 0.21, 0.05).
2. The Experiment 5 10D six-component mixture has weights (0.30, 0.25, 0.20, 0.15, 0.08, 0.02); overlapping components are not equivalent to six distinct density modes.
3. `SAChES_inspired` and `AdaptiveMixtureMH` were **custom comparators**, not faithful published SAChES or AIMM implementations.
4. Repulsive Dog Walker histories are exploratory, **not posterior samples**.
5. Corrected final MH draws and independently weighted importance draws are separate approaches.
6. Target-density call budgets are different from retained posterior sample sizes or wall-clock time.
7. Low R-hat alone does not establish mode coverage.
8. The exploratory model-fit/pivotal-quantity idea remained uncalibrated and is not an established goodness-of-fit test.

## Reproducibility expectations

For each scenario, method, replicate, and seed, log target specification and code version; equal target-evaluation budgets; acceptance; final retained draws; MCMC ESS and/or importance ESS; reference and estimated expectations; error and mode coverage. Do not select the primary endpoint after seeing results.

This repository is a record of attempts, including negative results, not a claim of a new generally superior algorithm.