# Historical R prototypes

This repository currently contains only [`flock_map_demo.R`](flock_map_demo.R). It depends on the separate `flock_map.R` source, which is **not yet committed**, so the demonstration is **not independently runnable from this checkout**.

The original October 2026 research archive contains 14 standalone R scripts covering FLOCK experiments 2–5, adaptive Metropolis comparisons, orthogonal pioneers, Dog Walker, and importance consolidation. See [the historical source inventory](../ARCHIVE_INVENTORY.md). None of those missing benchmark scripts should be inferred to be available from the directory listing.

The original scripts are research snapshots rather than a reusable library. In particular, their global R function definitions and output conventions may conflict if multiple scripts are sourced within one R session. No fresh end-to-end reproduction is claimed here.

See [reproducibility status](../docs/reproducibility.md).
