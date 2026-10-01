# novl_SOM

Self-organizing maps for time series, built from scratch in R and Rcpp: a
Kohonen SOM, a DTW-SOM, fourteen variants benchmarked on twelve data sets,
chaos, rings and tearing maps, and now a brainstorm of what to build next.

## Reports

| Report | What it holds |
|---|---|
| `som_next_steps.Rmd` | **Where to go next**: the lessons of the benchmarks, a verified catalogue of improvements for partial-discharge and movement data, a demonstrator of nine distances and six training ingredients on four map sizes with every measure, a detection test, and the precision ladder from fp32 to 256-bit MPFR |
| `som_vs_dtw_som.Rmd`, `som_metrics.Rmd`, `som_benchmark.Rmd`, `chaotic_units.Rmd`, `tearing_*.Rmd`, `growing_vs_classic.Rmd`, `som_variants_summary.Rmd` | the earlier reports of the project (kept in the author's working copy; this branch adds the next-steps report and its code) |

## som_next_steps.Rmd

Code in `next_steps/`:

| File | What it holds |
|---|---|
| `ns_core.cpp` | the online SOM as a template over the scalar type (fp32, fp64, 80-bit, 128-bit), nine distances in the contest and in the update (Euclidean, Manhattan, Chebyshev, Minkowski, cosine, correlation, DTW, soft-DTW, shape-based), the soft winner (deterministic annealing), conscience, the clipped pull, momentum; distances for judging; the topographic product and the ways along the map |
| `ns_mpfr.cpp` | the same loop in MPFR arbitrary precision (Euclidean and DTW), the reference of the precision ladder |
| `R/ns_data.R` | synthetic stand-ins (CBF, partial-discharge pulses, three-axis gestures) and the loader of the project's real data when `som_bench/data/datasets.rds` is present |
| `R/ns_maps.R` | lattices (hexagonal, rectangular, torus), the schedule, the seventeen variants, one training function |
| `R/ns_quality.R` | QE, EV, TE, KL, TP, used units, purity, ARI, NMI, trustworthiness, continuity, AUROC |
| `R/ns_bench.R` | the jobs, the parallel runs, the detection runs, the cache |
| `R/ns_precision.R` | the precision ladder |
| `R/ns_report.R` | chart and table style, map drawing |
| `tests/test_core.R` | checks of the compiled core against R implementations |
| `run_benchmark.R` | computes the cache outside the knit: `Rscript next_steps/run_benchmark.R [seeds] [passes] [cores]` |
| `ideas/ideas.json` | the brainstorm: ideas from seven angles, each verified by a skeptic, ranked into tiers |

Needs R with Rcpp, rmarkdown, plotly, DT, htmltools, jsonlite, digest, xfun,
Rmpfr, a C++17 compiler, and the GMP and MPFR libraries with headers
(`libmpfr-dev`, `libgmp-dev` on Debian and Ubuntu; `-lquadmath` comes with GCC).

Knit with `rmarkdown::render("som_next_steps.Rmd")`. The first knit trains
every variant on every data set and size with five seeds (about two hours of
CPU, half of it soft-DTW) and caches the results under `next_steps/cache`;
`run_benchmark.R` does the same on all cores. *Knit with parameters* changes
the seeds and passes, and the quick mode knits a small version for
development.
