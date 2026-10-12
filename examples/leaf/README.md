# Leaf light response

Run `Rscript examples/leaf/run.R` from the repository root, then
`Rscript examples/leaf/render.R output/leaf`. An optional first argument
to `run.R` selects a different output directory. Rendering without an input
argument uses the bundled posterior summaries from the manuscript.

The [data notes](data/README.md) describe the 16 raw records for leaf 1818,
an Eriophorum angustifolium leaf measured in Utqiagvik, Alaska. The script
applies the study's filters and removes the dark measurement, leaving 13
observations. Irradiance is centered and scaled to [-1, 1]. Net assimilation
is modeled with a concave rank-24 HSGP, half-width 2, a squared exponential
kernel, intercept prior N(0, 20^2), and independent Gaussian noise.

`prepare_leaf_example()` supplies the paper's parameters at full precision.
They are estimated by ordinary marginal likelihood at rank 32, then used in
both rank-24 posterior fits. To repeat the optimization:

```r
source("R/load.R")
source("examples/leaf/example.R")
estimate_leaf_parameters()
```

The constrained fit uses four chains with 1,500 warm-up and 3,000 retained
iterations, a dense metric, target acceptance 0.99, maximum tree depth 12,
and threshold tolerance 1e-8. The ordinary fit uses 20,000 independent
Gaussian posterior draws. Numerical maximization of the marginal density
supplies the constrained fit's starting point. First and second derivatives are transformed
back to irradiance units. The curvature figure multiplies the second
derivative by 10^6.

The output contains `curves.csv`, `application_summary.csv`,
`predictive_checks.csv`, and `fit/` with Stan draws and diagnostics.
Predictive checks compare the maximum absolute standardized residual and
the absolute sum of adjacent residual products in measurement-time order,
using one replicate dataset per posterior curve. The two application panels
and the supplementary curvature figure use the manuscript's shared theme.
