# Atmospheric mercury concentration

Run `Rscript examples/lidar/run.R` from the repository root, then
`Rscript examples/lidar/render.R output/lidar`. An optional first argument
to `run.R` selects a different output directory. Rendering without an input
argument uses the bundled posterior summaries from the manuscript.

The example uses all 221 SemiPar observations over 390–720 meters. It models
the negative log signal ratio, with distance rescaled to [0, 1] and response
centered and divided by its observed range. The model has rank 128,
half-width 2, a Matérn 7/2 kernel, a Gaussian intercept, and a noise standard
deviation whose logarithm is quadratic in distance. `lidar_settings()` lists
the fitted parameters at full precision.

The script uses four chains with 400 warm-up and 2,000 retained iterations,
target acceptance 0.9, and maximum tree depth 10. Expectation propagation at
129 derivative locations supplies the starting point and dense metric. The
metric stays fixed during warm-up while the step size adapts. Sampling uses
the continuum threshold with tolerance 1e-12 in standardized coefficient
units. The derivative is converted back to distance in meters and divided
by 1.6e-5 to obtain mercury concentration under the DIAL approximation.

The output contains `curves.csv`, `predictive_checks.csv`, and `fit/` with
Stan draws and diagnostics. The five predictive discrepancies are the mean
squared standardized residual, its maximum absolute value, absolute adjacent
residual correlation, and signed mean residuals in the two halves of the
observed distance range. One replicate dataset is drawn for each posterior
curve.

The parameters are estimated by ordinary underlying-GP marginal likelihood
and shared by the ordinary and constrained fits. To repeat that optimization:

```r
source("R/load.R")
source("examples/lidar/example.R")
estimate_lidar_parameters(read_lidar_observations("examples/lidar"))
```

The [data notes](data/README.md) identify the source and preserve its license
notice. Both panels use the same plotting settings as the leaf application.
