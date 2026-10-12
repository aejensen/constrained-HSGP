# Shape-constrained HSGP regression

Code and manuscript for *Positive, monotone, and convex Gaussian process
regression with Laplacian eigenfunctions*, by Andreas Kryger Jensen.
The implementation is self-contained.

The Stan models fit increasing, decreasing, convex, or concave regression
curves under Gaussian observations. They impose the derivative restriction
throughout the specified interval using the coefficient construction in the
paper. The R helpers prepare the ordinary Gaussian posterior, call Stan,
and reconstruct curves and derivatives.

## Requirements

Use R with `cmdstanr`, `posterior`, and `ggplot2`, a CmdStan installation,
and FLINT 3 with `pkg-config`. The C++ threshold calculations use FLINT's
interval arithmetic. On macOS, FLINT and pkg-config are available through
`brew install flint pkg-config`. On other systems, install FLINT 3 and make
its `flint.pc` file available to pkg-config.

```r
install.packages(c("posterior", "ggplot2"))
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()  # only if CmdStan is not already installed
```

See the [CmdStanR installation guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)
for compiler setup. These scripts were checked with CmdStan 2.37.0 and
FLINT 3.6.0 on macOS. Use build paths without spaces. The supplied compiler
flag `-ffp-contract=off` is needed for the interval calculations, and
`-ffast-math` or `-Ofast` must not be enabled.

## Fit your own data

Work from this directory. With numeric observation vectors `x` and `y`:

```r
source("R/load.R")
model <- prepare_hsgp(x, y, noise_sd = .1,
  kernel = "matern", nu = 3.5, amplitude = 1, lengthscale = .3,
  rank = 24L, interval = c(0, 1), center = .5, half_width = 2,
  intercept_mean = 0, intercept_sd = 2, direction = "increasing")
fit <- sample_hsgp(model, output_dir = "output/my_fit")
fit$diagnostics
fit$summary  # coefficient and sampling-coordinate R-hat and effective sample sizes
predict_hsgp(fit, seq(0, 1, length.out = 201), derivative = 1L)
ordinary_hsgp(model, seq(0, 1, length.out = 201), derivative = 1L)
```

Set `direction` to `"increasing"`, `"decreasing"`, `"convex"`, or `"concave"`.
These are separate restrictions. A convex or concave fit does not impose
monotonicity. `prepare_monotone_hsgp()` and `prepare_convex_hsgp()` are also
available. Both include a Gaussian intercept. Supply `noise_sd` as a scalar
or one standard deviation per observation, or use `noise_cov` for correlated
Gaussian errors. `kernel = "se"` selects the squared exponential kernel.
Covariance and noise parameters are supplied before sampling.

Monotonicity uses ranks 4–128 and requires the shape interval to lie inside
`center +/- half_width/2`. Convexity uses ranks 2–32 and requires the interval
to lie strictly inside `center +/- half_width`. These rank limits belong to
the supplied numerical routines. Derivatives are with respect to the `x`
units supplied to preparation. Returned intervals are pointwise latent
credible intervals, without observation noise.

`threshold_absolute_tolerance` controls the conservative numerical threshold
(default `1e-10`, in coefficient units). The C++ routines bound its maximum
over the interval. Stan uses an envelope gradient evaluated at a numerically
maximizing location. Inspect the usual sampling diagnostics before inference.
Build products go in `.build/` and fit output goes in the requested directory.
After editing a C++ header, call `build_hsgp(..., force_recompile = TRUE)`.

## Reproduce the applications

```sh
Rscript examples/lidar/run.R
Rscript examples/leaf/run.R
Rscript examples/lidar/render.R output/lidar
Rscript examples/leaf/render.R output/leaf
```

The scripts use the paper's model parameters and sampling settings, and save
posterior summaries and predictive checks. Each accepts an optional output
directory as its first argument. The [lidar notes](examples/lidar/README.md)
and [leaf notes](examples/leaf/README.md) explain the data, transformations,
and optional parameter re-estimation. New Monte Carlo runs will differ from
the recorded values by sampling variation.

To reproduce the published figures from the bundled compact summaries,
without fitting either model:

```sh
Rscript examples/lidar/render.R
Rscript examples/leaf/render.R
```

## Manuscript and files

Four manuscript versions are included:

- [Main paper without proofs](manuscript/main_only.pdf).
- [Standalone supplement with proofs](manuscript/supplement.pdf).
- [Paper with proofs in the supplement](manuscript/main.pdf).
- [Paper with proofs inline](manuscript/inline_proofs.pdf).

The [LaTeX sources](manuscript/main.tex), bibliography, and figures are
included. Keep these individual files for editing. From `manuscript/`, run
`python3 export.py` to create submission sources in `submission/`:

- `main.tex` contains the main paper without proofs.
- `supplement.tex` contains the supplementary material with proofs.

Each exported file embeds its bibliography and cross-references. Only its
figures remain separate. The exporter uses `pdflatex` and BibTeX in a temporary
directory to resolve the references, then removes compilation intermediates.
Add `--pdf` to retain PDFs, or `--target all` to also export the combined
layouts with proofs at the end and inline. `make -C manuscript submission`
exports the submission pair with PDFs, while `make -C manuscript` creates
all four layouts. `python3 export.py --check` checks source expansion without
compilation. The [manuscript instructions](manuscript/README.md) give details.

| Directory | Contents |
| --- | --- |
| `R/` | Gaussian preparation, Stan fitting, and prediction |
| `stan/` | Two Stan models and their C++ interfaces |
| `include/` | Numerical threshold calculations |
| `examples/` | Observations, application scripts, and compact plotting data |
| `manuscript/` | Four PDF versions, LaTeX sources, figures, and export script |
