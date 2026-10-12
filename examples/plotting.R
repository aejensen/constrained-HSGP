# Shared figure rendering from compact posterior summaries. Both applications
# use identical sizes, colors, line widths, and typography. No fitting occurs.
suppressPackageStartupMessages(library(ggplot2))
colors <- c(ordinary = "#BC7040", constrained = "#16649B")
style <- list(panel_title_pt = 11.5, axis_title_pt = 9.5,
              axis_text_pt = 8.5, legend_text_pt = 9,
              line_width = .7, ribbon_alpha = .16,
              point_size = 1.1, point_alpha = .75,
              point_color = "#343A40", width_in = 6.5, height_in = 3.1)
# Keep panel headings larger than axis titles and put the legend above the
# two side-by-side application panels. Dimensions below are in inches.
application_theme <- function() {
  theme_minimal(base_size = 10, base_family = "sans") +
    theme(legend.position = "top", legend.title = element_blank(),
          legend.text = element_text(size = style$legend_text_pt),
          axis.title = element_text(size = style$axis_title_pt, color = "#222222"),
          axis.text = element_text(size = style$axis_text_pt, color = "#555555"),
          strip.text = element_text(size = style$panel_title_pt, color = "#222222",
                                    margin = margin(4, 0, 6, 0)),
          panel.grid.minor = element_blank(),
          panel.spacing.x = grid::unit(1.2, "lines"),
          plot.margin = margin(6, 8, 6, 6))
}
# Use consistent ordinary/constrained colors with an application-specific label.
model_scales <- function(constrained_label) {
  model_labels <- c("Ordinary GP", constrained_label)
  list(scale_color_manual(values = colors, breaks = names(colors), labels = model_labels),
       scale_fill_manual(values = colors, breaks = names(colors), labels = model_labels))
}

# curves contains response and concentration columns from lidar_plot_summaries.
# observed supplies original negative-logratio measurements for the first panel.
render_lidar <- function(curves, observed, output) {
keep <- c("ordinary_hsgp", "global_monotone_hsgp")
curves <- curves[curves$method %in% keep, , drop = FALSE]
stopifnot(nrow(observed) == 221L, length(unique(curves$method)) == 2L)
panels <- c("atop(textstyle(bold('Negative log signal ratio')),scriptstyle('(dimensionless)'))",
            "atop(textstyle(bold('Mercury concentration')),scriptstyle(ng~m^{-3}))")
plot_data <- rbind(
  data.frame(range_m = curves$range_m, method = curves$method, target = panels[1L],
             mean = curves$mean, lower = curves$lower, upper = curves$upper),
  data.frame(range_m = curves$range_m, method = curves$method, target = panels[2L],
             mean = curves$mercury_ng_m3_mean, lower = curves$mercury_ng_m3_lower,
             upper = curves$mercury_ng_m3_upper))
plot_data$target <- factor(plot_data$target, levels = panels)
plot_data$method <- factor(plot_data$method, levels = keep,
                          labels = names(colors))
points <- data.frame(range_m = observed$range_m, value = observed$negative_logratio,
                     target = factor(panels[1L], levels = panels))
zero <- data.frame(target = factor(panels[2L], levels = panels), value = 0)
lidar <- ggplot(plot_data, aes(range_m, mean, color = method, fill = method)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = style$ribbon_alpha, color = NA) +
  geom_hline(data = zero, aes(yintercept = value), inherit.aes = FALSE,
             color = "#666666", linetype = "dashed", linewidth = .35) +
  geom_line(linewidth = style$line_width) +
  geom_point(data = points, aes(range_m, value), inherit.aes = FALSE,
             color = style$point_color, size = style$point_size, alpha = style$point_alpha) +
  facet_wrap(~target, ncol = 2L, scales = "free_y", labeller = label_parsed) +
  model_scales("Monotone GP") +
  labs(x = "Distance (m)", y = NULL, color = NULL, fill = NULL) +
  application_theme()
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(output, "lidar_monotone_fit.pdf"), lidar,
    width = style$width_in, height = style$height_in, bg = "white")
}

# leaf contains the three quantities from leaf_plot_summaries. obs supplies
# recorded Photo and PARi values. Save response/slope together and curvature
# separately, matching the main-paper and supplementary figures.
render_leaf <- function(leaf, obs, output) {
leaf$model <- factor(leaf$model, levels = c("Unconstrained", "Concave"), labels = names(colors))
curvature <- leaf[leaf$quantity == "Curvature (scaled by 10^6)", ]
leaf <- leaf[leaf$quantity %in% c("Net assimilation", "Marginal assimilation gain"), ]
labels <- c("Net assimilation" = "atop(textstyle(bold('Net assimilation')),scriptstyle(mu*mol~CO[2]~m^{-2}~s^{-1}))",
            "Marginal assimilation gain" = "atop(textstyle(bold('Marginal assimilation gain')),scriptstyle(mol~CO[2]~'per mol quanta'))")
leaf$quantity <- factor(leaf$quantity, levels = names(labels), labels = unname(labels))
points <- data.frame(q = obs$PARi, value = obs$Photo,
                     quantity = factor(unname(labels[1L]), levels = unname(labels)))
irradiance_label <- expression(paste("Incident irradiance (", mu, "mol quanta ", m^{-2}, " ", s^{-1}, ")"))
light <- ggplot(leaf, aes(q, mean, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = style$ribbon_alpha, color = NA) +
  geom_line(linewidth = style$line_width) +
  geom_point(data = points, aes(q, value), inherit.aes = FALSE,
             color = style$point_color, size = style$point_size, alpha = style$point_alpha) +
  facet_wrap(~quantity, ncol = 2L, scales = "free_y", labeller = label_parsed) +
  model_scales("Concave GP") +
  labs(x = irradiance_label, y = NULL, color = NULL, fill = NULL) +
  application_theme()
curvature_plot <- ggplot(curvature, aes(q, mean, color = model, fill = model)) +
  geom_hline(yintercept = 0, color = "#666666", linewidth = .35) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = style$ribbon_alpha, color = NA) +
  geom_line(linewidth = style$line_width) + model_scales("Concave GP") +
  labs(x = irradiance_label, y = expression("Second derivative" %*% 10^6),
       color = NULL, fill = NULL) + application_theme()
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(output, "light_response.pdf"), light,
    width = style$width_in, height = style$height_in, bg = "white")
  ggsave(file.path(output, "curvature.pdf"), curvature_plot,
    width = style$width_in, height = style$height_in, bg = "white")
}
