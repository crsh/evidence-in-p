library("ggplot2")
library("geomtextpath")

p_boundaries <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 1)

dat <- expand.grid(
  p = exp(seq(log(0.00005), log(1), length.out = 100))
  , n = exp(seq(log(3), log(10000), length.out = 100))
) |>
  transform(
    jab = 1 / (sqrt(pi / 2) * sqrt(n) * exp(-0.5 * qchisq(1 - p, df = 1)))
    , jab_approx = 1 / jab::jab_p(p, n)
  )

plot_settings <- list(
  scale_x_continuous(
    expand = expansion(0, 0)
    , breaks = c(5, 10, 20, 50, 100, 250, 500, 1000, 2500, 5000, 10000)
    , trans = "log"
    , name = bquote("Sample size" ~ italic(n))
  )
  , scale_y_continuous(
    expand = expansion(0, 0)
    , breaks = p_boundaries
    , labels = format(p_boundaries, scientific = FALSE, drop0trailing = TRUE)
    , trans = "log"
    , name = bquote(italic(p) * "-value")
  )
  , scale_fill_viridis_c(guide = "none")
  , theme_minimal(base_size = 16)
  , theme(
    axis.ticks.length = unit(5, "pt")
    , axis.ticks.x = element_line()
    , axis.ticks.y = element_line()
    , plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    , axis.title.x = element_text(margin = margin(t = 0.1, unit = "cm"))
    , axis.title.y = element_text(margin = margin(r = -0, unit = "cm"))
    , axis.text.x = element_text(angle = 30, hjust = 1)
  )
)

evidence_labels <- data.frame(
  n = c(17, 45, 75, 150, 350, 800, 1500, 3700)
  , p = c(0.00015, 0.0005, 0.00125, 0.003, 0.012, 0.055, 0.16, 0.4)
  , label = c(
    "Extreme", "Very strong", "Strong", "Moderate"
    , "Anecdotal", "Moderate", "Strong", "Very strong"
  )
  , angle = -c(20, 21, 21, 21, 22, 23, 24, 25)
)

to_reciprocal <- function(x) {
  ifelse(
    x > 1
    , as.character(round(x))
    , paste0("1/", round(1/x))
  )
}
 
ggplot(dat) +
  aes(x = n, y = p) +
  geom_raster(aes(fill = log(jab)), interpolate = TRUE) +
  # geom_hline(yintercept = p_boundaries, color = "white", alpha = 0.2) +
  geom_contour(
    aes(z = jab_approx)
    , color = "white"
    , linetype = "dotted"
    , breaks = c(1/30, 1/10, 1/3, 3, 10, 30, 100)
  ) +
  geom_textcontour(
    aes(z = jab, label = to_reciprocal(after_stat(level)))
    , color = "white"
    , breaks = c(1/30, 1/10, 1/3, 3, 10, 30, 100)
  ) +
  geom_text(
    aes(x = n, y = p, label = label, angle = angle)
    , data = evidence_labels
    , color = "white"
    , fontface = "bold"
    , size = 5
  ) +
  plot_settings
