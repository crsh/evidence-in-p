
pwappr <- function(p) {
  dplyr::case_when(
    p < 0.1 ~ -2*log(3*p),
    p < 0.5 ~ -4/3 * log(p) - log(16/9),
    p >= 0.5 ~ -0.5 * log(p)
  )
}

# exp(-0.5 * -2*log(3*p))               = 3p
# exp(-0.5 * -4/3 * log(p) - log(16/9)) = 9/16 * p^(2/3)
# exp(-0.5 * -0.5 * log(p))             = p^(1/4)

pwappr2 <- function(p) {
  dplyr::case_when(
    p < 0.1 ~ -sqrt(pi)*log(sqrt(2*exp(1))*p),
    p < 0.6 ~ -4/3 * log(p) - log(16/10),
    p >= 0.6 ~ -0.4 * log(p)
  )
}

exp(-0.5* -sqrt(pi)*log(sqrt(2*exp(1))*p)) = (pi - 1) * p^(sqrt(pi)/2)
exp(-0.5*-4/3 * log(p) - log(16/10))       = 5/8 * p^(2/3)
exp(-0.5*-0.4 * log(p))                    = p^(1/5)

# W = \begin{cases}
# -\sqrt{\pi} \log(\sqrt{2e}p) & \text{if } p < 0.1 \\
# -\frac{4}{3} \log(p) - \log(\frac{16}{10}) & \text{if } 0.1 \leq p < 0.6 \\
# -0.4 \log(p) & \text{if } p \geq 0.6
# \end{cases}

# \text{JAB}_{01} = \begin{cases}
# (\pi-1) & p^{9/10} & \sqrt{n} & \text{if } p < 0.1 \\
# 4/\pi & p^{2/3} & \sqrt{n} & \text{if } 0.1 \leq p < 0.6 \\
# & p^{1/5} & \sqrt{n} & \text{if } p \geq 0.6
# \end{cases}


dat <- tibble::tibble(
  p = seq(0, 1, length.out = 500)
  , W = qchisq(1-p, df = 1)
  , W1 = pwappr(p)
  , W2 = pwappr2(p)
) |>
  dplyr::mutate(
    piece = dplyr::case_when(
      p < 0.1 ~ "p < 0.1",
      p < 0.5 ~ "p < 0.5",
      p >= 0.5 ~ "p >= 0.5"
    )
  )

p1 <- dat |>
  ggplot() +
    aes(x = W) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_line(aes(y = W1, col = "3pn")) +
    geom_line(aes(y = W2, col = "2eppin")) +
    labs(x = "W", y = "~W") +
    coord_fixed()

p2 <- filter(dat, p < 0.1) |>
  ggplot() +
    aes(x = W) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_line(aes(y = W1, col = "3pn")) +
    geom_line(aes(y = W2, col = "2eppin")) +
    labs(x = "W", y = "~W") +
    coord_fixed()

p3 <- filter(dat, p < 0.5 & p > 0.1) |>
  ggplot() +
    aes(x = W) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_line(aes(y = W1, col = "3pn")) +
    geom_line(aes(y = W2, col = "2eppin")) +
    labs(x = "W", y = "~W") +
    coord_fixed()

p4 <- filter(dat, p > 0.5) |>
  ggplot() +
    aes(x = W) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_line(aes(y = W1, col = "3pn")) +
    geom_line(aes(y = W2, col = "2eppin")) +
    labs(x = "W", y = "~W") +
    coord_fixed()

(p1 + p2) / 
  (p3 + p4) +
  plot_layout(guides = "collect", axes = "collect")


p1 <- filter(dat, W < 1) |>
  ggplot() +
    aes(x = W) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_line(aes(y = W1, col = "3pn")) +
    geom_line(aes(y = W2, col = "2eppin")) +
    labs(x = "W", y = "~W") +
    coord_fixed()

p2 <- filter(dat, W > 2 & W < 3) |>
  ggplot() +
    aes(x = W) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_line(aes(y = W1, col = "3pn")) +
    geom_line(aes(y = W2, col = "2eppin")) +
    labs(x = "W", y = "~W") +
    coord_fixed()

p1 + p2 +
  plot_layout(guides = "collect", axes = "collect")


# rbind(
#   dat
#   , mutate(dat, piece = "all") |>
#     filter(W < 6)
# ) |>
#   ggplot() +
#   aes(x = 1 - p) +
#   geom_line(aes(y = W), linewidth = 2, col = "black") +
#   geom_line(aes(y = W1, col = "3pn")) +
#   geom_line(aes(y = W2, col = "2eppin")) +
#   scale_x_continuous(trans = "log") +
#   ggh4x::facet_grid2(~piece, scales = "free", independent = "all") +
#   theme_minimal(base_size = 16)

# # plot(
# #   1-p
# #   , qchisq(1-p, df = 1)
# #   , type = "l"
# #   , ylim = c(0, 10)
# #   , lwd = 2.5
# #   , ylab = "W"
# # )
# # lines(
# #   1-p
# #   , pwappr(p)
# #   , col = "skyblue"
# #   , lwd = 2.5
# # )
# # lines(
# #   1-p
# #   , pwappr2(p)
# #   , col = "red"
# #   , lwd = 2.5
# # )
