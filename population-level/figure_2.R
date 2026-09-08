# Hyphal tip orientation selection based on the model of Edelstein-Keshet & Ermentrout (1989)
# Tracks:
#  n(theta)   = density of tips in direction theta
#  rho(theta) = hyphal length in direction theta
# With branching angle (phi) <~16.8° anastomosis/removal of transverse tips outpaces branching and tips align

library(ggplot2) # For polar plot support

# Parameters (from Edelstein-Keshet & Ermentrout 1989 Fig. 4)
params <- list(
  alpha = 1.2,   # branch rate
  beta  = 5.0,   # anastomosis rate
  gamma = 0.2,   # death rate
  v     = 1.0    # growth rate
)

# Split orientation into N angles around circle
N      <- 180
theta  <- seq(0, 2 * pi, length.out = N + 1)
theta  <- theta[-length(theta)]
dtheta <- theta[2] - theta[1]

# Integration settings
t_end   <- 70
dt      <- 0.001
perturb <- 0.01 # 1% perturbation

# Collision kernel K(theta) proportional to |sin(theta)| (Eq. 17)
# Scaled so each row adds to 1 (Eq. 31b)
# Branches at right angles most likely to fuse
Kmat <- abs(sin(outer(theta, theta, "-")))
Kmat <- Kmat / sum(Kmat[1, ])

# Solve Eqs. 22a-b from a slightly perturbed equilibrium
run <- function(phi) {
  alpha <- params$alpha
  beta  <- params$beta
  gamma <- params$gamma
  v     <- params$v
  
  # Calculate equilibrium/uniform steady state (Eq. 22a-b with dn/dt = drho/dt = 0)
  rho_bar <- alpha / beta
  n_bar   <- gamma * alpha / (v * beta)
  
  shift <- as.integer(round(phi / dtheta)) %% N
  
  # Indices for n(theta - phi) and n(theta + phi)
  idx       <- seq_len(N)
  idx_minus <- ((idx - shift - 1) %% N) + 1
  idx_plus  <- ((idx + shift - 1) %% N) + 1
  
  rho <- rep(rho_bar, N)
  n   <- rep(n_bar, N) + perturb * n_bar * rnorm(N)
  
  for (step in seq_len(round(t_end / dt))) {
    # Branching: parent tip at theta disappears, two new tips at theta +/- phi (Eq. 12)
    branch_term <- alpha * (n[idx_minus] - n + n[idx_plus])
    
    # Anastomosis: tips fuse with branches weighted by kernel K (more likely when perpendicular)
    conv_rho    <- as.vector(Kmat %*% rho)
    cross_term  <- beta * n * conv_rho
    
    # Rate of change (Eqs. 22a-b)
    dn_dt   <- branch_term - cross_term
    drho_dt <- v * n - gamma * rho
    
    n   <- n   + dt * dn_dt
    rho <- rho + dt * drho_dt
    
    # Densities can't go negative
    n   <- pmax(n, 0)
    rho <- pmax(rho, 0)
  }
  n
}

# Run with small and large branching angles
phi_small <- 0.10
phi_large <- 0.40
set.seed(0)
results <- lapply(c(small = phi_small, large = phi_large), run)

# Bin tip densities for plot
bin_deg <- 10

bin_angles <- function(r) {
  theta_deg   <- theta * 180 / pi
  bin_breaks  <- seq(0, 360, by = bin_deg)
  bin_centers <- bin_breaks[-length(bin_breaks)] + bin_deg / 2
  bin_idx     <- cut(theta_deg, breaks = bin_breaks,
                     labels = FALSE, include.lowest = TRUE, right = FALSE)
  bin_idx     <- factor(bin_idx, levels = seq_along(bin_centers))
  binned_r    <- tapply(r, bin_idx, mean)
  data.frame(theta_deg = bin_centers, r = as.numeric(binned_r))
}

df <- do.call(rbind, Map(function(r, case) {
  d      <- bin_angles(r)
  d$r    <- 100 * d$r / sum(d$r, na.rm = TRUE)
  d$case <- case
  d
}, results, names(results)))

r_max    <- max(df$r, na.rm = TRUE)
uniform  <- 100 / (360 / bin_deg)   # expected % per bin if orientations were random
df$case  <- factor(df$case, levels = c("large", "small"))
y_breaks <- pretty(c(0, r_max), n = 4)
y_breaks <- c(y_breaks[y_breaks < r_max], r_max)

axis_deg <- seq(0, 315, by = 45)

# Small phi = alignment. Large phi = diffuse. Red dashed line = uniform distribution.
ggplot(df, aes(x = theta_deg, y = r)) +
  geom_col(width = bin_deg, fill = "blue",
           colour = "black", linewidth = 0.1) +
  geom_hline(yintercept = uniform, colour = "red",
             linetype = "22", linewidth = 0.4) +
  coord_radial(start = -pi / 2, expand = FALSE, clip = "off") +
  facet_wrap(~ case) +
  scale_x_continuous(
    limits = c(0, 360),
    breaks = axis_deg,
    labels = paste0(axis_deg, "°")
  ) +
  scale_y_sqrt(
    limits = c(0, r_max),
    breaks = y_breaks
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text       = element_blank(),
    axis.title       = element_blank(),
    axis.text.y      = element_blank(),
    axis.text.x      = element_text(size = 9, margin = margin(t = 8)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, colour = "black"),
    panel.spacing    = unit(2, "lines"),
    plot.margin      = margin(10, 10, 10, 10)
  )