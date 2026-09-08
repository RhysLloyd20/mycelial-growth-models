# Mycelial growth model based on Regalado et al. (1996)
# Tracks:
#  a = activator
#  h = inhibitor
#  s = resource
#  y = whether a hypha has grown there

# Parameters (from Regalado et al. 1996 Table 1)
params <- list(
  D_a   = 0.015,   # activator spread rate
  D_h   = 0.2,     # inhibitor spread rate
  D_s   = 0.01,    # resource spread rate
  u     = 0.13,    # activator decay
  v     = 0.12,    # inhibitor decay
  cc    = 0.004,   # activator production rate
  cp    = 0.004,   # inhibitor production rate
  k     = 0.04,    # saturation
  rho   = 0.055,   # baseline activator production
  eps   = 0.009,   # resource consumption rate
  d     = 0.0014,  # how strongly activator drives hyphal growth
  e     = 0.1,     # activator/hyphal decay
  f     = 9.0,     # activator/hyphal self-reinforcement
  astar = 1.8      # activator level needed to trigger hyphal growth
)

# Grid: 250x250 cells
N      <- 250
centre <- N %/% 2

# Integration settings
t_end <- 6000
dt    <- 0.5

# Resource: higher at top vs bottom
s_low  <- 1.1
s_high <- 1.5

# Works out for each grid cell how much of a chemical flows in or out from its 4 neighbours
# (From high to lower concentration cells)
laplacian <- function(field) {
  lap <- matrix(0, N, N)
  lap[2:(N-1), 2:(N-1)] <- (
    field[1:(N-2), 2:(N-1)] +
      field[3: N, 2:(N-1)] +
      field[2:(N-1), 1:(N-2)] +
      field[2:(N-1), 3:N] -
      4 * field[2:(N-1), 2:(N-1)]
  )
  lap
}

# Edge cells copy inner neighbours so nothing flows in/out
no_flux <- function(field) {
  field[1, ] <- field[2, ]; field[N, ] <- field[N-1, ]
  field[, 1] <- field[, 2]; field[, N] <- field[, N-1]
  field
}

run <- function() {
  
  D_a <- params$D_a; D_h <- params$D_h; D_s <- params$D_s
  u <- params$u; v <- params$v; k <- params$k
  cc <- params$cc; cp <- params$cp; rho <- params$rho
  eps <- params$eps; d <- params$d; e <- params$e
  f <- params$f; astar <- params$astar
  
  # Starting state (no mycelium)
  a <- matrix(0, N, N)
  h <- matrix(0.1, N, N)
  y <- matrix(0, N, N)
  s <- matrix(s_low, N, N)
  s[, (centre+1):N] <- s_high
  
  # Small random variation in production rates
  # (Otherwise pattern is symmetrical)
  cc_noise <- cc * (1 + 0.04 * (2 * matrix(runif(N*N), N, N) - 1))
  cp_noise <- cp * (1 + 0.04 * (2 * matrix(runif(N*N), N, N) - 1))
  
  # Starting spore at centre (circular blob of activator)
  idx <- (centre-3):(centre+3)
  dsq <- outer((-3:3)^2, (-3:3)^2, "+")
  a[idx, idx] <- ifelse(dsq <= 9, astar * 1.5 * exp(-dsq / 4), 0)
  y[idx, idx] <- ifelse(dsq <= 9, 0.1, 0)
  
  # Set steps to take snapshots
  snapshots  <- list()
  nsteps     <- round(t_end / dt)
  save_steps <- nsteps * c(1, 2, 3) / 3
  
  for (step in seq_len(nsteps)) {
    
    # How much each substance spreads to neighbours
    lap_a <- laplacian(a)
    lap_h <- laplacian(h)
    lap_s <- laplacian(s)
    
    # Activator promotes its own production fuelled by resource
    numer <- cc_noise * (a ^ 2) * s
    denom <- h * (1 + k * (a ^ 2)) # kept in check by inhibitor
    auto <- numer / denom   # Net activator produced
    
    # Rate of change (Eqs. 1a-c)
    da_dt <- D_a * lap_a + auto + rho * y - u * a    # spreading + growth - decay
    dh_dt <- D_h * lap_h + cp_noise * (a ^ 2) - v * h
    ds_dt <- D_s * lap_s - eps * s * y               # resource consumed where hyphae are
    
    # Hyphae form where activator >= threshold (Eq. 1d)
    active <- (a >= astar)
    dy_dt <- (d * a - e * y + y^2 / (1 + f * y^2)) * active
    
    # Time step
    a <- a + dt * da_dt
    h <- h + dt * dh_dt
    s <- s + dt * ds_dt
    y <- y + dt * dy_dt
    
    # Stop impossible negative values
    a <- pmax(a, 0)
    h <- pmax(h, 1e-6) # just above 0 to prevent any /0 attempt
    s <- pmax(s, 0)
    y <- pmax(y, 0)
    
    a <- no_flux(a)
    h <- no_flux(h)
    s <- no_flux(s)
    y <- no_flux(y)
    
    # Save snapshots at steps chosen above
    if (step %in% save_steps) {
      snapshots[[length(snapshots) + 1]] <- list(s = s, y = y)
    }
  }
  
  snapshots
}

# Run simulation
set.seed(123)
snapshots <- run()

# Plot snapshots
layout(matrix(1:4, nrow = 1), widths = c(1, 1, 1, 0.6))
cols <- rev(heat.colors(100))
zlim <- c(0, s_high)

par(mar = c(0.5, 0.5, 0.5, 0.5))
for (i in 1:3) {
  image(snapshots[[i]]$s, col = cols, zlim = zlim, axes = FALSE)
  overlay <- ifelse(snapshots[[i]]$y > 0.1, 1, NA)
  image(overlay, col = "black", axes = FALSE, add = TRUE)
}

# Colour bar and legend
plot.new()
yy <- seq(0.15, 0.9, length.out = 101)
for (j in 1:100) rect(0.1, yy[j], 0.4, yy[j + 1], col = cols[j], border = NA)
rect(0.1, 0.15, 0.4, 0.9, border = "black")
text(0.25, 0.95, "Resource", cex = 1.25)
text(0.45, 0.15, zlim[1], adj = c(0, 0.5))
text(0.45, 0.9, zlim[2], adj = c(0, 0.5))
rect(0.1, 0.03, 0.2, 0.07, col = "black")
text(0.25, 0.05, "Hyphae", adj = c(0, 0.5), cex = 1.25)