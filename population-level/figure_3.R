# Mycelial growth model based on Davidson (1998)
# Tracks:
#  m  = biomass
#  si = internal resource
#  se = external resource

library(deSolve)
library(fields)

# Parameters (from Davidson 1998 Fig. 5)
params <- list(
  c1 = 6.0,     # biomass growth rate
  c2 = 6.1,     # internal resource consumption rate
  c3 = 0.05,    # active uptake rate
  c4 = 1.0e-3,  # passive diffusion across cell wall
  c5 = 0.01,    # external-to-internal volume ratio
  c6 = 2.0,     # active uptake scaling
  k1 = 1.0,     # half-saturation constant (internal)
  k2 = 0.1,     # half-saturation constant (external)
  Di = 0.01,    # internal resource diffusion (translocation)
  Dm = 0.02     # biomass diffusion
)

# Grid: 3x3cm split into 51x51 cells
Lx <- 3.0
Ly <- 3.0
Nx <- 51
Ny <- 51
dx <- Lx / (Nx - 1)
dy <- Ly / (Ny - 1)
xgrid <- seq(0, Lx, length.out = Nx)
ygrid <- seq(0, Ly, length.out = Ny)
X <- outer(rep(1, Ny), xgrid)
Y <- outer(ygrid, rep(1, Nx))
Ncells <- Nx * Ny

# Integration settings (smaller atol/rtol = more accurate but slower)
t_end <- 400
atol  <- 1e-8     # absolute error
rtol  <- 1e-6     # relative error
lrw   <- 1200000  # memory for solver

# Distance from grid centre
r_centre <- sqrt((X - 1.5)^2 + (Y - 1.5)^2)

# Small circular biomass inoculum and internal resource at centre (radius 0.2)
m_init  <- ifelse(r_centre <= 0.2, 0.15, 0)
si_init <- ifelse(r_centre <= 0.2, 0.15, 0)

# Patch of external resource at centre (radius 0.25)
# Rest of grid has no resource
se_patch <- 0.12 # patch resource concentration
se_init  <- matrix(0, nrow = Ny, ncol = Nx)
se_init[r_centre <= 0.25] <- se_patch

# Computes div(a * grad(u)) : (describes how a substance 'u' spreads across the grid)
# Weighted by the concentration of another substance 'a' (Biomass moves faster where there is more internal resource)
div_a_grad_u <- function(a, u) {
  
  # Average of a between each pair of neighbouring cells
  a_x <- 0.5 * (a[, 2:Nx] + a[, 1:(Nx-1)])
  a_y <- 0.5 * (a[2:Ny, ] + a[1:(Ny-1), ])
  
  # How steeply u changes between neighbours
  du_dx_face <- (u[, 2:Nx] - u[, 1:(Nx-1)]) / dx
  du_dy_face <- (u[2:Ny, ] - u[1:(Ny-1), ]) / dy
  
  # Flow of u between neighbours (weighted by a)
  flux_x <- a_x * du_dx_face
  flux_y <- a_y * du_dy_face
  
  # Net flow into each cell
  div <- matrix(0, nrow = Ny, ncol = Nx)
  # x direction
  div[, 1] <- flux_x[, 1] / dx
  div[, 2:(Nx-1)] <- (flux_x[, 2:(Nx-1)] - flux_x[, 1:(Nx-2)]) / dx
  div[, Nx] <- -flux_x[, Nx-1] / dx
  # y direction
  div[1, ] <- div[1, ] + flux_y[1, ] / dy
  div[2:(Ny-1), ] <- div[2:(Ny-1), ] + (flux_y[2:(Ny-1), ] - flux_y[1:(Ny-2), ]) / dy
  div[Ny, ] <- div[Ny, ] - flux_y[Ny-1, ] / dy
  
  div
}

# Flatten three 2D grids into 1D vector (for ODE solver)
pack <- function(m, si, se) c(as.vector(m), as.vector(si), as.vector(se))

# Unpack output back into three 2D grids
unpack <- function(y) {
  m  <- matrix(y[1:Ncells], nrow = Ny, ncol = Nx)
  si <- matrix(y[(Ncells+1):(2*Ncells)], nrow = Ny, ncol = Nx)
  se <- matrix(y[(2*Ncells+1):(3*Ncells)], nrow = Ny, ncol = Nx)
  list(m = m, si = si, se = se)
}

# Rate of change for each substance (Eq. 2)
rhs <- function(t, y, run_params) {
  state <- unpack(y)
  
  # Densities can't go negative
  m  <- pmax(state$m, 0)
  si <- pmax(state$si, 0)
  se <- pmax(state$se, 0)
  
  c1 <- run_params$c1; c2 <- run_params$c2; c3 <- run_params$c3; c4 <- run_params$c4
  c5 <- run_params$c5; c6 <- run_params$c6; k1 <- run_params$k1; k2 <- run_params$k2
  Di <- run_params$Di; Dm <- run_params$Dm
  
  div_m  <- Dm * div_a_grad_u(si, m)    # Biomass spread proportional to internal resource
  div_si <- Di * div_a_grad_u(m,  si)   # Internal resource spread proportional to biomass (translocation)
  div_se <- 0                           # External resource doesn't spread
  
  # Biomass: autocatalytic growth fuelled by internal resource (limited by space)
  fm <- c1 * m^2 * ( si / (k1 + si) - m )
  
  # Internal resource: consumed by growth and gained from external uptake
  # Lost/gained by passive exchange across hyphal wall
  fi <- -c2 * m^2 * ( si / (k1 + si) ) +
    c3 * si * ( se / (k2 + se) ) -
    c4 * m * (si - se)
  
  # External resource: lost via uptake
  # Lost/gained by passive exchange across hyphal wall (scaled by c5)
  fe <- c5 * ( -c6 * m * ( se / (k2 + se) ) + c4 * m * (si - se) )
  
  list(pack(div_m + fm, div_si + fi, div_se + fe))
}

# Run from t=0 to t_end
run <- function(Di) {
  y0 <- pack(m_init, si_init, se_init)
  
  sol <- ode.2D(y = y0, times = c(0, t_end), func = rhs,
                parms = modifyList(params, list(Di = Di)),
                nspec = 3, dimens = c(Ny, Nx), method = "lsodes",
                lrw = lrw, atol = atol, rtol = rtol)
  unpack(pmax(sol[2, -1], 0))
}

# Run with and without translocation
with_Di    <- run(params$Di)
without_Di <- run(0)

# Plot: 3 columns (external resource, biomass, internal resource)
# 2 rows (without / with translocation) colourbar above each column
par(cex = 1.5)
cols <- tim.colors(256)

# Highest value on each colour scale
zlim_se <- c(0, se_patch)
zlim_m  <- c(0, ceiling(max(with_Di$m,  without_Di$m)  / 0.01) * 0.01)
zlim_si <- c(0, ceiling(max(with_Di$si, without_Di$si) / 0.01) * 0.01)
n_ticks <- 5

# Set up plot grid
layout(matrix(c(1, 0, 2, 0, 3,
                4, 0, 5, 0, 6,
                7, 0, 8, 0, 9), nrow = 3, byrow = TRUE),
       widths  = c(1, 0.12, 1, 0.12, 1),
       heights = c(0.22, 1, 1))
par(oma = c(1, 1, 1, 0.5))

# Draw colourbars
colourbar <- function(zlim, label) {
  
  at_vals <- seq(zlim[1], zlim[2], length.out = n_ticks)
  breaks  <- seq(zlim[1], zlim[2], length.out = length(cols) + 1)
  bar_mat <- matrix(seq_len(length(cols)), nrow = length(cols), ncol = 1)
  
  image(x = breaks, y = c(0, 1), z = bar_mat,
        col = cols,
        xaxt = "n", yaxt = "n",
        xlab = "", ylab = "",
        useRaster = TRUE)
  
  axis(1, at = at_vals, labels = format(at_vals), cex.axis = 1.2)
  box()
  title(main = label, line = 0.2, cex.main = 1.7)
}

# Draw heatmaps
plot_panel <- function(data, zlim, show_x = FALSE, show_y = FALSE) {
  image(xgrid, ygrid, t(data),
        col = cols,
        zlim = zlim,
        xaxt = "n", yaxt = "n",
        xlab = "", ylab = "",
        useRaster = TRUE)
  
  if (show_x) axis(1, at = c(0,1,2,3), cex.axis = 1.2)
  if (show_y) axis(2, at = c(0,1,2,3), cex.axis = 1.2)
  
  box()
}

# Row 1: colour bars
par(mar = c(1.2, 2.5, 2.2, 0.5))
colourbar(zlim_se, "External resource")
par(mar = c(1.2, 0.8, 2.2, 0.5))
colourbar(zlim_m, "Biomass")
par(mar = c(1.2, 0.8, 2.2, 0.5))
colourbar(zlim_si, "Internal resource")

# Row 2: without translocation (Di = 0)
par(mar = c(0.8, 2.5, 0.8, 0.5))
plot_panel(se_init, zlim_se, show_x = FALSE, show_y = TRUE)
par(mar = c(0.8, 0.8, 0.8, 0.5))
plot_panel(without_Di$m, zlim_m, show_x = FALSE, show_y = FALSE)
par(mar = c(0.8, 0.8, 0.8, 0.5))
plot_panel(without_Di$si, zlim_si, show_x = FALSE, show_y = FALSE)

# Row 3: with translocation (Di = 0.01)
par(mar = c(2.2, 2.5, 0.8, 0.5))
plot_panel(se_init, zlim_se, show_x = TRUE, show_y = TRUE)
par(mar = c(2.2, 0.8, 0.8, 0.5))
plot_panel(with_Di$m, zlim_m, show_x = TRUE, show_y = FALSE)
par(mar = c(2.2, 0.8, 0.8, 0.5))
plot_panel(with_Di$si, zlim_si, show_x = TRUE, show_y = FALSE)