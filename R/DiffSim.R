diffusion_SDT <- function(N, a, v, t0, z, sv, st0, sz=0, s=1, crit) {

  z <- z*a # convert relative starting point to absolute
  dt <- .001 # time step size

  sim_data <- matrix(0, nrow = N, ncol=3,
                     dimnames = list(NULL, c("RT","speeded_resp","delayed_resp"))
                     )

  NDT <- runif(N, t0 - st0/2, t0 + st0/2) # non-decision_times
  evidence <- rnorm(N, v, sv) # Sample SDT evidence strengths / drift rates
  drifts <- evidence* dt # Sample sampled evidence scale to instantaneous drift
  s <- sqrt(s^2 * dt) # scale drift coefficient to instantaneous s.d.

  if (sz > 0) {
    start_points <- runif(N, z-.5*sz, z+.5*sz)
  } else {
    start_points <- rep(z, N)
  }

  for (i in 1L:N){

    step <- 0
    pos <- start_points[i]
    v_inst <- drifts[i]

    while(!(pos > a | pos < 0 )) {
      step <- step + 1
      pos <- pos + rnorm(1, v_inst, s)
    }

    if (pos >= a){
      sim_data[i, "speeded_resp"] <- 1
    }

    sim_data[i, "RT"] <- step # Record number of steps
  }

  sim_data[, "RT"] <- NDT + sim_data[, "RT"]*dt # convert steps to RT

  above_crit_o <- evidence > crit[1]
  above_crit_n <- evidence > crit[2]
  said_old <- sim_data[, "speeded_resp"] == 1
  sim_data[above_crit_o & said_old, "delayed_resp"] <- 1
  sim_data[above_crit_n & !said_old, "delayed_resp"] <- 1

  return(sim_data)
}
