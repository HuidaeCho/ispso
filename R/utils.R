#' Create ISPSO control parameters.
#'
#' Create a named list of optional algorithm and execution control parameters
#' for \code{ispso()}. All arguments are optional; if not specified, sensible
#' defaults are used. Problem-dependent quantities (e.g., dimension and bound
#' ranges) are derived inside \code{ispso()} from \code{bounds}.
#'
#' @param S Integer, optional. Swarm size (number of particles). If \code{NULL},
#'   a default is chosen as a function of the problem dimension.
#' @param exclusion_factor Numeric, optional. Factor used in the stopping
#'   criterion based on exclusions per particle.
#' @param maxiter Integer, optional. Maximum number of iterations.
#' @param xeps Numeric, optional. Threshold for position-based convergence.
#' @param feps Numeric, optional. Threshold for objective-value convergence.
#' @param age Integer, optional. Particle age threshold used in nesting criteria.
#'
#' @param vmax_factor Numeric, optional. Factor for maximum particle velocity
#'   relative to the parameter-range span.
#' @param vmax0_factor Numeric, optional. Factor for maximum initial particle
#'   velocity relative to the diagonal span of the search space.
#' @param rprey_factor Numeric, optional. Factor for the prey radius relative to
#'   the diagonal span of the search space.
#' @param rspecies_factor Numeric, optional. Factor for the speciation radius
#'   relative to the diagonal span of the search space.
#' @param rnest_factor Numeric, optional. Factor for the nesting radius relative
#'   to the diagonal span of the search space.
#'
#' @param cluster Cluster object, optional. If provided, parallel evaluation is
#'   performed using this cluster. The cluster is not created or stopped by
#'   \code{ispso()}.
#'
#' @param .deterministic Logical, optional. Internal/debug parameter controlling
#'   deterministic execution.
#' @param .stop_after_solutions Integer, optional. Internal/debug parameter
#'   controlling early stopping after a given number of solutions.
#' @param .distance_to_solution Numeric, optional. Internal/debug distance
#'   threshold (fraction of diagonal span) used to identify solutions.
#' @param .plot_method Character, optional. Internal/debug plotting mode.
#' @param .plot_delay Numeric, optional. Internal/debug plotting delay.
#'
#' @return Named list of control parameters for \code{ispso()}.
#'
#' @export
ispso_control <- function(
  # ---- Core swarm / stopping parameters ----
  S = NULL,
  # Stopping criteria: Stop if the number of exclusions per particle since the
  # last minimum is greater than exclusion_factor * max sol iter / average sol
  # iter. The more difficult the problem is (i.e., high max sol iter / average
  # sol iter), the more iterations the algoritm requires to stop.
  exclusion_factor = 3,
  # Maximum iteration
  maxiter = 2000L,
  # Small positive number close to 0
  xeps = 0.001,
  feps = 0.0001,
  # Nesting criteria for global and local optima using particles' ages
  # (NEST_BY_AGE)
  age = 10L,

  # ---- Problem-scaled factors ----
  # Maximum particle velocity
  # vmax = (s$xmax - s$xmin) * vmax_factor
  vmax_factor = 0.1,
  # Maximum initial particle velocity
  # vmax0 = diagonal(s) * vmax0_factor
  vmax0_factor = 0.001,
  # Search radius for preys: One particle has two memories (i.e., x and pbest).
  # When two particles collide with each other within prey, one particle takes
  # more desirable x and pbest from the two particles' memories, and the other
  # particle is replaced with a quasi-random particle using scrambled Sobol'
  # sequences (PREY).
  rprey_factor = 0.0001,
  # Speciation radius: Li (2004) recommends 0.05*L<=rspecies<=0.1*L.
  rspecies_factor = 0.1,
  # Nesting radius
  rnest_factor = 0.01,

  # ---- Parallel optimization ----
  cluster = NULL,

  # ---- Internal / debugging parameters ----
  # Deterministic run?
  .deterministic = FALSE,
  # Stop if all the solutions are found! This is only for writing a paper, not
  # for real problems because the number of actual solutions is not known in
  # most cases.
  .stop_after_solutions = 0,
  #.stop_after_solutions = -1,
  # (0, 1]: Fraction of the diagonal span of the search space.
  .distance_to_solution = 0.01,
  # Plot method
  .plot_method = "",
  #.plot_method = "density",
  #.plot_method = "movement",
  #.plot_method = "density,species",
  #.plot_method = "movement,species",
  # Plotting delay in seconds
  .plot_delay = 0
) {
  list(
    S = S,
    exclusion_factor = exclusion_factor,
    maxiter = as.integer(maxiter),
    xeps = xeps,
    feps = feps,
    age = as.integer(age),

    vmax_factor = vmax_factor,
    vmax0_factor = vmax0_factor,
    rprey_factor = rprey_factor,
    rspecies_factor = rspecies_factor,
    rnest_factor = rnest_factor,

    cluster = cluster,

    .deterministic = .deterministic,
    .stop_after_solutions = .stop_after_solutions,
    .distance_to_solution = .distance_to_solution,
    .plot_method = .plot_method,
    .plot_delay = .plot_delay
  )
}

#' Compute the ISPSO swarm size.
#'
#' Computes the swarm size \code{S} from the number of parameters implied by
#' \code{bounds}. If \code{control$S} is provided, it is returned unchanged.
#' Otherwise, the default rule is used: \eqn{S = 10 + \lfloor 2\sqrt{D} \rfloor},
#' where \eqn{D} is the number of parameters.
#'
#' This helper is useful when external workflows (e.g., SWAT+ runs) need to
#' allocate per-worker resources before calling \code{\link{ispso}}.
#'
#' @param bounds Named list of length \eqn{D}. Each element is a numeric vector
#'   of length 2 giving the lower and upper bounds for one parameter.
#' @param control Named list, optional. Control parameters created by
#'   \code{\link{ispso_control}}. If \code{control$S} is not \code{NULL}, it
#'   overrides the default rule.
#'
#' @return Integer. Suggested swarm size \code{S}.
#'
#' @export
ispso_swarm_size <- function(bounds, control = NULL) {
  if (
    !is.null(control) &&
      !is.null(control$S) &&
      is.numeric(control$S) &&
      control$S > 1
  ) {
    return(as.integer(control$S))
  }

  as.integer(10L + floor(2 * sqrt(length(bounds))))
}
