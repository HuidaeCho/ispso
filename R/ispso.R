################################################################################
# Isolated-Speciation-Based Particle Swarm Optimization (ISPSO) extends
# Species-Based PSO (SPSO) for finding multiple global and local minima.
#
# Author: Huidae Cho, Ph.D. <grass4u@gmail.com>, Texas A&M University
#
# Requires: R <http://r-project.org> and R packages: fOptions, plotrix
#
# Available at: https://idea.isnew.info/ispso.html
#
# Cite this software as:
#   Cho, H., Kim, D., Olivera, F., Guikema, S. D., 2011. Enhanced Speciation in
#   Particle Swarm Optimization for Multi-Modal Problems. European Journal of
#   Operational Research 213 (1), 15-23.
#
# Isolated-Speciation-Based Particle Swarm Optimization (ISPSO)
# Copyright (C) 2008-2026, Huidae Cho <https://idea.isnew.info>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

#' Run the Isolated-Speciation-Based Particle Swarm Optimization (ISPSO)
#' algorithm
#'
#' @references
#' Cho, H., Kim, D., Olivera, F., Guikema, S. D., 2011. Enhanced Speciation in
#' Particle Swarm Optimization for Multi-Modal Problems. European Journal of
#' Operational Research 213 (1), 15–23.
#' <https://dot.org/10.1016/j.ejor.2011.02.026>.
#'
#' @param fn Function. Objective function to be minimized. Must accept a
#'   numeric vector of parameter values and return a single numeric value.
#' @param bounds Named list. Parameter bounds. Each element must be a
#'   numeric vector of length two specifying \code{c(min, max)} for that
#'   parameter. Parameter names are used internally to label results.
#' @param control Named list, optional. Algorithm and execution control
#'   parameters created by \code{ispso_control()}. If \code{NULL}, default
#'   settings are used.
#' @param init_pos Numeric matrix, optional. User-provided particle positions
#'   for a warm start. The first \code{min(control$S, nrow(init_pos))} particles
#'   are initialized from \code{init_pos} (rows correspond to particles).
#'   Remaining particles, if any, are initialized by the default initialization
#'   procedure. Must have \code{ncol(init_pos) == ndim}.
#' @param prev_pop Data frame, optional. Population object returned by a
#'   previous \code{ispso()} run for resuming the optimization. Only the last
#'   generation in \code{prev_pop} is used to resume the optimization.
#' @param prev_nests Data frame, optional. Nest object returned by a previous
#'   \code{ispso()} run for resuming the optimization.
#'
#' @export
ispso <- function(
  fn,
  bounds,
  control = ispso_control(),
  init_pos = c(),
  prev_pop = c(),
  prev_nests = c()
) {
  ##############################################################################
  # SUBROUTINES
  ##############################################################################

  mytryCatch <- function(expr, error) tryCatch(expr, error = function(e) error)
  plotswarm <- function(control) {
    control$.plot_method != "" &&
      control$.plot_method != "profile" &&
      control$.plot_method != "diversity" &&
      control$.plot_method != "mean_diversity" &&
      1
  }
  plotmethod <- function(control, method) {
    any(unlist(strsplit(control$.plot_method, ",")) == method)
  }
  mydist <- function(control, ...) {
    as.matrix(stats::dist(if (ndim == 1) c(...) else rbind(...)))
  }
  mydist2 <- function(...) sqrt(sum((...)^2))
  mynrow <- function(x) {
    if (is.null(x)) {
      0
    } else if (is.vector(x)) {
      1
    } else {
      nrow(x)
    }
  }
  myncol <- function(x) {
    if (is.null(x)) {
      0
    } else if (is.vector(x)) {
      length(x)
    } else {
      ncol(x)
    }
  }
  colmax <- function(x, ...) {
    ret <- c()
    for (i in 1:myncol(x)) {
      ret[i] <- max(x[, i], ...)
    }
    ret
  }
  colmin <- function(x, ...) {
    ret <- c()
    for (i in 1:myncol(x)) {
      ret[i] <- min(x[, i], ...)
    }
    ret
  }
  myround <- function(x, digits = 0) floor(x * 10^digits + 0.5) / 10^digits

  ##############################################################################
  # Evaluate function values.
  ##############################################################################
  evaluate_f <- function() {
    #-DEBUG---------------------------------------------------------------------
    if (plotswarm(control)) {
      if (ndim == 1) {
        .f <- c()
        for (.x in x) {
          .f <- c(.f, fn(.x))
        }
        if (plotmethod(control, "movement")) {
          graphics::plot(
            x,
            .f,
            col = 1:control$S,
            xlim = c(xmin[1], xmax[1]),
            ylim = c(min(.F), max(.F)),
            xlab = "x1",
            ylab = "f(x1)"
          )
        }
        #{DEBUG: We don't know true solutions in real problems.
        if (.have_sols) {
          .x <- fn(x[1, ], list(sol = TRUE))
          .fnests <- c()
          for (.x1 in .x) {
            .fnests <- c(.fnests, fn(.x1))
          }
          graphics::points(.x, .fnests, pch = 3, cex = 2, lwd = 2, col = "red")
          plotrix::draw.arc(
            .x,
            .fnests,
            .plot_rsol,
            0,
            2 * pi,
            col = "gray",
            lty = "dotted"
          )
        }
        #}
        if (!is.null(nests)) {
          .fnests <- c()
          for (.x in nests[, 1:ndim]) {
            .fnests <- c(.fnests, fn(.x))
          }
          graphics::points(
            nests[, 1:ndim],
            .fnests,
            pch = 4,
            cex = 2,
            lwd = 2
          )
          plotrix::draw.arc(
            nests[, 1:ndim],
            .fnests,
            rnest,
            0,
            2 * pi,
            col = "black",
            lty = "dotted"
          )
        }
        if (plotmethod(control, "density")) {
          graphics::points(x, .f, col = 1:control$S)
        }
      } else if (ndim >= 2) {
        #-----------------------------------------------------------------------
        if (plotmethod(control, "movement")) {
          graphics::plot(
            x[, control$.plot_x],
            col = 1:control$S,
            xlim = c(xmin[control$.plot_x[1]], xmax[control$.plot_x[1]]),
            ylim = c(xmin[control$.plot_x[2]], xmax[control$.plot_x[2]]),
            xlab = sprintf("x%d", control$.plot_x[1]),
            ylab = sprintf("x%d", control$.plot_x[2])
          )
          if (!is.null(prev_x)) {
            graphics::arrows(
              prev_x[, control$.plot_x[1]],
              prev_x[, control$.plot_x[2]],
              x[, control$.plot_x[1]],
              x[, control$.plot_x[2]],
              length = .1,
              col = 1:control$S
            )
          }
        }
        #{DEBUG: We don't know true solutions in real problems.
        if (.have_sols) {
          .x <- fn(x[1, ], list(sol = TRUE))
          graphics::points(
            matrix(.x[, control$.plot_x], nrow(.x), 2),
            pch = 3,
            cex = 2,
            lwd = 2,
            col = "red"
          )
          plotrix::draw.arc(
            .x[, control$.plot_x[1]],
            .x[, control$.plot_x[2]],
            .plot_rsol,
            0,
            2 * pi,
            col = "gray",
            lty = "dotted"
          )
        }
        #}
        if (!is.null(nests)) {
          graphics::points(
            matrix(nests[, control$.plot_x], nrow(nests), 2),
            pch = 4,
            cex = 2,
            lwd = 2
          )
          plotrix::draw.arc(
            nests[, control$.plot_x[1]],
            nests[, control$.plot_x[2]],
            rnest,
            0,
            2 * pi,
            col = "black",
            lty = "dotted"
          )
        }
        if (plotmethod(control, "density")) {
          graphics::points(x[, control$.plot_x], col = 1:control$S)
        }

        prev_x <<- x
        prev_f <<- f
      }
    }
    #---------------------------------------------------------------------------
    f <<- if (parallel) {
      unlist(parallel::clusterApplyLB(
        control$cluster,
        seq_len(control$S),
        function(i) {
          fn(
            x[i, ],
            list(
              worker_id = get(
                "worker_id",
                envir = .ispso_state,
                inherits = FALSE
              ),
              S = control$S,
              iter = iter,
              run = (iter - 1) * control$S + i
            )
          )
        }
      ))
    } else {
      c()
    }
    for (i in 1:control$S) {
      if (!parallel) {
        f[i] <<- fn(x[i, ])
      }
      if (
        f[i] < pbest[i, ndim + 1] ||
          (f[i] == Inf && pbest[i, ndim + 1] == Inf)
      ) {
        pbest[i, ] <<- c(x[i, ], f[i])
        if (
          f[i] < gbest[ndim + 1] ||
            (f[i] == Inf && gbest[ndim + 1] == Inf)
        ) {
          gbest <<- c(x[i, ], f[i])
          gb <<- i
        }
      }
    }
    evals <<- evals + control$S
    age <<- age + 1
    #-DEBUG---------------------------------------------------------------------

    if (ndim == 1) {
      if (plotmethod(control, "movement")) {
        if (!is.null(prev_f)) {
          graphics::arrows(
            prev_x[, 1],
            prev_f,
            x[, 1],
            f,
            length = .1,
            col = 1:control$S
          )
        }

        prev_x <<- x
        prev_f <<- f
      }
    }
    #---------------------------------------------------------------------------
  }

  ##############################################################################
  # Update particle velocities. (SPSO, Li, 2004)
  ##############################################################################
  update_v <- function() {
    ############################################################################
    # SPSO
    lbest <- matrix(nrow = control$S, ncol = ndim)
    l <- order(f)
    species <<- c()
    seed <<- c()
    isolated <- rep(1, control$S)
    for (i in 1:control$S) {
      if (is.null(seed)) {
        lbest[l[i], ] <- x[l[i], ]
        seed <<- l[i]
        species[seed] <<- seed
        next
      }
      n <- length(seed)
      found <- FALSE
      for (j in 1:n) {
        if (mydist2(x[seed[j], ] - x[l[i], ]) <= rspecies) {
          found <- TRUE
          isolated[c(l[i], seed[j])] <- 0
          species[l[i]] <<- seed[j]
          lbest[l[i], ] <- x[seed[j], ]
          break
        }
      }
      if (found) {
        next
      }
      n <- n + 1
      seed[n] <<- l[i]
      species[seed[n]] <<- seed[n]
      lbest[l[i], ] <- x[l[i], ]

      fseed <- f[l[i]]

      #{SPSO_NEIGHBOURING_SPECIES
      # A new species seed searches for superior particles within its
      # speciation radius that failed to form their own species, but
      # happened to belong to seeds with better fitness values.  This
      # behaviour allows superior particles to share their
      # information with neighbouring species having relatively poor
      # fitness values.
      for (j in (1:control$S)[-l[i]]) {
        if (
          f[j] < fseed &&
            mydist2(x[j, ] - x[l[i], ]) <= rspecies
        ) {
          lbest[l[i], ] <- x[j, ]
          fseed <- f[j]
        }
      }
      # End of SPSO_NEIGHBOURING_SPECIES}

      #{SPSO_NEIGHBOURING_PBESTS
      # A new species seed searches for superior pbests within its
      # speciation radius. This behaviour allows superior pbests to
      # share their information with neighbouring species having
      # relatively poor fitness values.
      for (j in (1:control$S)[-l[i]]) {
        if (
          pbest[j, ndim + 1] < fseed &&
            mydist2(pbest[j, 1:ndim] - x[l[i], ]) <= rspecies
        ) {
          lbest[l[i], ] <- pbest[j, 1:ndim]
          fseed <- pbest[j, ndim + 1]
        }
      }
      # End of SPSO_NEIGHBOURING_PBESTS}
    }

    #{SPSO_ISOLATED_SPECIES
    # Isolated particles form one species.
    if (any(isolated == 1)) {
      # Ignore isolated seeds
      seed <<- seed[isolated[seed] == 0]
      age[isolated == 1] <<- 1
      n <- length(seed) + 1
      tmp <- cbind(1:control$S, isolated, f)
      tmp <- tmp[order(tmp[, 3]), ]
      seed[n] <<- tmp[tmp[, 2] == 1, 1][1]
      species[seed[n]] <<- -seed[n]
      for (i in which(isolated == 1)) {
        lbest[i, ] <- x[seed[n], ]
        species[i] <<- -seed[n]
      }
    }
    # End of SPSO_ISOLATED_SPECIES}
    #-DEBUG---------------------------------------------------------------------
    if (plotswarm(control)) {
      if (ndim == 1) {
        .f <- c()
        for (.x in x[seed, ]) {
          .f <- c(.f, fn(.x))
        }
        graphics::points(x[seed, ], .f, col = seed, pch = 20)
        if (plotmethod(control, "species")) {
          plotrix::draw.arc(x[seed, ], .f, rspecies, 0, 2 * pi, col = seed)
        }
      } else if (ndim >= 2) {
        if (length(seed) == 1) {
          graphics::points(
            x[seed, control$.plot_x[1]],
            x[seed, control$.plot_x[2]],
            col = seed,
            pch = 20
          )
        } else {
          graphics::points(x[seed, control$.plot_x], col = seed, pch = 20)
        }
        if (plotmethod(control, "species")) {
          plotrix::draw.arc(
            x[seed, control$.plot_x[1]],
            x[seed, control$.plot_x[2]],
            rspecies,
            0,
            2 * pi,
            col = seed
          )
        }
      }
    }
    #---------------------------------------------------------------------------
    # Constriction PSO (Clerc and Kennedy, 2000)
    v <<- control$w *
      (v +
        control$c1 * t(stats::runif(ndim) * t(pbest[, 1:ndim] - x)) +
        control$c2 * t(stats::runif(ndim) * t(lbest[, 1:ndim] - x)))
    #{CHECK_NESTS
    # Check existing nests before flying to new points.
    if (!is.null(nests)) {
      for (i in seed) {
        if (
          any(
            mydist(control, x[i, ] + v[i, ], nests[, 1:ndim])[
              1,
              2:(nrow(nests) + 1)
            ] <
              2 * rnest
          )
        ) {
          num_exclusions <<- num_exclusions + 1
          # turbulence area: 2*rnest
          v[i, ] <<- v[i, ] + rspecies * stats::runif(ndim)
        }
      }
    }
    # End of CHECK_NESTS}

    v <<- t(pmax(t(v), -vmax))
    v <<- t(pmin(t(v), vmax))
    adjust_v()

    V <<- sqrt(rowSums(v^2))
    pop <<- rbind(pop, cbind(x, f, v = V, age))
  }

  #-----------------------------------------------------------------------------
  fly_away_and_substitute <- function(neighbours) {
    n <- sum(neighbours)
    x[neighbours, ] <<- new_x(n)
    v[neighbours, ] <<- new_v(n)
    age[neighbours] <<- 0
  }
  #-----------------------------------------------------------------------------

  ##############################################################################
  # Convergence check for SPSO
  ##############################################################################
  check_for_convergence <- function() {
    #{NEST_BY_AGE
    # Nesting criteria using particles' ages.
    for (i in seed) {
      if (age[i] < control$age) {
        next
      }

      halflife <- pop[control$S * (iter - 1:halflife.age) + i, ]

      if (
        exp(mean(log(
          (colmax(
            if (ndim == 1) {
              t(halflife[, 1:ndim])
            } else {
              halflife[, 1:ndim]
            }
          ) -
            colmin(
              if (ndim == 1) {
                t(halflife[, 1:ndim])
              } else {
                halflife[, 1:ndim]
              }
            )) /
            (xmax - xmin)
        ))) <=
          control$xeps &&
          stats::sd(halflife[, "f"]) <= control$feps
      ) {
        run <- evals - control$S + i
        nests <<- rbind(
          nests,
          c(x[i, ], f[i], V[i], age[i], run, evals),
          deparse.level = 0
        )
        num_exclusions_per_nest[nrow(nests)] <<- num_exclusions
        num_exclusions <<- 0
        fly_away_and_substitute(
          mydist(control, x[, 1:ndim])[i, ] <= rspecies
        )
        cat(sprintf(
          "\b\b\b\bf(x[%d,])=%g added at iter=%d, run=%d, evals=%d, nests=%d\n",
          i,
          f[i],
          iter,
          run,
          evals,
          nrow(nests)
        ))
      }
    }
    adjust_v()

    if (
      !is.null(nests) && !is.null(control$exclusion_factor) && num_exclusions
    ) {
      delta_sol_iters <- (nests[, ndim + 4] -
        c(0, nests[-nrow(nests), ndim + 4])) /
        control$S
      average_delta_sol_iter <- mean(delta_sol_iters)
      delta_curr_iter <- delta_sol_iters[nrow(nests)]

      func_difficulty <- max(delta_sol_iters) / average_delta_sol_iter

      cat(sprintf(
        "\b\b\b\b%d%%",
        min(
          100,
          myround(
            100 *
              (num_exclusions / control$S) /
              (control$exclusion_factor * func_difficulty)
          )
        )
      ))

      if (
        !control$dont_stop &&
          num_exclusions / control$S >
            control$exclusion_factor * func_difficulty
      ) {
        return(TRUE)
      }
    }

    if (!control$dont_stop) {
      if (control$.stop_after_solutions > 0) {
        return(control$.stop_after_solutions == mynrow(nests))
      }

      if (
        control$.stop_after_solutions < 0 &&
          .have_sols &&
          !is.null(nests) &&
          nrow(nests) == length(found)
      ) {
        return(TRUE)
      }
    }

    return(FALSE)
    # End of NEST_BY_AGE}
  }

  ##############################################################################
  # Update particle positions.
  ##############################################################################
  update_x <- function() {
    x <<- x + v

    #{PREY
    # Inferior particles are preyed by superior ones in neighbours. Note
    # that x and f do not correspond because x has been updated since f was
    # evaluated. Therefore, information sharing is based on the past
    # experiences (pbest and the previous position's f value)
    d <- mydist(control, x)
    n <- nrow(d)
    preyed <- rep(0, control$S)
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (preyed[j] || d[i, j] > rprey) {
          next
        }
        preyed[j] <- 1
        # share pbest infos
        if (pbest[i, ndim + 1] > pbest[j, ndim + 1]) {
          pbest[i, ] <<- pbest[j, ]
        }
        # share infos about the previous positions
        if (f[i] > f[j]) {
          x[i, ] <<- x[j, ]
          v[i, ] <<- v[j, ]
        }
        # new particles from Sobol' sequences
        x[j, ] <<- new_x()
        v[j, ] <<- new_v()
        pbest[j, ] <- c(x[j, ], BEST)
        age[j] <<- 0
      }
    }
    adjust_v()
    # End of PREY}

    #{CHECK_NESTS
    # Check existing nests before flying to new points.
    if (!is.null(nests)) {
      d <- mydist(control, nests[, 1:ndim], x)
      n <- nrow(nests)
      for (i in 1:n) {
        rnst <- rnest
        for (j in which(d[i, n + 1:control$S] <= rnst)) {
          if (j == gb) {
            gb <<- order(pbest[, ndim + 1])[2]
            gbest <<- pbest[gb, ]
          }
          num_exclusions <<- num_exclusions + 1
          x[j, ] <<- new_x()
          v[j, ] <<- new_v()
          pbest[j, ] <<- c(x[j, ], BEST)
          age[j] <<- 0
        }
      }
    }
    adjust_v()
    # End of CHECK_NESTS}
  }

  ##############################################################################
  # New particles' positions
  ##############################################################################
  new_x <- function(n = 1, seed = -1) {
    r <- if (seed >= 0) {
      randtoolbox::sobol(n, ndim, TRUE, 3, seed)
    } else {
      randtoolbox::sobol(n, ndim, FALSE, 3)
    }
    t(xmin + (xmax - xmin) * t(r))
  }

  ##############################################################################
  # New particles' velocities
  ##############################################################################
  new_v <- function(n = 1) {
    v <- c()
    for (i in 1:n) {
      r <- stats::runif(ndim)
      v <- rbind(v, vmax0 / sqrt(sum(r^2)) * r)
    }
    v
  }

  .runif.sobol <- function(
    n,
    dimension,
    init = TRUE,
    scrambling = 0,
    seed = 4711
  ) {
    matrix(stats::runif(n * dimension), n, dimension)
  }

  ##############################################################################
  # Adjust velocities
  ##############################################################################
  adjust_v <- function() {
    #{Confinement: random forth
    for (i in 1:control$S) {
      j <- which(x[i, ] + v[i, ] < xmin)
      k <- length(j)
      if (k) {
        v[i, j] <<- (x[i, j] - xmin[j]) * stats::runif(k)
      }
      j <- which(x[i, ] + v[i, ] > xmax)
      k <- length(j)
      if (k) {
        v[i, j] <<- (xmax[j] - x[i, j]) * stats::runif(k)
      }
    }
    # End of confinement}
  }

  ##############################################################################
  # VARIABLES
  ##############################################################################

  ndim <- length(bounds)
  xmin <- vapply(bounds, `[`, numeric(1), 1)
  xmax <- vapply(bounds, `[`, numeric(1), 2)

  diag_span <- sqrt(sum((xmax - xmin)^2))

  control$S <- ispso_swarm_size(bounds, control)

  vmax <- (xmax - xmin) * control$vmax_factor
  vmax0 <- diag_span * control$vmax0_factor
  rprey <- diag_span * control$rprey_factor
  rspecies <- diag_span * control$rspecies_factor
  rnest <- diag_span * control$rnest_factor

  parallel <- !is.null(control$cluster)
  if (parallel && length(control$cluster) > control$S) {
    warning(sprintf(
      "Cluster size (%d) exceeds swarm size (%d); extra workers will remain idle.",
      length(control$cluster),
      control$S
    ))
  }

  #-----------------------------------------------------------------------------
  # default values for debugging variables
  # Deterministic run?
  if (is.null(control$.deterministic)) {
    control$.deterministic <- FALSE
  }
  # Stop if all the solutions are found! This is only for writing a
  # paper, not for real problems because the number of actual solutions
  # is not known in most cases.
  if (is.null(control$.stop_after_solutions)) {
    control$.stop_after_solutions <- 0
  }
  # (0, 1]: Fraction of the diagonal span of the search space.
  if (is.null(control$.plot_distance_to_solution)) {
    control$.plot_distance_to_solution <- 0.05
  }
  # default plotting method for 1 and 2 dimensional problems: movement
  if (is.null(control$.plot_method)) {
    control$.plot_method <- "movement"
  }
  # default 2-d plot
  if (is.null(control$.plot_x)) {
    control$.plot_x <- 1:2
  }
  # no delay between plots
  if (is.null(control$.plot_delay)) {
    control$.plot_delay <- 0
  }
  # don't save plots
  if (is.null(control$.plot_save_prefix)) {
    control$.plot_save_prefix <- ""
  }
  #-----------------------------------------------------------------------------
  # Does the user provide the real solutions to fn()?
  if (
    any(
      mytryCatch(fn(rep(0, ndim), list(sol = TRUE)), "no_sols") == "no_sols"
    )
  ) {
    .have_sols <- FALSE
  } else {
    .have_sols <- TRUE
  }
  .plot_rsol <- sqrt(sum((xmax - xmin)^2)) *
    control$.plot_distance_to_solution

  if (control$S < 2) {
    stop("Swarm size (control$S) must be greater than 1.")
  }

  if (control$.plot_save_prefix == "") {
    .plot_save_format <- ""
  } else {
    .plot_save_format <- sprintf(
      "%s%%0%dd.png",
      control$.plot_save_prefix,
      floor(log10(control$maxiter) + 1)
    )
  }

  # Don't stop the algorithm until control$maxiter? FALSE by default
  if (is.null(control$dont_stop)) {
    control$dont_stop <- FALSE
  }

  # Constriction PSO (Clerc and Kennedy, 2000)
  if (is.null(control$c1)) {
    control$c1 <- 2.05
  }
  if (is.null(control$c2)) {
    control$c2 <- 2.05
  }
  if (is.null(control$w)) {
    control$w <- 2 /
      abs(
        2 -
          control$c1 -
          control$c2 -
          sqrt((control$c1 + control$c2)^2 - 4 * (control$c1 + control$c2))
      )
  }

  BEST <- Inf
  WORST <- -BEST

  # PSO
  if (control$.deterministic) {
    if (is.null(.ispso_state$seed_sobol)) {
      .ispso_state$seed_sobol <- 4711
    }
    if (!is.null(.ispso_state$seed_random)) {
      assign(".Random.seed", .ispso_state$seed_random, envir = .GlobalEnv)
    } else {
      set.seed(0)
      .ispso_state$seed_random <- .Random.seed
    }
  } else {
    .ispso_state$seed_sobol <- as.integer(stats::runif(1) * 100000)
    .ispso_state$seed_random <- .Random.seed
  }

  pop <- prev_pop
  nests <- prev_nests

  x <- if (is.null(pop)) {
    new_x(control$S, .ispso_state$seed_sobol)
  } else {
    pop[order(pop[, "f"])[1:control$S], 1:ndim]
  }
  if (!is.null(init_pos)) {
    n <- min(nrow(init_pos), control$S)
    x[1:n, ] <- init_pos[1:n, ]
  }
  v <- new_v(control$S)
  adjust_v()

  pbest <- matrix(nrow = control$S, ncol = ndim + 1)
  gbest <- c()
  gb <- 0
  pbest[, ndim + 1] <- gbest[ndim + 1] <- BEST
  prev_gbestf <- WORST
  f <- c()
  age <- rep(0, control$S)
  V <- c()
  halflife.age <- myround(0.5 * control$age)

  seed <- c()
  species <- c()
  #{DEBUG: We don't know true solutions in real problems.
  if (.have_sols) {
    found <- rep(0, nrow(fn(x[1, ], list(sol = TRUE))))
  }
  #}

  #-DEBUG-----------------------------------------------------------------------
  if (plotswarm(control)) {
    grDevices::palette(grDevices::rainbow(control$S))
    prev_x <- prev_f <- c()
    if (ndim == 1) {
      .X1 <- seq(
        xmin[1],
        xmax[1],
        (xmax[1] - xmin[1]) / 100
      )
      .F <- c()
      for (.x in .X1) {
        .F <- c(.F, fn(.x))
      }
      graphics::plot(
        .X1,
        .F,
        xlim = c(xmin[1], xmax[1]),
        type = "l",
        col = "lightgrey",
        xlab = "x1",
        ylab = "f(x1)"
      )
    } else if (ndim >= 2) {
      graphics::plot(
        2 * xmax[control$.plot_x[1]],
        2 * xmax[control$.plot_x[1]],
        xlim = c(xmin[control$.plot_x[1]], xmax[control$.plot_x[1]]),
        ylim = c(xmin[control$.plot_x[2]], xmax[control$.plot_x[2]]),
        xlab = sprintf("x%d", control$.plot_x[1]),
        ylab = sprintf("x%d", control$.plot_x[2])
      )
    }
  }
  #-----------------------------------------------------------------------------

  ##############################################################################
  # START HERE!
  ##############################################################################
  diversity <- mean_diversity <- c()
  evals <- iter <- 0
  num_exclusions_per_nest <- c()
  num_exclusions <- 0

  repeat {
    iter <- iter + 1
    diversity[iter] <- mean(sqrt(rowSums(t(t(x) - colMeans(x))^2)))
    mean_diversity[iter] <- mean(diversity[1:iter])
    if (plotmethod(control, "mean_diversity")) {
      graphics::plot(
        1:iter,
        mean_diversity,
        type = "b",
        pch = 20,
        xlab = "iters",
        ylab = "mean diversity"
      )
    } else if (plotmethod(control, "diversity")) {
      graphics::plot(
        1:iter,
        diversity,
        type = "b",
        pch = 20,
        xlab = "iters",
        ylab = "diversity"
      )
    }
    evaluate_f()
    update_v()
    if (
      check_for_convergence() || (control$maxiter && iter == control$maxiter)
    ) {
      break
    }
    update_x()
    #-DEBUG---------------------------------------------------------------------
    if (plotmethod(control, "profile")) {
      if (!is.null(nests)) {
        .n <- mynrow(nests)
        .evals <- c(0, nests[, ndim + 4], evals)
        .sols <- c(0:.n, .n)
      } else {
        .evals <- c(0, evals)
        .sols <- c(0, 0)
      }
      graphics::plot(
        .evals,
        .sols,
        type = "b",
        pch = 20,
        xlab = "evals",
        ylab = "sols"
      )
    }
    if (.plot_save_format != "") {
      if (!exists(".plots")) {
        .plots <- 0
      }
      .plots <- .plots + 1
      grDevices::savePlot(sprintf(.plot_save_format, .plots))
    }
    if (control$.plot_method != "" && control$.plot_delay > 0) {
      Sys.sleep(control$.plot_delay)
    }
    #---------------------------------------------------------------------------
  }

  #-DEBUG-----------------------------------------------------------------------
  if (plotmethod(control, "profile")) {
    if (!is.null(nests)) {
      .n <- mynrow(nests)
      .evals <- c(0, nests[, ndim + 4], evals)
      .sols <- c(0:.n, .n)
    } else {
      .evals <- c(0, evals)
      .sols <- c(0, 0)
    }
    graphics::plot(
      .evals,
      .sols,
      type = "b",
      pch = 20,
      xlab = "evals",
      ylab = "sols"
    )
  }
  if (plotswarm(control)) {
    if (ndim == 1) {
      graphics::lines(.X1, .F, col = "lightgrey")
      if (!is.null(nests)) {
        .fnests <- c()
        for (.x1 in nests[, 1:ndim]) {
          .fnests <- c(.fnests, fn(.x1))
        }
        graphics::points(nests[, 1:ndim], .fnests, pch = 4, cex = 2, lwd = 2)
      }
    } else if (ndim == 2) {
      .X1 <- seq(xmin[1], xmax[1], (xmax[1] - xmin[1]) / 50)
      .X2 <- seq(xmin[2], xmax[2], (xmax[2] - xmin[2]) / 50)
      .f <- c()
      for (.x2 in .X2) {
        for (.x1 in .X1) {
          .f <- c(.f, fn(c(.x1, .x2)))
        }
      }
      .F <- matrix(.f, length(.X1), length(.X2))
      graphics::contour(
        .X1,
        .X2,
        .F,
        xlim = c(xmin[1], xmax[1]),
        ylim = c(xmin[2], xmax[2]),
        col = "lightgrey",
        add = TRUE,
        nlevels = 50
      )
      if (!is.null(nests)) {
        graphics::points(
          matrix(nests[, 1:ndim], nrow(nests), ndim),
          pch = 4,
          cex = 2,
          lwd = 2
        )
      }
    }
  }
  if (.plot_save_format != "") {
    if (!exists(".plots")) {
      .plots <- 0
    }
    .plots <- .plots + 1
    grDevices::savePlot(sprintf(.plot_save_format, .plots))
  }
  #-----------------------------------------------------------------------------

  if (!is.null(nests)) {
    colnames(nests) <- c(
      paste(sep = "", "x", 1:ndim),
      "f",
      "v",
      "age",
      "run",
      "evals"
    )
  }
  colnames(pop) <- c(paste(sep = "", "x", 1:ndim), "f", "v", "age")

  invisible(list(
    iter = iter,
    evals = evals,
    nests = nests,
    pop = pop
  ))
}
