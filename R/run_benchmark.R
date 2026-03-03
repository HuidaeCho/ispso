mynrow <- function(x) {
  if (is.null(x)) {
    0
  } else if (is.vector(x)) {
    1
  } else {
    nrow(x)
  }
}

make_control <- function(parallel = FALSE) {
  control <- ispso_control(
    .plot_method = "movement"
  )

  if (isTRUE(parallel)) {
    control$cluster <- parallel::makeCluster(max(
      1L,
      parallel::detectCores() - 1L
    ))
    control$.owns_cluster <- TRUE
  } else {
    control$.owns_cluster <- FALSE
  }

  control
}

get_problem <- function(func_name) {
  fn <- get(func_name, mode = "function")

  if (func_name %in% c("f1", "f2", "f3", "f4")) {
    bounds <- list(x = c(0, 1))
  } else if (func_name %in% c("f5", "himmelblau")) {
    bounds <- list(x1 = c(-6, 6), x2 = c(-6, 6))
  } else if (func_name %in% c("f6", "rastrigin")) {
    bounds <- list(x1 = c(-1.5, 1.5), x2 = c(-1.5, 1.5))
  } else if (func_name %in% c("f7", "griewank")) {
    bounds <- list(x1 = c(-14, 14), x2 = c(-14, 14))
  } else if (func_name == "bimodal") {
    bounds <- list(x = c(-4, 8))
  } else if (func_name == "rosenbrock") {
    bounds <- list(x1 = c(-10, 10), x2 = c(-10, 10))
  } else if (func_name == "ackley") {
    bounds <- list(x1 = c(-32.768, 32.768), x2 = c(-32.768, 32.768))
  } else if (func_name == "levy5") {
    bounds <- list(x1 = c(-10, 10), x2 = c(-10, 10))
  } else if (func_name == "spherical") {
    bounds <- list(x1 = c(-5.12, 5.12), x2 = c(-5.12, 5.12))
  } else if (func_name == "quadric") {
    bounds <- list(x1 = c(-30, 30), x2 = c(-30, 30))
  } else {
    stop("Unknown func_name: ", func_name)
  }

  list(fn = fn, bounds = bounds)
}

run_benchmark <- function(func_name, parallel = FALSE) {
  prob <- get_problem(func_name)
  fn <- prob$fn
  bounds <- prob$bounds
  ndim <- length(bounds)

  control <- make_control(parallel = parallel)
  on.exit(
    {
      if (isTRUE(control$.owns_cluster) && !is.null(control$cluster)) {
        parallel::stopCluster(control$cluster)
      }
    },
    add = TRUE
  )

  # per-problem overrides (e.g., maxiter for griewank/f7)
  if (func_name %in% c("f7", "griewank")) {
    control$maxiter <- 3000L
  }

  ret <- ispso(fn, bounds, control = control)

  # True solutions (your function API expects list(sol=TRUE))
  sol <- fn(rep(0, ndim), list(sol = TRUE))
  nsols <- mynrow(sol)
  nnests <- mynrow(ret$nest)

  success <- FALSE
  if (nnests) {
    # derive diagonal span from bounds for the success check
    lower <- vapply(bounds, function(x) x[1L], numeric(1))
    upper <- vapply(bounds, function(x) x[2L], numeric(1))
    diag_span <- sqrt(sum((upper - lower)^2))

    rsol <- diag_span * control$.distance_to_solution

    too_far <- rep(0, nnests)
    error <- c()
    ndupnests <- rep(0, nsols)

    for (i in 1:nsols) {
      fmin <- fn(sol[i, ])
      for (j in 1:nnests) {
        d <- as.matrix(stats::dist(
          if (ndim == 1) {
            c(sol[i, 1], ret$nest[j, 1])
          } else {
            rbind(sol[i, ], ret$nest[j, 1:ndim])
          }
        ))[1, 2]

        if (d > rsol) {
          too_far[j] <- too_far[j] + 1
          next
        }

        error[j] <- abs(fmin - ret$nest[j, ndim + 1])
        ndupnests[i] <- ndupnests[i] + 1
      }
    }

    nest <- ret$nest[too_far < nsols, ]
    error <- error[too_far < nsols]
    max_error <- max(error)

    if (
      mynrow(ret$nest) == nsols &&
        sum(ndupnests == rep(1, nsols)) == nsols
    ) {
      cat(sprintf("\b\b\b\bGREAT!\n"))
      success <- TRUE
    }
  } else {
    nest <- ret$nest
    max_error <- NA
  }

  invisible(list(
    fn = fn,
    bounds = bounds,
    control = control,
    ret = ret,
    success = success,
    nest = nest,
    max_error = max_error
  ))
}
