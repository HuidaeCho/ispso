# Purpose: Provides an input file for various tests.
# Requires: ispso.R, funcs.R
# Variables: s, ret, res

source("ispso.R")
source("funcs.R")

diagonal <- function(s) sqrt(sum((s$xmax-s$xmin)^2))

set_parameters <- function(s){
	########################################################################
	# DEBUG
	# Deterministic run?
	s$.deterministic <- FALSE
	#s$.deterministic <- TRUE

	# Stop if all the solutions are found!  This is only for writing a
	# paper, not for real problems because the number of actual solutions
	# is not known in most cases.
	s$.stop_after_solutions <- 0
	#s$.stop_after_solutions <- -1

	# (0, 1]: Fraction of the diagonal span of the search space.
	s$.distance_to_solution <- 0.01

	# Plot method
	s$.plot_method <- "density"
	s$.plot_method <- "movement"
	#s$.plot_method <- sprintf("%s,species", s$.plot_method)
	s$.plot_delay <- 0
	########################################################################

	# Swarm size
	s$S <- 10 + floor(2*sqrt(s$D))

	# Maximum particle velocity
	s$vmax <- (s$xmax-s$xmin)*0.1

	# Maximum initial particle velocity
	s$vmax0 <- diagonal(s)*0.001

	# Stopping criteria: Stop if the number of exclusions per particle
	# since the last minimum is greater than exclusion_factor * max sol
	# iter / average sol iter. The more difficult the problem is (i.e.,
	# high max sol iter / average sol iter), the more iterations the
	# algoritm requires to stop.
	s$exclusion_factor <- 3
	# Maximum iteration
	s$maxiter <- 2000
	# Small positive number close to 0
	s$xeps <- 0.001
	s$feps <- 0.0001

	# Search radius for preys: One particle has two memories (i.e., x and
	# pbest).  When two particles collide with each other within prey, one
	# particle takes more desirable x and pbest from the two particles'
	# memories, and the other particle is replaced with a quasi-random
	# particle using scrambled Sobol' sequences (PREY).
	s$rprey <- diagonal(s)*0.0001

	# Nesting criteria for global and local optima using particles' ages
	# (NEST_BY_AGE).
	s$age <- 10

	# Speciation radius: Li (2004) recommends 0.05*L<=rspecies<=0.1*L.
	s$rspecies <- diagonal(s)*0.1

	# Nesting radius
	s$rnest <- diagonal(s)*0.01

	invisible(s)
}

test <- function(func.name){
	s <- list()
	if(func.name == "f1"){
		s$f <- eval(parse(text=func.name))
		s$D <- 1
		s$xmin <- 0
		s$xmax <- 1

		s <- set_parameters(s)

		ret <- ispso(s)
	}else
	if(func.name == "f2"){
		s$f <- eval(parse(text=func.name))
		s$D <- 1
		s$xmin <- 0
		s$xmax <- 1

		s <- set_parameters(s)

		ret <- ispso(s)
	}else
	if(func.name == "f3"){
		s$f <- eval(parse(text=func.name))
		s$D <- 1
		s$xmin <- 0
		s$xmax <- 1

		s <- set_parameters(s)

		ret <- ispso(s)
	}else
	if(func.name == "f4"){
		s$f <- eval(parse(text=func.name))
		s$D <- 1
		s$xmin <- 0
		s$xmax <- 1

		s <- set_parameters(s)

		ret <- ispso(s)
	}else
	if(func.name == "f5" || func.name == "himmelblau"){
		s$f <- eval(parse(text=func.name))
		s$D <- 2
		s$xmin <- rep(-6, s$D)
		s$xmax <- rep(6, s$D)

		s <- set_parameters(s)

		ret <- ispso(s)
	}else
	if(func.name == "f6" || func.name == "rastrigin"){
		s$f <- eval(parse(text=func.name))
		s$D <- 2
		s$xmin <- rep(-1.5, s$D)
		s$xmax <- rep(1.5, s$D)

		s <- set_parameters(s)

		ret <- ispso(s)
	}else
	if(func.name == "f7" || func.name == "griewank"){
		s$f <- eval(parse(text=func.name))
		s$D <- 2
		s$xmin <- rep(-14, s$D)
		s$xmax <- rep(14, s$D)

		s <- set_parameters(s)
		s$maxiter <- 3000

		ret <- ispso(s)
	}else
	if(func.name == "bimodal"){
		s$f <- eval(parse(text=func.name))
		s$D <- 1
		s$xmin <- -4
		s$xmax <- 8

		s <- set_parameters(s)

		ret <- ispso(s)
	}else
	if(func.name == "rosenbrock"){
		s$f <- eval(parse(text=func.name))
		s$D <- 2
		s$xmin <- rep(-10, s$D)
		s$xmax <- rep(10, s$D)

		s <- set_parameters(s)

		ret <- ispso(s)
	}else
	if(func.name == "ackley"){
		s$f <- eval(parse(text=func.name))
		s$D <- 2
		s$xmin <- rep(-32.768, s$D)
		s$xmax <- rep(32.768, s$D)

		s <- set_parameters(s)

		ret <- ispso(s)
	}else
	if(func.name == "levy5"){
		s$f <- eval(parse(text=func.name))
		s$D <- 2
		s$xmin <- rep(-10, s$D)
		s$xmax <- rep(10, s$D)

		s <- set_parameters(s)

		ret <- ispso(s)
	}else
	if(func.name == "spherical"){
		s$f <- eval(parse(text=func.name))
		s$D <- 2
		s$xmin <- rep(-5.12, s$D)
		s$xmax <- rep(5.12, s$D)

		s <- set_parameters(s)

		ret <- ispso(s)
	}else
	if(func.name == "quadric"){
		s$f <- eval(parse(text=func.name))
		s$D <- 2
		s$xmin <- rep(-30, s$D)
		s$xmax <- rep(30, s$D)

		s <- set_parameters(s)

		ret <- ispso(s)
	}

	res <- c()
	# True solutions
	sol <- s$f(rep(0, s$D), TRUE)
	# Number of solutions
	nsols <- mynrow(sol)
	# Number of nests
	nnests <- mynrow(ret$nest)
	# Success? FALSE by default
	success <- FALSE
	if(nnests){
		# Distance threshold
		rsol <- diagonal(s)*s$.distance_to_solution
		too_far <- rep(0, nnests)
		error <- c()
		# Number of duplicate nests
		ndupnests <- rep(0, nsols)
		for(i in 1:nsols){
			fmin <- s$f(sol[i,])
			for(j in 1:nnests){
				d <- as.matrix(dist(if(s$D == 1)
					c(sol[i, 1], ret$nest[j, 1]) else
					rbind(sol[i,],
					ret$nest[j, 1:s$D])))[1, 2]
				if(d > rsol){
					too_far[j] <- too_far[j] + 1
					next
				}
				error[j] <- abs(fmin-ret$nest[j, s$D+1])
				ndupnests[i] <- ndupnests[i] + 1
			}
		}
		# Good nests. There can be duplicate nests.
		nest <- ret$nest[too_far<nsols,]
		# Errors
		error <- error[too_far<nsols]
		max.error <- max(error)

		if(mynrow(ret$nest) == nsols &&
		   sum(ndupnests == rep(1, nsols)) == nsols){
			printf("\b\b\b\bGREAT!\n")
			success <- TRUE
		}
	}else{
		nest <- ret$nest
		max.error <- NA
	}

	invisible(list(
		s=s,
		ret=ret,
		success=success,
		nest=nest,
		max.error=max.error
	))
}
