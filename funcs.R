# Purpose: Defines test functions.
# Functions: f1, f2, f3, f4, himmelblau, rastrigin, griewank, bimodal,
# rosenbrock, ackley, levy5, spherical, quadric

f1 <- function(x, core=1, sol=FALSE){
	if(sol){
		# x = [0, 1]
		return(matrix(c(
			0.1,
			0.3,
			0.5,
			0.7,
			0.9
		), 5, 1, byrow=TRUE))
	}

	1 - sin(5*pi*x)^6
}

f2 <- function(x, core=1, sol=FALSE){
	if(sol){
		# x = [0, 1]
		return(matrix(c(
			0.1,
			0.299,
			0.499,
			0.698,
			0.898
		), 5, 1, byrow=TRUE))
	}

	1 - exp(-2*log(2)*((x-0.1)/0.8)^2)*sin(5*pi*x)^6
}

f3 <- function(x, core=1, sol=FALSE){
	if(sol){
		# x = [0, 1]
		return(matrix(c(
			0.08,
			0.247,
			0.451,
			0.681,
			0.934
		), 5, 1, byrow=TRUE))
	}

	1 - sin(5*pi*(x^(3/4)-0.05))^6
}

f4 <- function(x, core=1, sol=FALSE){
	if(sol){
		# x = [0, 1]
		return(matrix(c(
			0.08,
			0.246,
			0.449,
			0.679,
			0.93
		), 5, 1, byrow=TRUE))
	}

	1 - exp(-2*log(2)*((x-0.08)/0.854)^2)*sin(5*pi*(x^(3/4)-0.05))^6
}

f5 <- himmelblau <- function(x, core=1, sol=FALSE){
	if(sol){
		# x = [-6, 6]^2
		return(matrix(c(
			3.0, 2.0,
			-3.78, -3.28,
			3.58, -1.86,
			-2.815, 3.125
		), 4, 2, byrow=TRUE))
	}

	if(is.matrix(x))
		(x[,1]^2+x[,2]-11)^2+(x[,1]+x[,2]^2-7)^2
	else
		(x[1]^2+x[2]-11)^2+(x[1]+x[2]^2-7)^2
}

f6 <- rastrigin <- function(x, core=1, sol=FALSE){
	if(is.matrix(x))
		return(apply(x, 1, rastrigin))

	n <- length(x)
	if(sol){
		# x = [-1.5, 1.5]^n
		return(rastrigin.sols(rep(-1.5, n), rep(1.5, n)))
	}

	10*n+sum(x^2-10*cos(2*pi*x))
}
f6.sols <- rastrigin.sols <- function(xmin, xmax){
	n <- length(xmin)
	xmini <- ceiling(xmin)
	xmaxi <- floor(xmax)
	if(n == 1){
		return(matrix(xmini[n]:xmaxi[n], xmaxi[n]-xmini[n]+1, 1,
			byrow=TRUE))
	}

	x <- c()
	x1 <- rastrigin.sols(xmin[1:(n-1)], xmax[1:(n-1)])
	for(i in 1:nrow(x1)){
		for(x2 in xmini[n]:xmaxi[n])
			x <- rbind(x, c(x1[i,], x2))
	}
	return(x)
}

f7 <- griewank <- function(x, core=1, sol=FALSE){
	if(is.matrix(x))
		return(apply(x, 1, griewank))

	n <- length(x)
	if(sol){
		# x = [-14, 14]^n
		return(griewank.sols(rep(-14, n), rep(14, n)))
	}

	sum(x^2)/4000-prod(cos(x/sqrt(1:n)))+1
}
f7.sols <- griewank.sols <- function(xmin, xmax, odd=FALSE){
	n <- length(xmin)
	xmini <- ceiling(xmin[1:n]/(pi*sqrt(1:n)))
	xmaxi <- floor(xmax[1:n]/(pi*sqrt(1:n)))
	if(n == 1){
		x0 <- xmini[n] + if(odd) !mod(xmini[n], 2) else mod(xmini[n], 2)
		x <- pi*seq(x0, xmaxi[n], 2)
		return(matrix(x, length(x), 1, byrow=TRUE))
	}

	x <- c()
	x1 <- griewank.sols(xmin[1:(n-1)], xmax[1:(n-1)])
	x0 <- xmini[n] + if(odd) !mod(xmini[n], 2) else mod(xmini[n], 2)
	for(i in 1:nrow(x1)){
		for(x2 in pi*sqrt(n)*seq(x0, xmaxi[n], 2))
			x <- rbind(x, c(x1[i,], x2))
	}

	x1 <- griewank.sols(xmin[1:(n-1)], xmax[1:(n-1)], TRUE)
	x0 <- xmini[n] + if(odd) mod(xmini[n], 2) else !mod(xmini[n], 2)
	for(i in 1:nrow(x1)){
		for(x2 in pi*sqrt(n)*seq(x0, xmaxi[n], 2))
			x <- rbind(x, c(x1[i,], x2))
	}

	return(x)
}
f7.num_sols <- griewank.num_sols <- function(d, xmin, xmax){
	if(length(d) > 1){
		ret <- c()
		for(i in 1:length(d))
			ret[i] <- griewank.num_sols(d[i], xmin, xmax)
		return(ret)
	}

	num_evens <- function(d) floor(xmax[d]/(2*pi*sqrt(d)))-
		ceiling(xmin[d]/(2*pi*sqrt(d)))+1
	num_odds <- function(d) floor(xmax[d]/(2*pi*sqrt(d))+0.5)-
		ceiling(xmin[d]/(2*pi*sqrt(d))-0.5)
	num_modes <- function(d) prod(num_evens(1:d)+num_odds(1:d))

	if(d == 1){
		num_evens(d)
	}else{
		num_1 <- griewank.num_sols(d-1, xmin, xmax)
		num_1*num_evens(d)+(num_modes(d-1)-num_1)*num_odds(d)
	}
}

bimodal <- function(x, core=1, sol=FALSE){
	if(sol){
		# x = [-4, 8]
		return(matrix(c(
			0, 4
		), 2, 1, byrow=TRUE))
	}

	1 - (1/(sqrt(2*pi))*exp(-0.5*x^2) + 2/(sqrt(2*pi))*exp(-0.5*(2*x-8)^2))
}

rosenbrock <- function(x, core=1, sol=FALSE){
	n <- length(x)
	if(sol){
		# x = [-10, 10]^n
		# f(x) = 0
		return(matrix(rep(1, n), 1, n))
	}

	ret <- 0
	for(i in 1:(n-1))
		ret <- ret+100*(x[i+1]-x[i]^2)^2+(x[i]-1)^2
	ret
}


ackley <- function(x, core=1, sol=FALSE){
	n <- length(x)
	if(sol){
		# x = [-32.768, 32.768]^n
		# f(x) = 0
		return(matrix(rep(0, n), 1, n))
	}

	-20*exp(-0.2*sqrt(sum(x^2)/n))-exp(sum(cos(2*pi*x))/n)+20+exp(1)
}

levy5 <- function(x, core=1, sol=FALSE){
	if(sol){
		# x = [-10, 10]^n
		# 760 local minima
		# 1 global minimum f(x) = -176.1375 at x = (-1.3068, -1.4248)
		return(matrix(c(
			-1.3068, -1.4248
		), 1, 2, byrow=TRUE))
	}

	i <- 1:5
	sum(i*cos((i-1)*x[1]+i))*sum(i*cos((i+1)*x[2]+i))+
		(x[1]+1.42513)^2+(x[2]+0.80032)^2
}

spherical <- function(x, core=1, sol=FALSE){
	n <- length(x)
	if(sol){
		# x = [-5.12, 5.12]^n
		# one global minimum f(x) = 0 at x = (0, ..., 0)
		return(matrix(rep(0, n), 1, n))
	}

	sum(x^2)
}

quadric <- function(x, core=1, sol=FALSE){
	n <- length(x)
	if(sol){
		# f(x) = 0
		return(matrix(rep(0, n), 1, n))
	}

	ret <- 0
	for(i in 1:n)
		ret <- ret+sum(x[1:i])^2
	ret
}
