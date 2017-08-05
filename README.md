# ISPSO
Isolated-Speciation-based Particle Swarm Optimization <<https://idea.isnew.info/ispso.html>>

ispso.R implements the Isolated-Speciation-based Particle Swarm Optimization algorithm published in [Cho, Huidae, Kim, Dongkyun, Olivera, Francisco, Guikema, Seth D., 2011. Enhanced Speciation in Particle Swarm Optimization for Multi-Modal Problems. European Journal of Operational Research 213 (1), 15--23](http://www.sciencedirect.com/science/article/pii/S0377221711001810).

ISPSO is a multi-modal optimization algorithm that aims to discover global and local minima. This algorithm has successfully been used in climate change, storm tracking, hydrology, and hydraulics studies.

## Griewank Function

```
source("ispso.R")
source("funcs.R")

s <- list()
s$f <- griewank
s$D <- 2
s$xmin <- rep(-14, s$D)
s$xmax <- rep(14, s$D)
s$S <- 10 + floor(2*sqrt(s$D))
s$vmax <- (s$xmax-s$xmin)*0.1
s$vmax0 <- diagonal(s)*0.001
s$maxiter <- 2000
s$xeps <- 0.001
s$feps <- 0.0001
s$rprey <- diagonal(s)*0.0001
s$age <- 10
s$rspecies <- diagonal(s)*0.1
s$rnest <- diagonal(s)*0.01
s$.plot_distance_to_solution <- 0.01

ret <- ispso(s)
```

![Finding global and local minima in the Griewank function](griewank.gif "Finding global and local minima in the Griewank function")

## License

GNU General Public License Version 3
