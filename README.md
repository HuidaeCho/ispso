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

Copyright (C) 2008, Huidae Cho <<https://idea.isnew.info>>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <<http://www.gnu.org/licenses/>>.
