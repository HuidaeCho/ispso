# ISPSO
Isolated-Speciation-Based Particle Swarm Optimization <<https://idea.isnew.info/ispso.html>>

This R package implements the Isolated-Speciation-Based Particle Swarm Optimization algorithm published in [Cho, Huidae, Kim, Dongkyun, Olivera, Francisco, Guikema, Seth D., 2011. Enhanced Speciation in Particle Swarm Optimization for Multi-Modal Problems. European Journal of Operational Research 213 (1), 15--23](http://www.sciencedirect.com/science/article/pii/S0377221711001810).

ISPSO is a multi-modal optimization algorithm that aims to discover global and local minima. This algorithm has successfully been used in stochastic rainfall generation, climate change, storm tracking, hydrology, and hydraulics studies.

## Installation

```R
install.packages("remotes")
remotes::install_git("git@github.com:HuidaeCho/ispso.git")

# or if you want to build vignettes
remotes::install_git("git@github.com:HuidaeCho/ispso.git", build_vignettes = TRUE, dependencies = TRUE)

# read the vignette
vignette("ispso")
```

## Testing

```R
library(ispso)
source(system.file("benchmarks", "funcs.R", package = "ispso"))
source(system.file("benchmarks", "benchmark.R", package = "ispso"))
run_benchmark("griewank")
```

![Finding global and local minima in the Griewank function](vignettes/figures/griewank.gif "Finding global and local minima in the Griewank function")

## Acknowledgments

Development of the parallelization features in ISPSO (Parallel ISPSO or PISPSO) was supported in part by the [USGS](https://www.usgs.gov/) Water Resources Research Act 104(b) grant [NM_2023_Cho](https://water.usgs.gov/wrri/grant-details.php?ProjectID=2023NM163B&Type=Annual) through the New Mexico Water Resources Research Institute ([NM WRRI](https://nmwrri.nmsu.edu/)) under award GR0007017, as part of USGS Grant/Cooperative Agreement No. G21AP10635, along with an additional internal award from the NM WRRI. PISPSO is a general-purpose parallel optimization algorithm and was successfully applied to [SWAT+](https://swat.tamu.edu/software/plus/) optimization runs. This algorithm may be applied to other modeling frameworks using independent resources.

## License

Copyright (C) 2008-2026, Huidae Cho <<https://idea.isnew.info/>>

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
