

module ConsProb


using PyPlot, Roots, Optim, Debug
using FastGaussQuadrature: gausshermite
using Distributions: Normal, quantile

# setting up docile documentation
if VERSION < v"0.4.0-dev"
    using Docile
end


@document

"
Consumption Problems

This module has some different ways to solve a standard life-cycle consumption problem.
"

include("setup.jl")
include("funs.jl")
include("dchoice.jl")
include("plotting.jl")

export runall,
	   Model,
	   Param,
	   iidModel,
	   Euler!,
	   VFbi!,
	   EGM!,
	   plots,
	   solve!


end