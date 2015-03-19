

module ConsProb


using PyPlot, Roots, Optim
using FastGaussQuadrature: gausshermite

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
include("plotting.jl")

# export Param, Model

end