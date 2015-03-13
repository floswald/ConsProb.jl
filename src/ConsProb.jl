

module ConsProb

# using ApproXD
using PyPlot, Roots, Optim
using FastGaussQuadrature: gausshermite

include("setup.jl")
include("funs.jl")

export Param, Model, Model2

end