

module ConsProb


# setting up docile documentation
if VERSION < v"0.4.0-dev"
    using Docile
end


"""
Consumption Problems

This module has some different ways to solve a standard life-cycle consumption problem.
"""
ConsProb

include(joinpath("Models","Models.jl"))
include(joinpath("Standard","Standard.jl"))
include(joinpath("Dchoice","Dchoice.jl"))
include(joinpath("Plotting","Plotting.jl"))

export Dchoice, Standard, Models, Plotting

end