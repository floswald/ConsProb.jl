module Plotting

using PyPlot

	import ..Models:
		iidModel,Param,
		iidDModel,AR1Model,AR1Model_a,
	    cond,condvbound,
	    env, envvbound,
	    get_vbound,u




include("plots.jl")

	export plots

end