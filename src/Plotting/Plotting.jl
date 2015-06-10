module Plotting

using PyPlot

	import ..Models:
		iidModel,Param,
		iidDModel,AR1Model,AR1Model_a,
	    cond,condvbound,
	    env, envvbound,
	    get_vbound,u

	import ..Standard:
		runStd

	import ..Dchoice:
		runDchoice




include("plots.jl")

	export plots

end