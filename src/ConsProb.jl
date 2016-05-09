

VERSION >= v"0.4.0-dev+6521" && __precompile__(true)

module ConsProb

	# imports
	using FastGaussQuadrature: gausshermite
	using Distributions: Normal, quantile
	using PyPlot, Roots, Optim

	# export types
	export Param,
		   Model, 
		   iidModel, 
		   iidDModel, 
		   iidDebtModel, 
		   AR1Model, 
		   AR1Model_a,
		   Envelope,
		   Bfun


	# export methods
	export plots,
			EGM!,
			Euler!,
			VFbi!,
			v,
			vb

	# load files
	include(joinpath("Models","envelope.jl"))
	include(joinpath("Models","model-defs.jl"))
	include(joinpath("Models","utils.jl"))
	include(joinpath("Standard","solutions.jl"))
	include(joinpath("Dchoice","d-funs.jl"))
	include(joinpath("plots.jl"))

end