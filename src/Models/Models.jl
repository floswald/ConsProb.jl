

module Models

	using FastGaussQuadrature: gausshermite
	using Distributions: Normal, quantile

	include("envelope.jl")
	include("param.jl")
	include("utils.jl")
	include("model-defs.jl")

	# export Types
	export Param, 
		   Model, 
		   iidModel, 
		   iidDModel, 
		   iidDebtModel, 
		   AR1Model, 
		   AR1Model_a,
		   Envelope

    # export methods
    export linearapprox,
		   u,up,iup,
		   cond,condvbound,
		   env, envvbound,
		   set!, set_vbound!, get_vbound

end