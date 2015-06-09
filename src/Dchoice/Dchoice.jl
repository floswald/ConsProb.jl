


module Dchoice

using Roots: fzero

import ..Models:

	Param,
	iidDModel,
	linearapprox,u,up,iup,
	Envelope,
    cond,condvbound,
    env, envvbound,
    set!, set_vbound!,get_vbound

import ..Plotting:

	plots

include("d-funs.jl")

export run

end 