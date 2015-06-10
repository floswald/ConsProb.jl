


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

include("d-funs.jl")

export runDchoice

end 