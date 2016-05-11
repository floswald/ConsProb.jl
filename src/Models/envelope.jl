
# TODO: should this be specialized to hold n elements?
"""
# Bounded Function type

## Fields

* `v`: array of values on interior of domain
* `vbound`: value on lower boundary (at borrowing constraint)
"""
type Bfun
	v        :: Vector{Float64}		# values on interior of value function domain
	bound    :: Float64 				# value on the lower function boundary
end
v(x::Bfun) = x.v
b(x::Bfun) = x.bound
vb(x::Bfun) = vcat(x.bound,x.v)
function set!(x::Bfun,y::Vector{Float64}) 
	x.v = y
end
"Set all Bfuns equal to vector y"
function set!(x::Array{Bfun},y::Vector{Float64}) 
	for i in x
		set!(i,y)
	end
end
function set!(x::Bfun,y::LinSpace) 
	x.v = collect(y)
end
function set_bound!(x::Bfun,y::Float64) 
	x.bound = y
end
function set_bound!(x::Array{Bfun},y::Float64) 
	for i in x
		i.bound = y
	end
end



"""
# Envelope type

* A collection of `Bfun` types
* values conditional on discrete choice
* Envelopes over values conditional on discrete choice
* values of functions on lower boundary conditional on discrete choice
* Value of lower boundary on Envelope

"""
type Envelope
	cond        :: Dict{Int,Bfun} 	# dict of functions conditional on discrete choice
	env         :: Bfun				# actual envelope over discrete choices
	function Envelope()
		this = new()
		this.cond = [id => Bfun(Float64[],NaN) for id in 1:2 ]
		this.env = Bfun(Float64[],NaN)
		return this
	end
	function Envelope(cond,env)
		this = new()
		this.cond = copy(cond)
		this.env = copy(env)
		return this
	end
end

# function cond(e::Envelope,which::Int)
# 	e.cond[which]
# end
# function condvbound(e::Envelope,which::Int)
# 	vcat(e.cond_vbound[which],e.cond[which])
# end
# function env(e::Envelope)
# 	e.env
# end
# function envvbound(e::Envelope)
# 	vcat(e.env_vbound,e.env)
# end


# function cond(e::Dict{Int,Envelope},it::Int,which::Int)
# 	cond(e[it],which)
# end
# function condvbound(e::Dict{Int,Envelope},it::Int,which::Int)
# 	condvbound(e[it],which)
# end
# function env(e::Dict{Int,Envelope},it::Int)
# 	env(e[it])
# end
# function envvbound(e::Dict{Int,Envelope},it::Int)
# 	envvbound(e[it])
# end
# function set!(e::Dict{Int,Envelope},it::Int,v::Vector{Float64})
# 	e[it].env = deepcopy(v)
# end
# function set!(e::Dict{Int,Envelope},it::Int,which::Int,v::Vector{Float64})
# 	e[it].cond[which] = deepcopy(v)
# end
# function set_vbound!(e::Dict{Int,Envelope},it::Int,which::Int,v::Float64)
# 	e[it].cond_vbound[which] = v
# end
# function set_vbound!(e::Dict{Int,Envelope},it::Int,v::Float64)
# 	e[it].env_vbound = v
# end
# function get_vbound(e::Dict{Int,Envelope},it::Int,which::Int)
# 	e[it].cond_vbound[which]
# end
# function get_vbound(e::Dict{Int,Envelope},it::Int)
# 	e[it].env_vbound 
# end




