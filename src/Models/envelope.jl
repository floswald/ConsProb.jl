
"""
# Envelope type

* values conditional on discrete choice
* Envelopes over values conditional on discrete choice
* values of functions on lower boundary conditional on discrete choice
* Value of lower boundary on Envelope

"""
type Envelope
	cond        :: Dict{Int,Vector{Float64}} 	# dict of functions conditional on discrete choice
	cond_vbound :: Dict{Int,Float64} 	# dict of expected values on savings bound by discrete choice
	env         :: Vector{Float64}		# actual envelope over discrete choices
	env_vbound  :: Float64 				# expected value at savings bound on envelope
	function Envelope()
		this = new()
		this.cond = [id => Float64[] for id in 1:2 ]
		this.cond_vbound = [id => 0.0 for id in 1:2 ]
		this.env = Float64[]
		this.env_vbound = 0.0
		return this
	end
	function Envelope(cond,cond_vbound,env,vbound)
		this = new()
		this.cond = copy(cond)
		this.cond_vbound = copy(cond_vbound)
		this.env = copy(env)
		this.env_vbound = copy(vbound)
		return this
	end
end

function cond(e::Envelope,which::Int)
	e.cond[which]
end
function condvbound(e::Envelope,which::Int)
	vcat(e.cond_vbound[which],e.cond[which])
end
function env(e::Envelope)
	e.env
end
function envvbound(e::Envelope)
	vcat(e.env_vbound,e.env)
end


function cond(e::Dict{Int,Envelope},it::Int,which::Int)
	cond(e[it],which)
end
function condvbound(e::Dict{Int,Envelope},it::Int,which::Int)
	condvbound(e[it],which)
end
function env(e::Dict{Int,Envelope},it::Int)
	env(e[it])
end
function envvbound(e::Dict{Int,Envelope},it::Int)
	envvbound(e[it])
end
function set!(e::Dict{Int,Envelope},it::Int,v::Vector{Float64})
	e[it].env = deepcopy(v)
end
function set!(e::Dict{Int,Envelope},it::Int,which::Int,v::Vector{Float64})
	e[it].cond[which] = deepcopy(v)
end
function set_vbound!(e::Dict{Int,Envelope},it::Int,which::Int,v::Float64)
	e[it].cond_vbound[which] = v
end
function set_vbound!(e::Dict{Int,Envelope},it::Int,v::Float64)
	e[it].env_vbound = v
end
function get_vbound(e::Dict{Int,Envelope},it::Int,which::Int)
	e[it].cond_vbound[which]
end
function get_vbound(e::Dict{Int,Envelope},it::Int)
	e[it].env_vbound 
end


"""
Value Function type
"""
type Vfun
	v         :: Vector{Float64}		# values on interior of value function domain
	vbound    :: Float64 				# value on the lower function boundary
end
v(x::Vfun) = x.v
v0(x::Vfun) = vcat(vbound,v)


