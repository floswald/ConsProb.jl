


@doc """
Holds the user-set parameter values. 

**values to be set:**

* `gamma`: CRRA 
* `beta`: discount factor
* `R`: gross interest rate (i.e. 1+r )
* `na`: number of grid points for assets
* `ny`: number of grid points for income
* `nT`: number of time periods
* `a_low`: lower bound on assets
* `a_high`: upper bound on assets
* `mu`: unconditional mean of income (iid case)
* `sigma`: unconditional variance of income (iid case)
* `rho_z`: AR1 coefficient of income (AR1 case)
* `eps_z`: standard deviation of AR1 innovation (AR1 case)
* `dorefinements`: boolean switch of whether to filter out whiggles
* `alpha`: disutility of work

""" ->
type Param

	# CRRA
	gamma                 :: Float64
	neg_gamma             :: Float64
	oneminusgamma         :: Float64
	oneover_oneminusgamma :: Float64
	neg_oneover_gamma     :: Float64

	beta::Float64
	R::Float64

	# grids
	na :: Int 	# asset grid
	ny :: Int   # income grid (support points)
	nT :: Int   # maximal age
	a_high::Float64
	a_low::Float64
	nD:: Int # number of discrete choices

	# consumption floor
	cfloor :: Float64

	# disutility of work
	alpha :: Float64

	# iid income unertainty
	mu::Float64
	sigma::Float64

	# AR1 income uncertainty
	# y_t = z_t
	# ln z_t = rho_z * ln z_{t-1} + eps_z_{t}
	rho_z::Float64
	eps_z::Float64

	dorefinements::Bool
	printdebug::Bool

	# constructor for discrete choice model: log utility
	function Param(gamm::Float64)

		gamma                 = gamm
		neg_gamma             = (-1.0) * gamma
		oneminusgamma         = 1.0 - gamma
		oneover_oneminusgamma = 1.0 / oneminusgamma
		neg_oneover_gamma     = (-1.0) / gamma

		beta                  = 0.95
		R                     = 1.05

		na     = 100
		ny     = 50
		nT     = 6
		a_high = 50.0
		a_low  = 0.0
		nD     = 2

		cfloor = 0.001
		alpha = 2.3

		# iid income uncertainty params
		# mu = 10 	# mean income: 30K
		# sigma = 1  # sd income
		mu = 0 	# mean income: 30K
		sigma = 0.25  # sd income

		# AR1 income uncertainty
		# params from Ayiagari
		rho_z = 0.9
		eps_z = 1

		dorefine=true
		printdebug=false

		return new(gamma,neg_gamma,oneminusgamma,oneover_oneminusgamma,neg_oneover_gamma,beta,R,na,ny,nT,a_high,a_low,nD,cfloor,alpha,mu,sigma,rho_z,eps_z,dorefine,printdebug)
	end

	function Param(;mu=1)


		gamma                 = 2.0
		neg_gamma             = (-1.0) * gamma
		oneminusgamma         = 1.0 - gamma
		oneover_oneminusgamma = 1.0 / oneminusgamma
		neg_oneover_gamma     = (-1.0) / gamma

		beta                  = 0.95
		R                     = 1.05

		na = 200
		ny = 10
		nT = 8
		a_high = 300.0
		a_low  = 1e-6
		nD = 2

		cfloor = 0.001
		alpha = 1.0

		# iid income uncertainty params
		 # = 1 	# mean income: 30K
		sigma = 0.1  # sd income

		# AR1 income uncertainty
		# params from Ayiagari
		rho_z = 0.9
		eps_z = 1

		dorefine=false

		return new(gamma,neg_gamma,oneminusgamma,oneover_oneminusgamma,neg_oneover_gamma,beta,R,na,ny,nT,a_high,a_low,nD,cfloor,alpha,mu,sigma,rho_z,eps_z,dorefine)
	end
end

# Abstract Model type
# this is an overall model type, with potential subtypes

@doc """
Abstract Model type

**There are three model types:**

1. model with iid income uncertainty, use cash-on-hand m=y+a as state variable: V(a)
2. model with AR1 income uncertainty, use cash-on-hand m=y+a and y as state variables: V(m,y)
3. model with AR1 income uncertainty, use current assets a and y as state variables: V(a,y)
""" ->
abstract Model


# there are 3 different models:



@doc """
Model with iid income uncertainty after Deaton (1991)

uses cash-on-hand m=y+a as state variable
""" ->
type iidModel <: Model

	# computation grids
	avec::Vector{Float64}
	yvec::Vector{Float64}   # income support
	ywgt::Vector{Float64}   # income weights

	# intermediate objects (na,ny)
	m1::Array{Float64,2}	# matrix (na,ny)
	c1::Array{Float64,2}
	ev::Array{Float64,2}
	# intermediate objects (ny,1)
	m2::Vector{Float64} 
	c2::Vector{Float64} 

	# result objects
	C::Array{Float64,2} 	# consumption function on (na,nT)
	S::Array{Float64,2} 	# savings function on (na,nT)
	M::Array{Float64,2} 	# endogenous cash on hand on (na,nT)
	V::Array{Float64,2} 	# value function on (na,nT). Optional.
	Vzero::Array{Float64,1} 	# value function of saving zero

	toc::Float64   # elapsed time

	solver::ASCIIString

	@doc "Constructor for iid Model" ->
	function iidModel(p::Param,solver::String="EGM")

		avec          = scaleGrid(p.a_low,p.a_high,p.na,2)
		nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature

		# for y ~ N(mu,sigma), integrate y with 
		# http://en.wikipedia.org/wiki/Gauss-Hermite_quadrature
		yvec = sqrt(2.0) * p.sigma .* nodes .+ p.mu
		ywgt = weights .* pi^(-0.5)

		# precompute next period's cash on hand.
		m1 = Float64[avec[ia]*p.R + yvec[iy] for ia in 1:p.na, iy in 1:p.ny]

		# if you want a deterministic age profile in income, use income().
		# you would have to change the params of income() though.
		# m1 = Float64[avec[ia]*p.R + income(yvec[iy],it) for ia in 1:p.na, iy in 1:p.ny, it in 1:p.nT]

		c1 = zeros(p.na,p.ny)
		ev = zeros(p.na,p.ny)

		m2 = zeros(p.ny)
		c2 = zeros(p.ny)

		C = zeros(p.na,p.nT)
		S = zeros(p.na,p.nT)
		M = zeros(p.na,p.nT)
		V = zeros(p.na,p.nT)
		Vzero = zeros(p.nT)

		toc = 0.0


		return new(avec,yvec,ywgt,m1,c1,ev,m2,c2,C,S,M,V,Vzero,toc,solver)
	end
end

@doc """
Model with iid income uncertainty and unsecured debt

uses cash-on-hand m=y+a as state variable
""" ->
type iidDebtModel <: Model

	# computation grids
	bounds::Vector{Float64}
	avec::Matrix{Float64}
	yvec::Vector{Float64}   # income support
	ywgt::Vector{Float64}   # income weights

	# intermediate objects (na,ny)
	m1::Array{Float64,3}	# matrix (na,ny,nt)
	c1::Array{Float64,2}
	ev::Array{Float64,2}
	# intermediate objects (ny,1)
	m2::Vector{Float64} 
	c2::Vector{Float64} 

	# result objects
	C::Array{Float64,2} 	# consumption function on (na,nT)
	S::Array{Float64,2} 	# savings function on (na,nT)
	M::Array{Float64,2} 	# endogenous cash on hand on (na,nT)
	V::Array{Float64,2} 	# value function on (na,nT). Optional.
	Vzero::Array{Float64,1} 	# value function of saving zero
	
	dont::Array{Bool,3}	

	@doc "Constructor for iid Debt Model" ->
	function iidDebtModel(p::Param)

		nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature

		# for y ~ N(mu,sigma), integrate y with 
		# http://en.wikipedia.org/wiki/Gauss-Hermite_quadrature
		yvec = sqrt(2.0) * p.sigma .* nodes .+ p.mu
		ywgt = weights .* pi^(-0.5)

		# end of period assets a are subject to borrowing constraints
		# borrowing limits: natural debt limit approach.
		# you can borrow the smaller of 3 times average income or
		# the natural debt limit, which is the lowest income state discounted 
		# bounds = [(-1)*yvec[1]*(1-p.R^i) /(1- p.R) for i in (p.nT-2):-1:1]

		# can borrow only up to lowest income value
		bounds = [(-1)*yvec[1] for i in (p.nT-2):-1:1]
		# bounds[bounds .< (-5)*mean(yvec)] = (-5)*mean(yvec)
		# bounds = vcat(bounds,0.0)  #last period
		# bounds = Float64[(-5)*mean(yvec) for i in 1:(p.nT-2)]
		bounds = vcat(bounds,0.0,0.0)  #last and penultimate periods end of period assets must be positive

		avec = zeros(p.na,p.nT)
		for i=(p.nT):-1:1
			avec[:,(p.nT)-i+1] = linspace(bounds[(p.nT)-i+1],p.a_high,p.na)
		end
		# avec = linspace((-5)*mean(yvec),p.a_high,p.na)

		# precompute next period's cash on hand.
		m1 = Float64[avec[ia,it+1]*p.R + yvec[iy] for ia in 1:p.na, iy in 1:p.ny, it in 1:(p.nT-1)]

		# if you want a deterministic age profile in income, use income().
		# you would have to change the params of income() though.
		# m1 = Float64[avec[ia]*p.R + income(yvec[iy],it) for ia in 1:p.na, iy in 1:p.ny, it in 1:p.nT]

		c1 = zeros(p.na,p.ny)
		ev = zeros(p.na,p.ny)

		m2 = zeros(p.ny)
		c2 = zeros(p.ny)

		C = zeros(p.na,p.nT)
		S = zeros(p.na,p.nT)
		M = zeros(p.na,p.nT)
		V = zeros(p.na,p.nT)
		Vzero = zeros(p.nT)

		dont = falses(p.na,p.ny,p.nT)


		return new(bounds,avec,yvec,ywgt,m1,c1,ev,m2,c2,C,S,M,V,Vzero,dont)
	end
end

type Envelope
	cond       :: Dict{Int,Vector{Float64}} 	# dict of functions conditional on discrete choice
	# cond_vzero :: Dict{Int,Float64} 	# dict of expected value of saving zero by discrete choice
	# vzero      :: Float64 				# expected value of saving zero on envelope
	env        :: Vector{Float64}		# actual envelope over discrete choices
	function Envelope()
		this = new()
		this.cond = [id => Float64[] for id in 1:2 ]
		# this.cond_vzero = [id => 0.0 for id in 1:2 ]
		# this.vzero = 0.0
		this.env = Float64[]
		return this
	end
	function Envelope(cond,cond_vzero,vzero,env)
		this = new()
		this.cond = copy(cond)
		# this.cond_vzero = copy(cond_vzero)
		# this.vzero = copy(vzero)
		this.env = copy(env)
		return this
	end
end

function cond(e::Envelope,which::Int)
	e.cond[which]
end
function cond(e::Dict{Int,Envelope},it::Int,which::Int)
	cond(e[it],which)
end
function env(e::Dict{Int,Envelope},it::Int)
	e[it].env
end
function set!(e::Dict{Int,Envelope},it::Int,v::Vector{Float64})
	e[it].env = deepcopy(v)
end
function set!(e::Dict{Int,Envelope},it::Int,which::Int,v::Vector{Float64})
	e[it].cond[which] = deepcopy(v)
end
function setvzero!(e::Dict{Int,Envelope},it::Int,which::Int,ind::Union(Int64,UnitRange),v::Float64)
	e[it].cond[which][ind] = v
end
function setvzero!(e::Dict{Int,Envelope},it::Int,which::Int,ind::UnitRange,v::Vector{Float64})
	e[it].cond[which][ind] = v
end
function vzero(e::Dict{Int,Envelope},it::Int,which::Int)
	e[it].cond[which][1]
end
function vzero(e::Dict{Int,Envelope},it::Int)
	e[it].env[1]
end

@doc "get a matrix of conditional values by period" ->
function dmat(e::Dict{Int,Envelope},it)
	v = e[it].cond[1]
	for id in 2:length(e[it].cond)
		v = hcat(v,e[it].cond[id])
	end
	v
end

@doc "get a matrix of all values conditional on dchoice" ->
function mat(e::Dict{Int,Envelope},which)
	v = e[1].cond[which]
	for it in 2:(length(e)-1)
		v = hcat(v,e[it].cond[which])
	end
	v
end

@doc "get a matrix of all envelope values" ->
function mat(e::Dict{Int,Envelope})
	v = env(e,1)
	for it in 1:(length(e)-1)
		v = hcat(v,env(e,it))
	end
	v
end



@doc """
Binary Choice Model with iid income uncertainty

uses cash-on-hand m=y+a as state variable
""" ->
type iidDModel <: Model

	# nD is number of discrete choices: nD = 2

	# computation grids
	avec::Vector{Float64}
	yvec::Vector{Float64}   # income support
	ywgt::Vector{Float64}   # income weights

	# intermediate objects (na,ny,nD)
	m1::Dict{Int,Dict}	# a dict[it] for each period
	c1::Matrix{Float64}
	ev::Matrix{Float64}

	# result objects
	m :: Dict{Int,Envelope}
	v :: Dict{Int,Envelope}
	c :: Dict{Int,Envelope}

	dchoice::Dict{Int,Dict} # a dict[it] for each period

	# C::Array{Float64,3} 	# consumption function on (na,nT,nD)
	# S::Array{Float64,3} 	# savings function on (na,nT,nD)
	# M::Array{Float64,3} 	# endogenous cash on hand on (na,nT,nD)
	# V::Array{Float64,2} 	# value function on (na,nT). Optional.
	# v::Array{Float64,3} 	# value function on (na,nT,nD). 
	# vzero::Array{Float64,2} 	# value function of saving zero

	# toc::Float64   # elapsed time

	@doc "Constructor for iid Dchoice Model" ->
	function iidDModel(p::Param)

		# avec          = scaleGrid(0.0,p.a_high,p.na,2)
		avec          = linspace(0.0,p.a_high,p.na)
		# nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature
		nodes,weights = quadpoints(p.ny,0,1)  # from FastGaussQuadrature
		N = Normal(0,1)
		nodes = quantile(N,nodes)

		# for y ~ N(mu,sigma), integrate y with 
		# http://en.wikipedia.org/wiki/Gauss-Hermite_quadrature
		# yvec = sqrt(2.0) * p.sigma .* nodes .+ p.mu
		yvec = nodes * p.sigma
		# ywgt = weights .* pi^(-0.5)
		ywgt = weights

		# precompute next period's cash on hand.
		# (na,ny,nD)
		# iD = 1: no work
		# iD = 2: work
		m1 = [it => [id => Float64[avec[ia]*p.R + income(it,yvec[iy]) * (id-1) for ia in 1:p.na, iy in 1:p.ny  ] for id=1:p.nD] for it=1:p.nT]
		c1 = zeros(p.na,p.ny)
		ev = zeros(p.na,p.ny)

		# dicts
		# m = [it => Envelope([id => zeros(p.na) for id in 1:2],[id => 0.0 for id in 1:2], 0.0, zeros(p.na)) for it in 1:p.nT]
		m = [it => Envelope() for it in 1:p.nT]
		v = [it => Envelope() for it in 1:p.nT]
		c = [it => Envelope() for it in 1:p.nT]
		dchoice = [it => ["d" => zeros(Int,p.na), "Vzero" => 0.0, "threshold" => 0.0] for it=1:p.nT]

		return new(avec,yvec,ywgt,m1,c1,ev,m,v,c,dchoice)
	end
end

@doc """
Model with AR1 income uncertainty V(m,y)

uses cash-on-hand m=y+a and current income state y as state variables. Tracking y is necessary in order to compute the conditional expectation on y.
""" ->
type AR1Model <: Model

	# computation grids
	avec::Vector{Float64}
	zvec::Vector{Float64}   # shock support
	yvec::Vector{Float64}   # income support
	ywgt::Matrix{Float64}   # transition matrix for income

	# intermediate objects (na,ny)
	m1::Array{Float64,2}	# matrix (na,ny)
	c1::Array{Float64,2}
	ev::Array{Float64,2}
	# intermediate objects (ny,1)
	m2::Vector{Float64} 
	c2::Vector{Float64} 

	# result objects on (na,ny,nT)
	C::Array{Float64,3} 	
	S::Array{Float64,3} 	
	M::Array{Float64,3} 	
	V::Array{Float64,3} 	
	Vzero::Array{Float64,2} 	# value function of saving zero

	toc::Float64   # elapsed time

	@doc "Constructor for AR1 Model" ->
	function AR1Model(p::Param)

		avec = linspace(p.a_low,p.a_high,p.na)

		# get grid and transition matrix for z
		z,ywgt = rouwenhorst(p.rho_z,0.0,p.eps_z,p.ny)

		yvec = z + p.mu

		# precompute next period's cash on hand.
		m1 = p.R * Float64[avec[i] + yvec[j] for i=1:p.na, j=1:p.ny]
		c1 = zeros(p.na,p.ny)
		ev = zeros(p.na,p.ny)

		m2 = zeros(p.ny)
		c2 = zeros(p.ny)

		C = zeros(p.na,p.ny,p.nT)
		S = zeros(p.na,p.ny,p.nT)
		M = zeros(p.na,p.ny,p.nT)
		V = zeros(p.na,p.ny,p.nT)
		Vzero = zeros(p.ny,p.nT)

		toc = 0.0

		return new(avec,z,yvec,ywgt,m1,c1,ev,m2,c2,C,S,M,V,Vzero,toc)
	end
end


# 3) model with AR1 income uncertainty, use current assets a and y as state variables: V(a,y)
type AR1Model_a <: Model

	# computation grids
	avec::Vector{Float64}
	zvec::Vector{Float64}   # shock support
	yvec::Vector{Float64}   # income support
	ywgt::Matrix{Float64}   # transition matrix for income

	# result objects on (na,ny,nT)
	C::Array{Float64,3} 	
	S::Array{Float64,3} 	
	M::Array{Float64,3} 	
	V::Array{Float64,3} 	
	EV::Array{Float64,1} 	

	toc::Float64   # elapsed time

	@doc "Constructor for AR1 Model" ->
	function AR1Model_a(p::Param)

		avec = linspace(p.a_low,p.a_high,p.na)

		# get grid and transition matrix for z
		z,ywgt = rouwenhorst(p.rho_z,0.0,p.eps_z,p.ny)

		yvec = z + p.mu

		C = zeros(p.na,p.ny,p.nT)
		S = zeros(p.na,p.ny,p.nT)
		M = zeros(p.na,p.ny,p.nT)
		V = zeros(p.na,p.ny,p.nT)
		EV = zeros(p.na)

		toc = 0.0

		return new(avec,z,yvec,ywgt,C,S,M,V,EV,toc)
	end
end

@doc """
rouwenhorst AR1 approximation 


This is taken from [http://karenkopecky.net/RouwenhorstPaperFinal.pdf](Karen Kopecky's paper)

"""->
function rouwenhorst(rho::Float64,mu_eps,sigma_eps,n)
	q = (rho+1)/2
	nu = ((n-1)/(1-rho^2))^(1/2) * sigma_eps
	P = reshape([q,1-q,1-q,q],2,2)

	for i=2:n-1

		P = q * vcat(hcat(P , zeros(i,1)),zeros(1,i+1)) .+ (1-q).* vcat( hcat(zeros(i,1),P), zeros(1,i+1)) .+ 
		(1-q) .* vcat(zeros(1,i+1),hcat(P,zeros(i,1))) .+ q .*vcat(zeros(1,i+1),hcat(zeros(i,1),P))
		P[2:i,:] = P[2:i,:] ./ 2

	end

	z = linspace(mu_eps/(1-rho)-nu,mu_eps/(1-rho)+nu,n);
	return (z,P)
end

# asset grid scaling
function scaleGrid(lb::Float64,ub::Float64,n::Int,logorder::Int=1) 
	out = zeros(n)
	if logorder==1
		off = 1
		if lb<0 
			off = 1 - lb #  adjust in case of neg bound
		end
		out[1] = log(lb + off) 
		out[n] = log(ub + off) 
		out    = linspace(out[1],out[n],n)
		out    = exp(out) - off  
	elseif logorder==2
		off = 1
		if lb<0 
			off = 1 - lb #  adjust in case of neg bound
		end
		out[1] = log( log(lb + off) + off )
		out[n] = log( log(ub + off) + off )
		out    = linspace(out[1],out[n],n)
		out    = exp( exp(out) - off ) - off
	else
		error("supports only double log grid")
	end
end
