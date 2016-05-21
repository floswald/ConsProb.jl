

"""
Holds the user-set parameter values. 

**values to be set:**

* `gamma`: CRRA 
* `beta`: discount factor
* `R`: gross interest rate (i.e. 1+r )
* `na`: number of grid points for assets
* `ny`: number of grid points for income
* `nT`: number of time periods
* `a_lowT`: lower bound on assets in final period
* `a_low`: lower bound on assets
* `a_high`: upper bound on assets
* `mu`: unconditional mean of income (iid case)
* `sigma`: unconditional variance of income (iid case)
* `rho_z`: AR1 coefficient of income (AR1 case)
* `eps_z`: standard deviation of AR1 innovation (AR1 case)
* `dorefinements`: boolean switch of whether to filter out whiggles
* `alpha`: disutility of work

"""
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
	a_lowT::Float64
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
	function Param(gamm::Float64,alow=0.0)

		gamma                 = gamm
		neg_gamma             = (-1.0) * gamma
		oneminusgamma         = 1.0 - gamma
		oneover_oneminusgamma = 1.0 / oneminusgamma
		neg_oneover_gamma     = (-1.0) / gamma

		beta                  = 0.95
		R                     = 1.05

		na     = 200
		ny     = 50
		nT     = 25
		a_high = 50.0
		a_low  = alow
		a_lowT  = 0.0
		nD     = 2

		cfloor = 0.001
		alpha = 2.3

		# iid income uncertainty params
		# mu = 10 	# mean income: 30K
		# sigma = 1  # sd income
		mu = 0 	# mean income: 30K
		sigma = 0.20  # sd income

		# AR1 income uncertainty
		# params from Ayiagari
		rho_z = 0.9
		eps_z = 1

		dorefine=true
		printdebug=true

		return new(gamma,neg_gamma,oneminusgamma,oneover_oneminusgamma,neg_oneover_gamma,beta,R,na,ny,nT,a_high,a_low,a_lowT,nD,cfloor,alpha,mu,sigma,rho_z,eps_z,dorefine,printdebug)
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
		nT = 25
		a_high = 300.0
		a_lowT  = -25.0
		a_low  = 1e-6
		nD = 2

		cfloor = 0.001
		alpha = 1.0

		# iid income uncertainty params
		sigma = 0.1  # sd income

		# AR1 income uncertainty
		# params from Ayiagari
		rho_z = 0.9
		eps_z = 1

		dorefine=false

		return new(gamma,neg_gamma,oneminusgamma,oneover_oneminusgamma,neg_oneover_gamma,beta,R,na,ny,nT,a_high,a_low,a_lowT,nD,cfloor,alpha,mu,sigma,rho_z,eps_z,dorefine)
	end
end


# Abstract Model type
# this is an overall model type, with potential subtypes

"""
Abstract Model type

**There are three model types:**

1. model with iid income uncertainty, use cash-on-hand m=y+a as state variable: V(a)
2. model with AR1 income uncertainty, use cash-on-hand m=y+a and y as state variables: V(m,y)
3. model with AR1 income uncertainty, use current assets a and y as state variables: V(a,y)
"""
abstract Model




"""
Model with iid income uncertainty after Deaton (1991)

uses cash-on-hand m=y+a as state variable
"""
type iidModel <: Model

	# computation grids
	avec::Array{Vector{Float64}}	# different avec in each period possible
	yvec::Vector{Float64}   # income shock support
	ywgt::Vector{Float64}   # income weights
	incmat::Matrix{Float64}   # income levels

	# intermediate objects (na,ny)
	mnext::Array{Float64,3}	# matrix (na,ny,nt)
	cnext::Array{Float64,2}
	ev::Array{Float64,2}
	# intermediate objects (ny,1)
	m2::Vector{Float64} 
	c2::Vector{Float64} 

	# result objects
	# these are functions with a special value at their lower bound.
	# there is one for each period
	C::Array{Bfun} 	# consumption function
	S::Array{Bfun} 	# savings function
	M::Array{Bfun} 	# endogenous cash on hand
	V::Array{Bfun} 	# value function. b(V) gives expected value at lower bound
	# Vzero::Array{Float64,1} 	# expected value function of saving zero

	toc::Float64   # elapsed time
	it::Int  # current time period


	"""
	Constructor for iid Model
	"""
	function iidModel(p::Param)

		nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature

		# for y ~ N(mu,sigma), integrate y with 
		# http://en.wikipedia.org/wiki/Gauss-Hermite_quadrature
		yvec = sqrt(2.0) * p.sigma .* nodes .+ p.mu
		ywgt = weights .* pi^(-0.5)
		incmat = hcat([income(yvec,it,p) for it=1:p.nT]...)

		# precompute next period's cash on hand
		# first: what is lower bound on end-of-period assets in each period?
		if p.a_low < 0
			lowest_incomes = [income(yvec[1],it,p) for it=1:p.nT]
			bounds = blim(lowest_incomes,p.a_low,p.R,p.nT)
		else 
			bounds = [p.a_low for i=1:p.nT]
		end
		avec = [scaleGrid(bounds[i],p.a_high,p.na,2) for i=1:p.nT]

		# precompute next period's cash on hand.
		# mnext = Float64[avec[ia]*p.R + yvec[iy] for ia in 1:p.na, iy in 1:p.ny]

		# if you want a deterministic age profile in income, use income().
		# you would have to change the params of income() though.
		mnext = Float64[avec[it][ia]*p.R + income(yvec[iy],it,p) for ia in 1:p.na, iy in 1:p.ny, it in 1:p.nT]

		# mnext = zeros(p.na,p.ny)
		cnext = zeros(p.na,p.ny)
		ev = zeros(p.na,p.ny)

		m2 = zeros(p.ny)
		c2 = zeros(p.ny)

		C = [Bfun(zeros(p.na),p.cfloor) for i=1:p.nT]
		S = [Bfun(zeros(p.na),NaN) for i=1:p.nT]
		M = [Bfun(zeros(p.na),NaN) for i=1:p.nT]
		V = [Bfun(zeros(p.na),NaN) for i=1:p.nT]
		it = p.nT
		toc = 0.0


		return new(avec,yvec,ywgt,incmat,mnext,cnext,ev,m2,c2,C,S,M,V,it,toc)
	end
end


"""
	cashnext(a::Float64,y::Float64,p::Param)

Compute next periods cash in hand.
"""
function cashnext(a::Float64,y::Float64,p::Param)
	a*p.R +y
end

"""
	invcashnext(cashnext::Float64,y::Float64,p::Param)

Compute end-of-period asset corresponding to next periods cash in hand.
"""
function invcashnext(cashnext::Float64,y::Float64,p::Param)
	(cashnext - y) / p.R
end


"""
Binary Choice Model with iid income uncertainty

uses cash-on-hand m=y+a as state variable
"""
type iidDModel <: Model

	# nD is number of discrete choices: nD = 2

	# computation grids
	avec::Vector{Float64}
	yvec::Vector{Float64}   # income support
	ywgt::Vector{Float64}   # income weights

	# intermediate objects (na,ny,nD)
	mnext::Dict{Int,Dict}	# a dict[it] for each period
	cnext::Matrix{Float64}
	ev::Matrix{Float64}

	# result objects
	m :: Dict{Int,Envelope}
	v :: Dict{Int,Envelope}
	c :: Dict{Int,Envelope}

	# dchoice::Dict{Int,Dict} # a dict[it] for each period

	# C::Array{Float64,3} 	# consumption function on (na,nT,nD)
	# S::Array{Float64,3} 	# savings function on (na,nT,nD)
	# M::Array{Float64,3} 	# endogenous cash on hand on (na,nT,nD)
	# V::Array{Float64,2} 	# value function on (na,nT). Optional.
	# v::Array{Float64,3} 	# value function on (na,nT,nD). 
	# vzero::Array{Float64,2} 	# value function of saving zero

	# toc::Float64   # elapsed time

	"""
	Constructor for iid Dchoice Model
	"""
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
		mnext = [it => [id => Float64[avec[ia]*p.R + income(it,yvec[iy]) * (id-1) for ia in 1:p.na, iy in 1:p.ny  ] for id=1:p.nD] for it=1:p.nT]
		cnext = zeros(p.na,p.ny)
		ev = zeros(p.na,p.ny)

		# dicts
		# m = [it => Envelope([id => zeros(p.na) for id in 1:2],[id => 0.0 for id in 1:2], 0.0, zeros(p.na)) for it in 1:p.nT]
		m = [it => Envelope() for it in 1:p.nT]
		v = [it => Envelope() for it in 1:p.nT]
		c = [it => Envelope() for it in 1:p.nT]
		# dchoice = [it => ["d" => zeros(Int,p.na), "Vzero" => 0.0, "threshold" => 0.0] for it=1:p.nT]

		return new(avec,yvec,ywgt,mnext,cnext,ev,m,v,c)
	end
end

"""
Model with AR1 income uncertainty V(m,y)

uses cash-on-hand m=y+a and current income state y as state variables. Tracking y is necessary in order to compute the conditional expectation on y.
"""
type AR1Model <: Model

	# computation grids
	avec::Vector{Float64}
	zvec::Vector{Float64}   # shock support
	yvec::Vector{Float64}   # income support
	ywgt::Matrix{Float64}   # transition matrix for income

	# intermediate objects (na,ny)
	mnext::Array{Float64,2}	# matrix (na,ny)
	cnext::Array{Float64,2}
	ev::Array{Float64,2}
	# intermediate objects (ny,1)
	m2::Vector{Float64} 
	c2::Vector{Float64} 

	# result objects on (na,ny,nT)
	C::Array{Bfun,2} 	# consumption function
	S::Array{Bfun,2} 	# savings function
	M::Array{Bfun,2} 	# endogenous cash on hand
	V::Array{Bfun,2} 	# value function. b(V) gives expected value at lower bound

	toc::Float64   # elapsed time

	"""
	Constructor for AR1 Model
	"""
	function AR1Model(p::Param)

		avec = linspace(p.a_low,p.a_high,p.na)

		# get grid and transition matrix for z
		z,ywgt = rouwenhorst(p.rho_z,0.0,p.eps_z,p.ny)

		yvec = z + p.mu

		# precompute next period's cash on hand.
		mnext = Float64[p.R * avec[i] + yvec[j] for i=1:p.na, j=1:p.ny]
		cnext = zeros(p.na,p.ny)
		ev = zeros(p.na,p.ny)

		m2 = zeros(p.ny)
		c2 = zeros(p.ny)

		C = [Bfun(zeros(p.na),p.cfloor) for j=1:p.ny, i=1:p.nT]
		S = [Bfun(zeros(p.na),NaN)      for j=1:p.ny, i=1:p.nT]
		M = [Bfun(zeros(p.na),NaN)      for j=1:p.ny, i=1:p.nT]
		V = [Bfun(zeros(p.na),NaN)      for j=1:p.ny, i=1:p.nT]
		it = p.nT
		toc = 0.0

		return new(avec,z,yvec,ywgt,mnext,cnext,ev,m2,c2,C,S,M,V,toc)
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

	"""
	Constructor for AR1 Model
	"""
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

