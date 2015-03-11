

# param type 

type Param

	# CRRA
	gamma                 :: Float64
	neg_gamma             :: Float64
	oneminusgamma         :: Float64
	oneover_oneminusgamma :: Float64
	neg_oneover_gamma     :: Float64

	beta::Float64
	R::Float64

	# income unertainty
	mu::Float64
	sigma::Float64
	dist::ASCIIString

	# if income is AR1
	# rho::Float64

	# grids
	na :: Int 	# asset grid
	ny :: Int   # income grid (support points)
	nT :: Int   # maximal age
	a_high::Float64
	a_low::Float64


	function Param()

		gamma                 = 2.5
		neg_gamma             = (-1.0) * gamma
		oneminusgamma         = 1.0 - gamma
		oneover_oneminusgamma = 1.0 / oneminusgamma
		neg_oneover_gamma     = (-1.0) / gamma

		beta                  = 0.96
		R                     = 1.03

		mu = 30 	# mean income: 30K
		sigma = 3  # sd income
		dist = "normal"

		na = 100
		ny = 10
		nT = 8

		a_high = 100.0
		a_low  = 1.0

		return new(gamma,neg_gamma,oneminusgamma,oneover_oneminusgamma,neg_oneover_gamma,beta,R,mu,sigma,dist,na,ny,nT,a_high,a_low)
	end
end


type Model

	# computation grids
	avec::Vector{Float64}
	yvec::Vector{Float64}   # income support
	ywgt::Vector{Float64}   # income weights

	# intermediate objects (na,ny)
	m1::Array{Float64,2}	# matrix (na,ny)
	c1::Array{Float64,2}
	# intermediate objects (ny,1)
	m2::Vector{Float64} 
	c2::Vector{Float64} 

	# result objects
	C::Array{Float64,2} 	# consumption function on (na,nT)
	S::Array{Float64,2} 	# savings function on (na,nT)
	M::Array{Float64,2} 	# endogenous cash on hand on (na,nT)
	V::Array{Float64,2} 	# value function on (na,nT). Optional.

	toc::Float64   # elapsed time

	function Model(p::Param)

		avec          = linspace(p.a_low,p.a_high,p.na)
		nodes,weights = gausshermite(p.ny)  # from FastGaussQuadrature

		# for y ~ N(mu,sigma), integrate y with 
		# http://en.wikipedia.org/wiki/Gauss-Hermite_quadrature

		if p.dist == "lognormal"
			yvec = exp( sqrt(2.0) * log(p.sigma) .* nodes .+ log(p.mu) )
		else
			yvec = sqrt(2.0) * p.sigma .* nodes .+ p.mu
		end
		ywgt          = weights .* pi^(-0.5)

		# precompute next period's cash on hand.
		m1 = p.R * kron(avec,ones(1,p.ny)) .+ kron(ones(p.na,1),transpose(yvec))
		c1 = zeros(p.na,p.ny)

		m2 = zeros(p.ny)
		c2 = zeros(p.ny)

		C = zeros(p.na,p.nT)
		S = zeros(p.na,p.nT)
		M = zeros(p.na,p.nT)
		V = zeros(p.na,p.nT)

		toc = 0.0


		return new(avec,yvec,ywgt,m1,c1,m2,c2,C,S,M,V,toc)
	end
end


