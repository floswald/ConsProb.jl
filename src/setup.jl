

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

	# grids
	na :: Int 	# asset grid
	ny :: Int   # income grid (support points)
	nT :: Int   # maximal age
	a_high::Float64
	a_low::Float64

	# iid income unertainty
	mu::Float64
	sigma::Float64

	# AR1 income uncertainty
	# y_t = z_t
	# ln z_t = rho_z * ln z_{t-1} + eps_z_{t}
	rho_z::Float64
	eps_z::Float64
	AgeProf::Vector{Float64}

	function Param()

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

		# iid income uncertainty params
		mu = 100 	# mean income: 30K
		sigma = 10  # sd income

		# AR1 income uncertainty
		# params from Ayiagari
		rho_z = 0.9
		eps_z = 0.12

		return new(gamma,neg_gamma,oneminusgamma,oneover_oneminusgamma,neg_oneover_gamma,beta,R,na,ny,nT,a_high,a_low,mu,sigma,rho_z,eps_z)
	end
end

# consumption model with iid income uncertainty
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
		yvec = sqrt(2.0) * p.sigma .* nodes .+ p.mu
		ywgt = weights .* pi^(-0.5)

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


# consumption model with AR1 income uncertainty
type Model2

	# computation grids
	avec::Vector{Float64}
	yvec::Vector{Float64}   # income support
	ywgt::Matrix{Float64}   # transition matrix for income

	# intermediate objects (na,ny)
	m1::Array{Float64,2}	# matrix (na,ny)
	c1::Array{Float64,2}
	# intermediate objects (ny,1)
	m2::Vector{Float64} 
	c2::Vector{Float64} 

	# result objects on (na,ny,nT)
	C::Array{Float64,3} 	
	S::Array{Float64,3} 	
	M::Array{Float64,3} 	
	V::Array{Float64,3} 	
	EV::Array{Float64,3} 	

	toc::Float64   # elapsed time

	function Model2(p::Param)

		avec          = linspace(p.a_low,p.a_high,p.na)

		# get grid and transition matrix for z
		z,ywgt = rouwenhorst(p.rho_z,0.0,p.eps_z,p.ny)

		yvec = exp(z) + 20

		# precompute next period's cash on hand.
		m1 = p.R * kron(avec,ones(1,p.ny)) .+ kron(ones(p.na,1),transpose(yvec))
		c1 = zeros(p.na,p.ny)

		m2 = zeros(p.ny)
		c2 = zeros(p.ny)

		C = zeros(p.na,p.ny,p.nT)
		S = zeros(p.na,p.ny,p.nT)
		M = zeros(p.na,p.ny,p.nT)
		V = zeros(p.na,p.ny,p.nT)
		EV = zeros(p.na,p.ny,p.nT)

		toc = 0.0


		return new(avec,yvec,ywgt,m1,c1,m2,c2,C,S,M,V,EV,toc)
	end
end


# rouwenhorst AR1 approximation after 
# http://karenkopecky.net/RouwenhorstPaperFinal.pdf
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
