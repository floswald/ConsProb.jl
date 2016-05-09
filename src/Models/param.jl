
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
		a_lowT  = 0.0
		nD     = 2

		cfloor = 0.001
		alpha = 2.3

		# iid income uncertainty params
		mu = 10 	# mean income: 30K
		sigma = 1  # sd income
		# mu = 0 	# mean income: 30K
		# sigma = 0.25  # sd income

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
		nT = 8
		a_high = 300.0
		a_lowT  = 1e-6
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
