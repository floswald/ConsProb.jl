


# run all
function run()

	p = Param()

	# solve by standard euler equation root finding
	EEmod = Model(p)
	EEmod.toc = @elapsed Euler(EEmod,p)



	# solve maximizing the value function backward iteration
	VFmod = Model(p)
	VFmod.toc = @elapsed VFbi(VFmod,p)

	# solve by EGM
	EGMmod = Model(p)
	EGMmod.toc = @elapsed EGM(EGMmod,p)
	

	# plot results
	plots(EEmod,EGMmod,VFmod,p)

end




#utility functions

# utility
function u(x::Float64,p::Param)
	p.oneover_oneminusgamma * (x^p.oneminusgamma)
end
function u{T}(x::Array{T},p::Param)
	n = length(x)
	y = zeros(T,n)
	for i in 1:n
		y[i] = p.oneover_oneminusgamma * (x[i]^p.oneminusgamma)
	end
	y
end

# partial derivative of utility wrt c
function up(c::Float64,p::Param)
	c ^ (p.neg_gamma)
end
function up(c::Array{Float64,2},p::Param)
	n = length(c)
	x = zeros(size(c))
	for i in 1:n
		x[i] = up(c[i],p)
	end
	x
end
function up(c::Array{Float64},p::Param)
	n = length(c)
	x = zeros(n)
	for i in 1:n
		x[i] = up(c[i],p)
	end
	x
end

# inverse of partial derivative
function iup(u::Float64,p::Param)
	u ^ p.neg_oneover_gamma
end
function iup(u::Array{Float64},p::Param)
	n = length(u)
	x = zeros(n)
	for i in 1:n
		x[i] = iup(u[i],p)
	end
	return x
end



# endogenous grid method
function EGM(m::Model,p::Param)

	# final period: consume everything.
	m.M[:,p.nT] = linspace(p.a_low,p.a_high*4,p.na)
	m.C[:,p.nT] = linspace(p.a_low,p.a_high*4,p.na)

	# preceding periods
	for it in (p.nT-1):-1:1

		# interpolate optimal consumption from next period on all cash-on-hand states
		# using C[:,it+1] and M[:,it+1], find c(m,it)

		tmpx = [0.0, m.M[:,it+1] ] 
		tmpy = [0.0, m.C[:,it+1] ]
		for ia in 1:p.na
			for iy in 1:p.ny
				m.c1[ia+p.na*(iy-1)] = linearapprox(tmpx,tmpy,m.m1[ia+p.na*(iy-1)],1,p.na)
			end
		end

		# get expected marginal value of saving: RHS of euler equation
		# beta * R * E[ u'(c_{t+1}) ] 
		Eu = p.R * p.beta .* up(m.c1,p) * m.ywgt

		# get optimal consumption today from euler equation: invert marginal utility
		m.C[:,it] = iup(Eu,p)

		# get endogenous grid today
		m.M[:,it] = m.C[:,it] .+ m.avec


	end
end



# solving the euler equation
function Euler(m::Model,p::Param)

	# final period: consume everything.
	m.M[:,p.nT] = linspace(p.a_low,p.a_high*4,p.na)
	m.C[:,p.nT] = linspace(p.a_low,p.a_high*4,p.na)

	# preceding periods
	for it in (p.nT-1):-1:1

		for ia in 1:p.na

			cash = m.avec[ia] 	# current cash on hand

			# consumption equal to cash on hand
			res = EulerResid(cash,cash,m.C[:,it+1],p,m)

			if res < 0
				m.C[ia,it] = cash
			else
				m.C[ia,it] = fzero((x)->EulerResid(x,cash,m.C[:,it+1],p,m),cash/2,[1e-6,cash])
			end

		end


	end

end

# Euler Equation Residual
function EulerResid(c::Float64,cash::Float64,cplus::Vector{Float64},p::Param,m::Model)

	# given current c, what is next period's cash on hand
	m.m2 = p.R * (cash - c) .+ m.yvec  # (ny,1)

	# interpolate optimal consumption c(t+1), given c(t), on each y-state
	for iy in 1:p.ny
		m.c2[iy] = linearapprox([0,m.avec],[0,cplus],m.m2[iy],1,p.na)
	end

	# Expected marginal utility of consumption (RHS of euler equation)
	Euc = dot(p.R .* p.beta .* up(m.c2,p), m.ywgt) 	# (1,1)

	# residual
	c - iup(Euc,p)

end



# finding the maximum of the value function backward iteration
function VFbi(m::Model,p::Param)

	# final period: consume everything.
	m.C[:,p.nT] = linspace(p.a_low,p.a_high*4,p.na)
	m.V[:,p.nT] = u(m.C[:,p.nT],p)

	# preceding periods
	for it in (p.nT-1):-1:1

		for ia in 1:p.na

			cash = m.avec[ia] 	# current cash on hand

			# Brent's method for minimizing a function
			# withotu derivatives
			x = optimize((x)->VFobj(x,cash,m.V[:,it+1],m,p),1e-6,cash)
			m.V[ia,it] = -x.f_minimum
			m.C[ia,it] = x.minimum
		end
	end
end


function VFobj(c::Float64,cash::Float64,Vplus::Vector{Float64},m::Model,p::Param)

	# implied cash on hand tomorrow
	s = (cash - c) * p.R .+ m.yvec
	for iy in 1:p.ny
		m.m2[iy] = linearapprox(m.avec,Vplus,s[iy],1,p.na)
	end
	v = u(c,p) + p.beta * dot(m.m2,m.ywgt)
	return -v
end


function linearapprox(x::Vector{Float64},y::Vector{Float64},xi::Float64,lo::Int,hi::Int)
	r = 0.0
	n = length(x)

	# determining bounds 
	if xi == x[1]
		r = y[1] 
		return r
	elseif xi < x[1]
		# get linear approx below
		@inbounds r = y[1] + (y[2] - y[1]) * (xi - x[1])  / (x[2] - x[1])
		return r
	end
	if xi == x[n]
		r = y[n] 
		return (r,n)
	elseif xi > x[n]
		# get linear approx above
		@inbounds r = y[n] + (y[n] - y[n-1]) * (xi - x[n])  / (x[n] - x[n-1])
		return r
	end

	# if have to find interval
	if hi - lo > 1
		jinf = searchsortedlast(x,xi,lo,hi,Base.Forward)	# get rid
	# if not, lo is jinf
	else
		jinf = lo
	end
	@inbounds r = (y[jinf] * (x[jinf+1] - xi) + y[jinf+1] * (xi - x[jinf]) ) / (x[jinf+1] - x[jinf])
	return r
end





# plotting
function plots(EE::Model,EGM::Model,VF::Model,p::Param)

	figure()

	# plot consumption funcitons
	# ==========================

	# solve by maximizing the VF
	subplot(1,3,1)
	plot(VF.avec,VF.C[:,1:(p.nT-1)])
	xlabel("assets")
	ylabel("consumption")
	ylim([p.a_low,p.a_high])
	xlim([p.a_low,p.a_high])
	plot(ylim(),ylim())
	title("maximize the value function: $(round(VF.toc,5)) secs")

	subplot(1,3,2)
	plot(EE.avec,EE.C[:,1:(p.nT-1)])
	xlabel("assets")
	ylim([p.a_low,p.a_high])
	xlim([p.a_low,p.a_high])
	plot(ylim(),ylim())
	title("Solving the Euler Equation: $(round(EE.toc,5)) secs")

	# solve by EGM
	subplot(1,3,3)
	plot(EGM.M[:,1:(p.nT-1)],EGM.C[:,1:(p.nT-1)])
	xlabel("cash on hand")
	ylim([p.a_low,p.a_high])
	xlim([p.a_low,p.a_high])
	plot(ylim(),ylim())
	title("Endogenous Grid Method: $(round(EGM.toc,5)) secs")

end
