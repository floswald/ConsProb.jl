


# run all
function runall()

	# iid income model

	p = Param()

	# solve by standard euler equation root finding
	EEmod = iidModel(p)
	EEmod.toc = @elapsed Euler!(EEmod,p)



	# solve maximizing the value function backward iteration
	VFmod = iidModel(p)
	VFmod.toc = @elapsed VFbi!(VFmod,p)

	# solve by EGM
	EGMmod = iidModel(p)
	EGMmod.toc = @elapsed EGM!(EGMmod,p)

	# plot results
	plots(EEmod,EGMmod,VFmod,p)

	# AR1 income model
	# ================

	EGMmod = AR1Model(p)
	EGMmod.toc = @elapsed EGM!(EGMmod,p)

	VFmod = AR1Model(p)
	VFmod.toc = @elapsed VFbi!(VFmod,p)

	# does it matter whether I compute the model on 
	# current assets, given y, or
	# cash-on-hand, given y?
	VFmod_a = AR1Model_a(p)
	VFmod_a.toc = @elapsed VFbi!(VFmod_a,p)

	# plot results
	plots(EGMmod,VFmod,VFmod_a,p,1)  # plot period 1
end




#utility function

# utility without discrete choice
function u(x::Float64,p::Param)
	p.oneover_oneminusgamma * (x^p.oneminusgamma)
end
function u{T}(x::Array{T},p::Param)
	n = length(x)
	y = zeros(T,n)
	for i in 1:n
		y[i] = u(x[i],p)
	end
	y
end

# utility with discrete choice
function u(x::Float64,working::Bool,p::Param)
	p.oneover_oneminusgamma * (x^p.oneminusgamma) + p.alpha*working
end
function u{T}(x::Array{T},working::Bool,p::Param)
	n = length(x)
	y = zeros(T,n)
	for i in 1:n
		y[i] = u(x[i],working,p)
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



"lifecycle profile in income

returns income as function of deterministic age profile + shock. Assume that it=1 is age 20. Taken from Ikshkakov at https://github.com/fediskhakov/egdst
"
function income(it::Int,shock::Float64)
	age = it + 19
	exp( 1.5 + age*0.04 - 0.0004*(age^2) + shock)
end


# endogenous grid method
function EGM!(m::iidModel,p::Param)

	# final period: consume everything.
	m.M[:,p.nT] = linspace(p.a_low,p.a_high*4,p.na)
	m.C[:,p.nT] = linspace(p.a_low,p.a_high*4,p.na)
	m.C[m.C[:,p.nT].<p.cfloor,p.nT] = p.cfloor

	m.V[:,p.nT] = u(m.C[:,p.nT],p) + p.beta * 0.0

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

		# floor consumption
		m.C[m.C[:,it].<p.cfloor,it] = p.cfloor


		# get endogenous grid today
		m.M[:,it] = m.C[:,it] .+ m.avec

		# compute value function
		# ======================

		# expected value function (na,ny)
		fill!(m.ev,NaN)
		# dont: don't interpolate anything.
		if it==(p.nT-1)
			dont = trues(size(m.m1))
		else
			dont = m.m1 .< m.M[1,it+1]	# wherever potential next period's cash on hand (m.m1) is less than the lowest grid point of the endogenous grid next period (m.M), the agent will be credit constrained and will be saving zero (m.EV[1,it+1])
		end

		vv = m.V[:,it+1]
		tmpx = m.M[:,it+1]  
		for ia in 1:p.na
			for iy in 1:p.ny
				idx = ia+p.na*(iy-1)
				if dont[idx]
					m.ev[idx] = u(m.m1[idx],p) + p.beta * m.Vzero[it+1]
				else
					m.ev[idx] = linearapprox(tmpx,vv,m.m1[idx],1,p.na)
				end
			end
		end
		ev = m.ev * m.ywgt
		# if abs(m.avec[1]) > 1e-6
		# 	error("first element of avec is assumed to be zero: it's not!")
		# end
		m.Vzero[it] = ev[1] # save expected value of saving zero in first element.
		m.V[:,it]  = u(m.C[:,it],p) + p.beta * ev 
	end
end



# endogenous grid method for AR1 model
function EGM!(m::AR1Model,p::Param)

	# final period: consume everything.
	m.M[:,:,p.nT] = repmat(linspace(p.a_low,p.a_high*4,p.na),1,p.ny)
	m.C[:,:,p.nT] = repmat(linspace(p.a_low,p.a_high*4,p.na),1,p.ny)
	cc = m.C[:,:,p.nT]
	cc[cc.<p.cfloor] = p.cfloor
	m.C[:,:,p.nT] = cc

	m.V[:,:,p.nT] = u(m.C[:,:,p.nT],p) + p.beta * 0.0

	# preceding periods
	for it in (p.nT-1):-1:1

		# conditional on current income state
		for iy in 1:p.ny

			# interpolate optimal consumption from next period on all cash-on-hand states
			# using C[:,it+1] and M[:,it+1], find c(m,it)

			# next period's income index
			for iiy in 1:p.ny
				tmpx = [0.0, m.M[:,iiy,it+1] ] 
				tmpy = [0.0, m.C[:,iiy,it+1] ]
				for ia in 1:p.na
					m.c1[ia+p.na*(iiy-1)] = linearapprox(tmpx,tmpy,m.m1[ia+p.na*(iiy-1)],1,p.na)
					# m.c1[ia,iiy] = linearapprox(tmpx,tmpy,m.m1[ia,iiy],1,p.na)
				end
			end

			# get expected marginal value of saving: RHS of euler equation
			# beta * R * E[ u'(c_{t+1}) ] 
			Eu = p.R * p.beta .* m.ywgt[iy,:] * transpose(up(m.c1,p))

			# get optimal consumption today from euler equation: invert marginal utility
			m.C[:,iy,it] = iup(Eu,p)
			
			# floor consumption
			m.C[m.C[:,iy,it].<p.cfloor,iy,it] = p.cfloor

			# get endogenous grid today
			m.M[:,iy,it] = m.C[:,iy,it] .+ m.avec


			# compute value function
			# ======================

			# expected value function (na,ny)
			fill!(m.ev,NaN)
			# dont: don't interpolate anything.
			if it==(p.nT-1)
				# if next period is final period, don't have to worry about
				# next (i.e. period T+1) savings
				dont = trues(size(m.m1))
			else
				dont = m.m1 .< m.M[1,iy,it+1]	# wherever potential next period's cash on hand (m.m1) is less than the lowest grid point of the endogenous grid next period (m.M), the agent will be credit constrained and will be saving zero (here stored in m.Vzero[iy,it+1])
			end

			# for next period's income state
			for iiy in 1:p.ny
				tmpx = m.M[:,iiy,it+1]  
				vv   = m.V[:,iiy,it+1]
				for ia in 1:p.na
					idx = ia+p.na*(iiy-1)
					if dont[idx]
						m.ev[idx] = u(m.m1[idx],p) + p.beta * m.Vzero[iiy,it+1]
					else
						m.ev[idx] = linearapprox(tmpx,vv,m.m1[idx],1,p.na)
					end
				end
			end
			ev = transpose(m.ywgt[iy,:] * transpose( m.ev ))
			if abs(m.avec[1]) > 1e-6
				error("first element of avec is assumed to be zero: it's not!")
			end
			m.Vzero[iy,it] = ev[1] # save expected value of saving zero in first element.
			m.V[:,iy,it]  = u(m.C[:,iy,it],p) .+ p.beta * ev
		end  # current income
	end  # age
end



# solving the euler equation
function Euler!(m::iidModel,p::Param)

	# final period: consume everything.
	m.M[:,p.nT] = linspace(p.a_low,p.a_high*4,p.na)
	m.C[:,p.nT] = linspace(p.a_low,p.a_high*4,p.na)

	# preceding periods
	for it in (p.nT-1):-1:1

		for ia in 1:p.na

			cash = m.avec[ia] 	# current cash on hand

			# consumption equal to cash on hand
			res = EulerResid(cash,cash,m.C[:,it+1],p,m)

			# this is an implication of 
			# equation (6) in Deaton ECTA (1991):
			# u'(c_t) = max[ u'(m_t), beta * R * u'(c_{t+1}) ]
			# with constraint a_t >= 0, if
			# c_t = m_t => a_{t+1} = 0   (consume everything, possible constrained: wanted to consume more by borrowing, but could not.)
			# c_t < m_t => a_{t+1} > 0   (consume less than m and save some for tomorrow)
			# c_t > m_t => a_{t+1} < 0   (consume more than m by borrowing)
			# the residual function EulerResid returns
			# c - u'^(-1) [ beta * E[ u'(beta * R * E[ u'(c_{t+1}) ] )] ]
			# where c_{t+1} is implied by today's choice c.
			# if that difference is negative when c_t = m_t, this means that 
			# c_t < c_{t+1} or
			# u'(c_t) > u'(c_{t+1}), or
			# u'(m_t) > u'(c_{t+1}), and therefore
			# u'(c_t) = max[ u'(m_t), beta * R * u'(c_{t+1}) ] implies
			# that this consumer is borrowing constrained and consumes all cash in hand.
			if res < 0
				m.C[ia,it] = cash
			else
				m.C[ia,it] = fzero((x)->EulerResid(x,cash,m.C[:,it+1],p,m),cash/2,[1e-6,cash])
			end

		end


	end

end

# Euler Equation Residual
function EulerResid(c::Float64,cash::Float64,cplus::Vector{Float64},p::Param,m::iidModel)

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
function VFbi!(m::iidModel,p::Param)

	# final period: consume everything.
	m.C[:,p.nT] = linspace(p.a_low,p.a_high*4,p.na)
	m.V[:,p.nT] = u(m.C[:,p.nT],p)

	# preceding periods
	for it in (p.nT-1):-1:1

		for ia in 1:p.na

			cash = m.avec[ia] 	# current cash on hand

			# Brent's method for minimizing a function
			# withotu derivatives
			x = optimize((x)->VFobj(x,cash,m.V[:,it+1],m,p,m.ywgt),1e-7,cash)
			m.V[ia,it] = -x.f_minimum
			m.C[ia,it] = x.minimum
		end
	end
end
function VFobj(c::Float64,cash::Float64,Vplus::Array{Float64},m::iidModel,p::Param,integw::Vector{Float64})

	# implied cash on hand tomorrow
	s = (cash - c) * p.R .+ m.yvec
	for iy in 1:p.ny
		m.m2[iy] = linearapprox(m.avec,Vplus,s[iy],1,p.na)
	end
	v = u(c,p) + p.beta * dot(m.m2,integw)
	return -v
end



# finding the maximum of the value function backward iteration
# AR1 model

function VFbi!(m::AR1Model,p::Param)

	# final period: consume everything.
	m.C[:,:,p.nT] = repmat(linspace(p.a_low,p.a_high*4,p.na),1,p.ny)
	m.V[:,:,p.nT] = u(m.C[:,:,p.nT],p)

	# preceding periods
	for it in (p.nT-1):-1:1

		# compute conditional expected value function
		# m.EV[:,:,it] = transpose(m.ywgt * transpose(m.V[:,:,it+1]))

		# conditional on current income state
		for iy in 1:p.ny

			for ia in 1:p.na

				cash = m.avec[ia]  # take asset grid as current cash on hand (add random income to next period savings)

				# Brent's method for minimizing a function
				# without derivatives
				x = optimize((x)->VFobj(x,cash,m.V[:,:,it+1],m,p,m.ywgt[iy,:][:]),1e-7,cash)
				m.V[ia,iy,it] = -x.f_minimum
				m.C[ia,iy,it] = x.minimum
			end
		end
	end
end
function VFobj(c::Float64,cash::Float64,Vplus::Matrix{Float64},m::AR1Model,p::Param,integw::Vector{Float64})

	# implied cash on hand tomorrow
	s = (cash - c) * p.R + m.yvec
	# integrate out tomorrow's income uncertainty
	for iy in 1:p.ny
		m.m2[iy] = linearapprox(m.avec,Vplus[:,iy],s[iy],1,p.na)
	end
	v = u(c,p) + p.beta * dot(m.m2,integw)
	return -v
end




function VFbi!(m::AR1Model_a,p::Param)

	# final period: consume everything.
	m.C[:,:,p.nT] = repmat(linspace(p.a_low,p.a_high*4,p.na),1,p.ny)
	m.V[:,:,p.nT] = u(m.C[:,:,p.nT],p)

	# preceding periods
	for it in (p.nT-1):-1:1

		# conditional on current income state
		for iy in 1:p.ny

			# compute conditional expected value function
			m.EV =  m.V[:,:,it+1] * m.ywgt[iy,:][:]  # (1,ny) * (ny,na) = (1,na)

			for ia in 1:p.na

				cash = m.avec[ia] + m.yvec[iy] # current cash on hand

				# Brent's method for minimizing a function
				# withotu derivatives
				x = optimize((x)->VFobj(x,cash,m,p),1e-7,cash)
				m.V[ia,iy,it] = -x.f_minimum
				m.C[ia,iy,it] = x.minimum
			end
		end
	end
end
function VFobj(c::Float64,cash::Float64,m::AR1Model_a,p::Param)

	# implied savings tomorrow
	s = (cash - c) * p.R 
	# get the expected value function at savings s
	EV = linearapprox(m.avec,m.EV,s,1,p.na)
	v = u(c,p) + p.beta * EV
	return -v
end

function linearapprox(x::Vector{Float64},y::Vector{Float64},xi::Float64,lo::Int,hi::Int)
	r = 0.0
	n = length(x)
	@assert n==length(y)

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



