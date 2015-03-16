


# run all
function run()

	# iid income model

	p = Param()

	# solve by standard euler equation root finding
	EEmod = iidModel(p)
	EEmod.toc = @elapsed Euler(EEmod,p)



	# solve maximizing the value function backward iteration
	VFmod = iidModel(p)
	VFmod.toc = @elapsed VFbi(VFmod,p)

	# solve by EGM
	EGMmod = iidModel(p)
	EGMmod.toc = @elapsed EGM(EGMmod,p)

	# plot results
	plots(EEmod,EGMmod,VFmod,p)

	# AR1 income model
	# ================

	EGMmod = AR1Model(p)
	EGMmod.toc = @elapsed EGM(EGMmod,p)

	VFmod = AR1Model(p)
	VFmod.toc = @elapsed VFbi(VFmod,p)

	# does it matter whether I compute the model on 
	# current assets, given y, or
	# cash-on-hand, given y?
	VFmod_a = AR1Model_a(p)
	VFmod_a.toc = @elapsed VFbi(VFmod_a,p)

	# plot results
	plots(EGMmod,VFmod,VFmod_a,p,1)  # plot period 1
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
function EGM(m::iidModel,p::Param)

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



# endogenous grid method for AR1 model
function EGM(m::AR1Model,p::Param)

	# final period: consume everything.
	m.M[:,:,p.nT] = repmat(linspace(p.a_low,p.a_high*4,p.na),1,p.ny)
	m.C[:,:,p.nT] = repmat(linspace(p.a_low,p.a_high*4,p.na),1,p.ny)

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

			# get endogenous grid today
			m.M[:,iy,it] = m.C[:,iy,it] .+ m.avec

		end  # current income
	end  # age
end



# solving the euler equation
function Euler(m::iidModel,p::Param)

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
function VFbi(m::iidModel,p::Param)

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

function VFbi(m::AR1Model,p::Param)

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




function VFbi(m::AR1Model_a,p::Param)

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
function plots(EE::iidModel,EGM::iidModel,VF::iidModel,p::Param)

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
	suptitle("model with iid income uncertainty")

end


function plots(EGM::AR1Model,VF::AR1Model,VF_2::AR1Model_a,p::Param,it::Int)

	figure()

	# plot consumption funcitons
	# ==========================


	# solve by EGM
	subplot(1,3,1)
	plot(squeeze(EGM.M[:,:,it],3),squeeze(EGM.C[:,:,it],3))
	xlabel("cash on hand")
	ylim([p.a_low,p.a_high])
	xlim([p.a_low,p.a_high])
	plot(ylim(),ylim())
	title("Endogenous Grid\n Method: $(round(EGM.toc,5)) secs")

	subplot(1,3,2)
	# plot(VF.avec.+VF.yvec[iy],squeeze(VF.C[:,:,1],3))
	plot(VF.avec,squeeze(VF.C[:,:,it],3))
	# plot(VF.avec,squeeze(VF.C[:,iy,1:(p.nT-1)],2))
	xlabel("cash on hand")
	ylabel("consumption")
	ylim([p.a_low,p.a_high])
	xlim([p.a_low,p.a_high])
	plot(ylim(),ylim())
	title("max V computing expected\n cash-on-hand: $(round(VF.toc,5)) secs")

	subplot(1,3,3)
	plot([VF_2.avec[i] + VF_2.yvec[j] for i=1:p.na, j=1:p.ny],squeeze(VF_2.C[:,:,it],3))
	# plot(VF.avec,squeeze(VF.C[:,iy,1:(p.nT-1)],2))
	xlabel("cash on hand")
	ylabel("consumption")
	ylim([p.a_low,p.a_high])
	xlim([p.a_low,p.a_high])
	plot(ylim(),ylim())
	title("max V computing current\n cash-on-hand: $(round(VF_2.toc,5)) secs")

	# subplot(1,3,2)
	# plot(EE.avec,EE.C[:,iy:(p.nT-1)])
	# xlabel("assets")
	# ylim([p.a_low,p.a_high])
	# xlim([p.a_low,p.a_high])
	# plot(ylim(),ylim())
	# title("Solving the Euler Equation: $(round(EE.toc,5)) secs")


	suptitle("model with AR1 income, all y-states, period $it")

end
