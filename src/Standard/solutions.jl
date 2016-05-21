

function solve!(m::Model,p::Param)

	if m.solver == "EGM"
		EGM!(m,p)
	elseif m.solver == "VFbi"
		VFbi!(m,p)
	elseif m.solver == "Euler"
		Euler!(m,p)
	end
end

# run all
function runStd()

	# return a Dict of Dicts
	D = Dict{AbstractString,Any}()

	# iid model dict
	d = Dict{AbstractString,Model}()

	# iid income models

	p = Param(mu=10.0)

	# solve by standard euler equation root finding
	# d["Euler"] = iidModel(p)
	# # warm up the JIT
	# Euler!(d["Euler"],p)
	# # reset el
	# d["Euler"] = iidModel(p)
	# # measure time
	# d["Euler"].toc = @elapsed Euler!(d["Euler"],p)



	# solve maximizing the value function backward iteration
	d["VF"] = iidModel(p)
	VFbi!(d["VF"],p)
	d["VF"] = iidModel(p)
	d["VF"].toc = @elapsed VFbi!(d["VF"],p)

	# solve by EGM
	d["EGM"] = iidModel(p)
	EGM!(d["EGM"],p)
	d["EGM"] = iidModel(p)
	d["EGM"].toc = @elapsed EGM!(d["EGM"],p)

	D["iid"] = d	

	# plot results
	# plots(EE,EGM,VF,p)

	# AR1 income model
	# ================

	# d2 = Dict{AbstractString,Model}()
	# p = Param(mu=10.0)

	# d2["EGM"] = AR1Model(p)
	# d2["EGM"].toc = @elapsed EGM!(d2["EGM"],p)

	# d2["VF"] = AR1Model(p)
	# d2["VF"].toc = @elapsed VFbi!(d2["VF"],p)

	# # does it matter whether I compute the model on 
	# # current assets, given y, or
	# # cash-on-hand, given y?
	# d2["VF_a"] = AR1Model_a(p)
	# d2["VF_a"].toc = @elapsed VFbi!(d2["VF_a"],p)

	# D["AR1"] = d2

	# # plot results
	# # plots(EGMmod,VFmod,VFmod_a,p,1)  # plot period 1


	# # with debt
	# d3 = Dict{String,Model}()
	# p = Param()
	# d3["EGM"] = iidDebtModel(p)
	# EGM!(d3["EGM"],p)

	# D["iidDebt"] = d3

	return D
end



"""
	findLowestA(m::Model,a::Float64,it::Int)

Find lowest feasible end end-of-period asset in period `it`. 

# Examples
```julia
julia> findLowestA(m,-1.1,8)
-0.3
```
"""
function findLowestA(m::Model,cashn::Float64,p::Param,it::Int)
# take the lowest income state and check if positive consumption
# this assumes that ALL income states are feasible next period
	
	# tmpx = [0.0; m.M[:,it+1] ] 
	tmpx = vb(m.M[it+1])	# get next periods cash on hand withbound
	# you need to basically find a_low such that cons = 0
	# tmpx = [0.0; m.M[:,it+1] ] 
	tmpy = vb(m.C[it+1])	# get next periods consumption without bound
	# println("m = $tmpx")
	# println("c = $tmpy")

	# check lowest cash in hand
	# if implied consumption is negative
	if linearapprox(tmpx,tmpy,cashn) < 0
		# find cashnext such that c=p.cfloor
		cashn = p.cfloor + tmpx[1] - (tmpy[1] * (tmpx[2]-tmpx[1])) / (tmpy[2] - tmpy[1])
	end
	# return implied level of end-of-period asset that makes next period cash in hand = 0
	return invcashnext(cashn,income(m.yvec[1],it+1,p),p)	# find corresponding a level
end

"""
	findLowestA(mod::Model,m::Vector{Float64},c::Vector{Float64},a::Float64,it::Int)

Find lowest feasible end end-of-period asset in period `it`. 

# Examples
```julia
julia> findLowestA(m,-1.1,8)
-0.3
```
"""
function findLowestA(p::Param,m::Vector{Float64},c::Vector{Float64},cashn::Float64,ylow::Float64,it::Int)
# take the lowest income state and check if positive consumption
# this assumes that ALL income states are feasible next period
	
	# tmpx = [0.0; m.M[:,it+1] ] 
	# tmpx = vb(m.M[it+1])	# get next periods cash on hand withbound
	# you need to basically find a_low such that cons = 0
	# tmpx = [0.0; m.M[:,it+1] ] 
	# tmpy = vb(m.C[it+1])	# get next periods consumption without bound
	# println("m = $tmpx")
	# println("c = $tmpy")

	# check lowest cash in hand
	# if implied consumption is negative
	if linearapprox(m,c,cashn) < 0
		# find cashnext such that c=p.cfloor
		cashn = p.cfloor + m[1] - (c[1] * (m[2]-m[1])) / (c[2] - c[1])
	end
	# return implied level of end-of-period asset that makes next period cash in hand = 0
	return invcashnext(cashn,ylow,p)	# find corresponding a level
end

"""
	fillAvec!(m::iidModel,p::Param,it::Int)

Fill end-of-period asset vector avec with feasible sequence in period `it`.

# Examples
```julia
julia> fillAvec!(m,3)
```
"""
function fillAvec!(m::iidModel,p::Param,it::Int)
	a = 0.0
	# if lower bound of avec[it] < 0, need to find lowest feasible a.
	x = cashnext(m.avec[it][1],income(m.yvec[1],it+1,p),p)	
	if x < a
		# fill with grid from a to a_upper
		a = findLowestA(p,vb(m.M[it+1]),vb(m.C[it+1]),x,income(m.yvec[1],it+1,p),it)
		m.avec[it] = scaleGrid(a,p.a_high,p.na)
		# println("fillAvec scaling: $(m.avec[it])")
	else
		# don't do anything: all set
	end
end

"""
	fillAvec!(m::AR1Model,p::Param,iy::Int,it::Int)

Fill end-of-period asset vector avec with feasible sequence in period `it`.

# Examples
```julia
julia> fillAvec!(m,3)
```
"""
function fillAvec!(m::AR1Model,p::Param,iy::Int,it::Int)
	a = 0.0
	# if lower bound of avec[it] < 0, need to find lowest feasible a.
	x = cashnext(m.avec[it][1],income(m.yvec[1],it+1,p),p)	
	if x < a
		# fill with grid from a to a_upper
		a = findLowestA(p,vb(m.M[iy,it+1]),vb(m.C[iy,it+1]),x,income(m.yvec[1],it+1,p),it)
		m.avec[it] = scaleGrid(a,p.a_high,p.na)
		# println("fillAvec scaling: $(m.avec[it])")
	else
		# don't do anything: all set
	end
end

function fillCnext!(m::iidModel,p::Param,it::Int)
	idx = 0
	tmpx = vb(m.M[it+1])	# get next periods cash on hand with bound
	tmpy = vb(m.C[it+1])	# get next periods consumption 
	for ia in 1:p.na
		for iy in 1:p.ny
			idx = ia+p.na*(iy-1)
			midx = ia+p.na*(iy-1 + p.ny * (it+1-1))   # next period!
			# m.mnext[idx] = cashnext(m.avec[it][ia],income(m.yvec[iy],it+1,p),p)
			m.cnext[idx] = linearapprox(tmpx,tmpy,m.mnext[midx])
		end
	end
end

"""
	RHS(m::iidModel,p::Param,it::Int)

Compute the RHS of the Euler Equation.

# Details

Computes the right hand side of the Eurler Equation in the `iidModel` case. 
"""
function RHS(m::iidModel,p::Param,it::Int)
	p.R * p.beta .* up(m.cnext,p) * m.ywgt
end


"""
	EGM!(m::iidModel,p::Param)

Compute EGM solution for an `iidModel`.

"""
function EGM!(m::iidModel,p::Param)

	it = p.nT
	# initiate consumption function on arbitrary 2 (positive) points of m
	# final period: consume everything.
	# C is linear in last period.
	set!(m.M[it],[p.a_high/2;p.a_high])	
	set!(m.C[it],[p.a_high/2;p.a_high])	
	set!(m.V[it],u(v(m.C[it]),p) + p.beta * 0.0)   # future value in last period: 0.0

	# bounds
	set_bound!(m.M[it],m.avec[p.nT][1])	# here you decide whether one can die in debt or not
	# set_bound!(m.M[it],m.avec[it][1]*p.R)	# here you decide whether one can die in debt or not
	set_bound!(m.C[it],p.cfloor)
	set_bound!(m.V[it],0.0)

	# preceding periods
	for it in (p.nT-1):-1:1
		# fillAvec!(m,p,it)
		# get all future consumptions and compute RHS of euler equation
		fillCnext!(m,p,it)
		Eu = RHS(m,p,it)
		# get optimal consumption today from euler equation: invert marginal utility
		set!(m.C[it],iup(Eu,p))

		# get endogenous grid today
		set!(m.M[it],v(m.C[it]) + m.avec[it])
		# set bound on m: which m corresponds to zero consumption?
		set_bound!(m.M[it],m.avec[it][1])

		# compute value function
		computeVfun!(m,p,it)
	end
end

"""
	computeVfun!(m::iidModel,p::Param,it::Int)

Compute Value function in EGM step.

"""
function computeVfun!(m::iidModel,p::Param,it::Int)
	# expected value function (na,ny)
	fill!(m.ev,NaN)
	# dont: don't interpolate anything.
	if it==(p.nT-1)
		dont = trues(size(m.mnext))
	else
		dont = m.mnext .< v(m.M[it+1])[1]	# wherever potential next period's cash on hand (m.mnext) is less than the lowest grid point of the endogenous grid next period (m.M), the agent will be credit constrained and will be saving zero (with value b(m.V[it+1]))
	end

	vv = v(m.V[it+1])
	tmpx = v(m.M[it+1])
	for ia in 1:p.na
		for iy in 1:p.ny
			idx = ia+p.na*(iy-1)
			if dont[idx]
				m.ev[idx] = u(m.cnext[idx],p) + p.beta * b(m.V[it+1])
			else
				m.ev[idx] = linearapprox(tmpx,vv,m.mnext[idx])
			end
		end
	end
	ev = m.ev * m.ywgt
	# if abs(m.avec[1]) > 1e-6
	# 	error("first element of avec is assumed to be zero: it's not!")
	# end
	set_bound!(m.V[it],ev[1])
	set!(m.V[it],u(v(m.C[it]),p) + p.beta * ev)
end

# endogenous grid method for AR1 model
function EGM!(m::AR1Model,p::Param)

	# final period: consume everything.
	it = p.nT
	set!(m.M[:,it],[p.a_high/2;p.a_high])	
	set!(m.C[:,it],[p.a_high/2;p.a_high])	
	set!(m.V[:,it],u(v(m.C[it][1]),p) + p.beta * 0.0)   # future value in last period: 0.0

	# bounds
	set_bound!(m.M[:,it],0.0)	# here you decide whether one can die in debt or not
	set_bound!(m.C[:,it],p.cfloor)
	set_bound!(m.V[:,it],0.0)

	# preceding periods
	for it in (p.nT-1):-1:1

		# conditional on current income state
		for iy in 1:p.ny

			# interpolate optimal consumption from next period on all cash-on-hand states
			# using C[:,it+1] and M[:,it+1], find c(m,it)

			# next period's income index
			for iiy in 1:p.ny
				fillAvec!(m,p,iy,it)
				# find lowest feasible asset if allow for negative assets
				# findLowestA(p,vb(m.M[:,iiy,it+1]),vb(m.C[:,iiy,it+1]),cashn,ylow,it)
				tmpx = [0.0; m.M[:,iiy,it+1] ] 
				tmpy = [0.0; m.C[:,iiy,it+1] ]
				for ia in 1:p.na
					m.cnext[ia+p.na*(iiy-1)] = linearapprox(tmpx,tmpy,m.mnext[ia+p.na*(iiy-1)],1,p.na)
					# m.cnext[ia+p.na*(iiy-1)] = linearapprox(tmpx,tmpy,m.mnext[ia+p.na*(iy-1)],1,p.na)
				end
			end

			# get expected marginal value of saving: RHS of euler equation
			# beta * R * E[ u'(c_{t+1}) | y_t] 
			Eu = p.R * p.beta .* m.ywgt[iy,:] * transpose(up(m.cnext,p))

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
				dont = trues(size(m.mnext))
			else
				dont = m.mnext .< m.M[1,iy,it+1]	# wherever potential next period's cash on hand (m.mnext) is less than the lowest grid point of the endogenous grid next period (m.M), the agent will be credit constrained and will be saving zero (here stored in m.Vzero[iy,it+1])
			end

			# for next period's income state
			for iiy in 1:p.ny
				tmpx = m.M[:,iiy,it+1]  
				vv   = m.V[:,iiy,it+1]
				for ia in 1:p.na
					idx = ia+p.na*(iiy-1)
					# idx = ia+p.na*(iy-1)
					if dont[idx]
						m.ev[idx] = u(m.mnext[idx],p) + p.beta * m.Vzero[iiy,it+1]
					else
						m.ev[idx] = linearapprox(tmpx,vv,m.mnext[idx],1,p.na)
					end
				end
			end
			# ev = transpose(m.ywgt[iy,:] * transpose( m.ev ))
			ev = m.ev * m.ywgt[iy,:][:] 
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
	m.M[:,p.nT] = m.avec 
	m.C[:,p.nT] = m.avec
	m.C[m.C[:,p.nT].<p.cfloor,p.nT] = p.cfloor
	m.V[:,p.nT] = u(m.C[:,p.nT],p) + p.beta * 0.0

	# preceding periods
	for it in (p.nT-1):-1:1

		for ia in 1:p.na

			cash = m.avec[ia] 	# current cash on hand

			# consumption equal to cash on hand
			res = EulerResid(cash,cash,m.C[:,it+1],p,m,it)

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
			# u'(c_t) > u'(c_{t+1}), or (by c_t = m_t)
			# u'(m_t) > u'(c_{t+1}), and therefore
			# u'(c_t) = max[ u'(m_t), beta * R * u'(c_{t+1}) ] implies
			# that this consumer is borrowing constrained and consumes all cash in hand.
			if res < 0
				m.C[ia,it] = cash
			else
				# m.C[ia,it] = fzero((x)->EulerResid(x,cash,m.C[:,it+1],p,m,it),(cash + p.a_low-0.0001)/2,[p.a_low-0.0001,cash])
				# m.C[ia,it] = fzero((x)->EulerResid(x,cash,m.C[:,it+1],p,m,it),cash/2,[p.a_low,cash])
				m.C[ia,it] = fzero((x)->EulerResid(x,cash,m.C[:,it+1],p,m,it),[p.a_low,cash])
			end
			m.S[ia,it] = (cash - m.C[ia,it])*p.R

			# get expected value function
			EV = linearapprox(m.avec,m.V[:,it+1],m.S[ia,it].+ m.yvec)

			m.V[ia,it] = u(m.C[ia,it],p) + p.beta * dot(EV,m.ywgt)

		end
		ev



	end

end

# Euler Equation Residual
function EulerResid(c::Float64,cash::Float64,cplus::Vector{Float64},p::Param,m::iidModel,it::Int)

	# given current c, what is next period's cash on hand

	# if it == (p.nT-1)
	# 	# next period is last: no income!
	# 	Euc = p.R .* p.beta .* up(p.R * (cash - c),p)
	# else
		m.m2 = p.R * (cash - c) .+ m.yvec  # (ny,1)
		# interpolate optimal consumption c(t+1), given c(t), on each y-state
		m.c2 = linearapprox([0;m.avec],[0;cplus],m.m2)
		# for iy in 1:p.ny
		# 	# m.c2[iy] = linearapprox([0,m.avec],[0,cplus],m.m2[iy],1,p.na)
		# 	m.c2[iy] = linearapprox(m.avec,cplus,m.m2[iy])
		# end

		# Expected marginal utility of consumption (RHS of euler equation)
		Euc = dot(p.R .* p.beta .* up(m.c2,p), m.ywgt) 	# (1,1)
	# end

	# residual
	c - iup(Euc,p)

end



# finding the maximum of the value function backward iteration
function VFbi!(m::iidModel,p::Param)

	# final period: consume everything.
	m.C[:,p.nT] = m.avec
	m.C[m.C[:,p.nT].<p.cfloor,p.nT] = p.cfloor
	m.V[:,p.nT] = u(m.C[:,p.nT],p)

	# preceding periods
	for it in (p.nT-1):-1:1

		for ia in 1:p.na

			cash = m.avec[ia] 	# current cash on hand

			# Brent's method for minimizing a function
			# withotu derivatives
			x = optimize((x)->VFobj(x,cash,m.V[:,it+1],m,p,m.ywgt,it),p.a_low-100*eps(),cash)
			m.V[ia,it] = -x.f_minimum
			m.C[ia,it] = x.minimum
		end
	end
end
function VFobj(c::Float64,cash::Float64,Vplus::Array{Float64},m::iidModel,p::Param,integw::Vector{Float64},it::Int)

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
	x = m.avec
	x[x.<p.cfloor] = p.cfloor
	m.C[:,:,p.nT] = repmat(x,1,p.ny)
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
				x = optimize((x)->VFobj(x,cash,m.V[:,:,it+1],m,p,m.ywgt[iy,:][:]),p.a_low-100*eps(),cash)
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
	x = m.avec
	x[x.<p.cfloor] = p.cfloor
	m.C[:,:,p.nT] = repmat(x,1,p.ny)
	m.V[:,:,p.nT] = u(m.C[:,:,p.nT],p)

	# preceding periods
	for it in (p.nT-1):-1:1

		# conditional on current income state
		for iy in 1:p.ny

			# compute conditional expected value function
			m.EV =  m.V[:,:,it+1] * m.ywgt[iy,:][:]  # (na,ny) * (ny,1) = (na,	1)

			for ia in 1:p.na

				cash = m.avec[ia] + m.yvec[iy] # current cash on hand

				# Brent's method for minimizing a function
				# withotu derivatives
				x = optimize((x)->VFobj(x,cash,m,p),p.a_low-100*eps(),cash)
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
	EV = linearapprox(m.avec,m.EV,s)
	v = u(c,p) + p.beta * EV
	return -v
end



