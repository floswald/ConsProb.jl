


"iid income uncertainty disrete choice model
solved with EGM."
function EGM!(m::iidDModel,p::Param)

	@assert m.avec[1] == 0.0

	# loop over time periods
	for it in p.nT:-1:1

		# loop over dchoices
		for id in 1:p.nD

			working = convert(Bool,id - 1)

			# get out today's d-choice conditional optimal policies
			# from the m object
			m0 = m.dpolicy[it][id]["m"]
			c0 = m.dpolicy[it][id]["c"]
			v0 = m.dpolicy[it][id]["v"]

			if it==p.nT
				# final period: consume everyting.
				m0 = linspace(m.avec[1],p.a_high*4,p.na)
				c0 = linspace(m.avec[1],p.a_high*4,p.na)
				c0[c0.<p.cfloor] = p.cfloor
				v0 = u(c0,working,p) + p.beta * 0.0
			else
				# previous periods
				# interpolate optimal consumption from next period on all cash-on-hand states
				# using C[:,it+1] and M[:,it+1], find c(m,it)

				# precomputed next period's cash on hand on all income states
				mm1 = m.m1[it+1][id]["m1"]

				# next period's endogenous grid and cons function
				m1 = m.dpolicy[it+1][id]["m"]
				c1 = m.dpolicy[it+1][id]["c"]

				tmpx = [0.0, m1 ] 
				tmpy = [0.0, c1 ]
				for ia in 1:p.na
					for iy in 1:p.ny
						# fast:
						# m.c1[ia+p.na*(iy-1+p.ny*(it-1+p.nT*(id-1)))] = linearapprox(tmpx,tmpy,m.m1[ia+p.na*(iy-1+p.ny*(it-1+p.nT*(id-1)))],1,p.na)
						# slow:
						tmpc = linearapprox(tmpx,tmpy,mm1[ia+p.na*(iy-1)],1,p.na)

						# m.c1 is next period's cons at all possible income states
						m.c1[ia+p.na*(iy-1)] = tmpc > p.cfloor ? tmpc : p.cfloor
					end
				end

				# get expected marginal value of saving: RHS of euler equation
				# beta * R * E[ u'(c_{t+1}) ] 
				Eu = p.R * p.beta .* up(m.c1,p) * m.ywgt

				# get optimal consumption today from euler equation: invert marginal utility
				c0 = iup(Eu,p)

				# get endogenous grid today
				m0 = c0 .+ m.avec

				# compute value function
				# ======================

				vv = m.envelope[it+1]["v"]
				tmpx = m.envelope[it+1]["m"]  

				# expected value function (na,ny)
				fill!(m.ev,NaN)
				# dont: don't interpolate anything.
				if it==(p.nT-1)
					dont = trues(size(mm1))
				else
					dont = m1 .< m.envelope[it+1]["m"][1]	# wherever potential next period's cash on hand (m.m1) is less than the lowest grid point of the endogenous grid next period, the agent will be credit constrained and will be saving zero 
				end

				# again, compute EV directly for the credit constrained cases.
				# remember that EV is the envelope over 2 conditional vfuns

				for ia in 1:p.na
					for iy in 1:p.ny
						idx = ia+p.na*(iy-1)
						# if credit constrained and next period is last: will choose retire and zero savings!
						if dont[idx] && it==(p.nT-1)
							m.ev[idx] = u(m1[idx],false,p) + p.beta * m.envelope[it+1]["Vzero"]

						# if credit constrained and next period is not last: will choose work and zero savings!
						elseif dont[idx] && it!=(p.nT-1)
							m.ev[idx] = u(m1[idx],true,p) + p.beta * m.envelope[it+1]["Vzero"]
						else
						# if not credit constrained: will save!
							m.ev[idx] = linearapprox(tmpx,vv,mm1[idx],1,p.na)
						end
					end
				end
				# remember that mm1 is next period's cash on hand if TODAY discrete choice is id (e.g. work)
				ev = m.ev * m.ywgt
				if abs(m.avec[1]) < eps()
					m.dpolicy[it][id]["Vzero"] = ev[1]
				end
				v0  = u(c0,working,p) + p.beta * ev 

				# however, at this point there may be a problem in m.dpolicy[it][id]["v"] because of secondary kinks coming in from m.envelope[it+1]["v"]. Will now find wiggles where grid folds back onto itself.

				if p.dorefinements

					# alternative: use issorted(m0)

					# indices of non-ordered elements
					j = find(m0[2:end] .< m0[1:(end-1)])   # TODO implement as loop
					# j[1] is last ok element: M_b


					while length(j) > 0
						j1 = find(m0 .< m0[j[1]]) # which grid points are lower than that last ok one
						j1 = sum(j1.>j[1])  # number of indices greater than j[1], i.e. gridpoints that SHOULD be greater than m0[j[1]]

						# j1 is the number of points we need to correct

						# interpolate all points inside the fold using part of the grid we know is ok
						# get interpolated value from non-folded grid
						fold_on_ok = linearapprox(m0[1:j],v0[1:j],m0[(j+1):(j+j1)])

						# number of points where value is below true value on the fold
						j2 = sum(v0[(j+1):(j+j1)] .< fold_on_ok)
						# j2 is the number of points we want to drop

						j3 = find(m0 .> m0[j+j2+1])	 # indices of gridpoints greater than the point where values are lower than points on v
						j3 = sum(j3 .<= j[1]) # number of those lower than M_b
						#

						# delete inferior points
						deleteat!(m0,(j-j3+1):(j+j2))
						deleteat!(v0,(j-j3+1):(j+j2))
						deleteat!(c0,(j-j3+1):(j+j2))

						# search for next fold over region
						j = find(m0[2:end] .< m0[1:(end-1)])  
					end
					# 
				end # if refinements
			end  # if final period
		end  # loop over discrete choice

		# compute envelopes: compare discrete choice vfuns
		# ------------------------------------------------
		if it == p.nT
			# in last period it's optimal to retire.
			# and consume everything
			m.envelope[it]["m"] = linspace(m.avec[1],p.a_high*4,p.na)
			m.envelope[it]["c"] = linspace(m.avec[1],p.a_high*4,p.na)

		else

			# the endogenous grids in m.dpolicy[it][1]["m"] and 
			# m.dpolicy[it][1]["m"] are different.
			# we must compute an upper envelope over both value functions

			# in this example, there is only one intersection, and we
			# know this cannot lie in the credit-contrained region

			# 1) points from RETIRED grid that should appear in envelope
			# ------------------------------------------------

			# mark points of RETIRED grid that fall into the credit constraint
			# region of WORKING. 
			mark = m.dpolicy[it][1]["m"] .< m.dpolicy[it][2]["m"][1]

			# add indicators in credit constrained region
			# this is fedor's version:
			# mark1 = m.dpolicy[it][1]["v"][mark] .> u(max(p.cfloor,m.dpolicy[it][1]["m"][mark]),false,p) + p.beta * m.dpolicy[it+1][2]["Vzero"]
			# I think he meant that: where is the retirment value better than work
			mark1 = m.dpolicy[it][1]["v"][mark] .> u(max(p.cfloor,m.dpolicy[it][2]["m"][mark]),true,p) + p.beta * m.dpolicy[it+1][2]["Vzero"]

			# invert mark, i.e. where RETIRED is NOT in credit constrained region of work
			mark = !mark

			# where is retiring better than work outside of the credit constrained region?
			mark1 = [mark1, m.dpolicy[it][1]["v"][mark] .> linearapprox(m.dpolicy[it][2]["m"][mark],m.dpolicy[it][2]["v"][mark], m.dpolicy[it][1]["m"][mark])]

			# 2) points from WORKING grid that should appear in envelope
			# ------------------------------------------------

		end

	end
end


