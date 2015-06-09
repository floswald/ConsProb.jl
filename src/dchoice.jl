

"iid income uncertainty disrete choice model
solved with EGM."
function EGM!(m::iidDModel,p::Param)

	@assert m.avec[1] == 0.0

	# loop over time periods
	for it in p.nT:-1:1

		println("period: $it of $(p.nT)")

		# compute optimal consumption and values conditional on discrete choice
		d_EGM!(m,p,it)

		# compute envelopes over conditional functions
		d_env!(m,p,it)

	end

end


function d_EGM!(m::iidDModel,p::Param,it::Int)
	for id in 1:p.nD
		println("id: $id")
		working = convert(Bool,id - 1)

		if it==p.nT
			# final period: consume everyting.
			set!(m.m,it,id,vcat(0.0,p.a_high))
			set!(m.c,it,id,vcat(0.0,p.a_high))
			set!(m.v,it,id,vcat(0.0, NaN))
			# v0 = u(c0,working,p) + p.beta * 0.0
		else
			# previous periods
			# first elements: 0 and NaN
			push!(m.c[it].cond[id],0.0)
			push!(m.m[it].cond[id],0.0)
			push!(m.v[it].cond[id],NaN)  # element # 1 of value function is value of saving zero

			# precomputed next period's cash on hand on all income states
			# what's next period's cash on hand given you work/not today?
			mm1 = m.m1[it+1][id]

			# next period's endogenous grid and cons function
			# on the envelope!
			# this must include the lower bound on both m and c
			m1 = env(m.m,it+1)
			c1 = env(m.c,it+1) 

			# prepend next period's optimal policies with a zero to capture credit constraint
			for ia in 1:p.na
				for iy in 1:p.ny
					tmpc = linearapprox(m1,c1,mm1[ia+p.na*(iy-1)])
					m.c1[ia+p.na*(iy-1)] = max(tmpc,p.cfloor)
				end
			end

			# get expected marginal value of saving: RHS of euler equation
			# beta * R * E[ u'(c_{t+1}) ] 
			Eu = up(m.c1,p) * m.ywgt
			rhs = p.R * p.beta .* Eu

			# get optimal consumption today from euler equation: invert marginal utility
			append!(m.c[it].cond[id],iup(Eu,p))
			# get endogenous grid today
			append!(m.m[it].cond[id],cond(m.c,it,id)[2:end] .+ m.avec)


			# compute value function
			# ======================

			vv   = env(m.v,it+1)[2:end]
			tmpx = env(m.m,it+1)[2:end]

			# expected value function (na,ny)
			fill!(m.ev,NaN)
			# mask off values where we don't have to interpolate
			if it==(p.nT-1)
				mask = trues(size(mm1))
			else
				mask = mm1 .< env(m.m,it+1)[2]	# wherever potential next period's cash on hand (m.m1) is less than the second lowest grid point of the endogenous grid next period (lowest is 0), the agent will be credit constrained and will be saving zero 
				# in that region, can use analytic form of v(work)
			end

			# again, compute EV directly for the credit constrained cases.
			# remember that EV is the envelope over 2 conditional vfuns

			for ia in 1:p.na
				for iy in 1:p.ny
					idx = ia+p.na*(iy-1)
					# if credit constrained and next period is last: will choose retire and zero savings!
					if it==(p.nT-1)  # retired next period
						m.ev[idx] = u(max(mm1[idx],p.cfloor),false,p) + p.beta * vzero(m.v,it+1)
					else
						# if credit constrained and next period is not last: will choose work and zero savings!
						if mask[idx]
							m.ev[idx] = u(max(mm1[idx],p.cfloor),true,p) + p.beta * vzero(m.v,it+1)
						else
						# if not credit constrained: will save!
							m.ev[idx] = linearapprox(tmpx,vv,mm1[idx])
						end
					end
				end
			end
			# remember that mm1 is next period's cash on hand if TODAY discrete choice is id (e.g. work)
			ev = m.ev * m.ywgt
			# set vzero
			m.v[it].cond[id] = [ev[1]]
			# set rest of grid points
			append!(m.v[it].cond[id], u(cond(m.c,it,id)[2:end],working,p) + p.beta * ev)

			# however, at this point there may be a problem in m.dpolicy[it][id]["v"] because of secondary kinks coming in from m.envelope[it+1]["v"]. Will now find wiggles where grid folds back onto itself.

			if p.dorefinements

				refine_grids!(m,it,id,p)
				
			end # if refinements
		end  # if final period
	end  # loop over discrete choice
end

function d_env!(m::iidDModel, p::Param, it::Int)

	if it == p.nT
		# in last period it's optimal to retire.
		# and consume everything
		set!(m.m,it,vcat(0.0,p.a_high))
		set!(m.c,it,vcat(0.0,p.a_high))
		set!(m.v,it,vcat(0.0 , NaN))

	else

		# the endogenous grids in for each discrete choice are different.
		# we must compute an upper envelope over both value functions
		# in general, the endogenous grids are different for each discrete choice.

		# in this example, there is only one intersection, and we
		# know this cannot lie in the credit-contrained region, so 
		# we can use this knowledge. We will compute the analytic form 
		# of the value functions in the credit constrained region rather
		# than approximating them.

		# 1) points from RETIRED grid that should appear in envelope
		# ------------------------------------------------

		# `mask` covers up grid points where points of RETIRED grid 
		# fall into the credit constraint of working
		# in that region, we can use the analytic form of v(work)
		mask = cond(m.m,it,1) .< cond(m.m,it,2)[2]
		mask[1] = false 	# skip the first point (with special value vf(1)

		# mask1: indicator for where retirement is optimal
		# add indicators in credit constrained region
		# here is the retirment value better than work
		# v(work) is analytic here
		mask1 = vcat(false, cond(m.v,it,1)[mask] .> u(max(p.cfloor,cond(m.m,it,1)[mask]),true,p) + p.beta * cond(m.v,it,2)[1] )
		# remember that  cond(m.v,it,2)[1] is EV( saving zero )

		# invert mask, i.e. where RETIRED is NOT in credit constrained region of work
		# and thus we have to approximate v(work)
		mask = !mask
		mask[1] = false # skip the first point (with special value vf(1)

		# where is retiring better than work outside of the credit constrained region? must interpolate.
		mask1 = vcat( mask1, cond(m.v,it,1)[mask] .> linearapprox(cond(m.m,it,2)[2:end],cond(m.v,it,2)[2:end], cond(m.m,it,1)[mask]) )

		# 2) points from WORKING grid that should appear in envelope
		# ----------------------------------------------------------

		# where is WORKING better than RETIRE?
		# could use same approach as above. but credit constrained region for "retire"
		# is very small, so not worth the effort.
		# notice that we interpolate RETIRE on it's endo grid, so that we can get it's values at the endo grid of WORK - i.e. the points that v(work) is defined on!

		# take out first point with special value
		mask2 = vcat( false, cond(m.v,it,2)[2:end] .> linearapprox(cond(m.m,it,1)[2:end],cond(m.v,it,1)[2:end],cond(m.m,it,2)[2:end])  )

		
		# 3) find intersection points
		# ---------------------------

		ix1 = find(mask1)
		ix2 = find(mask2)   # left poitn of working grid before intersection
		if length(ix1) > 0
			i1 = ix1[1] -1
		end
		if length(ix2)>0 
			i2 = ix2[end]
			if i2==length(mask2)
				i2 -= 1
			end
		end

		# if intersection lies in credit constraint of working: retirement never optimal
		if length(ix1)==0
			isectm = Float64[]
			isectc = Float64[]
			isectv = Float64[]
		else
			if length(ix2) == 0
				# intersection in credit constraint of WORKING
				# use analytic forms
				func = (x)->u(x,true,p) + p.beta * vzero(m.v,it,2) - linearapprox(cond(m.m,it,1)[2:end],cond(m.v,it,1)[2:end],x)
				isectm = fzero(func,[p.cfloor,cond(m.m,it,2)[end]])
				# isectm = fzero(func,cond(m.m,it,2)[2])
				# isectm = fzero(func,m.dpolicy[it][2][1])
				# why do you use two points here?
				isectm = vcat(isectm, isectm + 100.0*eps())
				isectv = vcat(u(max(p.cfloor,isectm[1]),true,p) + p.beta*vzero(m.v,it,2), 
					          u(max(p.cfloor,isectm[2]),true,p) + p.beta*vzero(m.v,it,2))
			else
				# just get the linear segment connecting 
				a1 = (cond(m.v,it,1)[i1+1] - cond(m.v,it,1)[i1]) / (cond(m.m,it,1)[i1+1] - cond(m.m,it,1)[i1])
				a2 = (cond(m.v,it,2)[i2+1] - cond(m.v,it,2)[i2]) / (cond(m.m,it,2)[i2+1] - cond(m.m,it,2)[i2])
				b1 = cond(m.v,it,1)[i1] - a1*cond(m.m,it,1)[i1]
				b2 = cond(m.v,it,2)[i2] - a2*cond(m.m,it,2)[i2]

				isectm = vcat((b2-b1)/(a1-a2), 100.0*eps() + (b2-b1)/(a1-a2))

				isectv = vcat(a2*isectm[1] + b2, a1*isectm[2] + b1)

			end

			# consumption function from both rules
			isectc = vcat( linearapprox(cond(m.m,it,2),cond(m.c,it,2),isectm[1]), linearapprox(cond(m.m,it,1),cond(m.c,it,1),isectm[2]) )
			if p.printdebug
				println("isectv = $isectv")
				println("isectm = $isectm")
				println("isectc = $isectc")
			end
		end

		# combine points from both grids into envelope
		set!(m.m,it,vcat(0.0, cond(m.m,it,2)[mask2], isectm, cond(m.m,it,1)[mask1]))
		set!(m.c,it,vcat(0.0, cond(m.c,it,2)[mask2], isectc, cond(m.c,it,1)[mask1]))
		set!(m.v,it,vcat(vzero(m.v,it,2), cond(m.v,it,2)[mask2], isectv, cond(m.v,it,1)[mask1]))
	end
end


function refine_grids!(m::iidDModel,it::Int,id::Int,p::Param)
	j = find(cond(m.m,it,id)[2:end] .< cond(m.m,it,id)[1:(end-1)])
	# j[1] is last ok element: M_b

	if p.printdebug
		println("j = $j")
	end
	while length(j) > 0
		j1 = find(cond(m.m,it,id) .< cond(m.m,it,id)[j[1]]) # which grid points are lower than that last ok one
		j1 = sum(j1.>j[1])  # number of indices greater than j[1], i.e. gridpoints that SHOULD be greater than m0[j[1]]

		# j1 is the number of points we need to correct

		# interpolate all points inside the fold using part of the grid we know is ok
		# get interpolated value from non-folded grid
		folds_on_ok = linearapprox(cond(m.m,it,id)[1:j[1]],cond(m.v,it,id)[1:j[1]],cond(m.m,it,id)[(j[1]+1):(j[1]+j1)])

		# number of points where value is below true value on the fold
		j2 = sum(cond(m.v,it,id)[(j[1]+1):(j[1]+j1)] .< folds_on_ok)
		# j2 is the number of points we want to drop

		j3 = find(cond(m.m,it,id) .> cond(m.m,it,id)[j[1]+j2+1])	 # indices of gridpoints greater than the point where values are lower than points on v
		j3 = sum(j3 .<= j[1]) # number of those lower than M_b
		#

		if p.printdebug
			println("env2 (it=$it,id=$id), point $j: m=$(cond(m.m,it,id)[j[1]]), v=$(cond(m.v,it,id)[j[1]]).")
			println("# $j1 in fold over")
			println("# $j2 of points to be skipped after start of fold over")
			println("# $j3 of points to be skipped before start of fold over")
		end

		# delete inferior points
		deleteat!(m.m[it].cond[id],(j[1]-j3+1):(j[1]+j2))
		deleteat!(m.v[it].cond[id],(j[1]-j3+1):(j[1]+j2))
		deleteat!(m.c[it].cond[id],(j[1]-j3+1):(j[1]+j2))

		# search for next fold over region
		j = find(cond(m.m,it,id)[2:end] .< cond(m.m,it,id)[1:(end-1)])  
	end
end
