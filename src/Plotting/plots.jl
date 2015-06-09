


function plot_vf(m::iidModel,p::Param)

	# value function
	# 
	figure()
	subplot(1,2,1)
	plot(m.M[:,p.nT],m.V[:,p.nT],label="period $(p.nT)")
	for it in (p.nT-1):-1:1
		plot(m.M[:,it],m.V[:,it],label="period $it")
	end
	xlim([p.a_low,500])
	ylim([-0.7,0.0])
	title(L"Value functions: $V_t(M_t)$")
	xlabel(L"Cash on hand: $M_t$")
	grid("on")
	legend(loc="lower right",fancybox="true")

	# consumption function
	# 
	subplot(1,2,2)
	PyPlot.plot(m.M[:,p.nT],m.C[:,p.nT],label="period $(p.nT)")
	for it in (p.nT-1):-1:1
		PyPlot.plot(m.M[:,it],m.C[:,it],label="period $it")
	end
	ylim([p.a_low,100])
	xlim([p.a_low,100])
	ylim(xlim())
	PyPlot.plot(ylim(),ylim(),linestyle="--",linewidth=2.0,color="black")
	title(L"Consumption functions: $C_t(M_t)$")
	xlabel(L"Cash on hand: $M_t$")
	grid("on")
	# legend(loc="lower right",fancybox="true")

end


# plotting
function plots(d::Dict,p::Param)

	D = Any[]

    jet = ColorMap("jet")[:__call__] # get a color map function


	if haskey(d,"iid")

		# plot consumption functions
		# ==========================

		fig, axes = plt.subplots(1,3,figsize=(10,5))

	    currfig = 0
	    for (k,v) in d["iid"]
	    	currfig += 1
	    	ax = axes[currfig]
			if k=="EGM"
				ax[:plot]([0.0,v.M[:,1]],[0.0,v.C[:,1]],color=jet(0),lw=2,alpha=0.6,label="period 1")
				ax[:plot](v.M[:,p.nT],v.C[:,p.nT],color="black",lw=2,label="period $(p.nT)")
				for it in 2:(p.nT-1)
					x = [0.0,v.M[:,it]]
					ax[:plot](x,[0.0,v.C[:,it]],color=jet(it/(p.nT)),lw=2,alpha=0.6)
				end
			else
				ax[:plot](v.avec,v.C[:,1],color=jet(0),alpha=0.6,lw=2,label="period 1")
				ax[:plot](v.avec,v.C[:,p.nT],color="black",lw=2,label="period $(p.nT)")
				for it in 2:(p.nT-1)
					ax[:plot](v.avec,v.C[:,it],color=jet(it/(p.nT)),lw=2,alpha=0.6)
				end
			end
			if currfig==1
				ax[:set_ylabel]("consumption")
			end
			if k=="EGM"
				ax[:set_xlabel]("cash-on-hand")
			else
				ax[:set_xlabel]("assets")
			end
			ax[:grid]()
			ax[:set_ylim]([0,30])
			ax[:set_xlim]([0,30])
			# plot(ylim(),ylim(),color="black")
			ax[:set_title]("$k: $(round(v.toc,4)) secs")
			ax[:legend](loc="lower right")
		end
		fig[:suptitle]("iid Model Consumption functions by Age")
	

		# plot Value functions
		# ====================
		
		fig1, axes1 = plt.subplots(1,3,figsize=(10,5))

	    currfig = 0
	    for (k,v) in d["iid"]
	    	currfig += 1
	    	ax = axes1[currfig]
			if k=="EGM"
				ax[:plot](v.M[:,1],v.V[:,1],color=jet(0),lw=2,alpha=0.6,label="period 1")
				ax[:plot](v.M[:,p.nT],v.V[:,p.nT],color="black",lw=2,label="period $(p.nT)")
				for it in 2:(p.nT-1)
					x = v.M[:,it]
					ax[:plot](x,v.V[:,it],color=jet(it/(p.nT)),lw=2,alpha=0.6)
				end
			else
				ax[:plot](v.avec,v.V[:,1],color=jet(0),alpha=0.6,lw=2,label="period 1")
				ax[:plot](v.avec,v.V[:,p.nT],color="black",lw=2,label="period $(p.nT)")
				for it in 2:(p.nT-1)
					ax[:plot](v.avec,v.V[:,it],color=jet(it/(p.nT)),lw=2,alpha=0.6)
				end
			end
			if currfig==1
				ax[:set_ylabel]("Value")
			end
			if k=="EGM"
				ax[:set_xlabel]("cash-on-hand")
			else
				ax[:set_xlabel]("assets")
			end
			ax[:set_ylim]([-0.8,0])
			ax[:set_xlim]([0,300])
			# plot(ylim(),ylim(),color="black")
			ax[:grid]()
			ax[:set_title]("$k value functions")
			ax[:legend](loc="lower right")
		end
		fig1[:suptitle]("iid Model Value functions by Age")

		push!(D,fig,fig1)
	end

	if haskey(d,"AR1")

		# AR1 consumption functions in period 1

		it = 1

		fig2, axes2 = plt.subplots(1,3,figsize=(10,5))

	    currfig = 0
	    for (k,v) in d["AR1"]
	    	currfig += 1
	    	ax = axes2[currfig]
			if k=="EGM"
				ax[:plot]([0.0,v.M[:,1,it]],[0.0,v.C[:,1,it]],color=jet(0),lw=2,alpha=0.6,label="income 1")
				ax[:plot](v.M[:,p.ny,it],v.C[:,p.ny,it],color="black",lw=2,label="income $(p.ny)")
				for iy in 2:(p.ny-1)
					x = [0.0,v.M[:,iy,it]]
					ax[:plot](x,[0.0,v.C[:,iy,it]],color=jet(iy/(p.ny)),lw=2,alpha=0.6)
				end
			elseif k=="VF_a"
				ax[:plot](v.avec .+ v.yvec[1] ,v.C[:,1,it],color=jet(0),alpha=0.6,lw=2,label="income 1")
				ax[:plot](v.avec .+ v.yvec[p.ny] ,v.C[:,p.ny,it],color="black",lw=2,label="income $(p.ny)")
				for iy in 2:(p.ny-1)
					ax[:plot](v.avec .+ v.yvec[iy] ,v.C[:,iy,it],color=jet(iy/(p.ny)),lw=2,alpha=0.6)
				end
			else
				ax[:plot](v.avec,v.C[:,1,it],color=jet(0),alpha=0.6,lw=2,label="income 1")
				ax[:plot](v.avec,v.C[:,p.ny,it],color="black",lw=2,label="income $(p.ny)")
				for iy in 2:(p.ny-1)
					ax[:plot](v.avec,v.C[:,iy,it],color=jet(iy/(p.ny)),lw=2,alpha=0.6)
				end
			end
			if currfig==1
				ax[:set_ylabel]("consumption")
			end
			if k=="EGM"
				ax[:set_xlabel]("cash-on-hand")
			else
				ax[:set_xlabel]("assets")
			end
			ax[:grid]()
			ax[:set_ylim]([0,30])
			ax[:set_xlim]([0,60])
			# plot(ylim(),ylim(),color="black")
			ax[:set_title]("$k: $(round(v.toc,4)) secs")
			ax[:legend](loc="upper left")
		end
		fig2[:suptitle]("AR1 Model Consumption functions, age=$it")

		# AR1 value functions in period 1

		fig3, axes3 = plt.subplots(1,3,figsize=(10,5))

	    currfig = 0
	    for (k,v) in d["AR1"]
	    	currfig += 1
	    	ax = axes3[currfig]
			if k=="EGM"
				ax[:plot](v.M[:,1,it],v.V[:,1,it],color=jet(0),lw=2,alpha=0.6,label="income 1")
				ax[:plot](v.M[:,p.ny,it],v.V[:,p.ny,it],color="black",lw=2,label="income $(p.ny)")
				for iy in 2:(p.ny-1)
					x = v.M[:,iy,it]
					ax[:plot](x,v.V[:,iy,it],color=jet(iy/(p.ny)),lw=2,alpha=0.6)
				end
			elseif k=="VF_a"
				ax[:plot](v.avec .+ v.yvec[1] ,v.V[:,1,it],color=jet(0),alpha=0.6,lw=2,label="income 1")
				ax[:plot](v.avec .+ v.yvec[p.ny] ,v.V[:,p.ny,it],color="black",lw=2,label="income $(p.ny)")
				for iy in 2:(p.ny-1)
					ax[:plot](v.avec .+ v.yvec[iy] ,v.V[:,iy,it],color=jet(iy/(p.ny)),lw=2,alpha=0.6)
				end
			else
				ax[:plot](v.avec,v.V[:,1,it],color=jet(0),alpha=0.6,lw=2,label="period 1")
				ax[:plot](v.avec,v.V[:,p.ny,it],color="black",lw=2,label="period $(p.ny)")
				for iy in 2:(p.ny-1)
					ax[:plot](v.avec,v.V[:,iy,it],color=jet(iy/(p.ny)),lw=2,alpha=0.6)
				end
			end
			if currfig==1
				ax[:set_ylabel]("Value")
			end
			if k=="EGM"
				ax[:set_xlabel]("cash-on-hand")
			else
				ax[:set_xlabel]("assets")
			end
			ax[:grid]()
			ax[:set_ylim]([-1,-0.1])
			ax[:set_xlim]([0,300])
			# plot(ylim(),ylim(),color="black")
			ax[:set_title]("$k: $(round(v.toc,4)) secs")
			ax[:legend](loc="lower right")
		end
		fig3[:suptitle]("AR1 Model Value functions, age=$it")
		push!(D,fig2,fig3)

	end

	return D
end

	
function printplots()
	home = ENV["HOME"]
	r = runall()
	p = Param()
	f = plots(r,p)
	figure(f[1][:number])
	savefig("$home/Dropbox/public/ConsProb.jl/iidCons.png")
	figure(f[2][:number])
	savefig("$home/Dropbox/public/ConsProb.jl/iidVfun.png")
	figure(f[3][:number])
	savefig("$home/Dropbox/public/ConsProb.jl/AR1Cons.png")
	figure(f[4][:number])
	savefig("$home/Dropbox/public/ConsProb.jl/AR1Vfun.png")

	d = dchoice()
	figure(d[1][:number])
	savefig("$home/Dropbox/public/ConsProb.jl/Dchoice_condV.png")
	figure(d[2][:number])
	savefig("$home/Dropbox/public/ConsProb.jl/Dchoice_envC.png")
	figure(d[3][:number])
	savefig("$home/Dropbox/public/ConsProb.jl/Dchoice_envV.png")

end


function plots(EGM::AR1Model,VF::AR1Model,VF_2::AR1Model_a,p::Param,it::Int)

	figure("AR1cons",figsize=(9,6))

	# plot consumption funcitons
	# ==========================

	subplot(1,3,1)
	# plot(VF.avec.+VF.yvec[iy],squeeze(VF.C[:,:,1],3))
	plot(VF.avec,squeeze(VF.C[:,:,it],3))
	# plot(VF.avec,squeeze(VF.C[:,iy,1:(p.nT-1)],2))
	ylabel("consumption")
	xlabel("cash on hand")
	ylim([0,30])
	xlim([0,30])
	plot(ylim(),ylim())
	title("max V computing\n expected\n cash-on-hand: $(round(VF.toc,3)) secs")

	subplot(1,3,2)
	plot([VF_2.avec[i] + VF_2.yvec[j] for i=1:p.na, j=1:p.ny],squeeze(VF_2.C[:,:,it],3))
	# plot(VF.avec,squeeze(VF.C[:,iy,1:(p.nT-1)],2))
	xlabel("cash on hand")
	ylim([0,30])
	xlim([0,30])
	plot(ylim(),ylim())
	title("max V computing\n current cash-on-hand: $(round(VF_2.toc,3)) secs")

	# solve by EGM
	subplot(1,3,3)
	plot(squeeze(EGM.M[:,:,it],3),squeeze(EGM.C[:,:,it],3))
	xlabel("cash on hand")
	ylim([0,30])
	xlim([0,30])
	plot(ylim(),ylim())
	title("Endogenous Grid\n Method: $(round(EGM.toc,3)) secs")


	# value function
	figure("AR1vals",figsize=(9,6))
	ax1 = subplot(1,3,3)
	plot(squeeze(EGM.M[:,:,it],3),squeeze(EGM.V[:,:,it],3))
	xlabel("cash on hand")
	ylabel("Value")
	grid("on")
	title("Endogenous Grid\n Method")

	subplot(1,3,2,sharey=ax1)
	plot(VF.avec[2:end],squeeze(VF.V[2:end,:,it],3))
	# plot(VF.avec,squeeze(VF.C[:,iy,1:(p.nT-1)],2))
	xlabel("assets")
	title("max V expected\n cash-on-hand")
	grid("on")

	subplot(1,3,1,sharey=ax1)
	plot([VF_2.avec[i] + VF_2.yvec[j] for i=2:p.na, j=1:p.ny],squeeze(VF_2.V[2:end,:,it],3))
	xlabel("cash on hand")
	title("max V current\n cash-on-hand")
	grid("on")
end


function plots(m::iidDModel,p::Param)
	lwi = 1.5
	d_vf = plt.figure(figsize=(10,7))
	ax1 = plt.subplot(2,3,1)
	for j in 6:-1:1
		if j==6
			ax = ax1
			x = linspace(p.cfloor,cond(m.m,j,2)[2],50)
			ax[:plot](x,map(z -> u(z,true,p) + p.beta * get_vbound(m.v,j,2), x),label="working",color="black",lw=lwi)
			x = linspace(p.cfloor,cond(m.m,j,1)[2],50)
			ax[:plot](x,map(z -> u(z,false,p) + p.beta * get_vbound(m.v,j,1), x),label="retired",color="red",lw=lwi)
		else
			ax = plt.subplot(2,3,6-j+1,sharey=ax1)
			# working
			ax[:plot](cond(m.m,j,2),cond(m.v,j,2),color="black",lw=lwi)
			# analytic part
			x = linspace(p.cfloor,cond(m.m,j,2)[2],50)
			ax[:plot](x,map(z -> u(z,true,p) + p.beta * get_vbound(m.v,j,2), x),label="working",color="black",lw=lwi)

			# retired
			ax[:plot](cond(m.m,j,1),cond(m.v,j,1),label="retired",color="red",lw=lwi)
			# analytic part
			x = linspace(p.cfloor,cond(m.m,j,1)[2],50)
			ax[:plot](x,map(z -> u(z,false,p) + p.beta * get_vbound(m.v,j,1), x),color="red",lw=lwi)
		end
		ax[:set_title]("Period $j")
		ax[:set_xlim]([0,p.a_high])
		ax[:set_ylim]([-10,15])
		ax[:legend](loc="lower right")
		ax[:grid]()
	end
	d_vf[:suptitle]("Value functions for Retirement vs Working")

	# consumption function
	cons = plt.figure(figsize=(7,7))
    jet = ColorMap("jet")[:__call__] # get a color map function
    plot(env(m.m,1),env(m.c,1),color="blue",lw=lwi,label="t=1")
	for it in 2:(p.nT-1)
	    plot(env(m.m,it),env(m.c,it),color=jet(it/(p.nT)),lw=lwi,label="t=$it")
	end
    plot(env(m.m,p.nT),env(m.c,p.nT),color="black",lw=lwi,label="t=$(p.nT)")
	xlim((0,50))
	ylim(xlim())
	title("optimal consumption")
	legend(loc="upper left")
	xlabel("Cash on hand")
	ylabel("Consumption")
	grid()

	# envelope over value functions
	vf = plt.figure(figsize=(7,7))
	for it in 1:(p.nT)
	    plot(env(m.m,it),env(m.v,it),color=jet(it/(p.nT)),lw=lwi,label="t=$it")
		x = linspace(p.cfloor,env(m.m,it)[2],50)
		plot(x,map(z -> u(z,true,p) + p.beta * get_vbound(m.v,it), x),color=jet(it/(p.nT)),lw=lwi)
	end
	xlim((0,50))
	title("Envelopes over discrete Value functions")
	legend(loc="lower right")
	xlabel("Cash on hand")
	ylabel("Value")
	grid()
	return (d_vf,cons, vf)

	
end

