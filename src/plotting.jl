


function plot_vf(m::iidModel,p::Param)

	# value function
	# 
	figure()
	subplot(1,2,1)
	PyPlot.plot(m.M[:,p.nT],m.V[:,p.nT],label="period $(p.nT)")
	for it in (p.nT-1):-1:1
		PyPlot.plot(m.M[:,it],m.V[:,it],label="period $it")
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
function plots(EE::iidModel,EGM::iidModel,VF::iidModel,p::Param)

	figure("iidCons",figsize=(9,6))

	# plot consumption funcitons
	# ==========================

	# solve by maximizing the VF
	subplot(1,3,1)
	plot(VF.avec,VF.C[:,1:(p.nT-1)])
	xlabel("assets")
	ylabel("consumption")
	ylim([5,20])
	xlim([5,20])
	plot(ylim(),ylim())
	title("maximize\n value function:\n $(round(VF.toc,3)) secs")

	subplot(1,3,2)
	plot(EE.avec,EE.C[:,1:(p.nT-1)])
	xlabel("assets")
	ylim([5,20])
	xlim([5,20])
	plot(ylim(),ylim())
	title("Solving\n Euler Equation:\n $(round(EE.toc,5)) secs")

	# solve by EGM
	subplot(1,3,3)
	plot(EGM.M[:,1:(p.nT-1)],EGM.C[:,1:(p.nT-1)])
	xlabel("cash on hand")
	ylim([5,20])
	xlim([5,20])
	plot(ylim(),ylim())
	title("Endogenous Grid\n Method:$(round(EGM.toc,5)) secs")

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

	subplot(1,3,3)
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
