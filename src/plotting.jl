


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
	plot(m.M[:,p.nT],m.C[:,p.nT],label="period $(p.nT)")
	for it in (p.nT-1):-1:1
		plot(m.M[:,it],m.C[:,it],label="period $it")
	end
	ylim([p.a_low,100])
	xlim([p.a_low,100])
	ylim(xlim())
	plot(ylim(),ylim(),linestyle="--",linewidth=2.0,color="black")
	title(L"Consumption functions: $C_t(M_t)$")
	xlabel(L"Cash on hand: $M_t$")
	grid("on")
	# legend(loc="lower right",fancybox="true")

end

function plot(m::AR1Model,p::Param)
plot(m.V[:,:,1])
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
