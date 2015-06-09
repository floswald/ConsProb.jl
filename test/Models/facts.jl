
using ConsProb.Models

facts("Models.") do

	context("Param Type") do

		p = Param()
		@fact p.gamma => 2.0
		@fact p.mu => 1
		@fact u(1.0,p) => p.oneover_oneminusgamma
		@fact p.a_low => 1e-6
		@fact p.a_high => 300

		# change mu
		p = Param(mu=20.0)
		@fact p.mu => 20.0
		@fact p.gamma => 2.0
		@fact p.a_low => 1e-6
		@fact p.a_high => 300

		p = Param(1.0)
		@fact p.gamma => 1.0
		@fact p.mu => 0
		@fact p.a_low => 0.0
		@fact p.a_high => 50
	end

	context("Utilty Functions") do

		#log utility
		p = Param(1.0)
		@fact p.gamma => 1.0
		@fact p.mu => 0
		@fact u(1.0,p) => 0.0
		@fact up(1.0,p) => 1.0
		@fact iup(1.0,p) => 1.0

		x = rand()
		@fact u(x,p) => log(x)
		@fact up(x,p) => 1.0 / x
		@fact iup(x,p) => 1.0 / x
		# check utlity cost
		@fact u(1.0,false,p) => 0.0
		@fact u(1.0,true,p) => -p.alpha

		#crra utility
		p = Param()
		x = rand()
		@fact u(x,p) => p.oneover_oneminusgamma * (x^p.oneminusgamma)
		@fact up(x,p) => x^(-p.gamma)
		@fact iup(x,p) => x^(-1/p.gamma)
	end

	context("iidModel") do
		p = Param()
		m = iidModel(p)
		@fact size(m.ywgt) => (p.ny,)
		@fact m.avec[1] => roughly(1e-6)
		@fact m.avec[end] => roughly(300)
		@fact size(m.m1) => (p.na,p.ny)
	end

	context("iidDebtModel") do
		p = Param()
		m = iidDebtModel(p)
		@fact size(m.avec) => (p.na,p.nT)
		@fact size(m.ywgt) => (p.ny,)
		@fact size(m.m1) => (p.na,p.ny,p.nT-1)
		for it in 1:(p.nT-2)
			@fact m.avec[end,it] => roughly(300)
			@fact m.avec[1,it] => less_than(0)
		end
	end

	context("Binary Discrete Choice Model") do
		p = Param()
		m = iidDModel(p)
		@fact m.avec[1] => 0.0
		@fact length(m.m1) => p.nT
		@fact length(m.m1[1]) => p.nD
		@fact size(m.m1[1][1]) => (p.na,p.ny)
		@fact size(m.ywgt) => (p.ny,)
		@fact typeof(m.m) => Dict{Int,Envelope}
		@fact length(m.m) => p.nT
		@fact typeof(m.c) => Dict{Int,Envelope}
		@fact typeof(m.v) => Dict{Int,Envelope}
	end

	context("AR1 Model") do
		p = Param()
		m = AR1Model(p)
		@fact m.avec[1] => roughly(1e-6)
		@fact m.avec[end] => roughly(p.a_high)
		@fact size(m.m1) => (p.na,p.ny)
		@fact size(m.ywgt) => (p.ny,p.ny)
		@fact norm(sum(m.ywgt,2) .- ones(p.ny)) => roughly(0,atol=1e-12)

	end

end

facts("Envelope.") do

	context("constructors") do

		e = Envelope()
		@fact length(e.cond) => 2
		@fact length(e.cond_vbound) => 2
		@fact e.env_vbound => 0.0

		cc = [id=> rand(10) for id in 1:20]
		cv = [id=> rand() for id in 1:20]
		env = rand(15)
		envv = rand()
		e = Envelope(cc,cv,env,envv)
		@fact length(e.cond) => 20
		@fact length(e.cond_vbound) => 20

	end

	context("methods") do
		cc1 = [id=> rand(10) for id in 1:20]
		cv1 = [id=> rand() for id in 1:20]
		env1 = rand(15)
		envv1 = rand()
		e1 = Envelope(cc1,cv1,env1,envv1)

		cc2 = [id=> rand(10) for id in 1:20]
		cv2 = [id=> rand() for id in 1:20]
		env2 = rand(15)
		envv2 = rand()
		e2 = Envelope(cc2,cv2,env2,envv2)


		# single Envelope type
		@fact cond(e1,1) => cc1[1]
		@fact cond(e2,2) => cc2[2]
		@fact condvbound(e1,1) => vcat(cv1[1],cc1[1])
		@fact condvbound(e2,2) => vcat(cv2[2],cc2[2])
		@fact env(e1) => env1
		@fact env(e2) => env2
		@fact envvbound(e1) => vcat(envv1,env1 )
		@fact envvbound(e2) => vcat(envv2,env2 )

		# dict of envelopes
		ed = Dict{Int,Envelope}()
		ed[1] = e1
		ed[2] = e2
		@fact cond(ed,1,1) => cc1[1]
		@fact cond(ed,1,2) => cc1[2]
		@fact cond(ed,2,1) => cc2[1]
		@fact cond(ed,2,2) => cc2[2]

		@fact condvbound(ed,1,10) => vcat(cv1[10],cc1[10])
		@fact condvbound(ed,1,20) => vcat(cv1[20],cc1[20])
		@fact condvbound(ed,2,15) => vcat(cv2[15],cc2[15])
		@fact condvbound(ed,2,2) => vcat(cv2[2],cc2[2])

		@fact envvbound(ed,1) => vcat(envv1,env1)
		@fact envvbound(ed,2) => vcat(envv2,env2)

		x = rand(13)
		set!(ed,1,9,x)
		@fact cond(ed,1,9) => x
		set_vbound!(ed,1,9,x[3])
		@fact get_vbound(ed,1,9) => x[3]
		@fact condvbound(ed,1,9) => vcat(x[3],x)
		set_vbound!(ed,1,x[4])
		@fact get_vbound(ed,1) => x[4]
		@fact envvbound(ed,1) => vcat(x[4],env1)

	end





end


