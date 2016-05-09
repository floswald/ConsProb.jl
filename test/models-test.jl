
module Models_test 

using ConsProb, FactCheck

facts("Models.") do

	context("Param Type") do

		p = Param()
		@fact p.gamma => 2.0
		@fact p.mu => 1
		@fact ConsProb.u(1.0,p) => p.oneover_oneminusgamma
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
		@fact ConsProb.u(1.0,p) => 0.0
		@fact ConsProb.up(1.0,p) => 1.0
		@fact ConsProb.iup(1.0,p) => 1.0

		x = rand()
		@fact ConsProb.u(x,p) => log(x)
		@fact ConsProb.up(x,p) => 1.0 / x
		@fact ConsProb.iup(x,p) => 1.0 / x
		# check utlity cost
		@fact ConsProb.u(1.0,false,p) => 0.0
		@fact ConsProb.u(1.0,true,p) => -p.alpha

		#crra utility
		p = Param()
		x = rand()
		@fact ConsProb.u(x,p) => p.oneover_oneminusgamma * (x^p.oneminusgamma)
		@fact ConsProb.up(x,p) => x^(-p.gamma)
		@fact ConsProb.iup(x,p) => x^(-1/p.gamma)
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

facts("Bfun.") do

	context("constructor") do

		x = rand(3)
		l = 1.1
		b = Bfun(x,l)
		@fact ConsProb.v(b) --> x
		@fact ConsProb.vb(b) --> [l;x]

	end

	context("methods") do
		x = rand(3)
		l = -0.2
		b = Bfun(x,l)

		x2 = rand(3)
		ConsProb.set!(b,x2)

		@fact ConsProb.v(b) --> x2

		ConsProb.set_bound!(b,0.8)
		@fact ConsProb.b(b) --> 0.8

		@fact ConsProb.vb(b) --> [0.8;x2]
		
	end

end

end


