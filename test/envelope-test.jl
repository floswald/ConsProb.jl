
module Models_test 

using ConsProb, FactCheck

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

		y = rand(5)
		ly = rand()
		b2 = [b; Bfun(y,ly)]

		@fact ConsProb.b(b2[1]) --> 0.8
		@fact ConsProb.b(b2[2]) --> ly
		@fact ConsProb.vb(b2[1]) --> [0.8;x2]
		@fact ConsProb.vb(b2[2]) --> [ly;y]

	end

	context("methods for Arrays of Bfun") do

		b = [Bfun(zeros(j+(i-1)*3),j+(i-1)*3) for j=1:3, i=1:2]
		for j=1:3,i=1:2
			@fact ConsProb.v(b[j,i]) --> zeros(j+(i-1)*3)
			@fact ConsProb.b(b[j,i]) --> j+(i-1)*3
			@fact ConsProb.vb(b[j,i]) --> [j+(i-1)*3;zeros(j+(i-1)*3)]
		end

		x = rand(4)
		y = rand()
		ConsProb.set!(b,x)
		ConsProb.set_bound!(b,y)
		for j=1:3,i=1:2
			@fact ConsProb.v(b[j,i]) --> x
			@fact ConsProb.b(b[j,i]) --> y
			@fact ConsProb.vb(b[j,i]) --> [y;x]
		end




	end

end

# facts("Envelope.") do

# 	context("constructors") do

# 		e = Envelope()
# 		@fact length(e.cond) => 2
# 		@fact length(e.cond_vbound) => 2
# 		@fact e.env_vbound => 0.0

# 		cc = [id=> rand(10) for id in 1:20]
# 		cv = [id=> rand() for id in 1:20]
# 		env = rand(15)
# 		envv = rand()
# 		e = Envelope(cc,cv,env,envv)
# 		@fact length(e.cond) => 20
# 		@fact length(e.cond_vbound) => 20

# 	end

# 	context("methods") do
# 		cc1 = [id=> rand(10) for id in 1:20]
# 		cv1 = [id=> rand() for id in 1:20]
# 		env1 = rand(15)
# 		envv1 = rand()
# 		e1 = Envelope(cc1,cv1,env1,envv1)

# 		cc2 = [id=> rand(10) for id in 1:20]
# 		cv2 = [id=> rand() for id in 1:20]
# 		env2 = rand(15)
# 		envv2 = rand()
# 		e2 = Envelope(cc2,cv2,env2,envv2)


# 		# single Envelope type
# 		# @fact cond(e1,1) => cc1[1]
# 		# @fact cond(e2,2) => cc2[2]
# 		# @fact condvbound(e1,1) => vcat(cv1[1],cc1[1])
# 		# @fact condvbound(e2,2) => vcat(cv2[2],cc2[2])
# 		# @fact env(e1) => env1
# 		# @fact env(e2) => env2
# 		# @fact envvbound(e1) => vcat(envv1,env1 )
# 		# @fact envvbound(e2) => vcat(envv2,env2 )

# 		# dict of envelopes
# 		ed = Dict{Int,Envelope}()
# 		ed[1] = e1
# 		ed[2] = e2
# 		@fact cond(ed,1,1) => cc1[1]
# 		@fact cond(ed,1,2) => cc1[2]
# 		@fact cond(ed,2,1) => cc2[1]
# 		@fact cond(ed,2,2) => cc2[2]

# 		@fact condvbound(ed,1,10) => vcat(cv1[10],cc1[10])
# 		@fact condvbound(ed,1,20) => vcat(cv1[20],cc1[20])
# 		@fact condvbound(ed,2,15) => vcat(cv2[15],cc2[15])
# 		@fact condvbound(ed,2,2) => vcat(cv2[2],cc2[2])

# 		@fact envvbound(ed,1) => vcat(envv1,env1)
# 		@fact envvbound(ed,2) => vcat(envv2,env2)

# 		x = rand(13)
# 		ConsProb.set!(ed,1,9,x)
# 		@fact cond(ed,1,9) => x
# 		ConsProb.set_vbound!(ed,1,9,x[3])
# 		@fact ConsProb.get_vbound(ed,1,9) => x[3]
# 		@fact ConsProb.condvbound(ed,1,9) => vcat(x[3],x)
# 		ConsProb.set_vbound!(ed,1,x[4])
# 		@fact ConsProb.get_vbound(ed,1) => x[4]
# 		@fact ConsProb.envvbound(ed,1) => vcat(x[4],env1)

# 	end
# end



end


