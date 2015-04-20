

# testing funs

module test_funs

using ConsProb
using FactCheck

facts("test utility function") do

	p = Param()
	@fact ConsProb.u(1.0,p) => p.oneover_oneminusgamma

	x = rand()
	@fact ConsProb.up(x,p) => x^(-p.gamma)

	@fact ConsProb.iup(x,p) => x^(-1/p.gamma)
end


facts("check consumption function in period T") do

	context("iid Models") do

		p = Param()
		d = ["EGM" => iidModel(p,"EGM"), "VF" => iidModel(p,"VFbi"), "Euler" => iidModel(p,"Euler")]
		for (k,v) in d
			solve!(v,p)
		end

		# check consumption in period T

		for (k,v) in d
			x = v.avec
			x[x.<p.cfloor] = p.cfloor
			@fact maximum( abs(v.C[:,p.nT] .- x) ) => roughly(0.0) "wrong for $k"
		end
	end

	# context("AR1 Models") do

	# end


end

facts("check vectorized iup()") do

	context("iid Models") do

		p = Param()
		d = ["EGM" => iidModel(p,"EGM"), "VF" => iidModel(p,"VFbi"), "Euler" => iidModel(p,"Euler")]
		for (k,v) in d
			solve!(v,p)
		end

		# cons in period T

		for (k,v) in d
			x = v.avec
			x[x.<p.cfloor] = p.cfloor
			muc = p.beta * p.R * ConsProb.up(x,p)
			c0 = muc .^ (-1/p.gamma)
			@fact maxabs(c0 .- ConsProb.iup(muc,p)) => roughly(0.0) "iup() wrong for $k"
		end

	end


	# context("AR1 Models") do

	# end


end

end # module