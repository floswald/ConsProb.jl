

# testing funs

module test_funs

using ConsProb
using FactCheck

facts("test utility function") do

	p = Param()
	@fact ConsProb.u(1.0,p) => p.oneover_oneminusgamma

end



end