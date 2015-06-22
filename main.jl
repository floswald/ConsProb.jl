


home = ENV["HOME"]
cd("$home/git/EGM/ConsProb.jl/")

include("src/ConsProb.jl")
ConsProb.Plotting.doplots()
ConsProb.Plotting.printplots()


x = ConsProb.Dchoice.run()
include("test/runtests.jl")
x = ConsProb.dchoice();
p  = ConsProb.Param(1.0)
m  = ConsProb.iidDModel(p)
ConsProb.EGM!(m,p)
ConsProb.plots(m,p)

p  = ConsProb.Models.Param()
m  = ConsProb.Models.AR1Model(p)
ConsProb.Standard.EGM!(m,p)

p  = ConsProb.Models.Param(mu=100)
m  = ConsProb.Models.iidDebtModel(p)
ConsProb.Standard.EGM!(m,p)
ConsProb.Plotting.plot(m.M[:,1:7],m.V[:,1:7])
ConsProb.Plotting.xlim([-5,5])
ConsProb.Plotting.ylim([0,5])


# run all
r = ConsProb.Standard.runStd()

# plot all
p  = ConsProb.Models.Param()
f = ConsProb.Plotting.plots(r,p)
