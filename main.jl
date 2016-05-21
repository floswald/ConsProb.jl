

cd(Pkg.dir("ConsProb"))
include("src/ConsProb.jl")
include("test/runtests.jl")


p  = ConsProb.Param(2.0,-25.0)
m  = ConsProb.iidModel(p)
ConsProb.EGM!(m,p)
for i in 1:p.nT-1
	ConsProb.plot(ConsProb.v(m.M[i]),ConsProb.v(m.V[i]))
end
ConsProb.figure()
for i in 1:p.nT-1
	ConsProb.plot(ConsProb.vb(m.M[i]),ConsProb.vb(m.C[i]))
end
ConsProb.xlim(-3,3)
ConsProb.ylim(0,3)


ConsProb.doplots()
# ConsProb.printplots()


# x = ConsProb.run()
# include("test/runtests.jl")
# x = ConsProb.dchoice();
# p  = ConsProb.Param(1.0,-2.0)
# m  = ConsProb.iidDModel(p)
# ConsProb.EGM!(m,p)

# p  = ConsProb.Models.Param()
# m  = ConsProb.Models.AR1Model(p)
# ConsProb.Standard.EGM!(m,p)

# p  = ConsProb.Param(mu=100)
# m  = ConsProb.iidDebtModel(p)
# ConsProb.EGM!(m,p)
# ConsProb.plot(m.M[:,1:7],m.V[:,1:7])
# ConsProb.figure()
# ConsProb.plot(m.M[:,1:7],m.C[:,1:7])



# # run all
# r = ConsProb.Standard.runStd()

# # plot all
# p  = ConsProb.Models.Param()
# f = ConsProb.Plotting.plots(r,p)
