


home = ENV["HOME"]
cd("$home/git/EGM/ConsProb.jl/")

include("src/ConsProb.jl")


include("test/tfuns.jl")

p  = ConsProb.Param()
m  = ConsProb.iidDebtModel(p)
ConsProb.EGM!(m,p)
ConsProb.plot(m.M[:,1:7],m.V[:,1:7])
ConsProb.xlim([-5,100])

ConsProb.plot(m.M[:,1:7],m.C[:,1:7])
ConsProb.ylim([-1,2])
ConsProb.xlim([-1,2])



m  = ConsProb.iidModel(p,"Euler")
m  = ConsProb.iidModel(p,"Euler")
ConsProb.Euler!(m,p)
ConsProb.solve!(m,p)
ConsProb.EGM!(m,p)
ConsProb.plot_vf(m,p)

m2 = ConsProb.iidModel(p)


m  = ConsProb.AR1Model(p)
ConsProb.EGM(m,p)


ConsProb.VFbi(m,p)
ConsProb.VFbi_2(m2,p)



m3 = ConsProb.Model2(p)
ConsProb.VFbi(m,p)
ConsProb.EGM(m2,p)
ConsProb.EGM(m,p)
ConsProb.VFbi(m3,p)


# run all
r = ConsProb.runall()

# plot all
p  = ConsProb.Param()
f = ConsProb.plots(r,p)

ConsProb.figure(f[1][:number])
ConsProb.savefig("$home/Dropbox/public/ConsProb.jl/iidCons.png")


ConsProb.figure("AR1cons")
ConsProb.savefig("$home/Dropbox/public/ConsProb.jl/AR1cons.png")
ConsProb.figure("AR1vals")
ConsProb.savefig("$home/Dropbox/public/ConsProb.jl/AR1vals.png")

