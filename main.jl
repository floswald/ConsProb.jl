


home = ENV["HOME"]
cd("$home/git/EGM/ConsProb.jl/")

include("src/ConsProb.jl")

p  = ConsProb.Param()
m  = ConsProb.iidDModel(p)
ConsProb.EGM!(m,p)



m  = ConsProb.iidModel(p)
ConsProb.EGM!(m,p)
ConsProb.plot_vf(m,p)

m2 = ConsProb.Model2(p)


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
ConsProb.runall()
ConsProb.figure("iidCons")
ConsProb.savefig("$home/Dropbox/public/ConsProb.jl/iidCons.png")
ConsProb.figure("AR1cons")
ConsProb.savefig("$home/Dropbox/public/ConsProb.jl/AR1cons.png")
ConsProb.figure("AR1vals")
ConsProb.savefig("$home/Dropbox/public/ConsProb.jl/AR1vals.png")

