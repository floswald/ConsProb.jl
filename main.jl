


home = ENV["HOME"]
cd("$home/git/EGM/ConsProb.jl/")

include("src/ConsProb.jl")

p  = ConsProb.Param()
m  = ConsProb.Model(p)
m2 = ConsProb.Model2(p)
m3 = ConsProb.Model2(p)
ConsProb.VFbi(m,p)
ConsProb.EGM(m2,p)
ConsProb.EGM(m,p)
ConsProb.VFbi(m3,p)


# run all
ConsProb.run()
