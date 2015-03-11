


home = ENV["HOME"]
cd("$home/git/EGM/ConsProb.jl/")

include("src/ConsProb.jl")

# run all
ConsProb.run()
