

# Consumption / Savings / Discrete Choice Problems

1. Standard Lifecycle Savings
2. Lifecycle Savings with discrete Retirement choice
3. License


## Standard Lifecycle Savings

this repo has `julia` code that computes the solution to a finite time savings problem with iid or AR1 income uncertainty using 3 different methods:

1. value function backward iteration
2. euler equation backward iteration
3. endogenous grid method

## Lifecycle Savings with discrete Retirement choice

It also has julia implementation of the discrete choice model which is based on code from [Fedor Iskhakov's](https://github.com/fediskhakov/egdst) github. This implements the method in the [working paper by Iskhakov, Jorgensen, Rust and Schjerning](https://dl.dropboxusercontent.com/u/17240700/sync4web/dcegm.pdf)

## example

```julia
using ConsProb
r = ConsProb.runall()  # runs EGM, Euler and VFi savings models
p  = Param()
f = ConsProb.plots(r,p) # produces plots

ConsProb.dchoice()  # runs and plots EGM for discrete choice
```

## Output Savings Models

### iid income model, values at all ages:
[![iidcons](https://dl.dropboxusercontent.com/u/109115/ConsProb.jl/iidVfun.png)]()

### iid income model, consumption at all ages:
[![iidcons](https://dl.dropboxusercontent.com/u/109115/ConsProb.jl/iidCons.png)]()

### Value functions in AR1 income model, all y-states, age=1:
[![values](https://dl.dropboxusercontent.com/u/109115/ConsProb.jl/AR1Vfun.png)]()

### Consumption functions in AR1 income model, all y-states, age=1:
[![ar1cons](https://dl.dropboxusercontent.com/u/109115/ConsProb.jl/AR1Cons.png)]()


## Output Discrete Choice/Savings Model

### Conditinonal Value functions:
[![iidcons](https://dl.dropboxusercontent.com/u/109115/ConsProb.jl/Dchoice_condV.png)]()

### Envevelope over consumption functions, all ages:
[![iidcons](https://dl.dropboxusercontent.com/u/109115/ConsProb.jl/Dchoice_envC.png)]()

### Envevelope over Value functions, all ages:
[![iidcons](https://dl.dropboxusercontent.com/u/109115/ConsProb.jl/Dchoice_envV.png)]()

## License

Please observe the license when using this code in your work (see file `LICENSE`). Thank you.
