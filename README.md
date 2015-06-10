

# Consumption / Savings / Discrete Choice Problems

1. Installation
2. Standard Lifecycle Savings
3. Lifecycle Savings with discrete Retirement choice
4. License

## 1. Installation

```julia
# in Julia:
Pkg.clone("https://github.com/floswald/ConsProb.jl")
```

## 2. Standard Lifecycle Savings

this repo has `julia` code that computes the solution to a finite time savings problem with iid or AR1 income uncertainty using 3 different methods:

1. value function backward iteration
2. euler equation backward iteration
3. endogenous grid method

## 3. Lifecycle Savings with discrete Retirement choice

It also has julia implementation of the discrete choice model which is based on code from [Fedor Iskhakov's](https://github.com/fediskhakov/egdst) github. This implements the method in the [working paper by Iskhakov, Jorgensen, Rust and Schjerning](https://dl.dropboxusercontent.com/u/17240700/sync4web/dcegm.pdf)


----


### example

```julia
using ConsProb
Plotting.doplots()
Plotting.printplots(dir)
```

----

### Output Savings Models

#### iid income model, values at all ages:
[![iidcons](https://dl.dropboxusercontent.com/u/109115/ConsProb.jl/iidVfun.png)]()

#### iid income model, consumption at all ages:
[![iidcons](https://dl.dropboxusercontent.com/u/109115/ConsProb.jl/iidCons.png)]()

#### Value functions in AR1 income model, all y-states, age=1:
[![values](https://dl.dropboxusercontent.com/u/109115/ConsProb.jl/AR1Vfun.png)]()

#### Consumption functions in AR1 income model, all y-states, age=1:
[![ar1cons](https://dl.dropboxusercontent.com/u/109115/ConsProb.jl/AR1Cons.png)]()

----

### Output Discrete Choice/Savings Model

In this section, dashed lines in a graph stand for an analytic solution in that region (i.e. no approximation to be done).

#### Conditinonal Value functions:
[![iidcons](https://dl.dropboxusercontent.com/u/109115/ConsProb.jl/Dchoice_condV.png)]()

#### Envevelope over consumption functions, all ages:
[![iidcons](https://dl.dropboxusercontent.com/u/109115/ConsProb.jl/Dchoice_envC.png)]()

#### Envevelope over Value functions, all ages:
[![iidcons](https://dl.dropboxusercontent.com/u/109115/ConsProb.jl/Dchoice_envV.png)]()

----

## 4. License

Please observe the license when using this code in your work (see file `LICENSE`). Thank you.
