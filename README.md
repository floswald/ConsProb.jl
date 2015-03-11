

# Consumption / Savings problems in Economics

this repo has julia code that computes the solution to a finite time savings problem with iid income uncertainty using 3 different methods:

1. value function backward iteration
2. euler equation backward iteration
3. endogenous grid method

# example

this will produce a plot showing the resulting optimal consumption functions for each method.

```julia
using ConsProb
run()
```