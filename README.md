

# Consumption / Savings problems in Economics

this repo has `julia` code that computes the solution to a finite time savings problem with iid or AR1  income uncertainty using 3 different methods:

1. value function backward iteration
2. euler equation backward iteration
3. endogenous grid method

# example

this will produce a plot showing the resulting optimal consumption functions for each method and for each uncertainty scenario. Notice that the income uncertainty scenarios are not directly comparable, so you shouldn't expect the figures to be identical.

```julia
using ConsProb
ConsProb.run()
```