"""
# PLOT_MODEL(1) - Julia Function Documentation

## NAME
**plot_model** - Plot the solution of a solved SIR model using the `Plots.jl` library.

## SYNOPSIS
```
plot_model(solution)
```

## DESCRIPTION
The **plot_model** function takes an `ODESolution` object, which is the output of solving an SIR model using `DifferentialEquations.jl`, and generates a plot showing the dynamics of each compartment (Susceptible, Infected, and Recovered) over time. It uses the `plot` function from the `Plots.jl` library to visualize the data with time on the x-axis and the proportions of the compartments on the y-axis.

## PARAMETERS
**solution**  
: An `ODESolution` object that contains the time points and solution values for the SIR model compartments. The `solution` must be the output of a successfully solved SIR model using the `solve` function from `DifferentialEquations.jl`.

The solution object includes:

- **`t`**: A vector containing the time points at which the solution is evaluated.
- **`u`**: An array or vector of solution values corresponding to each time point. This should be a matrix with each compartment (Susceptible, Infected, Recovered) stored in separate columns.

## RETURN VALUE
The function generates a plot of the SIR model compartments (`Susceptible`, `Infected`, `Recovered`) against time. The plot is displayed using the `Plots.jl` library and includes labels and a title for clarity. The function does not return a value; it displays the plot directly.

## DIAGNOSTICS
**plot_model** displays a plot if the `solution` object is valid and contains the correct structure. The function will raise errors if:

- The input `solution` is not a valid `ODESolution` object.
- The `Plots.jl` library is not installed or imported correctly.
- The `solution` object does not contain the expected time points and values.

## NOTES
- Ensure that the `Plots` package is correctly installed and imported before using this function.
- The `solution` should contain at least three compartments (`Susceptible`, `Infected`, and `Recovered`) for proper plotting.

## SEE ALSO
- `plot` from the `Plots.jl` library for creating custom plots.
- `solve` from the `DifferentialEquations.jl` package for solving ODE problems.
- `ODESolution` for details on the solution object structure.
"""
function plot_model(solution)
    plot(solution, xlabel="Time", ylabel="Proportion", label=["Susceptible" "Infected" "Recovered"], title="SIR Model")
end

