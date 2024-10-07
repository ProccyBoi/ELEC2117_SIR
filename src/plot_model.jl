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


"""
# DATA_FITTING(1) - Julia Function Documentation

## NAME
**data_fitting** - Fit an SIR model to real-world data and visualize the comparison.

## SYNOPSIS
```
data_fitting(c, β, γ, N, S0, I0, R0, data)
```

## DESCRIPTION
The **data_fitting** function solves a system of ordinary differential equations (ODEs) using the SIR model with the given parameters and initial conditions. It then compares the model's infected population to real-world data, plotting both for visual inspection of the fit.

The function uses the `Force_of_infection!` model equations to simulate the spread of infection and plots the model's infected population over time alongside real data points. It is designed to provide a quick visual comparison to assess the model’s performance.

## PARAMETERS
**c**  
: Contact rate per day per individual. Affects how frequently susceptible individuals come into contact with infected individuals.

**β**  
: Transmission rate probability per contact. Determines the likelihood of disease transmission during each contact.

**γ**  
: Recovery rate. Represents the rate at which infected individuals recover and move into the recovered compartment.

**N**  
: Total population size. Used to scale initial proportions to absolute values for the ODE system.

**S0**  
: Initial proportion of the population that is susceptible.

**I0**  
: Initial proportion of the population that is infected.

**R0**  
: Initial proportion of the population that is recovered.

**data**  
: A vector containing the real-world infected data points to be used for comparison. The length of this vector determines the time span of the model simulation.

## RETURN VALUE
Returns a plot (`plt`) that displays:

- **Model Infected**: The number of infected individuals over time as predicted by the model.
- **Data**: The real-world infected data points for visual comparison.

The x-axis represents time in days, and the y-axis represents the number of infected individuals. The plot visually compares the model’s infected values with the given real data points.

## DIAGNOSTICS
**data_fitting** displays a plot of the model solution and real data if the inputs are valid. However, the function may raise errors if:

- The length of the `data` vector is not compatible with the time span of the model.
- The `DifferentialEquations.jl` or `Plots.jl` packages are not correctly installed.
- Invalid parameter values (e.g., negative population sizes or zero transmission rates) are provided.

## NOTES
- Ensure that both the `DifferentialEquations` and `Plots` packages are installed and correctly imported before using this function.
- Adjust the parameters (`c`, `β`, `γ`) to find the best fit for the real-world data. This can be used as part of a parameter optimization routine.

## SEE ALSO
- `solve` from the `DifferentialEquations.jl` package for solving ODE problems.
- `plot` and `scatter!` from the `Plots.jl` library for creating visualizations.
- `Force_of_infection!` for details on the model equations used.
"""
function data_fitting(c, β, γ, N, S0, I0, R0, data)
    # Define the time span based on the length of the data points (21 days)
    tspan = (0.0, length(data) - 1)  # Time span: 0 to 20

    # Create a corresponding time vector for the data points
    real_data_time = collect(tspan[1]:tspan[2])
    real_data_infected = data

    # Initial conditions
    u0 = [S0 * N, I0 * N, R0 * N]  # Scale initial proportions to absolute values
    # Define the ODE problem
    prob = ODEProblem(Force_of_infection!, u0, tspan, [c, β, γ, N])

    # Solve the problem with a dense solution for smooth plotting
    solution = solve(prob, saveat=0.1)

    # Extract the time points and infected values from the model solution
    model_time = solution.t
    model_infected = solution[2, :]  # Extract the infected individuals over time

    # Plot the model output vs. the real data
    plt = plot(model_time, model_infected, label="Model Infected", xlabel="Time (Days)", ylabel="Number Infected",
        title="Model vs. Real Data", lw=2, color=:blue)
    scatter!(real_data_time, real_data_infected, label="Data", marker=:circle, ms=3, color=:orange)

    return plt
end
