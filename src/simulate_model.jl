"""
# SIMULATE_MODEL(1) - Julia Function Documentation

## NAME
**simulate_model** - Simulate and visualize the SIR model dynamics, optionally comparing with real-world data.

## SYNOPSIS
```
simulate_model(model_type::Symbol, c, β, γ, herd_threshold, S0, I0, R0, tspan, N; real_data_time=[], real_data_infected=[])
```

## DESCRIPTION
The **simulate_model** function runs a specified SIR model based on the given parameters and initial conditions. It generates a plot showing the time evolution of the infected population. If real-world data points are provided, the function overlays the actual data onto the model’s output, allowing for a visual comparison between the simulated results and observed data.

This function is designed to provide an intuitive way to validate SIR models and assess their accuracy against real-world datasets.

## PARAMETERS
**model_type**  
: A `Symbol` representing the type of SIR model to run. Possible values include:
  - `:basic`: Basic SIR model.
  - `:force_of_infection`: SIR model with force of infection.
  - `:herd_immunity`: SIR model with herd immunity.

**c**  
: Contact rate per day per individual. Affects how frequently susceptible individuals come into contact with infected individuals.

**β**  
: Transmission rate probability per contact. Determines the likelihood of disease transmission during each contact.

**γ**  
: Recovery rate. Represents the rate at which infected individuals recover and move into the recovered compartment.

**herd_threshold**  
: Proportion of the population required to reach herd immunity. Only used if `model_type == :herd_immunity`.

**S0**  
: Initial proportion of the population that is susceptible.

**I0**  
: Initial proportion of the population that is infected.

**R0**  
: Initial proportion of the population that is recovered.

**tspan**  
: A tuple defining the start and end time of the simulation (e.g., `(0.0, 50.0)`).

**N**  
: Total population size. Used to scale initial proportions to absolute values for the ODE system.

**real_data_time**  
: (Optional) A vector of time points corresponding to the real-world data. If provided, the real-world data will be overlaid on the model’s output.

**real_data_infected**  
: (Optional) A vector of real-world data values corresponding to the number of infected individuals at the given `real_data_time` points.

## RETURN VALUE
Returns a `Plots.Plot` object (`plt`) that visualizes:

- **Model Infected**: The number of infected individuals over time as predicted by the model.
- **Real Data**: The real-world infected data points (if provided) for visual comparison.

The x-axis represents time in days, and the y-axis represents the number of infected individuals.

## DIAGNOSTICS
**simulate_model** generates a plot if all inputs are valid. However, the function may raise errors if:

- The `model_type` is not specified correctly.
- The `run_model` function encounters an error due to incorrect initial conditions or incompatible parameters.
- Invalid time span or parameter values are provided (e.g., negative values).

## NOTES
- Ensure that the `DifferentialEquations` and `Plots` packages are correctly installed and imported before using this function.
- To compare the model against real-world data, pass `real_data_time` and `real_data_infected` as vectors of equal length.

## SEE ALSO
- `run_model` for the main SIR model simulations.
- `plot` and `scatter!` from the `Plots.jl` library for customizing the visualizations.
"""
function simulate_model(model_type::Symbol, c, β, γ, herd_threshold, S0, I0, R0, tspan, N; real_data_time=[], real_data_infected=[])
  # Solve the problem
  solution = run_SIR(model_type, c, β, γ, herd_threshold, S0, I0, R0, tspan, N)

  # Extract model data
  time_points = solution.t
  infected = solution[2, :]

  # Create a smoother time vector for plotting
  fine_time_points = range(tspan[1], tspan[2], length=1000)
  infected_interp = [solution(t)[2] for t in fine_time_points]

  return plot_model(solution)
end
