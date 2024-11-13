# Define a struct to hold the SIRS model parameters with intervention
struct SIRSInterventionParameters
    β::Float64         # Transmission rate
    γ::Float64         # Recovery rate for mild cases
    δ::Float64         # Proportion progressing to severe infection
    γ_s::Float64       # Recovery rate for severe cases
    α::Float64         # Re-susceptibility rate after recovery
    N::Int             # Total population
    intervention_day::Int  # Day when intervention begins
    ε_i::Float64       # Efficacy of intervention (reduction in transmission rate)
    p_i::Float64       # Proportion of population complying with intervention
end

# Add a keyword-based constructor for SIRSInterventionParameters
function SIRSInterventionParameters(; β::Float64, γ::Float64, δ::Float64, γ_s::Float64, α::Float64, N::Int,
    intervention_day::Int, ε_i::Float64, p_i::Float64)
    return SIRSInterventionParameters(β, γ, δ, γ_s, α, N, intervention_day, ε_i, p_i)
end

# Define the SIRS model with intervention
"""
    intervention!(du, u, params, t)

Differential equations for the SIRS model with intervention, adjusting transmission rate after `intervention_day`.

# Arguments
- `du::Vector`: Derivative of state variables (output).
- `u::Vector`: Current state variables `[S, I, SevI, R]`.
- `params::SIRSInterventionParameters`: Model parameters including intervention.
- `t::Float64`: Time variable.

# Equations
The model dynamically adjusts transmission rate post-intervention day.
"""
function intervention!(du, u, params::SIRSInterventionParameters, t)
    S, I, SevI, R = u

    # Adjust β based on intervention day and compliance
    β_eff = t > params.intervention_day ? params.β * (1 - params.ε_i * params.p_i) : params.β

    # Force of infection, λ
    λ = β_eff * I / params.N

    # Differential equations, identical to SIRS_model! but with dynamic β
    du[1] = -λ * S + params.α * R                         # Susceptible (dS/dt)
    du[2] = λ * S - params.γ * I                          # Infectious (dI/dt)
    du[3] = params.δ * I * params.γ - params.γ_s * SevI   # Severe Infection (d(SevI)/dt)
    du[4] = (1 - params.δ) * params.γ * I + params.γ_s * SevI - params.α * R  # Recovered (dR/dt)
end

# Function to run the model with intervention parameters
"""
    run_intervention(params, S0, I0, SevI0, R0, tspan)

Runs the SIRS model with intervention.

# Arguments
- `params::SIRSInterventionParameters`: Struct containing model and intervention parameters.
- `S0::Float64`: Initial susceptible population count.
- `I0::Float64`: Initial infected population count.
- `SevI0::Float64`: Initial severe infection count.
- `R0::Float64`: Initial recovered count.
- `tspan::Tuple{Float64, Float64}`: Time span for the simulation.

# Returns
- `solution`: Solution object from DifferentialEquations.jl.
"""
function run_intervention(params::SIRSInterventionParameters, S0::Float64, I0::Float64, SevI0::Float64, R0::Float64, tspan::Tuple{Float64,Float64})
    u0 = [S0, I0, SevI0, R0]
    prob = ODEProblem(intervention!, u0, tspan, params)
    solution = solve(prob, Tsit5(), saveat=1.0)
    return solution
end

# Function to plot intervention results with absolute population values
"""
    plot_intervention(solution; title="SIRS Model with Intervention", ylabel="Population")

Plots the SIRS model solution with intervention in terms of absolute population.

# Arguments
- `solution`: Solution object from DifferentialEquations.jl.
- `title`: Title of the plot (optional).
- `ylabel`: Y-axis label (optional).

# Returns
- `plt`: Plot object.
"""
function plot_intervention(solution; title="SIRS Model with Intervention", ylabel="Population")
    N = solution.prob.p.N  # Retrieve total population from parameters

    plt = plot(
        solution.t, solution[1, :], label="Susceptible", xlabel="Time (days)", ylabel=ylabel, lw=2, color=:blue
    )
    plot!(plt, solution.t, solution[2, :], label="Infectious", lw=2, color=:orange)
    plot!(plt, solution.t, solution[3, :], label="Severe Illness", lw=2, color=:green)
    plot!(plt, solution.t, solution[4, :], label="Recovered", lw=2, color=:purple)
    title!(plt, title)
    return plt
end

# Simulation function
"""
    simulate_intervention()

Simulates and plots the SIRS model with intervention using default parameters.
"""
function simulate_intervention()
    # Define model parameters with intervention using the keyword-based constructor
    params = SIRSInterventionParameters(
        β=0.060 * 8,      # Transmission probability
        γ=1 / 7,           # Recovery rate for mild cases
        δ=0.2,             # Proportion progressing to severe infection
        γ_s=1 / 14,        # Recovery rate for severe cases
        α=1 / 30,          # Re-susceptibility rate after recovery
        N=6000,            # Total population
        intervention_day=30,  # Day on which intervention starts
        ε_i=0.3,           # Efficacy of intervention
        p_i=0.8            # Proportion complying with intervention
    )

    # Initial conditions (absolute population counts, converted to Float64)
    S0, I0, SevI0, R0 = Float64(params.N - 1), Float64(1), 0.0, 0.0
    tspan = (0.0, 180.0)  # Simulation for 300 days

    # Run model and plot results
    solution = run_intervention(params, S0, I0, SevI0, R0, tspan)
    plt = plot_intervention(solution)
    display(plt)
end

# Function to estimate the range of expected infections and severe cases
"""
    estimate_infection_severity_ranges(β_values, params, tspan)

Simulates the model over a range of β values and estimates the range and peak values of infected and severe cases.

# Arguments
- `β_values::Vector`: Range of transmission rates to evaluate.
- `params::SIRSInterventionParameters`: Base model parameters (excluding β, which varies).
- `tspan::Tuple{Float64, Float64}`: Time span for the simulation.

# Returns
- A tuple containing:
  - `β_peak_infections`: A dictionary mapping β values to their peak infections.
  - `β_peak_severe`: A dictionary mapping β values to their peak severe cases.
"""
function estimate_infection_severity_ranges()

    β_values = 0.00:0.001:0.50  # Example range of β values to explore

    # Define model parameters with intervention using the keyword-based constructor
    params = SIRSInterventionParameters(
        β=0.0365 * 8,      # Initial transmission probability (will be overridden in the loop)
        γ=1 / 7,           # Recovery rate for mild cases
        δ=0.2,             # Proportion progressing to severe infection
        γ_s=1 / 14,        # Recovery rate for severe cases
        α=1 / 30,          # Re-susceptibility rate after recovery
        N=6000,            # Total population
        intervention_day=30,  # Day on which intervention starts
        ε_i=0.3,           # Initial efficacy (will be overridden in the loop)
        p_i=0.816            # Proportion complying with intervention
    )

    # Define the time span for the simulation
    tspan = (0.0, 300.0)  # Example time span for 100 days

    # Default initial conditions (can be modified directly here)
    S0 = Float64(params.N - 1)  # Initial susceptible population
    I0 = Float64(1)             # Initial infected population
    SevI0 = Float64(0)          # Initial severe infection population
    R0 = Float64(0)             # Initial recovered population

    # Dictionaries to store peak values for infections and severe cases
    β_peak_infections = Dict{Float64,Float64}()
    β_peak_severe = Dict{Float64,Float64}()

    # Create an empty plot to display the results
    plt = plot(title="Infection and Severe Case Dynamics over β Range", xlabel="Time (days)", ylabel="Population")

    # Loop over each β value
    for β in β_values
        # Update the model parameters with the current β
        current_params = SIRSInterventionParameters(
            β=β,
            γ=params.γ,
            δ=params.δ,
            γ_s=params.γ_s,
            α=params.α,
            N=params.N,
            intervention_day=params.intervention_day,
            ε_i=params.ε_i,
            p_i=params.p_i
        )

        # Run the model with the current parameters
        solution = run_intervention(current_params, S0, I0, SevI0, R0, tspan)

        # Extract infected and severe infection values
        infected = solution[2, :]  # Infectious compartment over time
        severe_infected = solution[3, :]  # Severe infection compartment over time

        # Find the peak values
        peak_infected = maximum(infected)
        peak_severe = maximum(severe_infected)

        # Store the peaks in dictionaries
        β_peak_infections[β] = peak_infected
        β_peak_severe[β] = peak_severe

        # Plot the results for each β value
        plot!(plt, solution.t, infected, label="Infected (β=$β)", lw=2, color=:orange, linestyle=:solid)
        plot!(plt, solution.t, severe_infected, label="Severe (β=$β)", lw=2, color=:green, linestyle=:dash)
    end

    # Calculate and print the uncertainty range for peak infections and severe cases
    min_peak_infected = minimum(values(β_peak_infections))
    max_peak_infected = maximum(values(β_peak_infections))
    min_peak_severe = minimum(values(β_peak_severe))
    max_peak_severe = maximum(values(β_peak_severe))

    println("Range of peak infections: $min_peak_infected to $max_peak_infected")
    println("Range of peak severe cases: $min_peak_severe to $max_peak_severe")

    # Display the plot
    display(plot(plt, legend=false))

    return
end

# Function to estimate and plot trajectories for varying efficacy levels and β values
function plot_intervention_trajectories()
    β_values = [0.03, 0.04, 0.05, 0.06]           # Different β values to simulate
    efficacy_levels = [0.1, 0.2, 0.4, 0.5]         # Different efficacy levels (10%, 20%, 40%, 50%)
    tspan = (0.0, 300.0)                           # Time span for the simulation

    # Define model parameters with intervention
    params = SIRSInterventionParameters(
        β=0.0365 * 8,      # Initial transmission probability (will be overridden in the loop)
        γ=1 / 7,           # Recovery rate for mild cases
        δ=0.2,             # Proportion progressing to severe infection
        γ_s=1 / 14,        # Recovery rate for severe cases
        α=1 / 30,          # Re-susceptibility rate after recovery
        N=6000,            # Total population
        intervention_day=30,  # Day on which intervention starts
        ε_i=0.3,           # Initial efficacy (will be overridden in the loop)
        p_i=0.816            # Proportion complying with intervention
    )

    # Default initial conditions (can be modified directly here)
    S0 = Float64(params.N - 1)  # Initial susceptible population
    I0 = Float64(1)             # Initial infected population
    SevI0 = Float64(0)          # Initial severe infection population
    R0 = Float64(0)             # Initial recovered population

    # Initialize a 4x4 grid plot layout
    grid_plot = plot(layout=(4, 4), size=(1200, 1200), xlabel="Time (days)", ylabel="Population")

    # Loop over each β value (rows) and each efficacy level (columns)
    for (i, β) in enumerate(β_values)
        for (j, ε_i) in enumerate(efficacy_levels)
            # Update the model parameters with the current β and efficacy level
            current_params = SIRSInterventionParameters(
                β=β * 8,
                γ=params.γ,
                δ=params.δ,
                γ_s=params.γ_s,
                α=params.α,
                N=params.N,
                intervention_day=params.intervention_day,
                ε_i=ε_i,  # Adjust efficacy level
                p_i=params.p_i  # Proportion complying with intervention remains constant
            )

            # Run the model with the current parameters
            solution = run_intervention(current_params, S0, I0, SevI0, R0, tspan)

            # Extract compartment values over time
            susceptible = solution[1, :]
            infected = solution[2, :]
            severe_infected = solution[3, :]
            recovered = solution[4, :]

            # Plot in the corresponding subplot for each compartment
            plot!(grid_plot[i, j], solution.t, susceptible, label="Susceptible", lw=2, color=:blue)
            plot!(grid_plot[i, j], solution.t, infected, label="Infected", lw=2, color=:orange)
            plot!(grid_plot[i, j], solution.t, severe_infected, label="Severe", lw=2, color=:green)
            plot!(grid_plot[i, j], solution.t, recovered, label="Recovered", lw=2, color=:purple)
            # Set title for each subplot
            title!(grid_plot[i, j], "β=$β, ε_i=$ε_i")
        end
    end

    # Display the grid plot
    display(grid_plot)
end

function simulate_intervention_with_data_overlay()
    # Define model parameters and initial conditions
    β = 0.28
    γ = 1 / 7
    δ = 0.2
    γ_s = 1 / 14
    α = 1 / 30
    N = 6000
    intervention_day = 30
    ε_i = 0.3
    p_i = 0.5
    S0, I0, SevI0, R0 = Float64(N - 1), Float64(1), Float64(0), Float64(0)
    tspan = (0.0, 100.0)  # Full time span from 0 to 100 days

    # Observed data for infected and severe cases from day 30 to 80
    observed_infected = [
        155, 53, 67, 98, 130, 189, 92, 192, 145, 128, 68, 74, 126, 265, 154, 207,
        299, 273, 190, 152, 276, 408, 267, 462, 352, 385, 221, 420, 544, 329, 440,
        427, 369, 606, 416, 546, 475, 617, 593, 352, 337, 473, 673, 653, 523, 602,
        551, 686, 556, 600
    ]
    observed_severe = [
        22, 0, 15, 48, 38, 57, 9, 18, 20, 0, 41, 15, 35, 36, 27, 38, 24, 40, 34, 57,
        18, 29, 63, 66, 119, 76, 95, 28, 109, 136, 119, 104, 121, 93, 147, 129, 130,
        161, 133, 136, 138, 139, 181, 181, 218, 183, 167, 164, 219, 220
    ]

    # Set up model parameters
    params = SIRSInterventionParameters(
        β=β, γ=γ, δ=δ, γ_s=γ_s, α=α, N=N, intervention_day=intervention_day, ε_i=ε_i, p_i=p_i
    )

    # Run the model simulation
    solution = run_intervention(params, S0, I0, SevI0, R0, tspan)

    # Extract the full simulated infected and severe data for the time span 0 to 100
    model_infected = [solution(t)[2] for t in 0:100]  # Infectious compartment over time
    model_severe = [solution(t)[3] for t in 0:100]    # Severe infection compartment over time

    # Create the base plot for the model data over the full time span
    plt = plot(0:100, model_infected, label="Model Infected", lw=2, color=:blue)
    plot!(plt, 0:100, model_severe, label="Model Severe", lw=2, color=:red)

    # Overlay observed data for infected and severe cases from day 30 to 80
    plot!(plt, 30:79, observed_infected, seriestype=:scatter, label="Observed Infected", marker=:circle, color=:blue)
    plot!(plt, 30:79, observed_severe, seriestype=:scatter, label="Observed Severe", marker=:circle, color=:red)

    # Add labels and title
    xlabel!(plt, "Days")
    ylabel!(plt, "Population")
    title!(plt, "Model vs Observed Data for Infected and Severe Cases")

    # Display the combined plot
    display(plt)
end

# Error function to calculate the sum of squared differences
function calculate_total_error(solution, observed_infected, observed_severe, intervention_start)
    # Extract model predictions for infected and severe from day 30 onward
    model_infected = [solution(t)[2] for t in intervention_start:intervention_start+length(observed_infected)-1]
    model_severe = [solution(t)[3] for t in intervention_start:intervention_start+length(observed_severe)-1]

    # Calculate squared errors
    infected_error = sum((model_infected .- observed_infected) .^ 2)
    severe_error = sum((model_severe .- observed_severe) .^ 2)

    return infected_error + severe_error  # Total error
end

# Function to simulate the intervention with data overlay and find optimal compliance
function find_optimal_compliance()
    # Fixed model parameters
    β = 0.28
    β = 0.0361 * 8
    γ = 1 / 7
    δ = 0.2
    γ_s = 1 / 14
    α = 1 / 30
    N = 6000
    intervention_day = 30
    ε_i = 0.3
    tspan = (0.0, 100.0)  # Full time span from 0 to 100 days
    S0, I0, SevI0, R0 = Float64(N - 1), Float64(1), Float64(0), Float64(0)

    # Observed data for infected and severe cases from day 30 to 80
    observed_infected = [
        155, 53, 67, 98, 130, 189, 92, 192, 145, 128, 68, 74, 126, 265, 154, 207,
        299, 273, 190, 152, 276, 408, 267, 462, 352, 385, 221, 420, 544, 329, 440,
        427, 369, 606, 416, 546, 475, 617, 593, 352, 337, 473, 673, 653, 523, 602,
        551, 686, 556, 600
    ]
    observed_severe = [
        22, 0, 15, 48, 38, 57, 9, 18, 20, 0, 41, 15, 35, 36, 27, 38, 24, 40, 34, 57,
        18, 29, 63, 66, 119, 76, 95, 28, 109, 136, 119, 104, 121, 93, 147, 129, 130,
        161, 133, 136, 138, 139, 181, 181, 218, 183, 167, 164, 219, 220
    ]

    # Optimization function for compliance
    function objective_function(p_i)
        # Update parameters with current compliance level
        params = SIRSInterventionParameters(
            β=β, γ=γ, δ=δ, γ_s=γ_s, α=α, N=N, intervention_day=intervention_day, ε_i=ε_i, p_i=p_i
        )

        # Run the model simulation
        solution = run_intervention(params, S0, I0, SevI0, R0, tspan)

        # Calculate total error for this compliance level
        calculate_total_error(solution, observed_infected, observed_severe, intervention_day)
    end

    # Use optimization to find the best compliance level p_i (between 0 and 1)
    result = optimize(objective_function, 0.0, 1.0)
    best_compliance = Optim.minimizer(result)

    println("Optimal compliance level p_i: $best_compliance")

    # Run the simulation with the optimized compliance level and plot results
    params = SIRSInterventionParameters(
        β=β, γ=γ, δ=δ, γ_s=γ_s, α=α, N=N, intervention_day=intervention_day, ε_i=ε_i, p_i=best_compliance
    )
    solution = run_intervention(params, S0, I0, SevI0, R0, tspan)

    # Extract the full simulated infected and severe data for the time span 0 to 100
    model_infected = [solution(t)[2] for t in 0:100]
    model_severe = [solution(t)[3] for t in 0:100]

    # Create the base plot for the model data over the full time span
    plt = plot(0:100, model_infected, label="Model Infected", lw=2, color=:blue)
    plot!(plt, 0:100, model_severe, label="Model Severe", lw=2, color=:red)

    # Overlay observed data for infected and severe cases from day 30 to 80
    plot!(plt, 30:79, observed_infected, seriestype=:scatter, label="Observed Infected", marker=:circle, color=:blue)
    plot!(plt, 30:79, observed_severe, seriestype=:scatter, label="Observed Severe", marker=:circle, color=:red)

    # Add labels and title
    xlabel!(plt, "Days")
    ylabel!(plt, "Population")
    title!(plt, "Model vs Observed Data with Optimal Compliance Level")

    # Display the combined plot
    display(plt)

    return
end


function find_optimal_compliance_across_beta_range()
    # Fixed model parameters, excluding β and p_i, which vary
    γ = 1 / 7
    δ = 0.15
    γ_s = 1 / 14
    α = 1 / 30
    N = 6000
    intervention_day = 30
    ε_i = 0.3
    tspan = (0.0, 100.0)
    S0, I0, SevI0, R0 = Float64(N - 1), Float64(1), Float64(0), Float64(0)

    # Observed data for infected and severe cases from day 30 to 80
    observed_infected = [
        155, 53, 67, 98, 130, 189, 92, 192, 145, 128, 68, 74, 126, 265, 154, 207,
        299, 273, 190, 152, 276, 408, 267, 462, 352, 385, 221, 420, 544, 329, 440,
        427, 369, 606, 416, 546, 475, 617, 593, 352, 337, 473, 673, 653, 523, 602,
        551, 686, 556, 600
    ]
    observed_severe = [
        22, 0, 15, 48, 38, 57, 9, 18, 20, 0, 41, 15, 35, 36, 27, 38, 24, 40, 34, 57,
        18, 29, 63, 66, 119, 76, 95, 28, 109, 136, 119, 104, 121, 93, 147, 129, 130,
        161, 133, 136, 138, 139, 181, 181, 218, 183, 167, 164, 219, 220
    ]

    # Define a refined range for β and p_i
    β_values = 0.2:0.001:0.4
    p_i_values = 0.0:0.001:1.0

    # Dictionary to store the calculated error for each (β, p_i) pair
    error_map = Dict{Tuple{Float64, Float64}, Float64}()

    # Loop through each combination of β and p_i to calculate the error
    for β in β_values
        for p_i in p_i_values
            # Set up the parameters with current β and p_i
            params = SIRSInterventionParameters(
                β=β, γ=γ, δ=δ, γ_s=γ_s, α=α, N=N, intervention_day=intervention_day, ε_i=ε_i, p_i=p_i
            )

            # Run the model simulation
            solution = run_intervention(params, S0, I0, SevI0, R0, tspan)

            # Calculate the total error for the observed data and the model output
            error = calculate_total_error(solution, observed_infected, observed_severe, intervention_day)

            # Store the error for this (β, p_i) pair
            error_map[(β, p_i)] = error
        end
    end

    # Sort the (β, p_i) pairs by their error and select the top 5
    sorted_errors = sort(collect(error_map), by=x -> x[2])
    top_5_pairs = sorted_errors[1:5]

    # Print the top 5 best-fitting pairs
    println("Top 5 best-fitting (β, p_i) pairs with minimum error:")
    for ((β, p_i), error) in top_5_pairs
        β = β/8
        println("β = $β, p_i = $p_i, Error = $error")
    end

    # Plot the best-fitting results
    plt = plot(title="Top 5 Best-Fitting Model Results with Minimum Error", xlabel="Days", ylabel="Population")

    for ((β, p_i), _) in top_5_pairs
        # Set up the parameters for this best-fitting pair
        params = SIRSInterventionParameters(
            β=β, γ=γ, δ=δ, γ_s=γ_s, α=α, N=N, intervention_day=intervention_day, ε_i=ε_i, p_i=p_i
        )

        # Run the simulation for this pair
        solution = run_intervention(params, S0, I0, SevI0, R0, tspan)

        # Extract simulated data for infected and severe cases
        model_infected = [solution(t)[2] for t in 0:100]
        model_severe = [solution(t)[3] for t in 0:100]

        # Plot infected and severe cases for this pair
        β = β / 8
        plot!(plt, 0:100, model_infected, label="Model Infected (β=$β, p_i=$p_i)", lw=2, color=:blue)
        plot!(plt, 0:100, model_severe, label="Model Severe (β=$β, p_i=$p_i)", lw=2, color=:red)
    end

    # Overlay observed data for infected and severe cases from day 30 to 80
    plot!(plt, 30:79, observed_infected, seriestype=:scatter, label="Observed Infected", marker=:circle, color=:blue)
    plot!(plt, 30:79, observed_severe, seriestype=:scatter, label="Observed Severe", marker=:circle, color=:red)

    # Display the plot
    display(plt)

    return
end