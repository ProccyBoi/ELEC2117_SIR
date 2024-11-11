# Define a struct to hold the model parameters
struct SIRSParameters
    β::Float64      # Transmission rate
    γ::Float64      # Recovery rate for mild cases
    δ::Float64      # Proportion progressing to severe infection
    γ_s::Float64    # Recovery rate for severe cases
    α::Float64      # Re-susceptibility rate after recovery
    N::Int          # Total population
end

# Define a keyword-based constructor for SIRSParameters
function SIRSParameters(; β::Float64, γ::Float64, δ::Float64, γ_s::Float64, α::Float64, N::Int)
    return SIRSParameters(β, γ, δ, γ_s, α, N)
end

# Define the SIRS model with absolute population values
"""
    SIRS_model!(du, u, p, t)

Differential equation for the SIRS model with severe infection.

# Arguments
- `du::Vector`: Derivative of state variables (output).
- `u::Vector`: Current state variables `[S, I, SevI, R]`.
- `p::SIRSParameters`: Model parameters.
- `t::Float64`: Time variable (not used directly).

# Model Equations
The equations model the rates of change for Susceptible, Infectious, Severe Infection, and Recovered populations.
"""
function SIRS_model!(du, u, p::SIRSParameters, t)
    S, I, SevI, R = u
    β, γ, δ, γ_s, α, N = p.β, p.γ, p.δ, p.γ_s, p.α, p.N

    # Force of infection
    λ = β * I / N

    # Differential equations
    du[1] = -λ * S + α * R                         # Susceptible
    du[2] = λ * S - γ * I                          # Infectious
    du[3] = δ * I * γ - γ_s * SevI                 # Severe Infection
    du[4] = (1 - δ) * γ * I + γ_s * SevI - α * R   # Recovered
end

# Function to run the SIRS model
"""
    run_SIRS_model(params, S0, I0, SevI0, R0, tspan)

Runs the SIRS model simulation with the given parameters and initial conditions.

# Arguments
- `params::SIRSParameters`: Struct containing model parameters.
- `S0::Int`: Initial susceptible count.
- `I0::Int`: Initial infected count.
- `SevI0::Int`: Initial severe infection count.
- `R0::Int`: Initial recovered count.
- `tspan::Tuple{Float64, Float64}`: Time span for the simulation.

# Returns
- `solution`: Solution object from DifferentialEquations.jl.
"""
function run_SIRS_model(params::SIRSParameters, S0::Int, I0::Int, SevI0::Int, R0::Int, tspan::Tuple{Float64, Float64})
    # Initial conditions as absolute counts
    u0 = [S0, I0, SevI0, R0]

    # Define the ODE problem
    prob = ODEProblem(SIRS_model!, u0, tspan, params)

    # Solve the problem with daily intervals
    solution = solve(prob, Tsit5(), saveat=1.0)
    return solution
end

# Function to plot the results
"""
    plot_SIRS_model(solution; title="SIRS Model Results", ylabel="Population", color_scheme=:default)

Plots the solution of the SIRS model.

# Arguments
- `solution`: Solution object from DifferentialEquations.jl.
- `title`: Title of the plot (optional).
- `ylabel`: Y-axis label (optional).
- `color_scheme`: Color scheme for the plot (optional).

# Returns
- `plt`: Plot object.
"""
function plot_SIRS_model(solution; title="SIRS Model Results", ylabel="Population", color_scheme=:default)
    plt = plot(solution.t, solution[1, :], label="Susceptible", xlabel="Time (days)", ylabel=ylabel, lw=2, color=:blue)
    plot!(plt, solution.t, solution[2, :], label="Infectious", lw=2, color=:orange)
    plot!(plt, solution.t, solution[3, :], label="Severe Illness", lw=2, color=:green)
    plot!(plt, solution.t, solution[4, :], label="Recovered", lw=2, color=:purple)
    title!(plt, title)
    return plt
end

# Example usage function
"""
    simulate_SIRS()

Simulates and plots the SIRS model with default parameters.
"""
function simulate_SIRS()
    # Define default parameters using keyword arguments
    params = SIRSParameters(
        β = 0.28,             # Transmission rate
        γ = 1 / 7,            # Recovery rate for infectious individuals
        δ = 0.2,              # Proportion progressing to severe infection
        γ_s = 1 / 14,         # Recovery rate for severe cases
        α = 1 / 30,           # Re-susceptibility rate after recovery
        N = 6000              # Total population
    )

    # Initial conditions
    S0, I0, SevI0, R0 = params.N - 1, 1, 0, 0

    # Time span
    tspan = (0.0, 1000.0)  # Simulate for 1000 days

    # Run model and plot results
    solution = run_SIRS_model(params, S0, I0, SevI0, R0, tspan)
    plt = plot_SIRS_model(solution, title="SIRS Model with Population Dynamics", ylabel="Number of People")
    display(plt)
end

# Uncomment the line below to run an example simulation and plot the results
# simulate_SIRS()
