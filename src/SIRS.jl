# Define a struct to hold the model parameters for the SIRS model
struct SIRSParameters
    β::Float64      # Transmission rate: rate at which susceptible individuals become infected
    γ::Float64      # Recovery rate for mild cases: rate at which infectious individuals recover
    δ::Float64      # Proportion progressing to severe infection: fraction of infections that become severe
    γ_s::Float64    # Recovery rate for severe cases: rate at which severe cases recover
    α::Float64      # Re-susceptibility rate after recovery: rate at which recovered individuals become susceptible again
    N::Int          # Total population: size of the population being modeled
end

# Keyword-based constructor for SIRSParameters to initialize model parameters more flexibly
function SIRSParameters(; β::Float64, γ::Float64, δ::Float64, γ_s::Float64, α::Float64, N::Int)
    return SIRSParameters(β, γ, δ, γ_s, α, N)
end

# Differential equation function for the SIRS model including severe infections
"""
    SIRS_model!(du, u, p, t)

Defines the SIRS model equations with an additional severe infection compartment.

# Arguments
- `du::Vector`: Output vector for derivatives of state variables.
- `u::Vector`: Current state variables `[S, I, SevI, R]`.
- `p::SIRSParameters`: Model parameters.
- `t::Float64`: Time variable (included for ODE compatibility).

# Model Equations
Calculates the rate of change for each compartment: Susceptible, Infectious, Severe Infection, and Recovered.
"""
function SIRS_model!(du, u, p::SIRSParameters, t)
    # Unpack current state variables
    S, I, SevI, R = u
    # Unpack model parameters from the struct
    β, γ, δ, γ_s, α, N = p.β, p.γ, p.δ, p.γ_s, p.α, p.N

    # Force of infection (transmission rate scaled by infectious individuals)
    λ = β * I / N

    # Define the differential equations for each compartment
    du[1] = -λ * S + α * R                         # Susceptible (S)
    du[2] = λ * S - γ * I                          # Infectious (I)
    du[3] = δ * I * γ - γ_s * SevI                 # Severe Infection (SevI)
    du[4] = (1 - δ) * γ * I + γ_s * SevI - α * R   # Recovered (R)
end

# Function to run the SIRS model simulation with specified initial conditions and parameters
"""
    run_SIRS(params, S0, I0, SevI0, R0, tspan)

Runs the SIRS model simulation with given parameters and initial conditions.

# Arguments
- `params::SIRSParameters`: Struct containing model parameters.
- `S0::Int`: Initial count of susceptible individuals.
- `I0::Int`: Initial count of infectious individuals.
- `SevI0::Int`: Initial count of severe infection cases.
- `R0::Int`: Initial count of recovered individuals.
- `tspan::Tuple{Float64, Float64}`: Time span for the simulation (start and end times).

# Returns
- `solution`: Solution object from DifferentialEquations.jl for further analysis or plotting.
"""
function run_SIRS(params::SIRSParameters, S0::Int, I0::Int, SevI0::Int, R0::Int, tspan::Tuple{Float64,Float64})
    # Set initial conditions as a vector with compartment values
    u0 = [S0, I0, SevI0, R0]

    # Define the ODE problem with SIRS model equations
    prob = ODEProblem(SIRS_model!, u0, tspan, params)

    # Solve the ODE problem using a solver with daily intervals (saveat=1.0)
    solution = solve(prob, Tsit5(), saveat=1.0)
    return solution
end

# Function to plot results from the SIRS model simulation
"""
    plot_SIRS(solution; title="SIRS Model Results", ylabel="Population", color_scheme=:default)

Plots the solution of the SIRS model, showing the dynamics of each compartment.

# Arguments
- `solution`: Solution object from DifferentialEquations.jl.
- `title`: Title of the plot (optional).
- `ylabel`: Label for the Y-axis (optional).
- `color_scheme`: Color scheme for the plot (optional).

# Returns
- `plt`: Plot object.
"""
function plot_SIRS(solution; title="SIRS Model Results", ylabel="Population", color_scheme=:default)
    # Plot each compartment in the solution with appropriate labels and colors
    plt = plot(solution.t, solution[1, :], label="Susceptible", xlabel="Time (days)", ylabel=ylabel, lw=2, color=:blue)
    plot!(plt, solution.t, solution[2, :], label="Infectious", lw=2, color=:orange)
    plot!(plt, solution.t, solution[3, :], label="Severe Illness", lw=2, color=:green)
    plot!(plt, solution.t, solution[4, :], label="Recovered", lw=2, color=:purple)
    title!(plt, title)
    return plt
end

# Example function to simulate and plot the SIRS model with default parameters
"""
    simulate_SIRS()

Simulates and plots the SIRS model using default parameters.
"""
function simulate_SIRS()
    # Define default model parameters for the SIRS simulation
    params = SIRSParameters(
        β=0.28,             # Transmission rate
        γ=1 / 7,            # Recovery rate for infectious individuals
        δ=0.2,              # Fraction of infections progressing to severe illness
        γ_s=1 / 14,         # Recovery rate for severe cases
        α=1 / 30,           # Rate of becoming susceptible again after recovery
        N=6000              # Total population size
    )

    # Set initial conditions (mostly susceptible with one infectious case)
    S0, I0, SevI0, R0 = params.N - 1, 1, 0, 0

    # Time span for simulation (0 to 100 days)
    tspan = (0.0, 100.0)

    # Run the model and obtain solution
    solution = run_SIRS(params, S0, I0, SevI0, R0, tspan)

    # Plot the results with custom title and labels
    plt = plot_SIRS(solution, title="SIRS Model with Population Dynamics", ylabel="Number of People")
    display(plt)
end