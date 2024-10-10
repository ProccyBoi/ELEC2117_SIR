# Define the SIRS model equations with severe infection
function SIRS_model!(du, u, params, t)
    S, I, SevI, R = u  # Current values of S, I, SevI, and R compartments
    β, γ, δ, γ_s, α, N = params  # Parameters

    # Calculate force of infection: λ = β * I / N
    λ = β * I / N

    # SIRS model equations with severe infections
    du[1] = -λ * S + α * R                     # dS/dt
    du[2] = λ * S - γ * I                      # dI/dt
    du[3] = δ * I - γ_s * SevI                 # dSevI/dt
    du[4] = γ * (I - δ * I) + γ_s * SevI - α * R  # dR/dt
end

# Function to run the modified SIRS model
function run_SIRS_model(β, γ, δ, γ_s, α, N, S0, I0, SevI0, R0, tspan)
    # Initial conditions
    u0 = [S0, I0, SevI0, R0]  # Initial values for S, I, SevI, and R

    # Define the ODE problem
    prob = ODEProblem(SIRS_model!, u0, tspan, [β, γ, δ, γ_s, α, N])

    # Solve the problem with `saveat` ensuring one point per day
    solution = solve(prob, Tsit5(), saveat=collect(tspan[1]:1:tspan[2]))  # Generate solution at every day (integer)

    return solution
end


# Function to plot the SIRS model solution
function plot_SIRS(solution)
    plt = plot(solution, xlabel="Time (days)", ylabel="Number of People", label=["Susceptible" "Infectious" "Severe Illness" "Recovered"], lw=2)
    return plt
end

# Function to simulate the SIRS model using the same parameters as the Python code
function simulate_SIRS()
    # Parameters
    N = 6000          # Total population
    β = 0.25           # Infection rate (Beta * c)
    γ = 1 / 7           # Recovery rate for infected individuals
    ps = 0.2          # Proportion developing severe illness
    δ = ps * γ        # Transition rate to severe illness
    γ_s = 1 / 14        # Recovery rate for severe illness
    α = 1 / 30          # Rate of return to susceptibility

    # Initial conditions: [S, I, SevI, R]
    S0, I0, SevI0, R0 = 5999, 1, 0, 0

    # Time span
    tspan = (0.0, 160.0)

    # Run the model
    solution = run_SIRS_model(β, γ, δ, γ_s, α, N, S0, I0, SevI0, R0, tspan)

    # Plot the results
    plt = plot_SIRS(solution)
    display(plt)

    println("Simulation complete!")
end

# Function to plot the model against data
function data_fitting(β, γ, δ, γ_s, α, N, S0, I0, SevI0, R0, data_infected, data_severe)
    # Define the time span based on the length of the data points
    tspan = (0.0, length(data_infected) - 1)  # Time span: 0 to length of data - 1

    # Create a corresponding time vector for the data points
    real_data_time = collect(tspan[1]:tspan[2])
    real_data_infected = data_infected
    real_data_severe = data_severe

    # Initial conditions
    u0 = [S0 * N, I0 * N, SevI0 * N, R0 * N]  # Scale initial proportions to absolute values

    # Define the ODE problem
    prob = ODEProblem(SIRS_model!, u0, tspan, [β, γ, δ, γ_s, α, N])

    # Solve the problem with a dense solution for smooth plotting
    solution = solve(prob, Tsit5(), saveat=0.1)

    # Extract the time points and compartment values from the model solution
    model_time = solution.t
    model_infected = solution[2, :]  # Extract the infected individuals over time
    model_severe = solution[3, :]    # Extract the severe illness individuals over time

    # Plot the model output vs. the real data for both Infected and Severe Illness
    plt = plot(model_time, model_infected, label="Model Infected", xlabel="Time (Days)", ylabel="Number of Individuals",
        title="Model vs. Real Data", lw=2, color=:blue)
    scatter!(real_data_time, real_data_infected, label="Data Infected", marker=:circle, ms=3, color=:blue)

    plot!(model_time, model_severe, label="Model Severe Illness", lw=2, color=:purple)
    scatter!(real_data_time, real_data_severe, label="Data Severe Illness", marker=:diamond, ms=4, color=:purple)

    return plt
end

# data_infected = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 7, 20, 3, 29, 14, 11, 12, 16, 10, 58]
# data_severe = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 5]
# data_fitting(0.28, 1/7, 0.2 * (1/7), 1/14, 1/30, 6000, 5999/6000, 1/6000, 0, 0, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 7, 20, 3, 29, 14, 11, 12, 16, 10, 58], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 5])
