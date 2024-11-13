using DifferentialEquations, Plots

# Differential Equations for the SIR Model with Intervention
function intervention!(du, u, params::SIRSInterventionParameters, t)
    S, I, SevI, R = u
    β_eff = t >= params.intervention_day ? params.β * (1 - params.ε_i * params.p_i) : params.β
    λ = β_eff * I / params.N  # Force of infection

    # Differential equations, identical to SIRS_model! but with dynamic β
    du[1] = -λ * S + params.α * R                         # Susceptible (dS/dt)
    du[2] = λ * S - params.γ * I                          # Infectious (dI/dt)
    du[3] = params.δ * I * params.γ - params.γ_s * SevI   # Severe Infection (d(SevI)/dt)
    du[4] = (1 - params.δ) * params.γ * I + params.γ_s * SevI - params.α * R  # Recovered (dR/dt)
end

# Function to Run Simulation for the SIR Model with Intervention
function run_intervention(params::SIRSInterventionParameters, S0, I0, SevI0, R0, tspan)
    u0 = [S0, I0, SevI0, R0]
    prob = ODEProblem(intervention!, u0, tspan, params)
    solution = solve(prob, Tsit5(), saveat=1.0)
    return solution
end

# Function 1: Analyse Pathogen Variant Similarity
function analyse_variant_similarity(observed_town1_infected, observed_town1_severe, observed_town2_infected, observed_town2_severe)
    β_values = 0.1:0.005:0.4
    δ_values = 0.15:0.005:0.25
    best_error_town1, best_error_town2 = Inf, Inf
    best_β1, best_δ1, best_β2, best_δ2 = 0.0, 0.0, 0.0, 0.0

    for β in β_values, δ in δ_values
        # Calibrate for Town 1 with both infected and severely infected data
        error_town1 = calculate_sir_error(β, δ, observed_town1_infected, observed_town1_severe)
        if error_town1 < best_error_town1
            best_error_town1 = error_town1
            best_β1, best_δ1 = β, δ
        end

        # Calibrate for Town 2 with both infected and severely infected data
        error_town2 = calculate_sir_error(β, δ, observed_town2_infected, observed_town2_severe)
        if error_town2 < best_error_town2
            best_error_town2 = error_town2
            best_β2, best_δ2 = β, δ
        end
    end

    println("Best β for Town 1: $best_β1, Best δ for Town 1: $best_δ1")
    println("Best β for Town 2: $best_β2, Best δ for Town 2: $best_δ2")

    # Plot Town 2 results
    solution2 = run_simulation(0, best_β2, best_δ2, 0.3, 0.8)
    plt2 = plot(title="Town 2 Fit", xlabel="Days", ylabel="Population")
    model_infected_town2 = [solution2(t)[2] for t in 0:100]
    model_severe_town2 = [solution2(t)[3] for t in 0:100]
    plot!(plt2, 0:100, model_infected_town2, label="Model Infected", lw=2, color=:blue)
    plot!(plt2, 0:100, model_severe_town2, label="Model Severe", lw=2, color=:red)
    plot!(plt2, 30:79, observed_town2_infected, seriestype=:scatter, label="Observed Infected", marker=:circle, color=:blue)
    plot!(plt2, 30:79, observed_town2_severe, seriestype=:scatter, label="Observed Severe", marker=:circle, color=:red)
    display(plt2)

    β_difference = abs(best_β1 - best_β2)
    δ_difference = abs(best_δ1 - best_δ2)
    tolerance = 0.05
    return β_difference < tolerance && δ_difference < tolerance
end

# Helper to Calculate SIR Error for Variant Similarity
function calculate_sir_error(β, δ, observed_infected, observed_severe)
    γ, γ_s, α = 1 / 7, 1 / 14, 1 / 30
    N, S0, I0, SevI0, R0 = 10000, 9999.0, 1.0, 0.0, 0.0
    tspan = (0.0, 100.0)

    params = SIRSInterventionParameters(β=β, γ=γ, δ=δ, γ_s=γ_s, α=α, N=N, intervention_day=37, ε_i=0.3, p_i=0.8)
    solution = run_intervention(params, S0, I0, SevI0, R0, tspan)

    # Compute error for infected and severe cases separately, then sum
    model_infected = [solution(t)[2] for t in 27:length(observed_infected)+26]
    model_severe = [solution(t)[3] for t in 27:length(observed_severe)+26]
    infected_error = sum((observed_infected .- model_infected) .^ 2)
    severe_error = sum((observed_severe .- model_severe) .^ 2)

    return infected_error + severe_error
end

# Function 2: Estimate Infection Entry Day
function estimate_infection_entry_day(observed_infected, observed_severe)
    β, δ = 0.295, 0.15
    best_day, min_error = -1, Inf
    ε_i, p_i = 0.3, 0.8

    for entry_day in 1:27
        solution = run_simulation(entry_day, β, δ, ε_i, p_i)
        error = calculate_total_error(solution, observed_infected, observed_severe, entry_day)

        if error < min_error
            min_error = error
            best_day = entry_day
        end
    end

    println("Estimated infection entry day: $best_day with error $min_error")
    solution_best = run_simulation(best_day, β, δ, ε_i, p_i)

    # Plotting results
    plt = plot(title="Estimated Entry Day Fit", xlabel="Days", ylabel="Population")
    model_infected = [solution_best(t)[2] for t in 0:100]
    model_severe = [solution_best(t)[3] for t in 0:100]
    plot!(plt, 0:100, model_infected, label="Model Infected", lw=2, color=:blue)
    plot!(plt, 0:100, model_severe, label="Model Severe", lw=2, color=:red)
    display(plt)

    return best_day, min_error
end

# Function 3: Evaluate Hypothetical Early Intervention
function evaluate_intervention_hypothetical(observed_infected, observed_severe, entry_day; threshold_value=1e7)
    β, δ = 0.259, 0.15
    ε_i, p_i = 0.3, 0.8
    intervention_day = entry_day

    solution = run_simulation(intervention_day, β, δ, ε_i, p_i)
    error = calculate_total_error(solution, observed_infected, observed_severe, intervention_day)

    println("With early intervention on day $intervention_day, error compared to observed data is $error")

    # Plotting results
    plt = plot(title="Hypothetical Early Intervention Fit", xlabel="Days", ylabel="Population")
    model_infected = [solution(t)[2] for t in 0:100]
    model_severe = [solution(t)[3] for t in 0:100]
    plot!(plt, 0:100, model_infected, label="Model Infected", lw=2, color=:blue)
    plot!(plt, 0:100, model_severe, label="Model Severe", lw=2, color=:red)
    display(plt)

    return error < threshold_value
end

# Helper to Calculate SIR Error for Variant Similarity
function calculate_sir_error(β, δ, observed_data)
    γ, γ_s, α = 1 / 7, 1 / 14, 1 / 30
    N, S0, I0, SevI0, R0 = 10000, 9999.0, 1.0, 0.0, 0.0
    tspan = (0.0, 100.0)

    params = SIRSInterventionParameters(β=β, γ=γ, δ=δ, γ_s=γ_s, α=α, N=N, intervention_day=37, ε_i=0.3, p_i=0.8)
    solution = run_intervention(params, S0, I0, SevI0, R0, tspan)

    model_infected = [solution(t)[2] for t in 27:length(observed_data)+26]
    error = sum((observed_data .- model_infected) .^ 2)
    return error
end

# Helper to Calculate Total Error Against Observed Data
function calculate_total_error(solution, observed_infected, observed_severe, entry_day)
    model_infected = [solution(t)[2] for t in entry_day:(entry_day+length(observed_infected)-1)]
    model_severe = [solution(t)[3] for t in entry_day:(entry_day+length(observed_severe)-1)]

    infected_error = sum((observed_infected .- model_infected) .^ 2)
    severe_error = sum((observed_severe .- model_severe) .^ 2)
    total_error = infected_error + severe_error
    return total_error
end

# Simulation Function for Town Modelling with Interventions
function run_simulation(entry_day, β, δ, ε_i, p_i)
    γ, γ_s, α = 1 / 7, 1 / 14, 1 / 30
    N = 10000
    tspan = (0.0, 100.0)

    # Initial conditions: assume 1 initial case entering on `entry_day`
    S0, I0, SevI0, R0 = N - 1, 1.0, 0.0, 0.0
    params = SIRSInterventionParameters(β=β, γ=γ, δ=δ, γ_s=γ_s, α=α, N=N, intervention_day=entry_day, ε_i=ε_i, p_i=p_i)

    # Solve the differential equations with intervention dynamics
    solution = run_intervention(params, S0, I0, SevI0, R0, tspan)
    return solution
end
