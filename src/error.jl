# Calculate the error between model output and real data
"""
    calculate_error(model_output, data)

Computes the sum of squared errors (SSE) between the model output and real data.

# Arguments
- `model_output::Vector`: The output from the model simulation.
- `data::Vector`: Observed data values.

# Returns
- `Float64`: The sum of squared errors.
"""
function calculate_error(model_output::Vector, data::Vector)
    # Ensure that the lengths of the two vectors are the same
    @assert length(model_output) == length(data) "Model output and data must have the same length."

    # Calculate the sum of squared errors
    return sum((model_output .- data) .^ 2)
end

# Define the SIRS model equations with severe infection
"""
    SIRS_model!(du, u, p, t)

Differential equation for the SIRS model with severe infection.

# Arguments
- `du::Vector`: Derivative of state variables (output).
- `u::Vector`: Current state variables `[S, I, SevI, R]`.
- `p::Vector`: Parameters `[β, γ, δ, γ_s, α, N]`.
- `t::Float64`: Time variable.

# Model Equations
The equations model the rates of change for Susceptible, Infectious, Severe Infection, and Recovered populations.
"""
function SIRS_model!(du, u, p, t)
    S, I, SevI, R = u
    β, γ, δ, γ_s, α, N = p

    # Force of infection
    λ = β * I / N

    # Differential equations
    du[1] = -λ * S + α * R                         # dS/dt: Susceptible
    du[2] = λ * S - γ * I                          # dI/dt: Infectious
    du[3] = δ * I * γ - γ_s * SevI                 # d(SevI)/dt: Severe Infection
    du[4] = (1 - δ) * γ * I + γ_s * SevI - α * R   # dR/dt: Recovered
end

# Optimizes the β parameter by minimizing the error between model and observed data
"""
    optimise_beta(γ, δ, γ_s, α, N, S0, I0, SevI0, R0, data, β_values)

Finds the optimal value of β that minimizes the error between model output and observed data.

# Arguments
- `γ, δ, γ_s, α, N`: Model parameters.
- `S0, I0, SevI0, R0`: Initial conditions for state variables.
- `data::Vector`: Observed data values to fit the model against.
- `β_values::Vector`: Range of β values to evaluate.

# Returns
- `Tuple`: The best β and the minimum error.
"""
function optimise_beta(γ, δ, γ_s, α, N, S0, I0, SevI0, R0, data, β_values)
    min_error = Inf
    best_beta = 0.0
    tspan = (0.0, length(data) - 1)

    for β in β_values
        solution = run_SIRS_model(β, γ, δ, γ_s, α, N, S0 * N, I0 * N, SevI0 * N, R0 * N, tspan)
        model_output = solution[2, :]  # Infectious compartment

        current_error = calculate_error(model_output, data)
        if current_error < min_error
            min_error = current_error
            best_beta = β
        end
    end

    return best_beta, min_error
end

# Simulates error across a range of β values and compares model output with real data
"""
    simulate_error(β_values, c, γ, ps, γ_s, α, S0, I0, SevI0, R0, tspan, N; kwargs...)

Simulates and plots error for a range of β values by comparing model output with real data.

# Arguments
- `β_values::Vector`: Transmission rates to evaluate.
- `c::Float64`: Scaling factor for β values on the plot.
- Other model parameters and initial conditions.
- `real_data_time::Vector`: Time points of real data.
- `real_data_infected::Vector`: Observed infected data.
- `real_data_severe_time::Vector`: Time points for severe cases.
- `real_data_severe_infected::Vector`: Observed severe infection data.

# Returns
- `Plot`: Combined error plots.
"""
function simulate_error(β_values::AbstractVector, c::Float64, γ::Float64, ps::Float64, γ_s::Float64, α::Float64,
    S0::Float64, I0::Float64, SevI0::Float64, R0::Float64, tspan::Tuple{Float64,Float64}, N::Int;
    real_data_time::AbstractVector{Int} = [], real_data_infected::AbstractVector{Int} = [],
    real_data_severe_time::AbstractVector{Int} = [], real_data_severe_infected::AbstractVector{Int} = [])

    errors, errors_severe = Float64[], Float64[]

    for β in β_values
        solution = run_SIRS_model(β, γ, ps * γ, γ_s, α, N, S0 * N, I0 * N, SevI0 * N, R0 * N, tspan)
        @assert solution.retcode == :Success "Simulation failed for β = $β."

        model_infected = [solution(t)[2] for t in real_data_time]
        model_severe_infected = [solution(t)[3] for t in real_data_severe_time]

        push!(errors, calculate_error(model_infected, real_data_infected))
        push!(errors_severe, calculate_error(model_severe_infected, real_data_severe_infected))
    end

    # Plot the error as a function of β
    plt1 = plot(β_values / c, errors, xlabel="Transmission rate (β)", ylabel="Error",
        title="Error vs. β for SIRS Model (Infected)", lw=2, marker=:circle)
    plt2 = plot(β_values / c, errors_severe, xlabel="Transmission rate (β)", ylabel="Error",
        title="Error vs. β for SIRS Model (Severe Infected)", lw=2, marker=:circle)
    plt3 = plot(β_values / c, errors .+ errors_severe, xlabel="Transmission rate (β)", ylabel="Error",
        title="Error vs. β for SIRS Model (Combined)", lw=2, marker=:circle)

    combined_plot = plot(plt1, plt2, plt3, layout=(3, 1), size=(800, 600))
    return combined_plot
end

# Optimizes both β and p_s by minimizing error for infected and severe cases
"""
    optimize_ps_and_beta(β_values, c, γ, ps_range, γ_s, α, S0, I0, SevI0, R0, tspan, N; kwargs...)

Optimizes both β and p_s by finding the combination that minimizes error.

# Arguments
- Similar to `simulate_error`, with additional parameters `ps_range` and optimization outputs.

# Returns
- Tuple containing combined plots and optimized values.
"""
function optimize_ps_and_beta(β_values::AbstractVector, c::Float64, γ::Float64, ps_range::AbstractVector, γ_s::Float64, α::Float64,
    S0::Float64, I0::Float64, SevI0::Float64, R0::Float64, tspan::Tuple{Float64,Float64}, N::Int;
    real_data_time::AbstractVector{<:Number} = [], real_data_infected::AbstractVector{<:Number} = [],
    real_data_severe_time::AbstractVector{<:Number} = [], real_data_severe_infected::AbstractVector{<:Number} = [])

    best_ps, best_beta, min_total_error = 0.0, 0.0, Inf
    errors_ps, best_errors_infected, best_errors_severe, best_combined_errors = Float64[], Float64[], Float64[], Float64[]

    for ps in ps_range
        errors_infected, errors_severe = Float64[], Float64[]

        for β in β_values
            solution = run_SIRS_model(β, γ, ps * γ, γ_s, α, N, S0 * N, I0 * N, SevI0 * N, R0 * N, tspan)
            model_infected = [solution(t)[2] for t in real_data_time]
            model_severe_infected = [solution(t)[3] for t in real_data_severe_time]

            err_infected = calculate_error(model_infected, real_data_infected)
            err_severe = calculate_error(model_severe_infected, real_data_severe_infected)
            combined_error = err_infected + err_severe

            push!(errors_infected, err_infected)
            push!(errors_severe, err_severe)

            if combined_error < min_total_error
                min_total_error, best_ps, best_beta = combined_error, ps, β
                best_errors_infected, best_errors_severe = deepcopy(errors_infected), deepcopy(errors_severe)
                best_combined_errors = [i + j for (i, j) in zip(errors_infected, errors_severe)]
            end
        end

        push!(errors_ps, sum(errors_infected) + sum(errors_severe))
    end

    plt_ps = plot(ps_range, errors_ps, xlabel="Proportion Severe (p_s)", ylabel="Total Error",
        title="Total Error vs. p_s for SIRS Model", lw=2, marker=:circle)
    println("Optimized p_s: ", best_ps, ", Optimized β: ", best_beta, " with minimum total error: ", min_total_error)

    combined_plot = plot(plot(β_values / c, best_errors_infected, xlabel="Transmission rate (β)", ylabel="Error",
        title="Error vs. β (Infected)", lw=2, marker=:circle),
        plot(β_values / c, best_errors_severe, xlabel="Transmission rate (β)", ylabel="Error",
        title="Error vs. β (Severe Infected)", lw=2, marker=:circle),
        plot(β_values / c, best_combined_errors, xlabel="Transmission rate (β)", ylabel="Error",
        title="Error vs. β (Combined)", lw=2, marker=:circle), layout=(3, 1), size=(800, 600))

    return combined_plot, plt_ps, best_ps, best_beta, min_total_error
end
