"""
    SSE_error(model_output, data)

Computes the sum of squared errors (SSE) between model output and observed data.

# Arguments
- `model_output::Vector`: Model-simulated data.
- `data::Vector`: Observed real-world data.

# Returns
- Sum of squared errors as a Float64 value.
"""
function SSE_error(model_output::Vector, data::Vector)
    @assert length(model_output) == length(data) "Model output and data must have the same length."
    return sum((model_output .- data) .^ 2)
end

"""
    SIRS_model!(du, u, p, t)

Defines the SIRS model equations with an additional severe infection compartment.

# Arguments
- `du::Vector`: Derivative output vector for state variables.
- `u::Vector`: Current state variables `[S, I, SevI, R]`.
- `p::Vector`: Model parameters `[β, γ, δ, γ_s, α, N]`.
- `t::Float64`: Time variable (for compatibility, unused directly).

# Model Equations
Calculates the rate of change for Susceptible, Infected, Severe Infected, and Recovered compartments.
"""
# function SIRS_model!(du, u, p, t)
#     S, I, SevI, R = u
#     β, γ, δ, γ_s, α, N = p
#     λ = β * I / N

#     du[1] = -λ * S + α * R                     # dS/dt
#     du[2] = λ * S - γ * I                      # dI/dt
#     du[3] = δ * I * γ - γ_s * SevI             # dSevI/dt
#     du[4] = (1 - δ) * γ * I + γ_s * SevI - α * R # dR/dt
# end

"""
    simulate_error(β_values, c, γ, ps, γ_s, α, S0, I0, SevI0, R0, tspan, N; kwargs...)

Simulates and plots error across a range of β values, comparing model output with observed data.

# Arguments
- `β_values::Vector`: Range of β values.
- `c`: Scaling factor for β in plots.
- Model parameters and initial state values.
- Real data for time and infection counts (optional).

# Returns
- Combined error plots for regular and severe infections.
"""
function simulate_error(β_values::AbstractVector, c::Float64, γ::Float64, ps::Float64, γ_s::Float64, α::Float64,
    S0::Float64, I0::Float64, SevI0::Float64, R0::Float64, tspan::Tuple{Float64,Float64}, N::Int;
    real_data_time::AbstractVector{Int}=[], real_data_infected::AbstractVector{Int}=[],
    real_data_severe_time::AbstractVector{Int}=[], real_data_severe_infected::AbstractVector{Int}=[])
    errors, errors_severe = Float64[], Float64[]

    for β in β_values
        solution = run_SIRS_model(β, γ, ps * γ, γ_s, α, N, S0 * N, I0 * N, SevI0 * N, R0 * N, tspan)
        model_infected = [solution(t)[2] for t in real_data_time]
        model_severe_infected = [solution(t)[3] for t in real_data_severe_time]
        push!(errors, SSE_error(model_infected, real_data_infected))
        push!(errors_severe, SSE_error(model_severe_infected, real_data_severe_infected))
    end

    plt1 = plot(β_values / c, errors, xlabel="Transmission rate (β)", ylabel="Error", title="Infected Error", lw=2)
    plt2 = plot(β_values / c, errors_severe, xlabel="Transmission rate (β)", ylabel="Error", title="Severe Error", lw=2)
    combined_plot = plot(plt1, plt2, layout=(2, 1), size=(800, 600))
    return combined_plot
end

"""
    optimise_beta(γ, δ, γ_s, α, N, S0, I0, SevI0, R0, data, β_values)

Finds the optimal transmission rate β by minimizing the error between model output and observed data.

# Arguments
- `γ, δ, γ_s, α, N`: Model parameters.
- `S0, I0, SevI0, R0`: Initial state values.
- `data::Vector`: Observed data values.
- `β_values::Vector`: Range of β values for optimization.

# Returns
- A tuple with the optimal β and the minimum error.
"""
function optimise_beta(γ, δ, γ_s, α, N, S0, I0, SevI0, R0, data, β_values)
    min_error, best_beta = Inf, 0.0
    tspan = (0.0, Float64(length(data) - 1))

    for β in β_values
        # Construct SIRSParameters struct with the current β value
        params = SIRSParameters(β, γ, δ, γ_s, α, N)

        # Scale initial conditions to population size
        solution = run_SIRS(params, Int64(S0 * N), Int64(I0 * N), Int64(SevI0 * N), Int64(R0 * N), tspan)

        # Calculate the error between model output and observed data
        current_error = SSE_error(solution[2, :], data)

        # Update the best β if the current error is lower
        if current_error < min_error
            min_error, best_beta = current_error, β
        end
    end
    return best_beta, min_error
end

"""
    optimise_ps_and_beta(β_values, c, γ, ps_range, γ_s, α, S0, I0, SevI0, R0, tspan, N; kwargs...)

Optimises both β and p_s to minimize error in infected and severe cases.

# Arguments
- `β_values`, `ps_range`: Transmission rates and proportion severe values to test.
- Model parameters and initial states.
- Real data and time points for comparison.

# Returns
- Combined plot of error across β and p_s, with optimized values.
"""
function optimise_ps_and_beta(β_values::AbstractVector, c::Float64, γ::Float64, ps_range::AbstractVector, γ_s::Float64, α::Float64,
    S0::Float64, I0::Float64, SevI0::Float64, R0::Float64, tspan::Tuple{Float64,Float64}, N::Int;
    real_data_time::AbstractVector{Int}=[], real_data_infected::AbstractVector{Int}=[],
    real_data_severe_time::AbstractVector{Int}=[], real_data_severe_infected::AbstractVector{Int}=[])

    best_ps, best_beta, min_total_error = 0.0, 0.0, Inf
    errors_ps = Float64[]

    for ps in ps_range
        for β in β_values
            solution = run_SIRS_model(β, γ, ps * γ, γ_s, α, N, S0 * N, I0 * N, SevI0 * N, R0 * N, tspan)
            model_infected = [solution(t)[2] for t in real_data_time]
            model_severe_infected = [solution(t)[3] for t in real_data_severe_time]

            combined_error = SSE_error(model_infected, real_data_infected) + SSE_error(model_severe_infected, real_data_severe_infected)
            if combined_error < min_total_error
                min_total_error, best_ps, best_beta = combined_error, ps, β
            end
        end
        push!(errors_ps, min_total_error)
    end

    plot(ps_range, errors_ps, xlabel="Proportion Severe (p_s)", ylabel="Total Error", title="Error vs. p_s", lw=2)
end
