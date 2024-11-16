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
    simulate_error(β_values, c, γ, δ, γ_s, α, S0, I0, SevI0, R0, tspan, N; kwargs...)

Simulates and plots error across a range of β values, comparing model output with observed data.

# Arguments
- `β_values::Vector`: Range of β values.
- `c`: Scaling factor for β in plots.
- Model parameters and initial state values.
- Real data for time and infection counts (optional).

# Returns
- Combined error plots for regular and severe infections.
"""
function simulate_error(β_values::AbstractVector, c::Float64, γ::Float64, δ::Float64, γ_s::Float64, α::Float64,
    S0::Float64, I0::Float64, SevI0::Float64, R0::Float64, tspan::Tuple{Float64,Float64}, N::Int;
    real_data_time::AbstractVector{Int}, real_data_infected::AbstractVector{Int},
    real_data_severe_time::AbstractVector{Int}, real_data_severe_infected::AbstractVector{Int})

    errors, errors_severe = Float64[], Float64[]

    for β in β_values
        params = SIRSParameters(β, γ, δ, γ_s, α, N)
        solution = run_SIRS(params, round(Int, S0 * N), round(Int, I0 * N), round(Int, SevI0 * N), round(Int, R0 * N), tspan)
        model_infected = [solution(t)[2] for t in real_data_time]
        model_severe_infected = [solution(t)[3] for t in real_data_severe_time]

        push!(errors, SSE_error(model_infected, real_data_infected))
        push!(errors_severe, SSE_error(model_severe_infected, real_data_severe_infected))
    end

    if !isempty(errors) && !isempty(errors_severe)
        plt1 = plot(β_values / c, errors, xlabel="Transmission rate (β)", ylabel="Error",
            title="Error vs. β (Infected)", lw=2, marker=:circle)
        plt2 = plot(β_values / c, errors_severe, xlabel="Transmission rate (β)", ylabel="Error",
            title="Error vs. β (Severe Infected)", lw=2, marker=:circle)
        return plot(plt1, plt2, layout=(2, 1), size=(800, 600))
    else
        println("No successful simulations.")
        return nothing
    end
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
function optimise_beta(γ::Float64, δ::Float64, γ_s::Float64, α::Float64, N::Int, S0::Float64, I0::Float64, SevI0::Float64,
    R0::Float64, data::AbstractVector, β_values::AbstractVector)

    min_error, best_beta = Inf, 0.0
    tspan = (0.0, Float64(length(data) - 1))  # Ensure tspan is Tuple{Float64, Float64}

    for β in β_values
        params = SIRSParameters(β, γ, δ, γ_s, α, N)
        solution = run_SIRS(params, round(Int, S0 * N), round(Int, I0 * N), round(Int, SevI0 * N), round(Int, R0 * N), tspan)

        model_output = [solution(t)[2] for t in 0:length(data)-1]
        current_error = SSE_error(model_output, data)

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
function optimise_ps_and_beta(β_values::AbstractVector, c::Float64, γ::Float64, ps_range::AbstractVector, δ::Float64,
    γ_s::Float64, α::Float64, S0::Float64, I0::Float64, SevI0::Float64, R0::Float64, tspan::Tuple{Float64,Float64}, N::Int;
    real_data_time::AbstractVector{Int}, real_data_infected::AbstractVector{Int},
    real_data_severe_time::AbstractVector{Int}, real_data_severe_infected::AbstractVector{Int})

    min_error, best_ps, best_beta = Inf, 0.0, 0.0

    for ps in ps_range
        for β in β_values
            params = SIRSParameters(β, γ, δ, γ_s, α, N)
            solution = run_SIRS(params, round(Int, S0 * N), round(Int, I0 * N), round(Int, SevI0 * N), round(Int, R0 * N), tspan)

            model_infected = [solution(t)[2] for t in real_data_time]
            model_severe_infected = [solution(t)[3] for t in real_data_severe_time]

            combined_error = SSE_error(model_infected, real_data_infected) + SSE_error(model_severe_infected, real_data_severe_infected)

            if combined_error < min_error
                min_error, best_ps, best_beta = combined_error, ps, β
            end
        end
    end

    println("Optimal β: $best_beta, Optimal p_s: $best_ps with error: $min_error")
    return best_ps, best_beta, min_error
end


"""
    run_error(action; kwargs...)

Runs the specified error calculation or optimization function.

# Arguments
- `action::String`: Choose from 'simulate_error', 'optimise_beta', or 'optimise_ps_and_beta'.
- Additional keyword arguments for each function.

# Returns
- The result of the specified function.
"""

real_data_time = [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
real_data_infected = [11, 7, 20, 3, 29, 14, 11, 12, 16, 10, 58, 34, 26, 29, 51, 55]
real_data_severe_time = [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
real_data_severe_infected = [0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 5, 5, 5, 2, 9, 4]

function run_error(action::String;
    β_values=0.2:0.001:0.3, ps_range=0.15:0.01:0.25, δ=0.2,
    c=8.0, γ=1 / 7, γ_s=1 / 14, α=1 / 30, S0=5999/6000, I0=1/6000, SevI0=0.0, R0=0.0,
    tspan=(0.0, 50.0), N=6000,
    real_data_time=real_data_time,
    real_data_infected=real_data_infected,
    real_data_severe_time=real_data_severe_time,
    real_data_severe_infected=real_data_severe_infected)

    # Ensure observed data is provided
    if isempty(real_data_time) || isempty(real_data_infected)
        error("Observed data for `real_data_time` and `real_data_infected` must be provided.")
    end

    if action == "simulate_error"
        println("Running simulate_error...")
        return simulate_error(β_values, c, γ, δ, γ_s, α, S0, I0, SevI0, R0, tspan, N;
            real_data_time=real_data_time,
            real_data_infected=real_data_infected,
            real_data_severe_time=real_data_severe_time,
            real_data_severe_infected=real_data_severe_infected)

    elseif action == "optimise_beta"
        println("Running optimise_beta...")
        best_beta, min_error = optimise_beta(γ, δ, γ_s, α, N, S0, I0, SevI0, R0, real_data_infected, β_values)
        best_beta /= 8
        println("Optimal β: $best_beta with error: $min_error")
        return

    elseif action == "optimise_ps_and_beta"
        println("Running optimise_ps_and_beta...")
        best_ps, best_beta, min_error = optimise_ps_and_beta(
            β_values, c, γ, ps_range, δ, γ_s, α, S0, I0, SevI0, R0, tspan, N;
            real_data_time=real_data_time,
            real_data_infected=real_data_infected,
            real_data_severe_time=real_data_severe_time,
            real_data_severe_infected=real_data_severe_infected)
        best_beta /= 8
        println("Optimal β: $best_beta, Optimal p_s: $best_ps with error: $min_error")
        return

    else
        error("Invalid action: $action. Choose 'simulate_error', 'optimise_beta', or 'optimise_ps_and_beta'.")
    end
end

# run_error("simulate_error")
# run_error("optimise_beta")
# run_error("optimise_ps_and_beta")