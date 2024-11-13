# Function to optimize β, δ, and infection start day for best fit to both infected and severely infected data
function optimise_town2()
    # Fixed parameters
    γ, γ_s, α = 1 / 7, 1 / 14, 1 / 30
    N = 10000
    ε_i, p_i = 0.3, 0.85
    tspan = (0.0, 100.0)

    # Observed data for Town 2
    observed_infected = [
        21, 29, 25, 30, 28, 34, 28, 54, 57, 92, 73, 80, 109, 102, 128, 135, 163, 150, 211,
        196, 233, 247, 283, 286, 332, 371, 390, 404, 467, 529, 598, 641, 704, 702, 788,
        856, 854, 955, 995, 1065, 1106, 1159, 1217, 1269, 1298, 1328, 1339, 1383, 1431,
        1422, 1414, 1485, 1464, 1480
    ]
    observed_severe = [
        3, 3, 4, 7, 3, 8, 7, 5, 9, 13, 15, 3, 20, 13, 11, 20, 16, 11, 15, 18, 27, 24, 28,
        36, 41, 35, 41, 55, 63, 66, 72, 80, 90, 104, 109, 115, 127, 135, 147, 162, 163,
        186, 194, 200, 216, 223, 241, 249, 258, 275, 277, 299, 302, 300
    ]

    # Refined range for β and δ values
    β_values = 0.2:0.001:0.4
    δ_values = 0.15:0.01:0.25

    # Store the best parameters and minimum error
    min_error = Inf
    best_beta, best_delta, best_start_day = 0.0, 0.0, 1

    # Initial conditions
    S0, I0 = N - 1, 1.0
    SevI0, R0 = 0.0, 0.0

    # Loop over each possible start day (1 to 26) and adjust observed_days
    for start_day in 1:26
        adjusted_observed_days = start_day:(start_day+length(observed_infected)-1)
        adjusted_intervention_day = 36 + start_day - 1

        # Loop over each combination of β and δ
        for β in β_values
            for δ in δ_values
                # Set up parameters for current β, δ, and intervention day
                params = SIRSInterventionParameters(
                    β=β, γ=γ, δ=δ, γ_s=γ_s, α=α, N=N, intervention_day=adjusted_intervention_day, ε_i=ε_i, p_i=p_i
                )

                # Solve with `saveat` aligned to adjusted observed days for direct comparison
                solution = solve(ODEProblem(intervention!, [S0, I0, SevI0, R0], tspan, params), DP8(), saveat=adjusted_observed_days)

                # Extract model predictions on the adjusted observed days
                model_infected = [solution(t)[2] for t in adjusted_observed_days]
                model_severe = [solution(t)[3] for t in adjusted_observed_days]

                # Calculate error with emphasis on severe cases (tune weights if needed)
                infected_error = sum((observed_infected .- model_infected) .^ 2)
                severe_error = sum((observed_severe .- model_severe) .^ 2)
                total_error = infected_error + severe_error

                # Update best parameters if current error is lowest
                if total_error < min_error
                    min_error = total_error
                    best_beta, best_delta, best_start_day = β, δ, start_day
                end
            end
        end
    end

    # Output best-fit β, δ, start day, and associated error
    best_beta = best_beta / 8
    println("Optimal β: $best_beta, Optimal δ: $best_delta, Optimal start day: $best_start_day with error: $min_error")

    # Rerun model with best β, δ, and start day for visualization
    best_beta = best_beta * 8
    params = SIRSInterventionParameters(
        β=best_beta, γ=γ, δ=best_delta, γ_s=γ_s, α=α, N=N, intervention_day=(36 + best_start_day - 1), ε_i=ε_i, p_i=p_i
    )
    best_solution = solve(ODEProblem(intervention!, [S0, I0, SevI0, R0], tspan, params), DP8(), saveat=1.0)

    # Plot model predictions vs observed data
    best_beta = best_beta / 8
    plt = plot(title="Best-Fit Model", xlabel="Days", ylabel="Population")
    plot!(plt, 0:100, [best_solution(t)[2] for t in 0:100], label="Model Infected (β=$best_beta, δ=$best_delta)", lw=2, color=:blue)
    plot!(plt, 0:100, [best_solution(t)[3] for t in 0:100], label="Model Severe (β=$best_beta, δ=$best_delta)", lw=2, color=:red)
    scatter!(plt, best_start_day:(best_start_day+length(observed_infected)-1), observed_infected, label="Observed Infected", color=:blue, marker=:circle)
    scatter!(plt, best_start_day:(best_start_day+length(observed_severe)-1), observed_severe, label="Observed Severe", color=:red, marker=:circle)

    display(plt)
end
