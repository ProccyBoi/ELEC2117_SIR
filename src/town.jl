# Helper to calculate model error for Town 2 data
function error_town2(solution, observed_infected, observed_severe, start_day)
    model_inf = [solution(t)[2] for t in start_day:(start_day+length(observed_infected)-1)]
    model_sev = [solution(t)[3] for t in start_day:(start_day+length(observed_severe)-1)]
    inf_error = sum((model_inf .- observed_infected) .^ 2)
    sev_error = sum((model_sev .- observed_severe) .^ 2)
    return inf_error + sev_error
end

# Analysis function for Town 2 using Town 1 parameters to assess variant similarity
function analyse_town2()
    # Define model parameters based on best fit from Town 1
    β_town1, p_i_town1 = 0.0365*8, 0.816
    params = SIRSInterventionParameters(
        β=β_town1, γ=1 / 7, δ=0.2, γ_s=1 / 14, α=1 / 30, N=10000, intervention_day=37, ε_i=0.3, p_i=p_i_town1
    )
    S0, I0, SevI0, R0 = 9999.0, 1.0, 0.0, 0.0
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

    # Run simulation and calculate error to assess similarity
    solution = run_intervention(params, S0, I0, SevI0, R0, tspan)
    error = error_town2(solution, observed_infected, observed_severe, 28)

    println("Total error using Town 1 parameters in Town 2: $error")

    return
end

# Estimate the day the infection entered Town 2 based on model fit
function estimate_entry_day()
    β_town1, p_i_town1 = 0.0365 * 8, 0.816
    ε_i = 0.3
    min_error = Inf
    best_day = -1  # Placeholder for the best day
    best_solution = nothing  # Placeholder for the best solution

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

    for day in 1:27  # Test days from 1 to 27 to find the best entry day
        # Parameters for this simulation
        params = SIRSInterventionParameters(
            β=β_town1, γ=1 / 7, δ=0.2, γ_s=1 / 14, α=1 / 30, N=10000, intervention_day=37, ε_i=ε_i, p_i=p_i_town1
        )
        # Start with 1 initial case on the candidate entry day
        S0, I0, SevI0, R0 = 10000.0, 0.0, 0.0, 0.0
        tspan = (0.0, 100.0)

        # Introduce 1 infected person on the entry day
        solution = run_intervention(params, S0 - 1, 1.0, SevI0, R0, tspan)

        # Calculate error from the candidate entry day onwards
        error = error_town2(solution, observed_infected, observed_severe, day)

        # Update best day if the current error is lower than the minimum recorded error
        if error < min_error
            min_error = error
            best_day = day
            best_solution = solution
        end
    end

    println("Estimated entry day in Town 2: Day $best_day with error $min_error")

    # Plot the observed and model data for the best fit
    plt = plot(title="Estimated Entry Day Model vs Observed Data", xlabel="Days", ylabel="Population")
    model_infected = [best_solution(t)[2] for t in 0:100]
    model_severe = [best_solution(t)[3] for t in 0:100]

    # Plot the model predictions for infected and severe cases
    plot!(0:100, model_infected, label="Model Infected (Day $best_day)", lw=2, color=:blue)
    plot!(0:100, model_severe, label="Model Severe (Day $best_day)", lw=2, color=:red)

    # Overlay observed data for infected and severe cases from day 28 onward
    plot!(28:(28+length(observed_infected)-1), observed_infected, seriestype=:scatter, label="Observed Infected", marker=:circle, color=:blue)
    plot!(28:(28+length(observed_severe)-1), observed_severe, seriestype=:scatter, label="Observed Severe", marker=:circle, color=:red)

    display(plt)

    return best_day, min_error
end

# Hypothetical early intervention in Town 2 (e.g., on day 15)
function test_prevent_outbreak()
    β_town1, p_i_town1 = 0.2918, 0.8
    params = SIRSInterventionParameters(
        β=β_town1, 
        γ=1 / 7, 
        δ=0.2, 
        γ_s=1 / 14, 
        α=1 / 30, 
        N=10000, 
        intervention_day=15, 
        ε_i=0.3, 
        p_i=p_i_town1
    )
    S0, I0, SevI0, R0 = 9999.0, 1.0, 0.0, 0.0
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

    # Run model with early intervention
    solution = run_intervention(params, S0, I0, SevI0, R0, tspan)

    # Plot results
    plt = plot(title="Hypothetical Early Intervention in Town 2", xlabel="Days", ylabel="Population")
    model_inf = [solution(t)[2] for t in 0:100]
    model_sev = [solution(t)[3] for t in 0:100]
    plot!(plt, 0:100, model_inf, label="Model Infected (Early Intervention)", lw=2, color=:blue)
    plot!(plt, 0:100, model_sev, label="Model Severe (Early Intervention)", lw=2, color=:red)
    plot!(plt, 28:(28+length(observed_infected)-1), observed_infected, seriestype=:scatter, label="Observed Infected", marker=:circle, color=:blue)
    plot!(plt, 28:(28+length(observed_severe)-1), observed_severe, seriestype=:scatter, label="Observed Severe", marker=:circle, color=:red)
    display(plt)

    println("Hypothetical early intervention plotted for assessment.")
end


# Weighted error function to give more importance to severe cases if needed
function calculate_error_town2(model_infected, model_severe, observed_infected, observed_severe; weight=1.0)
    infected_error = sum((model_infected .- observed_infected) .^ 2)
    severe_error = sum((model_severe .- observed_severe) .^ 2) * weight
    return infected_error + severe_error
end

# Updated optimisation function with grid search for both β, δ, and intervention_day
function optimise_town2()
    # Define fixed parameters
    γ, γ_s, α = 1 / 7, 1 / 14, 1 / 30
    N, S0, I0, SevI0, R0 = 10000, 9999.0, 1.0, 0.0, 0.0
    ε_i, p_i = 0.3, 0.0  # `p_i` fixed at 0.8 for intervention
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
    observed_days = 27:80

    # Define the refined range of β and δ values
    β_values = 0.2:0.01:0.3
    δ_values = 0.05:0.01:0.25

    # Initialize variables to store the best parameters and minimum error
    min_error = Inf
    best_beta, best_delta, best_day = 0.0, 0.0, 0

    # Loop through each possible starting day of infection (1 to 28)
    for start_day in 1:28
        # Loop through each combination of β and δ
        for β in β_values
            for δ in δ_values
                params = SIRSInterventionParameters(
                    β=β, γ=γ, δ=δ, γ_s=γ_s, α=α, N=N, intervention_day=start_day, ε_i=ε_i, p_i=p_i
                )

                # Solve the ODE with `Tsit5()` solver and finer `saveat=0.5` for better precision
                solution = solve(ODEProblem(intervention!, [S0, I0, SevI0, R0], tspan, params), Tsit5(), saveat=0.1)

                # Calculate the model predictions only for observed data points
                model_infected = [solution(t)[2] for t in observed_days]
                model_severe = [solution(t)[3] for t in observed_days]

                # Calculate Mean Absolute Error (MAE) between model predictions and observed data
                infected_error = sum(abs, observed_infected .- model_infected)
                severe_error = sum(abs, observed_severe .- model_severe)
                total_error = infected_error + severe_error

                # Check if the current error is the minimum
                if total_error < min_error
                    min_error = total_error
                    best_beta, best_delta, best_day = β, δ, start_day
                end
            end
        end
    end

    # Output the best β, δ, and intervention day values and the associated minimum error
    println("Optimal start day: $best_day, Optimal β: $best_beta, Optimal δ: $best_delta with error: $min_error")

    # Plot observed vs model data for the best-fit parameters with a higher precision solver
    params = SIRSInterventionParameters(
        β=best_beta, γ=γ, δ=best_delta, γ_s=γ_s, α=α, N=N, intervention_day=best_day, ε_i=ε_i, p_i=p_i
    )
    best_solution = solve(ODEProblem(intervention!, [S0, I0, SevI0, R0], tspan, params), DP8(), saveat=0.1)

    # Plot model predictions vs observed data
    plt = plot(title="Best-Fit Model vs Observed Data for Optimized β, δ, and Start Day", xlabel="Days", ylabel="Population")
    plot!(plt, 0:0.1:100, [best_solution(t)[2] for t in 0:0.1:100], label="Model Infected (β=$best_beta, δ=$best_delta)", lw=2, color=:blue)
    plot!(plt, 0:0.1:100, [best_solution(t)[3] for t in 0:0.1:100], label="Model Severe (β=$best_beta, δ=$best_delta)", lw=2, color=:red)
    scatter!(plt, observed_days, observed_infected, label="Observed Infected", color=:blue, marker=:circle)
    scatter!(plt, observed_days, observed_severe, label="Observed Severe", color=:red, marker=:circle)

    display(plt)
end