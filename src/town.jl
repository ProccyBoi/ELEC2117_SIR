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
    β_town1, p_i_town1 = 0.036, 0.77
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

    return solution, error  # Returns solution to visually assess variant similarity
end

# Estimate day infection entered Town 2 based on model fit
function estimate_entry_day()
    β_town1, p_i_town1 = 0.28, 0.77
    ε_i = 0.3
    errors = []

    for day in 10:27
        params = SIRSInterventionParameters(
            β=β_town1, γ=1 / 7, δ=0.2, γ_s=1 / 14, α=1 / 30, N=10000, intervention_day=37, ε_i=ε_i, p_i=p_i_town1
        )
        S0, I0, SevI0, R0 = 9999.0, 1.0, 0.0, 0.0
        tspan = (0.0, 100.0)

        solution = run_intervention(params, S0, I0, SevI0, R0, tspan)
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

        error = error_town2(solution, observed_infected, observed_severe, day)
        push!(errors, (day, error))
    end

    best_day, best_error = findmin(errors)
    println("Estimated entry day in Town 2: Day $best_day with error $best_error")
    return
end

# Hypothetical early intervention in Town 2 (e.g., on day 15)
function test_prevent_outbreak()
    β_town1, p_i_town1 = 0.31, 0.77
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

# Function to calculate combined error for infected and severe cases
function calculate_error_town2(model_infected, model_severe, observed_infected, observed_severe)
    infected_error = sum((model_infected .- observed_infected) .^ 2)
    severe_error = sum((model_severe .- observed_severe) .^ 2)
    return infected_error + severe_error
end

# Function to optimize β for the best fit to observed Town 2 data
function optimise_beta_town2()
    # Define parameters that do not vary
    γ, δ, γ_s, α = 1 / 7, 0.2, 1 / 14, 1 / 30
    N, S0, I0, SevI0, R0 = 10000, 9999.0, 1.0, 0.0, 0.0
    p_i, ε_i = 0.77, 0.3
    tspan = (0.0, 100.0)
    intervention_day = 37

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

    # Define a range of β values to test
    β_values = 0.001:0.0001:0.5
    min_error = Inf
    best_beta = 0.0

    # Iterate over β values to find the one that minimizes the error
    for β in β_values
        params = SIRSInterventionParameters(
            β=β, γ=γ, δ=δ, γ_s=γ_s, α=α, N=N, intervention_day=intervention_day, ε_i=ε_i, p_i=p_i
        )
        solution = run_intervention(params, S0, I0, SevI0, R0, tspan)

        # Extract model predictions for infected and severe compartments
        model_infected = [solution(t)[2] for t in 28:(28+length(observed_infected)-1)]
        model_severe = [solution(t)[3] for t in 28:(28+length(observed_severe)-1)]

        # Calculate error between model predictions and observed data
        current_error = calculate_error_town2(model_infected, model_severe, observed_infected, observed_severe)

        # Update minimum error and best β if the current error is lower
        if current_error < min_error
            min_error = current_error
            best_beta = β
        end
    end

    println("Optimal β for best fit: $best_beta with error: $min_error")
    return
end