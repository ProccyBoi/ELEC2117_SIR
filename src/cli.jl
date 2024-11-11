# Data for Town 1 and Town 2
const DATA_SETS = Dict(
    "Town 1" => Dict(
        "N" => 6000,
        "observed_infected_days" => 15:80,
        "observed_infected" => [
            11, 7, 20, 3, 29, 14, 11, 12, 16, 10, 58, 34, 26, 29, 51, 55, 155, 53, 67, 98, 130,
            189, 92, 192, 145, 128, 68, 74, 126, 265, 154, 207, 299, 273, 190, 152, 276, 408,
            267, 462, 352, 385, 221, 420, 544, 329, 440, 427, 369, 606, 416, 546, 475, 617,
            593, 352, 337, 473, 673, 653, 523, 602, 551, 686, 556, 600
        ],
        "observed_severe_days" => 21:80,
        "observed_severe" => [
            0, 0, 1, 2, 5, 5, 5, 2, 9, 4, 22, 0, 15, 48, 38, 57, 9, 18, 20, 0, 41, 15, 35, 36,
            27, 38, 24, 40, 34, 57, 18, 29, 63, 66, 119, 76, 95, 28, 109, 136, 119, 104, 121,
            93, 147, 129, 130, 161, 133, 136, 138, 139, 181, 181, 218, 183, 167, 164, 219, 220
        ]
    ),
    "Town 2" => Dict(
        "N" => 10000,
        "observed_infected_days" => 27:80,
        "observed_infected" => [
            21, 29, 25, 30, 28, 34, 28, 54, 57, 92, 73, 80, 109, 102, 128, 135, 163, 150, 211,
            196, 233, 247, 283, 286, 332, 371, 390, 404, 467, 529, 598, 641, 704, 702, 788,
            856, 854, 955, 995, 1065, 1106, 1159, 1217, 1269, 1298, 1328, 1339, 1383, 1431,
            1422, 1414, 1485, 1464, 1480
        ],
        "observed_severe_days" => 27:80,
        "observed_severe" => [
            3, 4, 7, 3, 8, 7, 5, 9, 13, 15, 3, 20, 13, 11, 20, 16, 11, 15, 18, 27, 24, 28, 36,
            41, 35, 41, 55, 63, 66, 72, 80, 90, 104, 109, 115, 127, 135, 147, 162, 163, 186,
            194, 200, 216, 223, 241, 249, 258, 275, 277, 299, 302, 300
        ]
    )
)

# Function to select the data set based on user choice
function select_data_set()
    println("\nChoose the data set to model:")
    println("1. Town 1")
    println("2. Town 2")
    data_choice = readline()

    if data_choice == "1"
        return DATA_SETS["Town 1"]
    elseif data_choice == "2"
        return DATA_SETS["Town 2"]
    else
        println("Invalid choice. Defaulting to Town 1 data.")
        return DATA_SETS["Town 1"]
    end
end

# Main function to interact with the user and plot the SIRS model with or without intervention
function main()
    println("Welcome to the SIRS Model Visualizer!")

    # Step 1: Select the data set
    dataset_choice = select_data_set()
    N = dataset_choice["N"]
    observed_infected_days = dataset_choice["observed_infected_days"]
    observed_infected = dataset_choice["observed_infected"]
    observed_severe_days = dataset_choice["observed_severe_days"]
    observed_severe = dataset_choice["observed_severe"]

    # Step 2: Choose whether to include intervention or not
    println("\nSelect model type:")
    println("1. SIRS Model without Intervention")
    println("2. SIRS Model with Intervention")
    model_choice = readline()

    # Define default parameters
    β, γ, δ, γ_s, α = 0.28, 1 / 7, 0.2, 1 / 14, 1 / 30
    S0, I0, SevI0, R0 = N - 1, 1, 0, 0  # Use Int64 for non-intervention case
    tspan = (0.0, 180.0)

    # Run SIRS Model without intervention
    if model_choice == "1"
        params = SIRSParameters(β=β, γ=γ, δ=δ, γ_s=γ_s, α=α, N=N)
        solution = run_SIRS(params, S0, I0, SevI0, R0, tspan)
        plt = plot_SIRS(solution, title="SIRS Model without Intervention", ylabel="Population")

        # Run SIRS Model with intervention
    elseif model_choice == "2"
        intervention_day = 30
        ε_i = 0.3
        p_i = 0.8

        # Define intervention parameters
        params = SIRSInterventionParameters(β=β, γ=γ, δ=δ, γ_s=γ_s, α=α, N=N, intervention_day=intervention_day, ε_i=ε_i, p_i=p_i)
        # Use Float64 for intervention
        solution = run_intervention(params, Float64(S0), Float64(I0), Float64(SevI0), Float64(R0), tspan)
        plt = plot_intervention(solution, title="SIRS Model with Intervention", ylabel="Population")

        # Option to overlay observed data
        println("Would you like to overlay the observed data? (y/n)")
        overlay_data = readline()
        if overlay_data == "y"
            scatter!(plt, observed_infected_days, observed_infected, label="Observed Infected", marker=:circle, color=:red)
            scatter!(plt, observed_severe_days, observed_severe, label="Observed Severe", marker=:diamond, color=:green)
        end
    else
        println("Invalid choice.")
        return
    end

    # Display the plot
    display(plt)
    println("\nThank you for using the SIRS Model Visualizer!")
end
