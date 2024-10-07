using InteractiveUtils

function parse_input(prompt::String, suggested_value::Float64)
    # Display the prompt and suggested value, and get user input
    print(prompt)
    input = readline()

    # If the user presses Enter without input, return the suggested value
    if isempty(input)
        return suggested_value
    else
        # Try parsing the input to a Float64; fallback to suggested if parsing fails
        return tryparse(Float64, input) !== nothing ? parse(Float64, input) : suggested_value
    end
end

function parse_input_int(prompt::String, suggested_value::Int)
    print(prompt)
    input = readline()
    if isempty(input)
        return suggested_value
    else
        return tryparse(Int, input) !== nothing ? parse(Int, input) : suggested_value
    end
end

function use_SIR()
    println("Welcome to the SIR Model Simulation CLI!")
    println("This CLI allows you to simulate different variations of the SIR model: Basic and with Force of Infection.")
    println()

    # Step 1: Select the SIR Model Type
    println("Please choose the SIR model type you'd like to simulate:")
    println("1. Basic SIR Model")
    println("2. SIR Model with Force of Infection")
    model_choice = readline()

    # Validate user input
    if model_choice != "1" && model_choice != "2"
        println("Invalid choice! Please enter 1 or 2.")
        return
    end

    # Set the model type based on user selection
    model_type = model_choice == "1" ? :basic : :force_of_infection

    # Step 2: Explain the parameters for the SIR model
    println("\nNow let's define the parameters for the SIR model.")
    println("For each parameter, I'll provide a brief explanation and ask for your input.")
    println()

    # Transmission rate (β)
    println("Transmission rate (β): This is the probability of transmission per contact.")
    β = parse_input("Enter the transmission rate (Suggested: 0.03): ", 0.03)

    # Recovery rate (γ)
    println("Recovery rate (γ): The rate at which infected individuals recover and move into the recovered state.")
    γ = parse_input("Enter the recovery rate (Suggested: 0.1): ", 0.1)

    # Contact rate (c) - Only used for the SIR model with force of infection
    c = 0.0
    if model_type == :force_of_infection
        println("Contact rate (c): The number of contacts per day per person.")
        c = parse_input("Enter the contact rate (Suggested: 10.0): ", 10.0)
    end

    # Population size (N)
    println("Population size (N): The total number of individuals in the population.")
    N = parse_input_int("Enter the population size (Suggested: 5000): ", 5000)

    # Initial conditions: Proportions of Susceptible, Infected, and Recovered individuals
    println("\nInitial Conditions:")
    println("Proportion of susceptible individuals (S0):")
    S0 = parse_input("Enter the initial proportion of susceptible individuals (Suggested: 0.9998): ", 0.9998)

    println("Proportion of infected individuals (I0):")
    I0 = parse_input("Enter the initial proportion of infected individuals (Suggested: 0.0002): ", 0.0002)

    println("Proportion of recovered individuals (R0):")
    R0 = parse_input("Enter the initial proportion of recovered individuals (Suggested: 0.0): ", 0.0)

    # Validate that the sum of S0, I0, and R0 equals 1
    if !(S0 + I0 + R0 ≈ 1.0)
        println("Error: The sum of initial proportions (S0, I0, and R0) must equal 1.")
        return
    end

    # Time span for the simulation
    println("\nTime Span for the Simulation:")
    t_start = parse_input("Enter the start time (Suggested: 0.0): ", 0.0)
    t_end = parse_input("Enter the end time (Suggested: 180.0): ", 180.0)
    tspan = (t_start, t_end)

    # Step 3: Run the selected model
    println("\nRunning the selected SIR model...")

    solution = run_model(model_type, c, β, γ, 0.0, S0, I0, R0, tspan, N)

    # Step 4: Plot the results
    println("Simulation complete! Plotting the results...")
    display(plot_model(solution))

    println("\nSimulation finished successfully. Thank you for using the SIR Model CLI!")
end
