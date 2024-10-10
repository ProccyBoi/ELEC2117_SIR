# Define the custom error function for Least Squares
function calculate_error(model_output, data)
    # Ensure that the lengths of the two vectors are the same
    if length(model_output) != length(data)
        error("Model output and data must have the same length.")
    end

    # Calculate the sum of squared errors
    return sum((model_output .- data).^2)
end
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

# Updated optimize_beta function
function optimize_beta(γ, δ, γ_s, α, N, S0, I0, SevI0, R0, data, β_values)
    # Store the minimum error and best β value
    min_error = Inf
    best_beta = 0.0

    # Time span for the model based on the data length (match exactly)
    tspan = (0.0, length(data) - 1)

    # Iterate over each β value and calculate the error
    for β in β_values
        # Run the SIRS model for the current β value using the given time span
        solution = run_SIRS_model(β, γ, δ, γ_s, α, N, S0 * N, I0 * N, SevI0 * N, R0 * N, tspan)

        # Extract the number of infected individuals (I compartment) at the same time points as the data
        model_output = solution[2, :]  # This should match the length of `data`

        # Calculate the error using the defined function
        current_error = calculate_error(model_output, data)

        # Update the minimum error and best β if a lower error is found
        if current_error < min_error
            min_error = current_error
            best_beta = β
        end
    end

    return best_beta, min_error
end

# optimize_beta(1/7, 0.2 * (1/7), 1/14, 1/30, 6000, 5999/6000, 1/6000, 0, 0, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 7, 20, 3, 29, 14, 11, 12, 16, 10, 58], 0.01:0.01:0.1)
# optimize_beta(1/7, 0.2 * (1/7), 1/14, 1/30, 6000, 5999/6000, 1/6000, 0, 0, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 7, 20, 3, 29, 14, 11, 12, 16, 10, 58], 0.1:0.01:0.5)