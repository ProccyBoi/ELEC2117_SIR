# Define the custom error function for Least Squares
function calculate_error(model_output, data)
    # Ensure that the lengths of the two vectors are the same
    if length(model_output) != length(data)
        error("Model output and data must have the same length.")
    end

    # Calculate the sum of squared errors
    return sum((model_output .- data).^2)
end

function optimize_beta(c, γ, N, S0, I0, R0, data, β_values)
    # Store the minimum error and best β value
    min_error = Inf
    best_beta = 0.0

    # Iterate over each β value and calculate the error
    for β in β_values
        solution = run_model(:force_of_infection, c, β, γ, 0.0, S0, I0, R0, (0.0, length(data) - 1), N)
        model_output = solution[2, :]  # Extract the number of infected individuals

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
