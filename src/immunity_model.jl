"""
    reproduction_number(params::Vector)::Float64

Calculates the basic reproduction number R₀ of the disease using a parameter vector.

# Arguments
- `params::Vector`: A vector where:
  - `params[1]` is the contact rate (c),
  - `params[2]` is the probability of transmission per contact (β),
  - `params[3]` is the recovery rate (γ).

# Returns
- The basic reproduction number R₀, representing the expected number of secondary infections produced by one infected individual in a completely susceptible population.
"""
function reproduction_number(params::Vector)::Float64
    return params[1] * params[2] / params[3]
end

"""
    reproduction_number(c::Float64, β::Float64, γ::Float64)::Float64

Calculates the basic reproduction number R₀ of the disease given individual parameters.

# Arguments
- `c::Float64`: Contact rate (average number of contacts per individual per time unit).
- `β::Float64`: Transmission probability per contact.
- `γ::Float64`: Recovery rate (rate at which infected individuals recover).

# Returns
- The basic reproduction number R₀.
"""
function reproduction_number(c::Float64, β::Float64, γ::Float64)::Float64
    return (c * β) / γ
end

"""
    plot_herd_immunity_threshold()

Plots the herd immunity threshold p_c = 1 - (1 / R₀) as a function of R₀.

# Description
Creates a plot showing how the herd immunity threshold (pₙ) varies with different values of the basic reproduction number (R₀). The herd immunity threshold represents the proportion of the population that must be immune to prevent the disease from spreading.

# Note
Uses LaTeX syntax for axis labels and title in the plot for better readability.
"""
function plot_herd_immunity_threshold()
    R0 = 0.1:0.1:30  # Start from 0.1 to avoid division by zero
    pc = 1 .- 1 ./ R0  # Herd immunity threshold formula

    plot(R0, pc,
        xlabel="R₀",
        ylabel="pₙ",
        title="Herd Immunity Threshold vs R₀",
        label="pₙ",
        lw=2)
end

"""
    herd_immunity_threshold(c::Float64, β::Float64, γ::Float64)

Calculates and prints the herd immunity threshold for the given parameters.

# Arguments
- `c::Float64`: Contact rate.
- `β::Float64`: Transmission probability per contact.
- `γ::Float64`: Recovery rate.

# Description
Calculates the herd immunity threshold using the formula p_c = 1 - (1 / R₀), where R₀ is the basic reproduction number, and prints the result.
"""
function herd_immunity_threshold(c::Float64, β::Float64, γ::Float64)
    R0 = reproduction_number(c, β, γ)
    pc = 1 - 1 / R0
    println("Herd immunity threshold (pₙ): ", pc)
end