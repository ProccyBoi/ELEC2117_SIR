module SIR_Model

# Import necessary packages
using Plots
using DifferentialEquations
using LaTeXStrings
using Revise
using Statistics
using Optim
using InteractiveUtils

# SIR.jl
export Basic!, Force_of_infection!, Herd_immunity!, run_SIR

# solve_model.jl
export solve_model

# plot_model.jl
export plot_model

# simulate_model.jl
export simulate_model

# immunity_model.jl
export reproduction_number, plot_herd_immunity_threshold, herd_immunity_threshold

# error.jl
export SSE_error, simulate_error, optimise_beta, optimise_ps_and_beta, run_error

# SIRS.jl
export SIRSParameters, SIRS_model!, run_SIRS, plot_SIRS, simulate_SIRS

# cli.jl
export parse_input, parse_input_int, select_model_type, select_data_set, 
    get_town_data, get_parameters, select_action, main

# intervention.jl
export SIRSInterventionParameters, intervention!, run_intervention, 
    plot_intervention, simulate_intervention, 
    estimate_infection_severity_ranges, plot_intervention_trajectories, 
    simulate_intervention_with_data_overlay, calculate_total_error, 
    find_optimal_compliance, find_optimal_compliance_across_beta_range

# town.jl
export optimise_town2

# Include submodules
include("SIR.jl")
include("solve_model.jl")
include("plot_model.jl")
include("simulate_model.jl")
include("immunity_model.jl")
include("error.jl")
include("SIRS.jl")
include("cli.jl")
include("intervention.jl")
include("town.jl")

end
