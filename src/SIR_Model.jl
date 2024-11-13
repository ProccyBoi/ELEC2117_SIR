module SIR_Model

# Import necessary packages
using Plots
using DifferentialEquations
using LaTeXStrings
using Revise
using Statistics
using Optim
using InteractiveUtils

# Core SIR model functionality exports
export SIRS_model!, run_SIRS_model, plot_SIRS_model, simulate_SIRS

# Intervention-related functionality exports
export intervention!, run_intervention, plot_intervention, simulate_intervention,
    estimate_infection_severity_ranges, plot_intervention_trajectories, 
    simulate_intervention_with_data_overlay, calculate_total_error, find_optimal_compliance,
    find_optimal_compliance_across_beta_range

# Analysis and metrics functionality exports
export reproduction_number, plot_herd_immunity_threshold, Herd_immunity!,
    herd_immunity_threshold, Force_of_infection!, Basic!

# Error calculation and data fitting functionality exports
export calculate_error, data_fitting, optimise_beta, simulate_error,
    optimize_ps, optimize_ps_and_beta

export test_prevent_outbreak, estimate_entry_day, analyse_town2, error_town2, 
    calculate_error_town2, optimise_town2

# Miscellaneous or CLI functionality exports
export define_model, solve_model, plot_model, simulate_model,
    run_model, parse_input, parse_input_int, use_SIR, use_SIRS,
    plot_SIRS


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
export SSE_error, simulate_error, optimise_beta, optimise_ps_and_beta

# SIRS.jl
export SIRSParameters, SIRS_model!, run_SIRS, plot_SIRS, simulate_SIRS

# cli.jl
export parse_input, parse_input_int, select_model_type, select_data_set, 
    get_town_data, get_parameters, select_action, main

# TLT5.jl
export analyse_variant_similarity, estimate_infection_entry_day, evaluate_intervention_hypothetical,
    calculate_sir_error, calculate_total_error, run_simulation, run_intervention, plot_town

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
include("TLT5.jl")

end
