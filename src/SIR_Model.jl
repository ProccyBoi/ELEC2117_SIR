module SIR_Model

# Import necessary packages
using Plots
using DifferentialEquations
using LaTeXStrings
using Revise
using Statistics
using Optim

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
    calculate_error_town2, optimise_beta_town2

# Miscellaneous or CLI functionality exports
export define_model, solve_model, plot_model, simulate_model,
    run_model, parse_input, parse_input_int, use_SIR, use_SIRS,
    plot_SIRS

# Include submodules
include("define_model.jl")
include("solve_model.jl")
include("plot_model.jl")
include("simulate_model.jl")
include("immunity_model.jl")
include("error.jl")
include("cli.jl")
include("SIRS.jl")
include("intervention.jl")
include("town.jl")

end
