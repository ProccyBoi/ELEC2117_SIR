module SIR_Model

using Plots
using DifferentialEquations
using LaTeXStrings
using Revise

export define_model, solve_model, plot_model, simulate_model, reproduction_number, plot_herd_immunity_threshold, run_model, Basic!, Force_of_infection!, Herd_immunity!, data_fitting, calculate_error, run_SIRS_model, SIRS_model, use_SIR, parse_input, parse_input_int

include("define_model.jl")
include("solve_model.jl")
include("plot_model.jl")
include("simulate_model.jl")
include("immunity_model.jl")
include("error.jl")
include("cli.jl")

end
