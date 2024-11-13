# ELEC2117 Infectious Disease Modeling Project Documentation

## Project Overview
The ELEC2117 Infectious Disease Modeling Project uses mathematical models to simulate an infectious disease outbreak in a small town. The project’s primary goal is to aid the Department of Health in understanding outbreak dynamics, estimating infection rates, and assessing intervention efficacy. Implemented in Julia, the project incorporates SIR and SIRS epidemiological models with public health interventions to analyze and control infection spread.

## Project Files

This project is organized into several Julia scripts, each designed for specific modeling and analysis tasks:

- **SIR.jl**: Defines the core SIR model equations and functions to simulate infection spread in a closed population without reinfection or intervention.
  
- **SIRS.jl**: Extends the SIR model by incorporating a resusceptibility factor, allowing recovered individuals to become susceptible again after a set immunity period.
  
- **intervention.jl**: Integrates public health interventions within the SIR/SIRS framework, modifying transmission dynamics to account for interventions that reduce the transmission rate.
  
- **simulate_model.jl**: Runs simulations with variable parameter values, such as transmission probability (β), recovery rate, and intervention efficacy, to explore a range of intervention outcomes.
  
- **plot_model.jl**: Visualizes model predictions alongside observed data, enabling comparative analysis of infection trends and intervention effectiveness.
  
- **solve_model.jl**: A high-level script that manages model setup, equation solving, and parameter fitting for the SIRS model, including optimizations for tuning transmission rates.
  
- **town.jl**: Analyzes infection spread in a hypothetical second town, considering variant characteristics, estimated day of infection entry, and intervention strategies.
  
- **error.jl** and **immunity_model.jl**: Support sensitivity analysis and error tracking by managing error calculations and time-based immunity loss, helping refine model accuracy as parameters are adjusted.

## Setup

### Prerequisites
- **Julia 1.7** or newer
- Julia packages: `DifferentialEquations` and `Plots`