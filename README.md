### Project Name: Epidemic Modeling in Julia
## Description
This project is a comprehensive Julia-based implementation for simulating, analyzing, and visualizing epidemic models. The repository includes several scripts that define different components of epidemic models such as SIR, immunity, and plotting modules. Each script serves a specific role in creating a modular, extendable, and high-performance simulation environment for studying epidemiological dynamics.

## Repository Structure
The repository contains the following Julia scripts, each serving a unique purpose in the project:

cli.jl: Command-line interface script for executing different simulations and models. It provides an entry point for running the entire program, managing inputs, and setting various parameters for simulations.

define_model.jl: Contains the core functions and structures to define different types of epidemiological models. This includes the SIR model and other extensions for incorporating immunity and other factors.

error.jl: Handles error reporting and debugging utilities to ensure the robustness of the simulations. This script manages custom error types and messages specific to the models.

immunity_model.jl: Extends the base SIR model by incorporating immunity dynamics. It includes functions to simulate the development and loss of immunity within a population.

plot_model.jl: Provides various visualization utilities using plotting libraries. It allows the user to generate graphs and charts for the models, illustrating the dynamics of infection, recovery, and other variables over time.

simulate_model.jl: Includes core simulation routines for running the defined models over specified time steps. This script is the main engine driving the model simulations.

SIR_Model.jl: Implements the classic SIR (Susceptible-Infectious-Recovered) model, detailing the infection, recovery, and transmission dynamics of a basic epidemic model.

solve_model.jl: Contains numerical solvers and optimizers used to find solutions to the differential equations defining the epidemic models.

## Installation
To use this project, ensure that you have Julia installed on your system. Clone this repository, navigate to the project folder, and run the following command to install any necessary dependencies: