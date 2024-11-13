"""
# BASIC!(1) - Julia Function Documentation

## NAME
**Basic!** - Defines the basic SIR (Susceptible-Infected-Recovered) model equations for simulating infectious disease dynamics.

## SYNOPSIS
```
Basic!(du, u, params, t)
```

## DESCRIPTION
The **Basic!** function represents the fundamental Susceptible-Infected-Recovered (SIR) model for simulating the spread of infectious diseases. It defines the dynamics of three compartments—Susceptible (`S`), Infected (`I`), and Recovered (`R`)—using a set of ordinary differential equations (ODEs). The function calculates the rate of change for each compartment at a given time `t` based on the current values of `S`, `I`, and `R` as well as model parameters.

This function is typically used as part of an `ODEProblem` definition in the `DifferentialEquations.jl` package to simulate the progression of an epidemic in a closed population.

## PARAMETERS
**du**  
: A vector that stores the calculated derivatives (`dS/dt`, `dI/dt`, and `dR/dt`). This is modified in place.

**u**  
: A vector containing the current values of the compartments `S`, `I`, and `R` at time `t`. Each element of `u` represents:
  - `u[1]`: Susceptible individuals (`S`).
  - `u[2]`: Infected individuals (`I`).
  - `u[3]`: Recovered individuals (`R`).

**params**  
: A tuple or vector containing the model parameters, defined as:
  - **`β`**: Transmission rate probability per contact.
  - **`γ`**: Recovery rate, representing the rate at which infected individuals recover.
  - **`N`**: Total population size.

**t**  
: The current time point in the simulation. This is typically passed automatically by the `DifferentialEquations.jl` solver and does not need to be manually set.

## RETURN VALUE
The function modifies the `du` vector in place to update the derivatives for each compartment:

- **`du[1]`**: Change in the number of susceptible individuals (`dS/dt`).
- **`du[2]`**: Change in the number of infected individuals (`dI/dt`).
- **`du[3]`**: Change in the number of recovered individuals (`dR/dt`).

The derivatives are calculated using the standard SIR model equations:

- **Susceptible**:  
    dS/dt = -λS = βIS/N

- **Infected**:  
    dI/dt = λS - γI = cIS/N - γI

- **Recovered**:  
    dR/dt = γI

## DIAGNOSTICS
The **Basic!** function should not raise errors if provided with valid inputs. However, potential issues may arise if:

- The `u` vector does not contain exactly three elements (`S`, `I`, `R`).
- The `params` vector does not contain the expected parameters (`β`, `γ`, `N`).
- The initial conditions or parameter values are not positive, which may result in unrealistic or unstable model outputs.

## NOTES
- This function is designed to be used with `DifferentialEquations.jl` and should be defined before creating an `ODEProblem`.
- Ensure that the initial conditions for `S`, `I`, and `R` add up to the total population `N` for a closed population model.

## SEE ALSO
- `Force_of_infection!` for a variation of the SIR model that includes force of infection.
- `Herd_immunity!` for a model with a herd immunity threshold.
- `solve` from the `DifferentialEquations.jl` package for solving ODE problems.
"""
function Basic!(du, u, params, t)
    S, I, R = u         # Current values of S, I, R
    β, γ, N = params    # Parameters: transmission rate, recovery rate, and population size

    # Basic SIR model equations
    du[1] = -β * S * I / N         # dS/dt
    du[2] = β * S * I / N - γ * I  # dI/dt
    du[3] = γ * I                  # dR/dt
end

"""
# FORCE_OF_INFECTION!(1) - Julia Function Documentation

## NAME
**Force_of_infection!** - Defines the SIR (Susceptible-Infected-Recovered) model with an explicit force of infection term.

## SYNOPSIS
```
Force_of_infection!(du, u, params, t)
```

## DESCRIPTION
The **Force_of_infection!** function represents a variant of the standard Susceptible-Infected-Recovered (SIR) model that explicitly incorporates the force of infection (λ). This force of infection is calculated as a product of the contact rate, transmission rate, and the proportion of infected individuals. It defines the dynamics of three compartments—Susceptible (`S`), Infected (`I`), and Recovered (`R`)—using a set of ordinary differential equations (ODEs). The function computes the rate of change for each compartment at a given time `t` based on the current values of `S`, `I`, and `R` as well as model parameters.

This function is typically used as part of an `ODEProblem` definition in the `DifferentialEquations.jl` package to simulate the progression of an epidemic when a constant contact rate is explicitly considered.

## PARAMETERS
**du**  
: A vector that stores the calculated derivatives (`dS/dt`, `dI/dt`, and `dR/dt`). This is modified in place.

**u**  
: A vector containing the current values of the compartments `S`, `I`, and `R` at time `t`. Each element of `u` represents:
  - `u[1]`: Susceptible individuals (`S`).
  - `u[2]`: Infected individuals (`I`).
  - `u[3]`: Recovered individuals (`R`).

**params**  
: A tuple or vector containing the model parameters, defined as:
  - **`c`**: Contact rate per day per individual. Represents the average number of contacts a susceptible individual has per day.
  - **`β`**: Transmission rate probability per contact.
  - **`γ`**: Recovery rate, representing the rate at which infected individuals recover.
  - **`N`**: Total population size.

**t**  
: The current time point in the simulation. This is typically passed automatically by the `DifferentialEquations.jl` solver and does not need to be manually set.

## RETURN VALUE
The function modifies the `du` vector in place to update the derivatives for each compartment:

- **`du[1]`**: Change in the number of susceptible individuals (`dS/dt`).
- **`du[2]`**: Change in the number of infected individuals (`dI/dt`).
- **`du[3]`**: Change in the number of recovered individuals (`dR/dt`).

The derivatives are calculated using the SIR model equations with an explicit force of infection term:

- **Force of Infection**:  
    λ = cβI/N

- If **herd immunity is reached** (i.e., `I + R >= herd_threshold * N`):
    λ = 0.0

- **Susceptible**:  
    dS/dt = -λS = -cβIS/N

- **Infected**:  
    dI/dt = λS - γI = cβIS/N - γI

- **Recovered**:
    dR/dt = γI

## DIAGNOSTICS
The **Force_of_infection!** function should not raise errors if provided with valid inputs. However, potential issues may arise if:

- The `u` vector does not contain exactly three elements (`S`, `I`, `R`).
- The `params` vector does not contain the expected parameters (`c`, `β`, `γ`, `N`).
- The initial conditions or parameter values are not positive, which may result in unrealistic or unstable model outputs.

## NOTES
- This function is designed to be used with `DifferentialEquations.jl` and should be defined before creating an `ODEProblem`.
- Ensure that the initial conditions for `S`, `I`, and `R` add up to the total population `N` for a closed population model.

## SEE ALSO
- `Basic!` for the standard SIR model without force of infection.
- `Herd_immunity!` for a model with a herd immunity threshold.
- `solve` from the `DifferentialEquations.jl` package for solving ODE problems.
"""
function Force_of_infection!(du, u, params, t)
    S, I, R = u
    c, β, γ, N = params     # Include `N` in the parameter list

    # Calculate force of infection: λ = (c * β * I) / N
    λ = c * β * I / N

    # SIR equations with force of infection
    du[1] = -λ * S                 # dS/dt
    du[2] = λ * S - γ * I          # dI/dt
    du[3] = γ * I                  # dR/dt
end

"""
# HERD_IMMUNITY!(1) - Julia Function Documentation

## NAME
**Herd_immunity!** - Defines the SIR (Susceptible-Infected-Recovered) model with a herd immunity threshold.

## SYNOPSIS
```
Herd_immunity!(du, u, params, t)
```

## DESCRIPTION
The **Herd_immunity!** function represents a modified version of the SIR (Susceptible-Infected-Recovered) model that includes a herd immunity threshold. The function calculates the rate of change for each compartment at a given time `t` based on the current values of `S`, `I`, and `R` and uses a threshold to determine when herd immunity is achieved.

Once the sum of the Infected (`I`) and Recovered (`R`) individuals reaches a specified proportion of the population (`herd_threshold`), the force of infection is set to zero, effectively halting further spread of the disease. This simulates the effect of herd immunity in a population, where a sufficiently large number of individuals are immune, preventing the disease from spreading.

This function is typically used as part of an `ODEProblem` definition in the `DifferentialEquations.jl` package to simulate the progression of an epidemic in populations where herd immunity can be achieved.

## PARAMETERS
**du**  
: A vector that stores the calculated derivatives (`dS/dt`, `dI/dt`, and `dR/dt`). This is modified in place.

**u**  
: A vector containing the current values of the compartments `S`, `I`, and `R` at time `t`. Each element of `u` represents:
  - `u[1]`: Susceptible individuals (`S`).
  - `u[2]`: Infected individuals (`I`).
  - `u[3]`: Recovered individuals (`R`).

**params**  
: A tuple or vector containing the model parameters, defined as:
  - **`c`**: Contact rate per day per individual. Represents the average number of contacts a susceptible individual has per day.
  - **`β`**: Transmission rate probability per contact.
  - **`γ`**: Recovery rate, representing the rate at which infected individuals recover.
  - **`herd_threshold`**: Proportion of the population that must be immune (either recovered or currently infected) to achieve herd immunity (e.g., `0.6` for 60%).
  - **`N`**: Total population size.

**t**  
: The current time point in the simulation. This is typically passed automatically by the `DifferentialEquations.jl` solver and does not need to be manually set.

## RETURN VALUE
The function modifies the `du` vector in place to update the derivatives for each compartment:

- **`du[1]`**: Change in the number of susceptible individuals (`dS/dt`).
- **`du[2]`**: Change in the number of infected individuals (`dI/dt`).
- **`du[3]`**: Change in the number of recovered individuals (`dR/dt`).

The derivatives are calculated using the SIR model equations with a herd immunity adjustment:

- **Force of Infection**:  
    λ = cβI/N

- If **herd immunity is reached** (i.e., `I + R >= herd_threshold * N`):
    λ = 0.0

- **Susceptible**:  
    dS/dt = -λS = -cβIS/N

- **Infected**:  
    dI/dt = λS - γI = cβIS/N - γI

- **Recovered**:
    dR/dt = γI

## DIAGNOSTICS
The **Herd_immunity!** function should not raise errors if provided with valid inputs. However, potential issues may arise if:

- The `u` vector does not contain exactly three elements (`S`, `I`, `R`).
- The `params` vector does not contain the expected parameters (`c`, `β`, `γ`, `herd_threshold`, `N`).
- The initial conditions or parameter values are not positive, which may result in unrealistic or unstable model outputs.
- The `herd_threshold` parameter is not within the range of 0 to 1.

## NOTES
- This function is designed to be used with `DifferentialEquations.jl` and should be defined before creating an `ODEProblem`.
- Ensure that the initial conditions for `S`, `I`, and `R` add up to the total population `N` for a closed population model.
- Choose the `herd_threshold` parameter carefully based on the specific epidemiological characteristics of the disease being modeled.

## SEE ALSO
- `Basic!` for the standard SIR model without force of infection or herd immunity.
- `Force_of_infection!` for a model with a constant force of infection term.
- `solve` from the `DifferentialEquations.jl` package for solving ODE problems.
"""
function Herd_immunity!(du, u, params, t)
    S, I, R = u
    c, β, γ, herd_threshold, N = params

    # Calculate force of infection: λ = (c * β * I) / N
    λ = c * β * I / N

    # Adjust λ based on herd immunity threshold
    if I + R >= herd_threshold * N      # Check if R reaches herd immunity proportion (e.g., 0.6 * N)
        λ = 0.0                         # Set force of infection to zero once herd immunity is reached
    end

    # SIR equations with herd immunity threshold
    du[1] = -λ * S                  # dS/dt
    du[2] = λ * S - γ * I           # dI/dt
    du[3] = γ * I                   # dR/dt
end

"""
# run_SIR(1) - Julia Function Documentation

## NAME
**run_SIR** - Run different types of SIR (Susceptible-Infected-Recovered) models based on the specified parameters.

## SYNOPSIS
```
run_SIR(model_type::Symbol, c, β, γ, herd_threshold, S0, I0, R0, tspan, N)
```

## DESCRIPTION
The **run_SIR** function simulates various versions of the SIR model based on the specified `model_type`. It dynamically selects between different models, such as the basic SIR model, the SIR model with a force of infection, or the SIR model with herd immunity, and solves the chosen model using the `DifferentialEquations.jl` package.

The function initializes the model with appropriate parameters and initial conditions, constructs an `ODEProblem`, and solves it over the specified time span (`tspan`). This flexibility allows users to experiment with different epidemiological scenarios without modifying the core simulation logic.

## PARAMETERS
**model_type**  
: A `Symbol` that determines which type of SIR model to run. Possible values include:
  - **`:basic`**: Runs the basic SIR model.
  - **`:force_of_infection`**: Runs the SIR model with an explicit force of infection term.
  - **`:herd_immunity`**: Runs the SIR model with a herd immunity threshold.

**c**  
: Contact rate per day per individual. Used in models that consider force of infection or herd immunity. Represents the average number of contacts a susceptible individual has per day.

**β**  
: Transmission rate probability per contact. Determines the likelihood of disease transmission during each contact.

**γ**  
: Recovery rate. Represents the rate at which infected individuals recover and move into the recovered compartment.

**herd_threshold**  
: Proportion of the population required to reach herd immunity. Only used if `model_type == :herd_immunity`. For example, a value of `0.6` means 60% of the population needs to be recovered or immune for herd immunity to take effect.

**S0**  
: Initial proportion of the population that is susceptible. This value is scaled to the total population size `N`.

**I0**  
: Initial proportion of the population that is infected. This value is scaled to the total population size `N`.

**R0**  
: Initial proportion of the population that is recovered. This value is scaled to the total population size `N`.

**tspan**  
: A tuple defining the start and end time of the simulation (e.g., `(0.0, 50.0)`). The simulation will run over this time span.

**N**  
: Total population size. This is used to scale the initial proportions (`S0`, `I0`, `R0`) to absolute values in the ODE system.

## RETURN VALUE
Returns an `ODESolution` object containing the solution to the chosen SIR model. The `ODESolution` object includes:

- **`t`**: A vector containing the time points at which the solution is evaluated.
- **`u`**: A matrix of solution values corresponding to each time point. Each row represents a different compartment (`S`, `I`, and `R`) in the model.

## DIAGNOSTICS
**run_SIR** returns an `ODESolution` object if the problem is solved successfully. However, potential issues may arise if:

- **Invalid model type**: If `model_type` is not one of `:basic`, `:force_of_infection`, or `:herd_immunity`, the function raises an error: `Unknown model type: model_type`.
- **Parameter Mismatch**: If the number of elements in `params` does not match the required number for a given model, the function may raise an error during the construction of the `ODEProblem`.
- **Negative Initial Conditions**: If `S0`, `I0`, or `R0` are negative, the function will produce unrealistic results and may trigger numerical errors.

## NOTES
- Ensure that the `DifferentialEquations.jl` package is correctly installed and imported before using this function.
- The initial conditions (`S0`, `I0`, and `R0`) should add up to 1 (i.e., `S0 + I0 + R0 = 1`) to ensure a consistent population proportion.
- Choose a time span (`tspan`) that captures the dynamics of the epidemic to avoid missing important behaviors.

## SEE ALSO
- `Basic!` for the standard SIR model.
- `Force_of_infection!` for the SIR model with force of infection.
- `Herd_immunity!` for the SIR model with a herd immunity threshold.
- `solve` from the `DifferentialEquations.jl` package for solving ODE problems.
"""
function run_SIR(model_type::Symbol, c, β, γ, herd_threshold, S0, I0, R0, tspan, N)
    # Initial conditions
    u0 = [S0 * N, I0 * N, R0 * N]  # Scale initial proportions to match total population size `N`

    # Define the problem based on the model type
    if model_type == :basic
        # Basic SIR model (ignores `c` and `herd_threshold`)
        prob = ODEProblem(Basic!, u0, tspan, [β, γ, N])
    elseif model_type == :force_of_infection
        # SIR model with force of infection (uses `c`, but no herd immunity threshold)
        prob = ODEProblem(Force_of_infection!, u0, tspan, [c, β, γ, N])
    elseif model_type == :herd_immunity
        # SIR model with herd immunity (uses `c`, `herd_threshold`, and `N`)
        prob = ODEProblem(Herd_immunity!, u0, tspan, [c, β, γ, herd_threshold, N])
    else
        error("Unknown model type: $model_type")
    end

    # Solve the problem and return the solution
    return solve(prob)
end