using SIR_Model
using Test

# Parameters for testing
c = 10              # Contact rate per day per person
β = 0.03            # Transmission rate probability per contact
γ = 0.1             # Recovery rate (10 days) 1/10
N = 5000            # Population size (not directly used in proportions)


# Initial conditions (proportions of the population)
S0 = 0.9998             # 4999 of the population is initially susceptible
I0 = 0.0002             # 1 is initially infected
R0 = 0                  # 0% recovered

# Time span for simulation
tspan = (0.0, 180.0)

@testset "SIR Model Tests" begin
    # Test the basic SIR model
    @testset "Basic SIR Model" begin
        solution = run_model(:basic, c, β, γ, 0.0, S0, I0, R0, tspan, N)  # Herd threshold set to 0.0 as a placeholder

        # Check that the total population S + I + R remains constant at every time step
        for (S, I, R) in solution.u
            @test isapprox(S + I + R, N, atol=1e-6)
        end

        # Check the initial values match the given initial conditions
        @test isapprox(solution.u[1][1], S0 * N, atol=1e-6)  # Check initial susceptible
        @test isapprox(solution.u[1][2], I0 * N, atol=1e-6)  # Check initial infected
        @test isapprox(solution.u[1][3], R0 * N, atol=1e-6)  # Check initial recovered
    end

    # Test the SIR model with force of infection
    @testset "SIR Model with Force of Infection" begin
        solution = run_model(:force_of_infection, c, β, γ, 0.0, S0, I0, R0, tspan, N)  # Herd threshold set to 0.0 as a placeholder

        # Check that the total population S + I + R remains constant at every time step
        for (S, I, R) in solution.u
            @test isapprox(S + I + R, N, atol=1e-6)
        end

        # Calculate the force of infection λ = (c * β * I) / N and check it
        for (S, I, R) in solution.u
            λ = (c * β * I) / N
            @test λ >= 0.0  # Ensure λ is always non-negative
        end

        # Check the initial force of infection is correct using the absolute number of infected individuals
        initial_infected_abs = I0 * N  # Scale I0 to the absolute value
        initial_λ = (c * β * initial_infected_abs) / N
        calculated_λ = (c * β * solution.u[1][2]) / N
        @test isapprox(initial_λ, calculated_λ, atol=1e-6)
    end


    # Additional Edge Case Tests
    @testset "Edge Case Tests" begin
        # No initial infections
        solution_no_infection = run_model(:basic, c, β, γ, 0.0, S0, 0.0, R0, tspan, N)
        @test all(x -> isapprox(x[2], 0.0, atol=1e-6), solution_no_infection.u)

        # All recovered initially
        solution_all_recovered = run_model(:basic, c, β, γ, 0.0, 0.0, 0.0, 1.0, tspan, N)
        @test all(x -> isapprox(x[3], 1.0 * N, atol=1e-6), solution_all_recovered.u)
    end

    # Parameter Robustness Tests
    @testset "Parameter Robustness Tests" begin
        # High transmission rate (β = 1.0)
        solution_high_transmission = run_model(:basic, c, 1.0, γ, 0.0, S0, I0, R0, tspan, N)
        @test all(x -> x[2] >= 0.0, solution_high_transmission.u)

        # Zero recovery rate (γ = 0.0)
        solution_no_recovery = run_model(:basic, c, β, 0.0, 0.0, S0, I0, R0, tspan, N)
        @test all(x -> x[3] == 0.0, solution_no_recovery.u)
    end

    # Timestep Validation Tests
    @testset "Timestep Validation" begin
        # Short simulation time
        short_tspan = (0.0, 10.0)
        solution_short = run_model(:basic, c, β, γ, 0.0, S0, I0, R0, short_tspan, N)
        @test length(solution_short.t) > 0  # Ensure solution is not empty

        # Long simulation time
        long_tspan = (0.0, 500.0)
        solution_long = run_model(:basic, c, β, γ, 0.0, S0, I0, R0, long_tspan, N)
        @test length(solution_long.t) > 0  # Ensure solution is not empty
    end

    # Visual Verification of the Models (Optional)
    @testset "Plot Solutions for Visual Verification" begin
        # Generate and display plots for visual inspection (optional)
        plot_model(run_model(:basic, c, β, γ, 0.0, S0, I0, R0, tspan, N))
        plot_model(run_model(:force_of_infection, c, β, γ, 0.0, S0, I0, R0, tspan, N))
        true
    end
end
