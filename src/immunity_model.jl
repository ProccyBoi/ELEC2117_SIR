function reproduction_number(params::Vector)::Float64
    return params[1] * params[2] / params[3]
end

function reproduction_number(c::Float64, b::Float64, y::Float64)::Float64
    return (c * b) / y
end

function plot_herd_immunity_threshold()
    R0 = 0:0.1:30
    pc = 1 .- 1 ./ R0

    plot(R0, pc, xlabel=L"R_0", ylabel=L"p_c", title=L"Plot of $p_c = 1 - \frac{1}{R_0}$", label=L"p_c", lw=2)
end

function herd_immunity_threshold(c::Float64, b::Float64, y::Float64)
    println(1 - 1 / (reproduction_number(c, b, y)))
end