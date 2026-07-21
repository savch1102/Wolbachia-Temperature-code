## Libraries ----

using CairoMakie
using DataFrames
using DifferentialEquations
using DiffEqCallbacks

## Temperature parameters ----

function zeta(T)
    1 / (1 + exp(1.5 * (T - 30.5)))
end

function phi(T)
    1 - (1 / (1.5 + exp(0.6 * (-T + 32))))
end

function sigma(T)
    0.035 / (1 + exp(-1.2 * (T - 31)))
end

## Model ----

function wolb(u, p, t)
    # Parameters and conditions

    r_u, r_w, μ_u, μ_w, p_u, p_w, p_h, γ, br, id, K = p
    S, I, R, S_u, I_u, S_w, I_w = u

    if id == "A"
        temp = @. (1 * sin(2π * (t/365) - 2)) + (5 * sin(4π * (t/365) - 0.5)) + 25
        # temp = 25 + 5 * sin((4π * (t - 0) / 365))
    elseif id == "B"
        temp = @. (1 * sin(2π * (t/365) - 0)) + (5 * sin(4π * (t/365) - 0.5)) + 25
        # temp = 25 + 5 * sin((4π * (t - 60) / 365))
    else
        temp = @.  5 * sin(2π * (t/365) - 1) + 25
        # temp = 25 + 5 * sin((4π * (t - 120) / 365))
    end

    N_h = S + I + R
    N_v = S_u + I_u + S_w + I_w

    β_h = br * p_h * (I / N_h)
    β_u = br * p_u * (I_u / N_v)
    β_w = br * p_w * (I_w / N_v)

    freq_w = (S_w + I_w) / (S_u + S_w + I_u + I_w)

    b_u = (r_u * S_u * (1 - phi(temp) * freq_w)) + (r_w * S_w * (1 - zeta(temp)))
    b_w =  r_w * S_w * zeta(temp)

    # ODEs

    dS_u =  (b_u * (1 - N_v / K)) - (S_u * (μ_u + β_h))
    dI_u =  (S_u * β_h) - (I_u * μ_u)

    dS_w =  (b_w * (1 - N_v / K)) - (S_w * (μ_w + sigma(temp) + β_h))
    dI_w =  (S_w * β_h) - (I_w * (μ_w + sigma(temp)))

    dS   = -(β_u + β_w) * S
    dI   =  ((β_u + β_w) * S) - (γ * I)
    dR   =  γ * I

    return [dS, dI, dR, dS_u, dI_u, dS_w, dI_w]
end

## Parameters & Initial conditions ----

params = [
    0.30, # r_u
    0.25, # r_w
    0.05, # μ_u
    0.06, # μ_w
    0.30, # p_u
    0.10, # p_w
    0.40, # p_h
    0.07, # γ
    0.60, # br
    "A",  # id
    3.e4, # K
]

ic = [
    100,          # S
    0,            # I
    0,            # R
    2.5e3 * 1.00, # S_u
    2.5e3 * 0.00, # I_u
    0,            # S_w
    0,            # I_w
]

time = (0, 365)

## Simulation

id = ["A", "B", "C"]

days = [5, 7, 10, 14]
interventions = 500:500:3000

f = Figure(size = (1000, 300), fontsize = 12)
hm = nothing
positions = 1:length(days)
for k in eachindex(id)
    p = copy(params)
    p[end-1] = id[k]
    r = CartesianIndices((1,3))[k][1]
    c = CartesianIndices((1,3))[k][2]
    ax = Axis(f[r,c],
              title = "$(id[k])) Seasonal Regime $(id[k])",
              titlealign = :left,
              xlabel = "Release frequency (days)",
              ylabel = "Intervention (N)",
              xticks = (positions, string.(days)), 
              yticks = (interventions, string.(interventions))
        )

    Z = zeros(length(days), length(interventions))

    for i in eachindex(days)
        for j in eachindex(interventions)
            function add_wolb!(integrator)
                integrator.u[6] += interventions[j]
            end
            cb = PeriodicCallback(add_wolb!, days[i])

            res = DataFrame(solve(ODEProblem(wolb, ic, time, p),
                                  Vern9(), callback = cb, saveat = 0.1))
            rename!(res, [:time, :S, :I, :R, :S_u, :I_u, :S_w, :I_w])
            res.prop_w = (res.S_w .+ res.I_w) ./ (res.S_u .+ res.I_u .+ res.S_w .+ res.I_w)
            idx = findfirst(x -> x .> 0.60, res.prop_w)

            if !isnothing(idx)
                Z[i, j] = res.time[idx] 
            else
                Z[i, j] = NaN 
            end
        end
    end
    hm = heatmap!(ax, positions, interventions, Z, colormap = :thermal, colorrange = (0, 365), interpolate = false)
    xlims!(ax, 0.5, length(days) + 0.5)
end

Colorbar(f[1, 4], hm, label = "Time (days)")

f
