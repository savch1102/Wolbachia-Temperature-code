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

    r_u, r_w, μ_u, μ_w, p_u, p_w, p_h, γ, br, temp, K = p
    S, I, R, S_u, I_u, S_w, I_w = u

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
    30.0, # temp
    3.e4, # K
]

ic = [
    300,          # S
    0,            # I
    0,            # R
    2.5e3,        # S_u
    0,# 2.5e3 * 0.01, # I_u
    0,            # S_w
    0,            # I_w
]

time = (0, 365)

## Simulation

temps = [25, 29, 28, 30]
id = ["A)", "C)", "B)", "D)"]

f = Figure(size = (900, 600), fontsize = 12)

l1 = nothing
l2 = nothing

for k in eachindex(temps)
    p = copy(params)
    p[end-1] = temps[k]

    function add_wolb!(integrator)
        integrator.u[6] += 1500
    end

    cb_wolb = PeriodicCallback(add_wolb!, 7)

    function add_inf!(integrator)
        integrator.u[2] += 3
    end

    cb_inf = PresetTimeCallback(90, add_inf!)

    call = CallbackSet(cb_wolb, cb_inf)

    res = DataFrame(solve(ODEProblem(wolb, ic, time, p),
                            Rodas5P(), callback = call, saveat = 0.1))
    rename!(res, [:time, :S, :I, :R, :S_u, :I_u, :S_w, :I_w])
    res.prop_w = (res.S_w .+ res.I_w) ./ (res.S_u .+ res.I_u .+ res.S_w .+ res.I_w)

    r = CartesianIndices((2,2))[k][1]
    c = CartesianIndices((2,2))[k][2]

    ax1 = Axis(f[r,c],
               title = "$(id[k]) Temperature: $(temps[k]) °C",
               titlealign = :left,
               xlabel = "Time (days)",
               ylabel = "Proportion",
               limits = (nothing, (nothing, 1.1)),
               yticks = 0:0.2:1.0
        )
    l1 = lines!(ax1, res.time, res.prop_w, color = :indianred, linewidth = 2)

    ax2 = Axis(f[r,c], 
               yaxisposition = :right,
               ylabel = "Dengue Cases (N)",
               ygridvisible = false,
               limits = (nothing, (nothing, 110)),
               yticks = 0:20:100)
    hidespines!(ax2)
    hidexdecorations!(ax2)
    l2 = lines!(ax2, res.time, res.I, color = :midnightblue, linestyle = :dash, linewidth = 2)
    println(maximum(res.I))
end
colgap!(f.layout, 50)
Legend(f[3, 1:2], [l1, l2], ["Proportion of Wolbachia", "Dengue cases"], "Population", orientation = :horizontal)
f
