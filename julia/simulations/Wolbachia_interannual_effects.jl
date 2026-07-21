## Libraries ----

using CairoMakie
using DataFrames
using DifferentialEquations
using DiffEqCallbacks
using Distributions
using Statistics
using Random

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

    r_u, r_w, μ_u, μ_w, p_u, p_w, p_h, γ, br, K, λ, σ, ω, id = p
    S, I, R, S_u, I_u, S_w, I_w = u

    if id == "A"
        temp = @. (1 * sin(2π * (t/365) - 2)) + 
                (5 * sin(4π * (t/365) - 0.5)) + (25 + (λ * t))       # Change 1: Base line change
    elseif id == "B"
        temp = @. (1 + (σ * t)) * ((1 * sin(2π * (t/365) - 2)) + 
                            (5 * sin(4π * (t/365) - 0.5))) + 25      # Change 2: Intensity
    else
        temp = @. (1 * sin(2π * (t/365) - (2 + (ω * t)))) + 
           (5 * sin(4π * (t/365) - (0.5 + (ω * t)))) + 25            # Change 3: Phase displacement
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
    0.30,       # r_u
    0.25,       # r_w
    0.05,       # μ_u
    0.06,       # μ_w
    0.30,       # p_u
    0.10,       # p_w
    0.40,       # p_h
    0.07,       # γ
    0.60,       # br
    3.e4,       # K
    1.00,       # λ
    1.00,       # σ
    1.00,       # ω
    "A"         # id
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

time = (0, 365 * 5)

## Simulation ----

Random.seed!(1234)

id = ["A", "B", "C"]
λ_vals = 0:8.0e-5:0.008
σ_vals = 0:2.7e-6:0.00027
ω_vals = 0:1.5e-5:0.0015

df = DataFrame(
    time = Float64[], 
    S = Float64[], 
    I = Float64[], 
    R = Float64[],
    S_u = Float64[],
    I_u = Float64[],
    S_w = Float64[],
    I_w = Float64[],
    prop_w = Float64[],
    func = String[]
)

for k in eachindex(id)
    if id[k] == "A"
        for i in eachindex(λ_vals)
            p = copy(params)
            p[end] = id[k]
            p[end-3] = λ_vals[i]

            function add_wolb!(integrator)
                integrator.u[6] += 1500
            end
            cb = PeriodicCallback(add_wolb!, 7)

            res = DataFrame(solve(ODEProblem(wolb, ic, time, p), Vern9(), callback = cb, saveat = 0.1))
            rename!(res, [:time, :S, :I, :R, :S_u, :I_u, :S_w, :I_w])
            res.prop_w = (res.S_w .+ res.I_w) ./ (res.S_u .+ res.I_u .+ res.S_w .+ res.I_w)
            res.func = fill(id[k], nrow(res))

            df = vcat(df, res)
        end
    elseif id[k] == "B"
        for i in eachindex(σ_vals)
            p = copy(params)
            p[end] = id[k]
            p[end-2] = σ_vals[i]

            function add_wolb!(integrator)
                integrator.u[6] += 1500
            end
            cb = PeriodicCallback(add_wolb!, 7)

            res = DataFrame(solve(ODEProblem(wolb, ic, time, p), Vern9(), callback = cb, saveat = 0.1))
            rename!(res, [:time, :S, :I, :R, :S_u, :I_u, :S_w, :I_w])
            res.prop_w = (res.S_w .+ res.I_w) ./ (res.S_u .+ res.I_u .+ res.S_w .+ res.I_w)
            res.func = fill(id[k], nrow(res))
            
            df = vcat(df, res)
        end
    else
        for i in eachindex(ω_vals)
            p = copy(params)
            p[end] = id[k]
            p[end-1] = ω_vals[i]

            function add_wolb!(integrator)
                integrator.u[6] += 1500
            end
            cb = PeriodicCallback(add_wolb!, 7)

            res = DataFrame(solve(ODEProblem(wolb, ic, time, p), Vern9(), callback = cb, saveat = 0.1))
            rename!(res, [:time, :S, :I, :R, :S_u, :I_u, :S_w, :I_w])
            res.prop_w = (res.S_w .+ res.I_w) ./ (res.S_u .+ res.I_u .+ res.S_w .+ res.I_w)
            res.func = fill(id[k], nrow(res))
            
            df = vcat(df, res)
        end
    end
end

df_stats = combine(groupby(df, [:func, :time]),
                   :prop_w => mean => :mean,
                   :prop_w => (x -> quantile(x, 0.025)) => :lower,
                   :prop_w => (x -> quantile(x, 0.975)) => :upper)
df_stats.time = df_stats.time ./ 365

function recode(x::String)
    if x == "A"
        return "Base line change"
    elseif x == "B"
        return "Intensity change"
    else
        return "Phase displacement change"
    end
end

df_stats.newname = recode.(df_stats.func)

## Plot ----

colors0 = ["indianred", "seagreen", "steelblue"]

f = Figure(size = (1000, 300), fontsize = 12)
funcs = ["A", "B", "C"]
labs = ["A)", "B)", "C)"]

for i in eachindex(funcs)
    r = CartesianIndices((1,3))[i][1]
    c = CartesianIndices((1,3))[i][2]

    df_plot = subset(df_stats, :func => ByRow(x -> x == funcs[i]))

    ax = Axis(f[r,c],
          title = "$(labs[i]) $(unique(df_plot.newname)[1])",
          titlealign = :left,
          xlabel = "Time (years)",
          ylabel = "Proportion", 
          xticks = 0:1:df_stats.time[end], 
          yticks = 0:0.2:1)

    lines!(ax, df_plot.time, df_plot.mean, label = unique(df_plot.newname), color = colors0[i])
    band!(ax, df_plot.time, df_plot.lower, df_plot.upper, color = (colors0[i], 0.2))
    # axislegend(ax, position = :rb; labelsize = 10)
end

f
