## Libraries ----

using CairoMakie

## Functions ----

function zeta(T)
    1 / (1 + exp(2.0 * (T - 30.8)))
end

function phi(T)
    0.95 / (1 + exp(0.8 * (T - 32.5)))
end

function sigma(T)
    0.01 / (1 + exp(-1.2 * (T - 31)))
end

function zeta(T)
    1 / (1 + exp(1.5 * (T - 30.5)))
end

function phi(T)
    1 - (1 / (1.5 + exp(0.6 * (-T + 32))))
end

function sigma(T)
    0.035 / (1 + exp(-1.2 * (T - 31)))
end

## Plot ----

temp = 20:0.01:40
zeta0 = zeta.(temp)
phi0 = phi.(temp)
sigma0 = sigma.(temp)

# c_purple = "#9558b2"
# c_green  = "#389826"
# c_red    = "#cb3c33"

f = Figure(size = (700, 300))

ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Value", yticks = 0:0.2:1)
lines!(ax, temp, zeta0, color = :red, linewidth = 2, label = "ζ")
lines!(ax, temp, phi0, color = :blue, linewidth = 2, label = "ϕ")
axislegend(ax, position = :rt)

ax2 = Axis(f[1,2], xlabel = "Temperature (°C)", ylabel = "Value", yticks = 0:0.01:0.04)
lines!(ax2, temp, sigma0, color = :green, linewidth = 2, label = "σ")
axislegend(ax2, position = :lt)

f
