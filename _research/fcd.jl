using DifferentialEquations, Plots

τ = 2.0
# The system below assumes repression of Z by Y is strong compared to activation by X.
 # See https://www.sciencedirect.com/science/article/pii/S1097276509008594#fig3
function fcd!(du,u,h,p,t)
    x,y,z = u
    hist_y = h(p, t-τ; idxs=2) # history function, see DDE tutorial for info
    β₁, β₂, α₁, α₂ = p
    du[1] = dx = 0 # Constant input X.
    du[2] = dy = β₁*x - α₁*y
    du[3] = dz = (β₂*x)/hist_y - α₂*z
end


β₁ = 0.60
β₂ = 0.80
α₁ = 0.20
α₂ = 0.50

X₀ = 1.
Y₀ = β₁*X₀/α₁
Z₀ = 0

p = [β₁, β₂, α₁, α₂, τ]
u₀ = [X₀, Y₀, Z₀]

h(p,t) = [X₀, Y₀, Z₀]
h(p,t;idxs=2) = 1.0

tspan = (0.0,10.0)
prob = DDEProblem(fcd!, u₀, h, tspan, p; constant_lag=[τ])

sol = solve(prob)
plot(sol)