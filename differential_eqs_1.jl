#=
DifferentialEquations.jlを用いた覚え書き
=#
using DifferentialEquations
using Plots

function radius_Al(du,u,p,t)
    k,l = p
    du[1] = u[2]
    du[2] = -2*u[2]/t-(k^2-l*(l+1)/t^2)*u[1]
end

function main(l)
        #initial condition
        A = 1.0
        dA = 0.0
        ini_vec_Al = [A,dA]
        #parameters
        k = 1.0
        l = 0.0
        p = (k,l)
        #time interval
        tspan = (10e-8,20.0)

        prob = ODEProblem(radius_Al,ini_vec_Al,tspan,p)
        sol = solve(prob)

    plot(sol,
    title="Bessel function (l=$l)",
    linewidth = 3,
    lc = ["blue" "red"], 
    xaxis = "r",
    label=["Aₗ" "dAₗ"],
    xlims=(0,10),
    ylims=(-0.5,1)
    )
end

main(1.0)