cd(@__DIR__)

println("Working directory: ", pwd())

include("../src/components.jl")
include("../src/ThermoProps.jl")
include("../src/optimizer.jl")

using NLopt, Clapeyron
components = groups_from_smiles(["CC", "CCC", "CCCC", "CCCCC", "CCCCCC", "CCCCCCC","CCCCCCCC", "CCCCCCCCC", "CCCCCCCCCC"])
model = load_model(components)

toestimate = [    
    Dict( #epsilon CH2
        :param => :epsilon,
        :indices => (2,2),
        :lower => 200.,
        :upper => 500.,
        :guess => 470.
    ),
    Dict( #epsilon CH3
        :param => :epsilon,
        :indices => (1,1),
        :lower => 200.,
        :upper => 400.,
        :guess => 250.
    ),
    Dict( #sigma CH2
        :param => :sigma,
        :indices => (2,2),
        :recombine => true,
        :factor => 1e-10,
        :lower => 2.,
        :upper => 5.,
        :guess => 4.8
    ),    
    Dict( #sigma CH3
        :param => :sigma,
        :indices => (1,1),
        :recombine => true,
        :factor => 1e-10,
        :lower => 2.,
        :upper => 5.,
        :guess => 4.
    ),
    Dict( #Sk CH2
        :param => :shapefactor,
        :indices => (2,2),
        :lower => 0.15,
        :upper => 1.,
        :guess => 0.2
    ),
    Dict( #Sk CH3
        :param => :shapefactor,
        :indices => (1,1),
        :lower => 0.1,
        :upper => 1.,
        :guess => 0.5
    ),
    Dict( #lambda_r CH2
        :param => :lambda_r,
        :indices => (2,2),
        :recombine => true,
        :lower => 8.,
        :upper => 30.,
        :guess => 20.
    ),
    Dict( #lambda_r CH3
        :param => :lambda_r,
        :indices => (1,1),
        :recombine => true,
        :lower => 8.,
        :upper => 30.,
        :guess => 10.
    ),
    Dict( #epsilon CH3-CH2
        :param => :epsilon,
        :indices => (1,2),
        :lower => 200.,
        :upper => 500.,
        :guess => 350.
    ),

];

println("Loading estimator...")
data(f) = joinpath(@__DIR__, "Data", f)

estimator, objective, initial, upper, lower = Estimation(model, toestimate, [
    (1., data("ethane_rhol.csv")), 
    (1., data("propane_rhol.csv")), 
    (1., data("butane_rhol.csv")), 
    (1., data("pentane_rhol.csv")),
    (1., data("hexane_rhol.csv")),
    (1., data("heptane_rhol.csv")),
    (1., data("octane_rhol.csv")), 
    (1., data("nonane_rhol.csv")),  
    (1., data("decane_rhol.csv")),
    (1., data("ethane_sat_p.csv")), 
    (1., data("propane_sat_p.csv")), 
    (1., data("butane_sat_p.csv")),
    (1., data("pentane_sat_p.csv")),
    (1., data("hexane_sat_p.csv")),
    (1., data("heptane_sat_p.csv")),
    (1., data("octane_sat_p.csv")), 
    (1., data("nonane_sat_p.csv")),  
    (1., data("decane_sat_p.csv")),
    (1., data("ethane_sat_rhol.csv")), 
    (1., data("propane_sat_rhol.csv")), 
    (1., data("butane_sat_rhol.csv")),
    (1., data("pentane_sat_rhol.csv")),
    (1., data("heptane_sat_rhol.csv")),
    (1., data("hexane_sat_rhol.csv")),
    (1., data("octane_sat_rhol.csv")), 
    (1., data("nonane_sat_rhol.csv")),  
    (1., data("decane_sat_rhol.csv"))
], [:vrmodel])

n = length(initial)
f0 = objective(initial)
println(
"Numeber parameters: $n
Initial guess: $initial
Lower bounds: $lower
Upper bounds: $upper
objective: $(round(f0; digits = 5))"
)

params_lit = [473.39, 256.77,  4.8801, 4.0772, 0.22932, 0.57255, 19.871, 15.050, 350.77]
f_lit = objective(params_lit)
println("Objective lit: $f_lit")

model_opt = optimizer(estimator,objective,initial,upper,lower)
# export_model(model_opt)

