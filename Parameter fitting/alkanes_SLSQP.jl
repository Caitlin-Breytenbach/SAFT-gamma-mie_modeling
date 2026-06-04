cd(@__DIR__)

println("Working directory: ", pwd())

include("../src/components.jl")
include("../src/ThermoProps.jl")
include("../src/optimizer.jl")

using NLopt, Clapeyron, CSV, DataFrames
components = groups_from_smiles(["C"^i for i in 2:10])
model = SAFTgammaMie(components)

toestimate = [  
    Dict( #epsilon CH3
        :param => :epsilon,
        :indices => (1,1),
        # :recombine => true,
        :lower => 200.,
        :upper => 400.,
        :guess => 300.
    ),  

    Dict( #epsilon CH2
        :param => :epsilon,
        :indices => (2,2),
        # :recombine => true,
        :lower => 200.,
        :upper => 500.,
        :guess => 300.
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

    Dict( #sigma CH2
        :param => :sigma,
        :indices => (2,2),
        :recombine => true,
        :factor => 1e-10,
        :lower => 2.,
        :upper => 5.,
        :guess => 4.
    ), 
    
    Dict( #Sk CH3
        :param => :shapefactor,
        :indices => (1,1),
        :lower => 0.1,
        :upper => 1.,
        :guess => 0.5
    ),

    Dict( #Sk CH2
        :param => :shapefactor,
        :indices => (2,2),
        :lower => 0.15,
        :upper => 1.,
        :guess => 0.5
    ),

    Dict( #lambda_r CH3
        :param => :lambda_r,
        :indices => (1,1),
        :recombine => true,
        :lower => 8.,
        :upper => 30.,
        :guess => 15.
    ),

    Dict( #lambda_r CH2
        :param => :lambda_r,
        :indices => (2,2),
        :recombine => true,
        :lower => 8.,
        :upper => 30.,
        :guess => 15.
    ),

    Dict( #epsilon CH3-CH2
        :param => :epsilon,
        :indices => (1,2),
        :lower => 200.,
        :upper => 500.,
        :guess => 300.
    ),

];

println("Loading estimator...")
data(f) = joinpath(@__DIR__, "Data", f)

w1 = w2 = w3 = 1.
w4 = w5 = 0.4

estimator, objective, initial, upper, lower = Estimation(model, toestimate, [
    (w1, data("ethane_rhol.csv")), 
    (w1, data("propane_rhol.csv")), 
    (w1, data("butane_rhol.csv")), 
    (w1, data("pentane_rhol.csv")),
    (w1, data("hexane_rhol.csv")),
    (w1, data("heptane_rhol.csv")),
    (w1, data("octane_rhol.csv")), 
    (w1, data("nonane_rhol.csv")),  
    (w1, data("decane_rhol.csv")),
    (w2, data("ethane_sat_p.csv")), 
    (w2, data("propane_sat_p.csv")), 
    (w2, data("butane_sat_p.csv")),
    (w2, data("pentane_sat_p.csv")),
    (w2, data("hexane_sat_p.csv")),
    (w2, data("heptane_sat_p.csv")),
    (w2, data("octane_sat_p.csv")), 
    (w2, data("nonane_sat_p.csv")),  
    (w2, data("decane_sat_p.csv")),
    (w3, data("ethane_sat_rhol.csv")), 
    (w3, data("propane_sat_rhol.csv")), 
    (w3, data("butane_sat_rhol.csv")),
    (w3, data("pentane_sat_rhol.csv")),
    (w3, data("heptane_sat_rhol.csv")),
    (w3, data("hexane_sat_rhol.csv")),
    (w3, data("octane_sat_rhol.csv")), 
    (w3, data("nonane_sat_rhol.csv")),  
    (w3, data("decane_sat_rhol.csv")),
    (w4, data("ethane_sat_rhov.csv")), 
    (w4, data("propane_sat_rhov.csv")), 
    (w4, data("butane_sat_rhov.csv")),
    (w4, data("pentane_sat_rhov.csv")),
    (w4, data("heptane_sat_rhov.csv")),
    (w4, data("hexane_sat_rhov.csv")),
    (w4, data("octane_sat_rhov.csv")), 
    (w4, data("nonane_sat_rhov.csv")),  
    (w4, data("decane_sat_rhov.csv")),
    (w5, data("ethane_Cp.csv")), 
    (w5, data("propane_Cp.csv")), 
    (w5, data("butane_Cp.csv")),
    (w5, data("pentane_Cp.csv")),
    (w5, data("heptane_Cp.csv")),
    (w5, data("hexane_Cp.csv")),
    (w5, data("octane_Cp.csv")), 
    (w5, data("nonane_Cp.csv")),  
    (w5, data("decane_Cp.csv"))
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

params_lit = [ 256.77, 473.39, 4.0772, 4.8801, 0.57255, 0.22932, 15.050, 19.871,  350.77]
f_lit = objective(params_lit)
println("Objective lit: $f_lit")

model_opt = optimizer(estimator,objective,initial,upper,lower, 0)
export_model(model_opt)

final_epsilon  = model_opt.params.epsilon.values
final_sigma    = model_opt.params.sigma.values
final_sk       = model_opt.params.shapefactor.values
final_lambda_r = model_opt.params.lambda_r.values

println(
"epsilon: $final_epsilon
sigma: $(final_sigma*1e10)
shape factor: $final_sk
lambda_r: $final_lambda_r")