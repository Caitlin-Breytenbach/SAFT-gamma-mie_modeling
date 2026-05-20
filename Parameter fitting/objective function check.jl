cd(@__DIR__)

println("Working directory: ", pwd())

include("../src/components.jl")
include("../src/ThermoProps.jl")
include("../src/optimizer.jl")
include("../src/AAD_calc.jl")

using NLopt, Clapeyron

components = groups_from_smiles(["CC", "CCC", "CCCC", "CCCCC", "CCCCCC", "CCCCCCC","CCCCCCCC", "CCCCCCCCC", "CCCCCCCCCC"])
comp_names = ["ethane", "propane", "butane", "pentane", "hexane", "heptane", "octane", "nonane", "decane"]

model = SAFTgammaMie(components)
const param_dir = "../Results/Parameters/"

toestimate = [    
    Dict( #epsilon CH2
        :param => :epsilon,
        :indices => (2,2),
        # :recombine => true,
        :lower => 200.,
        :upper => 500.,
        :guess => 470.
    ),
    Dict( #epsilon CH3
        :param => :epsilon,
        :indices => (1,1),
        # :recombine => true,
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

params_lit = [473.39, 256.77,  4.8801, 4.0772, 0.22932, 0.57255, 19.871, 15.050, 350.77]
params_fitted = [336.73, 250.53, 4.3018, 4.0163, 0.31237, 0.59143, 13.456, 14.791, 294.96]
# params_intial = [470, 250, 4.8, 4.0, 0.2, 0.5, 20.0, 10.0, 350]

mse_lit = objective(params_lit)
# mse_initial = objective(params_intial)
mse_fitted = objective(params_fitted)

# println(
# "literature: $(round(mse_lit, digits=10))
# initial: $(round(mse_initial, digits=10))
# fitted: $(round(mse_fitted, digits=10))")

# mse_lit = Clapeyron.objective_function(estimator, params_lit)
# mse_initial = Clapeyron.objective_function(estimator, params_intial)
# mse_fitted = Clapeyron.objective_function(estimator, params_fitted)

# println(
#     "litterature: $mse_lit
#     initial: $mse_initial
#     fitted: $mse_fitted"
# )

function objective_per_file(data_file, params_lit, params_fitted, comps)
    total_mse_lit = 0.0
    total_mse_fitted = 0.0

    total_AAD_lit = 0.0
    total_AAD_fitted = 0.0

    i = 1
    for f in data_file
        estimator_single, objective_single, initial, upper, lower = Estimation(model, toestimate, [data(f)],[:vrmodel])

        model_lit = SAFTgammaMie(comps[i])
        model_fitted = SAFTgammaMie(comps[i]; userlocations = [
            param_dir * "singledata_SAFTgammaMie_1st.csv",
            param_dir * "pairdata_SAFTgammaMie_1st.csv",
            # param_dir * "assocdata_SAFTgammaMie.csv",
        ])

        mse_lit = objective_single(params_lit)
        mse_fitted = objective_single(params_fitted)

        exp = CSV.read(data(f), DataFrame; header=3)


        if occursin("sat_p", f)
            calc_lit = sat_envelope(model_lit, exp.T).p
            calc_fitted = sat_envelope(model_fitted, exp.T).p

            AAD_lit, per_AAD_lit = AAD(calc_lit, exp.out_P)
            AAD_fitted, per_AAD_fitted = AAD(calc_fitted, exp.out_P)            

        elseif  occursin("sat_rhol", f)
            calc_lit =sat_envelope(model_lit, exp.T).rhol
            calc_fitted =sat_envelope(model_fitted, exp.T).rhol

            AAD_lit, per_AAD_lit = AAD(calc_lit, exp.out_rhol)
            AAD_fitted, per_AAD_fitted = AAD(calc_fitted, exp.out_rhol)            
        
        elseif occursin("rhol", f)
            calc_lit = rhol_curve(model_lit, exp).rho_vals
            calc_fitted = rhol_curve(model_fitted, exp).rho_vals

            AAD_lit, per_AAD_lit = AAD(calc_lit, exp.out_rhol)
            AAD_fitted, per_AAD_fitted = AAD(calc_fitted, exp.out_rhol)
        end

        println(
        "file name: $f
        literature: mse = $(round(mse_lit, digits=5)) AAD = $(round(per_AAD_lit, digits=5)) 
        fitted: mse = $(round(mse_fitted, digits=5)) AAD = $(round(per_AAD_fitted, digits=5))")

        total_mse_lit += mse_lit
        total_mse_fitted += mse_fitted

        total_AAD_lit += per_AAD_lit
        total_AAD_fitted += per_AAD_fitted
        i += 1
    end
    n = length(comps)
    AAD_lit_avg = total_AAD_lit/n 
    AAD_fitted_avg = total_AAD_fitted/n
    return total_mse_lit, total_mse_fitted, AAD_lit_avg, AAD_fitted_avg
end

rhol_data = [    
    ("ethane_rhol.csv"), 
    ("propane_rhol.csv"), 
    ("butane_rhol.csv"), 
    ("pentane_rhol.csv"),
    ("hexane_rhol.csv"),
    ("heptane_rhol.csv"),
    ("octane_rhol.csv"), 
    ("nonane_rhol.csv"),  
    ("decane_rhol.csv")
]

rhol_obj_lit, rhol_obj_fitted, rhol_AAD_lit, rhol_AAD_fitted = objective_per_file(rhol_data, params_lit, params_fitted, components)

sat_p_data = [   
    ("ethane_sat_p.csv"), 
    ("propane_sat_p.csv"), 
    ("butane_sat_p.csv"),
    ("pentane_sat_p.csv"),
    ("hexane_sat_p.csv"),
    ("heptane_sat_p.csv"),
    ("octane_sat_p.csv"), 
    ("nonane_sat_p.csv"),  
    ("decane_sat_p.csv")
]

sat_p_obj_lit, sat_p_obj_fitted, sat_p_AAD_lit, sat_p_AAD_fitted = objective_per_file(sat_p_data, params_lit, params_fitted, components)

sat_rhol_data = [    
    ("ethane_sat_rhol.csv"), 
    ("propane_sat_rhol.csv"), 
    ("butane_sat_rhol.csv"),
    ("pentane_sat_rhol.csv"),
    ("hexane_sat_rhol.csv"),
    ("heptane_sat_rhol.csv"),
    ("octane_sat_rhol.csv"), 
    ("nonane_sat_rhol.csv"),  
    ("decane_sat_rhol.csv")
]

sat_rhol_obj_lit, sat_rhol_obj_fitted, sat_rhol_AAD_lit, sat_rhol_AAD_fitted = objective_per_file(sat_rhol_data, params_lit, params_fitted, components)

obj_lit = rhol_obj_lit + sat_p_obj_lit + sat_rhol_obj_lit
obj_fitted = rhol_obj_fitted + sat_p_obj_fitted +sat_rhol_obj_fitted

AAD_lit = (rhol_AAD_lit + sat_p_AAD_lit + sat_rhol_AAD_lit)/3
AAD_fitted = (rhol_AAD_fitted + sat_p_AAD_fitted + sat_rhol_AAD_fitted)/3


println(
"Objective function value
Liquid density: literature = $(round(rhol_obj_lit, digits=10)), fitted = $(round(rhol_obj_fitted, digits=10))
Saturated density: literature = $(round(sat_rhol_obj_lit, digits=10)), fitted = $(round(sat_rhol_obj_fitted, digits=10))
Saturated Pressure: literature = $(round(sat_p_obj_lit, digits=10)), fitted = $(round(sat_p_obj_fitted, digits=10))")

println(
"Average AAD
Liquid density: literature = $(round(rhol_AAD_lit, digits=10)), fitted = $(round(rhol_AAD_fitted, digits=10))
Saturated density: literature = $(round(sat_rhol_AAD_lit, digits=10)), fitted = $(round(sat_rhol_AAD_fitted, digits=10))
Saturated Pressure: literature = $(round(sat_p_AAD_lit, digits=10)), fitted = $(round(sat_p_AAD_fitted, digits=10))")

println(
"Total objective 
Literature: $(round(obj_lit, digits=5))
Fitted: $(round(obj_fitted, digits=5))
Average AAD
Literature: $(round(AAD_lit, digits=5))
Fitted: $(round(AAD_fitted, digits=5))")