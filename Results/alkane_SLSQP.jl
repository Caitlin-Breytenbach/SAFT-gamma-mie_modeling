cd(@__DIR__)
include("../src/PlotTheme.jl")
include("../src/ThermoProps.jl")
include("../src/components.jl")
using .PlotTheme
using CairoMakie
using GCIdentifier, ChemicalIdentifiers, Clapeyron, CSV, DataFrames

apply_theme!() 

const param_dir = "../Results/Parameters_testing/"
const data_dir = "../Parameter fitting/Data/"
components = groups_from_smiles(["CC", "CCC","CCCC", "CCCCC", "CCCCCC", "CCCCCCC", "CCCCCCCC", "CCCCCCCCC", "CCCCCCCCCC"])
comp_names = ["ethane", "propane", "butane", "pentane", "hexane", "heptane", "octane", "nonane", "decane"]
N = length(components)

for i in 1:N 
    name = comp_names[i]

    model_lit   = load_model(components[i])
    model_fitted = load_model(components[i]; userlocations = [
        param_dir * "singledata_SAFTgammaMie.csv",
        param_dir * "pairdata_SAFTgammaMie.csv",
        # param_dir * "assocdata_SAFTgammaMie.csv",
    ])
    
    Tc_lit,   pc_lit,   vc_lit   = crit_pure(model_lit)
    Tc_fitted, pc_fitted, vc_fitted = crit_pure(model_fitted)
    
    crits = Dict(
        "Papaioannou" => (Tc_lit,   pc_lit),
        "fitted" => (Tc_fitted, pc_fitted),
    )
    crits_full = Dict(
        "Papaioannou" => (Tc_lit,   pc_lit,   vc_lit),
        "fitted" => (Tc_fitted, pc_fitted, vc_fitted),
    )   
    
    P_exp    = CSV.read(data_dir * "$(name)_sat_p.csv",    DataFrame; header=3)
    rhol_sat_exp = CSV.read(data_dir * "$(name)_sat_rhol.csv", DataFrame; header=3)
    # rhol_exp = CSV.read(data_dir * "$(name)_rhol.csv", DataFrame; header=3)

    curves = Dict(
        "Papaioannou" => sat_envelope(model_lit, vcat(P_exp.T, Tc_lit)),
        "fitted" => sat_envelope(model_fitted, vcat(P_exp.T, Tc_fitted)),
    )

    fig_p = plot_saturation_pressure(curves, crits;
        exp_T    = P_exp.T,
        exp_p    = P_exp.out_P,
    )
    save("vapour_pressure_$(name).png", fig_p, px_per_unit=3)

    curves = Dict(
        "Papaioannou" => sat_envelope(model_lit, LinRange(minimum(rhol_sat_exp.T), Tc_lit, 200)),
        "fitted" => sat_envelope(model_fitted, LinRange(minimum(rhol_sat_exp.T), Tc_fitted, 200)),
    )

    fig_vle = plot_VLE_envelope(curves, crits_full;
        exp_rhol_T = rhol_sat_exp.T,
        exp_rhol   = rhol_sat_exp.out_rhol,
    )
    save("$(name)_VLE_envelope.png", fig_vle; px_per_unit=3)

    println("$name complete")
end