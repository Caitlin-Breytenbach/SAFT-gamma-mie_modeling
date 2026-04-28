include("../src/PlotTheme.jl")
include("../src/ThermoProps.jl")
using .PlotTheme
using CairoMakie
using GCIdentifier, ChemicalIdentifiers, Clapeyron, CSV, DataFrames
 
apply_theme!() 
println("Set up complete")

ethane = get_groups_from_smiles("CC", SAFTgammaMieGroups)       #T_triple = 90.35K, 0.9*T_c = 274.79K
propane = get_groups_from_smiles("CCC", SAFTgammaMieGroups)     #T_triple = 85.51K, 0.9*T_c = 332.90K
# butane = get_groups_from_smiles("CCCC", SAFTgammaMieGroups)
pentane = get_groups_from_smiles("CCCCC", SAFTgammaMieGroups)   #T_triple = 143.478K, 0.9*T_c = 422.69K
# hexane = get_groups_from_smiles("CCCCCC", SAFTgammaMieGroups)
# heptane = get_groups_from_smiles("CCCCCCC", SAFTgammaMieGroups) 
octane = get_groups_from_smiles("CCCCCCCC", SAFTgammaMieGroups)  #T_triple = 216.418K, 0.9*T_c = 569.57K
nonane = get_groups_from_smiles("CCCCCCCCC", SAFTgammaMieGroups) #T_triple = 219.68K, 0.9*T_c = 535.09K
decane = get_groups_from_smiles("CCCCCCCCCC", SAFTgammaMieGroups) #T_triple = 243.536K, 0.9*T_c = 555.93K
components = [ethane, propane, pentane, octane, nonane, decane]
comp_names = ["ethane", "propane", "pentane", "octane", "nonane", "decane"]
N = length(components)
 
const PARAM_DIR = "C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Results/Parameters_testing/"
const DATA_DIR  = "C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/"
for i in 1:N
    name = comp_names[i]

    model_B   = load_model(components[i])
    model_fitted = load_model(components[i]; userlocations = [
        PARAM_DIR * "singledata_SAFTgammaMie.csv",
        PARAM_DIR * "pairdata_SAFTgammaMie.csv",
        PARAM_DIR * "assocdata_SAFTgammaMie.csv",
    ])
    println("Models set up complete")
    
    Tc_B,   pc_B,   vc_B   = crit_pure(model_B)
    Tc_fitted, pc_fitted, vc_fitted = crit_pure(model_fitted)
    
    crits = Dict(
        "Papaioannou" => (Tc_B,   pc_B),
        "fitted" => (Tc_fitted, pc_fitted),
    )
    crits_full = Dict(
        "Papaioannou" => (Tc_B,   pc_B,   vc_B),
        "fitted" => (Tc_fitted, pc_fitted, vc_fitted),
    )
    println("critical points calculated")

    P_exp    = CSV.read(DATA_DIR * "$(name)_sat_p.csv",    DataFrame; header=3)
    rhol_sat_exp = CSV.read(DATA_DIR * "$(name)_sat_rhol.csv", DataFrame; header=3)
    rhol_exp = CSV.read(DATA_DIR * "$(name)_sat_rhol.csv", DataFrame; header=3)
    # H_exp    = CSV.read(DATA_DIR * "$(name)_enthalpy.csv", DataFrame; header=3)
    println("Experimental data loaded")

    # curves = Dict(
    #     "Papaioannou" => sat_envelope(model_B, vcat(P_exp.T, Tc_B)),
    #     "fitted" => sat_envelope(model_fitted, vcat(P_exp.T, Tc_fitted)),
    # )

    # fig_p = plot_saturation_pressure(curves, crits;
    #     exp_T    = P_exp.T,
    #     exp_p    = P_exp.out_P,
    # )
    # save("vapour_pressure_$(name).png", fig_p, px_per_unit=3)

    curves = Dict(
        "Papaioannou" => sat_envelope(model_B, LinRange(minimum(rhol_sat_exp.T), Tc_B, 200)),
        "fitted" => sat_envelope(model_fitted, LinRange(minimum(rhol_sat_exp.T), Tc_B, 200)),
    )

    fig_vle = plot_VLE_envelope(curves, crits_full;
        exp_rhol_T = rhol_sat_exp.T,
        exp_rhol   = rhol_sat_exp.out_rhol,
    )
    save("$(name)_VLE_envelope.png", fig_vle; px_per_unit=3)

    println("$name complete")


end