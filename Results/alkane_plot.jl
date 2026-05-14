include("../src/PlotTheme.jl")
include("../src/ThermoProps.jl")
using .PlotTheme
using CairoMakie
using GCIdentifier, ChemicalIdentifiers, Clapeyron, CSV, DataFrames
 
apply_theme!() 
println("Set up complete")

components = groups_from_smiles(["CC", "CCC", "CCCCC", "CCCCCCCC", "CCCCCCCCC", "CCCCCCCCCC"])
comp_names = ["ethane", "propane", "pentane", "octane", "nonane", "decane"]
N = length(components)
 
const PARAM_DIR = "C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Results/Parameters_testing/"
const DATA_DIR  = "C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/"
for i in 1:N
    name = comp_names[i]

    model_lit   = load_model(components[i])
    model_fitted = load_model(components[i]; userlocations = [
        PARAM_DIR * "singledata_SAFTgammaMie_alkane_SLSQP.csv",
        PARAM_DIR * "pairdata_SAFTgammaMie_alkane_SLSQP.csv",
        PARAM_DIR * "assocdata_SAFTgammaMie_alkane_SLSQP.csv",
    ])
    println("Models set up complete")
    
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
    println("critical points calculated")

    P_exp    = CSV.read(DATA_DIR * "$(name)_sat_p.csv",    DataFrame; header=3)
    rhol_sat_exp = CSV.read(DATA_DIR * "$(name)_sat_rhol.csv", DataFrame; header=3)
    rhol_exp = CSV.read(DATA_DIR * "$(name)_sat_rhol.csv", DataFrame; header=3)
    # H_exp    = CSV.read(DATA_DIR * "$(name)_enthalpy.csv", DataFrame; header=3)
    println("Experimental data loaded")

    # curves = Dict(
    #     "Papaioannou" => sat_envelope(model_lit, vcat(P_exp.T, Tc_lit)),
    #     "fitted" => sat_envelope(model_fitted, vcat(P_exp.T, Tc_fitted)),
    # )

    # fig_p = plot_saturation_pressure(curves, crits;
    #     exp_T    = P_exp.T,
    #     exp_p    = P_exp.out_P,
    # )
    # save("vapour_pressure_$(name).png", fig_p, px_per_unit=3)

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