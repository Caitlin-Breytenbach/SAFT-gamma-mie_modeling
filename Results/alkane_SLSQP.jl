cd(@__DIR__)
include("../src/PlotTheme.jl")
include("../src/ThermoProps.jl")
include("../src/components.jl")
include("../src/AAD_calc.jl")
using .PlotTheme
using CairoMakie
using GCIdentifier, ChemicalIdentifiers, Clapeyron, CSV, DataFrames

apply_theme!() 

const param_dir = "../Results/Parameters/"
const data_dir = "../Parameter fitting/Data/"
components = groups_from_smiles(["CC", "CCC","CCCC", "CCCCC", "CCCCCC", "CCCCCCC", "CCCCCCCC", "CCCCCCCCC", "CCCCCCCCCC"])
comp_names = ["ethane", "propane", "butane", "pentane", "hexane", "heptane", "octane", "nonane", "decane"]
N = length(components)

for i in 1:N 
    name = comp_names[i]

    model_lit   = load_model(components[i])
    model_fitted = load_model(components[i]; userlocations = [
        param_dir * "singledata_SAFTgammaMie_1st.csv",
        param_dir * "pairdata_SAFTgammaMie_1st.csv",
        # param_dir * "assocdata_SAFTgammaMie.csv",
    ])
        # model_weighted = load_model(components[i]; userlocations = [
        # param_dir * "singledata_SAFTgammaMie_2nd.csv",
        # param_dir * "pairdata_SAFTgammaMie_2nd.csv",
        # param_dir * "assocdata_SAFTgammaMie.csv",
    # ])
    
    # Tc_lit,   pc_lit,   vc_lit   = crit_pure(model_lit)
    # Tc_fitted, pc_fitted, vc_fitted = crit_pure(model_fitted)
    # Tc_weighted, pc_weighted, vc_weighted = crit_pure(model_weighted)
    
    # crits = Dict(
    #     "Papaioannou" => (Tc_lit,   pc_lit,   vc_lit),
    #     "fitted" => (Tc_fitted, pc_fitted, vc_fitted),
    #     "weighted" => (Tc_weighted, pc_weighted, vc_weighted)
    # )   
    
    P_exp    = CSV.read(data_dir * "$(name)_sat_p.csv",    DataFrame; header=3)
    rhol_sat_exp = CSV.read(data_dir * "$(name)_sat_rhol.csv", DataFrame; header=3)
    rhol_exp = CSV.read(data_dir * "$(name)_rhol.csv", DataFrame; header=3)

    # curves = Dict(
    #     "Papaioannou" => sat_envelope(model_lit, vcat(P_exp.T, Tc_lit)),
    #     "fitted" => sat_envelope(model_fitted, vcat(P_exp.T, Tc_fitted)),
    #     "weighted" => sat_envelope(model_weighted, vcat(P_exp.T, Tc_weighted))
    # )

    # fig_p = plot_saturation_pressure(curves, crits;
    #     exp_T    = P_exp.T,
    #     exp_p    = P_exp.out_P,
    # )
    # save("vapour_pressure_$(name).png", fig_p, px_per_unit=3)

    # curves = Dict(
    #     "Papaioannou" => sat_envelope(model_lit, LinRange(minimum(rhol_sat_exp.T), Tc_lit, 200)),
    #     "fitted" => sat_envelope(model_fitted, LinRange(minimum(rhol_sat_exp.T), Tc_fitted, 200)),
    #     "weighted" => sat_envelope(model_weighted, LinRange(minimum(rhol_sat_exp.T), Tc_weighted, 200))
    # )

    # fig_vle = plot_VLE_envelope(curves, crits;
    #     exp_rhol_T = rhol_sat_exp.T,
    #     exp_rhol   = rhol_sat_exp.out_rhol,
    # )
    # save("$(name)_VLE_envelope.png", fig_vle; px_per_unit=3)


    # rhol_curves = Dict(
    #     "Papaioannou" => rhol_curve(model_lit,    rhol_exp),
    #     "fitted"      => rhol_curve(model_fitted, rhol_exp),
    #     "weighted"      => rhol_curve(model_weighted, rhol_exp),
    # )

    # fig_rhol = plot_rhol(rhol_curves;
    #     exp_T   = rhol_exp.T,
    #     exp_p   = rhol_exp.p,
    #     exp_rho = rhol_exp.out_rhol,
    # )

    # save("$(name)_rhol_vs_pressure.png", fig_rhol; px_per_unit=3)

    rho_sat_lit = 1 ./(sat_envelope(model_lit, P_exp.T).vl)
    rho_sat_fitted = 1 ./(sat_envelope(model_fitted, P_exp.T).vl)
    rho_sat = rhol_sat_exp.out_rhol

    println("rho sat done")

    AAD_rho_lit = AAD(rho_sat_lit, rho_sat)
    AAD_rho_fitted = AAD(rho_sat_fitted, rho_sat)

    p_sat_lit = sat_envelope(model_lit, P_exp.T).p
    p_sat_fitted = sat_envelope(model_fitted, P_exp.T).p
    p_sat = P_exp.out_P 

    println("p sat done")

    AAD_p_lit = AAD(p_sat_lit, p_sat)
    AAD_p_fitted = AAD(p_sat_fitted, p_sat)

    rhol_lit = rhol_curve(model_lit, rhol_exp)
    rhol_fitted = rhol_curve(model_fitted, rhol_exp)
    rhol_exp_vals = rhol_exp.out_rhol

    AAD_rhol_lit = AAD(rhol_lit.rho_vals, rhol_exp_vals)
    AAD_rhol_fitted = AAD(rhol_fitted.rho_vals, rhol_exp_vals)

    println("rhol done")

    results = DataFrame(
        Property = ["Saturated Pressure", "Saturated liquid density", "liquid density"],
        AAD_lit = [AAD_p_lit, AAD_rho_lit, AAD_rhol_lit],
        AAD_fit = [AAD_p_fitted, AAD_rho_fitted, AAD_rhol_fitted]
    )

    CSV.write("alkane_ADD.csv", results)
end
