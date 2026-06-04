cd(@__DIR__)
include("../../src/PlotTheme.jl")
include("../../src/ThermoProps.jl")
include("../../src/components.jl")
include("../../src/AAD_calc.jl")
using .PlotTheme
using CairoMakie
using GCIdentifier, ChemicalIdentifiers, Clapeyron, CSV, DataFrames, Statistics

apply_theme!() 

const param_dir = "../Results/Parameters/"
const data_dir = "../Parameter fitting/Data/"
components = groups_from_smiles(["C"^i for i in 2:10])
comp_names = ["ethane", "propane", "butane", "pentane", "hexane", "heptane", "octane", "nonane", "decane"]
N = length(components)

AAD_results = DataFrame(
    component = String[],
    property = String[],
    AAD_lit = Float64[],
    AAD_fit = Float64[],
    percent_AAD_lit = Float64[],
    percent_AAD_fit = Float64[]
)

for i in 1:N 
    name = comp_names[i]

    model_lit   = SAFTgammaMie(components[i])
    model_fitted = SAFTgammaMie(components[i]; userlocations = [
        param_dir * "singledata_SAFTgammaMie_6th.csv",
        param_dir * "pairdata_SAFTgammaMie_6th.csv",
        # param_dir * "assocdata_SAFTgammaMie.csv",
    ])
        # model_weighted = SAFTgammaMie(components[i]; userlocations = [
        # param_dir * "singledata_SAFTgammaMie_3rd.csv",
        # param_dir * "pairdata_SAFTgammaMie_3rd.csv",
        # param_dir * "assocdata_SAFTgammaMie.csv",
    # ])
    
    Tc_lit,   pc_lit,   vc_lit   = crit_pure(model_lit)
    Tc_fitted, pc_fitted, vc_fitted = crit_pure(model_fitted)
    # Tc_weighted, pc_weighted, vc_weighted = crit_pure(model_weighted)
    
    crits = Dict(
        "Papaioannou" => (Tc_lit,   pc_lit,   vc_lit),
        "fitted" => (Tc_fitted, pc_fitted, vc_fitted),
        # "weighted" => (Tc_weighted, pc_weighted, vc_weighted)
    )   
    
    P_exp    = CSV.read(data_dir * "$(name)_sat_p.csv",    DataFrame; header=3)
    rhol_sat_exp = CSV.read(data_dir * "$(name)_sat_rhol.csv", DataFrame; header=3)
    rhov_sat_exp = CSV.read(data_dir * "$(name)_sat_rhov.csv", DataFrame; header=3)
    rhol_exp = CSV.read(data_dir * "$(name)_rhol.csv", DataFrame; header=3)

    # curves = Dict(
    #     "Papaioannou" => sat_envelope(model_lit, vcat(P_exp.T, Tc_lit)),
    #     "fitted" => sat_envelope(model_fitted, vcat(P_exp.T, Tc_fitted)),
    #     # "weighted" => sat_envelope(model_weighted, vcat(P_exp.T, Tc_weighted))
    # )

    # fig_p = plot_saturation_pressure(curves, crits;
    #     exp_T    = P_exp.T,
    #     exp_p    = P_exp.out_p,
    # )
    # save("vapour_pressure_$(name).png", fig_p, px_per_unit=3)

    # curves = Dict(
    #     "Papaioannou" => sat_envelope(model_lit, LinRange(minimum(rhol_sat_exp.T), Tc_lit, 200)),
    #     "fitted" => sat_envelope(model_fitted, LinRange(minimum(rhol_sat_exp.T), Tc_fitted, 200)),
    #     # "weighted" => sat_envelope(model_weighted, LinRange(minimum(rhol_sat_exp.T), Tc_weighted, 200))
    # )

    # fig_vle = plot_VLE_envelope(curves, crits;
    #     exp_rhol_T = rhol_sat_exp.T,
    #     exp_rhol   = rhol_sat_exp.out_rhol,
    #     exp_rhov_T = rhov_sat_exp.T,
    #     exp_rhov   = rhov_sat_exp.out_rhov
    # )
    # save("$(name)_VLE_envelope.png", fig_vle; px_per_unit=3)


    # rhol_curves = Dict(
    #     "Papaioannou" => rhol_curve(model_lit,    rhol_exp),
    #     "fitted"      => rhol_curve(model_fitted, rhol_exp),
    #     # "weighted"      => rhol_curve(model_weighted, rhol_exp),
    # )

    # fig_rhol = plot_rhol(rhol_curves;
    #     exp_T   = rhol_exp.T,
    #     exp_p   = rhol_exp.p,
    #     exp_rho = rhol_exp.out_rhol,
    # )

    # save("$(name)_rhol_vs_pressure.png", fig_rhol; px_per_unit=3)

    rhol_sat_lit = (sat_envelope(model_lit, rhol_sat_exp.T).rhol)
    rhol_sat_fitted = (sat_envelope(model_fitted, rhol_sat_exp.T).rhol)
    rhol_sat = rhol_sat_exp.out_rhol

    AAD_rhol_sat_lit = AAD(rhol_sat_lit, rhol_sat)
    AAD_rhol_sat_fitted = AAD(rhol_sat_fitted, rhol_sat)

    rhov_sat_lit = (sat_envelope(model_lit, rhov_sat_exp.T).rhov)
    rhov_sat_fitted = (sat_envelope(model_fitted, rhov_sat_exp.T).rhov)
    rhov_sat = rhov_sat_exp.out_rhov

    AAD_rhov_sat_lit = AAD(rhov_sat_lit, rhov_sat)
    AAD_rhov_sat_fitted = AAD(rhov_sat_fitted, rhov_sat)

    p_sat_lit = sat_envelope(model_lit, P_exp.T).p
    p_sat_fitted = sat_envelope(model_fitted, P_exp.T).p
    p_sat = P_exp.out_p 

    AAD_p_lit = AAD(p_sat_lit, p_sat)
    AAD_p_fitted = AAD(p_sat_fitted, p_sat)

    rhol_lit = rhol_curve(model_lit, rhol_exp)
    rhol_fitted = rhol_curve(model_fitted, rhol_exp)
    rhol_exp_vals = rhol_exp.out_rhol

    AAD_rhol_lit = AAD(rhol_lit.rho_vals, rhol_exp_vals)
    AAD_rhol_fitted = AAD(rhol_fitted.rho_vals, rhol_exp_vals)


    for (prop, aad_lit, aad_fit, paad_lit, paad_fit) in [
        ("Saturated Pressure",       AAD_p_lit[1],    AAD_p_fitted[1],    AAD_p_lit[2],    AAD_p_fitted[2]),
        ("Saturated liquid density", AAD_rhol_sat_lit[1],  AAD_rhol_sat_fitted[1],  AAD_rhol_sat_lit[2],  AAD_rhol_sat_fitted[2]),
        ("Saturated vapour density", AAD_rhov_sat_lit[1], AAD_rhov_sat_fitted[1], AAD_rhov_sat_lit[2], AAD_rhov_sat_fitted[2]),
        ("Liquid density",           AAD_rhol_lit[1], AAD_rhol_fitted[1], AAD_rhol_lit[2], AAD_rhol_fitted[2]),
    ]
        push!(AAD_results, (name, prop, aad_lit, aad_fit, paad_lit, paad_fit))
    end
end

CSV.write("alkane_AAD.csv", AAD_results)

summary = combine(groupby(AAD_results, :property),
    :AAD_lit        => mean => :AAD_lit,
    :AAD_fit        => mean => :AAD_fit,
    :percent_AAD_lit => mean => :percent_AAD_lit,
    :percent_AAD_fit => mean => :percent_AAD_fit,
)
insertcols!(summary, 1, :component => fill("Average", nrow(summary)))

CSV.write("alkane_AAD.csv", vcat(AAD_results, summary); append=false)
