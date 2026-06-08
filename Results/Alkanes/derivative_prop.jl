cd(@__DIR__)
include("../../src/PlotTheme.jl")
include("../../src/ThermoProps.jl")
include("../../src/components.jl")
include("../../src/AAD_calc.jl")
using .PlotTheme
using Clapeyron, CairoMakie, CSV, DataFrames, Statistics
using GCIdentifier, ChemicalIdentifiers

const src_dir   = joinpath(@__DIR__, "../../src/")
const data_dir  = joinpath(@__DIR__, "../../Parameter fitting/Data/")
const param_dir = joinpath(@__DIR__, "../../Results/Parameters/")

components = groups_from_smiles(["C"^i for i in 2:10])
comp_names = ["ethane", "propane", "butane", "pentane",  "hexane", "heptane", "octane", "nonane", "decane"]
N = length(comp_names)

AAD_results = DataFrame(
    component = String[],
    property = String[],
    percent_AAD_lit = Float64[],
    percent_AAD_fit = Float64[]
)

for i in 1:N
    name = comp_names[i]
    props_exp = CSV.read(data_dir * "$(name)_derivative_prop.csv", DataFrame; header=1)

    model_fitted = SAFTgammaMie(components[i]; userlocations = [
        param_dir * "singledata_SAFTgammaMie_5th.csv",
        param_dir * "pairdata_SAFTgammaMie_5th.csv"], idealmodel=JobackIdeal
    )
    # model_fitted_4th = SAFTgammaMie(components[i]; userlocations = [
    #     param_dir * "singledata_SAFTgammaMie.csv",
    #     param_dir * "pairdata_SAFTgammaMie.csv"], idealmodel=JobackIdeal
    # )
    model_lit = SAFTgammaMie(components[i]; idealmodel=JobackIdeal)

    Cp_curves = Dict(
        "lit"      => Cp_plot(model_lit, props_exp),
        "fitted"   => Cp_plot(model_fitted, props_exp),
        # "new" => Cp_plot(model_fitted_4th, props_exp)
    )

    fig_Cp = plot_Cp(Cp_curves;
        exp_T   = props_exp.T,
        exp_p   = props_exp.p,
        exp_Cp = props_exp.Cp,
    )

    save("$(name)_Cp.png", fig_Cp; px_per_unit=3)

    mask = props_exp.Phase .!= "supercritical"
    props_filtered = props_exp[mask, :]

    Cv_curves = Dict(
        "lit"      => Cv_plot(model_lit, props_filtered),
        "fitted"   => Cv_plot(model_fitted, props_filtered),
        # "new" => Cv_plot(model_fitted_4th, props_exp)
    )

    fig_Cv = plot_Cv_isothermal(Cv_curves;
        exp_T   = props_exp.T,
        exp_p   = props_exp.p,
        exp_Cv = props_exp.Cv,
    )

    save("$(name)_Cv_isothermal.png", fig_Cv; px_per_unit=3)

    fig_Cv = plot_Cv_isobaric(Cv_curves;
        exp_T   = props_filtered.T,
        exp_p   = props_filtered.p,
        exp_Cv = props_filtered.Cv,
    )

    save("$(name)_Cv_isobaric.png", fig_Cv; px_per_unit=3)

    u_curves = Dict(
        "lit"      => u_plot(model_lit, props_exp),
        "fitted"   => u_plot(model_fitted, props_exp),
        # "new" => u_plot(model_fitted_4th, props_exp)
    )

    fig_u = plot_u(u_curves;
        exp_T   = props_exp.T,
        exp_p   = props_exp.p,
        exp_u = props_exp.u,
    )

    save("$(name)_u.png", fig_u; px_per_unit=3)

    Cp_lit = (Cp_plot(model_lit, props_exp).Cp_vals)
    Cp_fitted = (Cp_plot(model_fitted, props_exp).Cp_vals)
    Cp = props_exp.Cp

    AAD_Cp_lit = AAD(Cp_lit, Cp)
    AAD_Cp_fitted = AAD(Cp_fitted, Cp)

    Cv_lit = (Cv_plot(model_lit, props_filtered).Cv_vals)
    Cv_fitted = (Cv_plot(model_fitted, props_filtered).Cv_vals)
    Cv = props_filtered.Cv

    AAD_Cv_lit = AAD(Cv_lit, Cv)
    AAD_Cv_fitted = AAD(Cv_fitted, Cv)

    u_lit = u_plot(model_lit, props_exp).u_vals
    u_fitted = u_plot(model_fitted, props_exp).u_vals
    u = props_exp.u

    AAD_u_lit = AAD(u_lit, u)
    AAD_u_fitted = AAD(u_fitted, u)


    for (prop, paad_lit, paad_fit) in [
        ("Speed of sound",          AAD_u_lit[2],    AAD_u_fitted[2]),
        ("Isobaric heat capacity",  AAD_Cp_lit[2],  AAD_Cp_fitted[2]),
        ("Isochoric heat capacity", AAD_Cv_lit[2], AAD_Cv_fitted[2]),
    ]
        push!(AAD_results, (name, prop, paad_lit, paad_fit))
    end
end

CSV.write("alkane_AAD.csv", AAD_results)

summary = combine(groupby(AAD_results, :property),
    :percent_AAD_lit => mean => :percent_AAD_lit,
    :percent_AAD_fit => mean => :percent_AAD_fit,
)
insertcols!(summary, 1, :component => fill("Average", nrow(summary)))

CSV.write("alkane_AAD_derive.csv", vcat(AAD_results, summary); append=false)
