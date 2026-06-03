cd(@__DIR__)
include("../../src/PlotTheme.jl")
include("../../src/ThermoProps.jl")
include("../../src/components.jl")
include("../../src/AAD_calc.jl")
using .PlotTheme
using Clapeyron, CairoMakie, CSV, DataFrames
using GCIdentifier, ChemicalIdentifiers

const src_dir   = joinpath(@__DIR__, "../../src/")
const data_dir  = joinpath(@__DIR__, "../../Parameter fitting/Data/")
const param_dir = joinpath(@__DIR__, "../../Results/Parameters/")

components = groups_from_smiles(["CC", "CCC", "CCCC", "CCCCC", "CCCCCC", "CCCCCCC", "CCCCCCCC", "CCCCCCCCC", "CCCCCCCCCC"])
comp_names = ["ethane", "propane", "butane", "pentane",  "hexane", "heptane", "octane", "nonane", "decane"]
N = length(comp_names)

for i in 1:N
    name = comp_names[i]
    props_exp = CSV.read(data_dir * "$(name)_derivative_prop.csv", DataFrame; header=1)

    model_fitted = SAFTgammaMie(components[i]; userlocations = [
        param_dir * "singledata_SAFTgammaMie_1st.csv",
        param_dir * "pairdata_SAFTgammaMie_1st.csv"], idealmodel=JobackIdeal
    )
    model_fitted_4th = SAFTgammaMie(components[i]; userlocations = [
        param_dir * "singledata_SAFTgammaMie.csv",
        param_dir * "pairdata_SAFTgammaMie.csv"], idealmodel=JobackIdeal
    )
    model_lit = SAFTgammaMie(components[i]; idealmodel=JobackIdeal)

    Cp_curves = Dict(
        "lit"      => Cp_plot(model_lit, props_exp),
        "fitted"   => Cp_plot(model_fitted, props_exp),
        "new" => Cp_plot(model_fitted_4th, props_exp)
    )

    fig_Cp = plot_Cp(Cp_curves;
        exp_T   = props_exp.T,
        exp_p   = props_exp.p,
        exp_Cp = props_exp.Cp,
    )

    save("$(name)_Cp.png", fig_Cp; px_per_unit=3)

    Cv_curves = Dict(
        "lit"      => Cv_plot(model_lit, props_exp),
        "fitted"   => Cv_plot(model_fitted, props_exp),
        "new" => Cv_plot(model_fitted_4th, props_exp)
    )

    # fig_Cv = plot_Cv_isothermal(Cv_curves;
    #     exp_T   = props_exp.T,
    #     exp_p   = props_exp.p,
    #     exp_Cv = props_exp.Cv,
    # )

    # save("$(name)_Cv_isothermal.png", fig_Cv; px_per_unit=3)

    fig_Cv = plot_Cv_isobaric(Cv_curves;
        exp_T   = props_exp.T,
        exp_p   = props_exp.p,
        exp_Cv = props_exp.Cv,
    )

    save("$(name)_Cv_isobaric.png", fig_Cv; px_per_unit=3)

    u_curves = Dict(
        "lit"      => u_plot(model_lit, props_exp),
        "fitted"   => u_plot(model_fitted, props_exp),
        "new" => u_plot(model_fitted_4th, props_exp)
    )

    fig_u = plot_u(u_curves;
        exp_T   = props_exp.T,
        exp_p   = props_exp.p,
        exp_u = props_exp.u,
    )

    save("$(name)_u.png", fig_u; px_per_unit=3)
end
