cd(@__DIR__)
include("../src/PlotTheme.jl")
include("../src/ThermoProps.jl")
include("../src/components.jl")
include("../src/AAD_calc.jl")
using .PlotTheme
using Clapeyron, CairoMakie, CSV, DataFrames
using GCIdentifier, ChemicalIdentifiers

const param_dir = "../Results/Parameters/"
const data_dir = "../Parameter fitting/Data/"

props_exp = CSV.read(data_dir * "decane_derivative_prop.csv", DataFrame; header=1)

components = groups_from_smiles(["CCCCCCCCCC"])
model_fitted = SAFTgammaMie(components; userlocations = [
        param_dir * "singledata_SAFTgammaMie_1st.csv",
        param_dir * "pairdata_SAFTgammaMie_1st.csv"], idealmodel=JobackIdeal
)
model_lit = SAFTgammaMie(components; idealmodel=JobackIdeal)

Cp_curves = Dict(
    "lit"      => Cp_plot(model_lit, props_exp),
    "fitted"   => Cp_plot(model_fitted, props_exp),
)

fig_Cp = plot_Cp(Cp_curves;
    exp_T   = props_exp.T,
    exp_p   = props_exp.p,
    exp_Cp = props_exp.Cp,
)

save("decane_Cp.png", fig_Cp; px_per_unit=3)


Cv_curves = Dict(
    "lit"      => Cv_plot(model_lit, props_exp),
    "fitted"   => Cv_plot(model_fitted, props_exp),
)

fig_Cv = plot_Cv(Cv_curves;
    exp_T   = props_exp.T,
    exp_p   = props_exp.p,
    exp_Cv = props_exp.Cv,
)

save("decane_Cv.png", fig_Cv; px_per_unit=3)

u_curves = Dict(
    "lit"      => u_plot(model_lit, props_exp),
    "fitted"   => u_plot(model_fitted, props_exp),
)

fig_u = plot_u(u_curves;
    exp_T   = props_exp.T,
    exp_p   = props_exp.p,
    exp_u = props_exp.u,
)

save("decane_u.png", fig_u; px_per_unit=3)