include("../src/PlotTheme.jl")
include("../src/ThermoProps.jl")
using .PlotTheme
using CairoMakie
using GCIdentifier, ChemicalIdentifiers, Clapeyron, CSV, DataFrames
 
apply_theme!() 
println("Set up complete")

cyclohexanol = get_groups_from_name("cyclohexanol", SAFTgammaMieGroups; check=true)
 
const PARAM_DIR = "C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Results/Parameters_testing/"
const DATA_DIR  = "C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/"
 
model_B   = load_model(cyclohexanol)
model_ECA = load_model(cyclohexanol; userlocations = [
    PARAM_DIR * "singledata_SAFTgammaMie_cCHOH_100ECA.csv",
    PARAM_DIR * "pairdata_SAFTgammaMie_cCHOH-cCH2_100ECA.csv",
    PARAM_DIR * "assocdata_SAFTgammaMie_cCHOH_100ECA.csv",
])
println("Models set up complete")
 
Tc_B,   pc_B,   vc_B   = crit_pure(model_B)
Tc_ECA, pc_ECA, vc_ECA = crit_pure(model_ECA)
 
curves = Dict(
    "Bernet" => sat_envelope(model_B,   LinRange(300.0, Tc_B,   100)),
    "ECA100" => sat_envelope(model_ECA, LinRange(300.0, Tc_ECA, 100)),
)
crits = Dict(
    "Bernet" => (Tc_B,   pc_B),
    "ECA100" => (Tc_ECA, pc_ECA),
)
crits_full = Dict(
    "Bernet" => (Tc_B,   pc_B,   vc_B),
    "ECA100" => (Tc_ECA, pc_ECA, vc_ECA),
)
println("critical points calculated")

P_exp    = CSV.read(DATA_DIR * "cyclohexanol_sat_p.csv",    DataFrame; header=3)
rhol_exp = CSV.read(DATA_DIR * "cyclohexanol_sat_rhol.csv", DataFrame; header=3)
rhov_exp = CSV.read(DATA_DIR * "cyclohexanol_sat_rhov.csv", DataFrame; header=3)
H_exp    = CSV.read(DATA_DIR * "cyclohexanol_enthalpy.csv", DataFrame; header=3)
println("Experimental data loaded")


fig_p = plot_saturation_pressure(curves, crits;
    exp_T    = P_exp.T,
    exp_p    = P_exp.out_p,
    exp_crit = (647.1, 4303.912e3),
    Tlims    = (300, 700),
)
save("vapour_pressure.png", fig_p, px_per_unit=3)
