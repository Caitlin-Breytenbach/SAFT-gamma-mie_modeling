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
components = groups_from_smiles(["C"^15, "C"^20, "C"^25, "C"^30])
comp_names = ["pentadecane", "eicosane", "pentacosane", "triacontane"]
N = length(components)

AAD_results = DataFrame(
    component = String[],
    property = String[],
    percent_AAD_lit = Float64[],
    percent_AAD_fit = Float64[],
)

all_curves = Dict{String, Any}()

# for i in 1:N
#     name = comp_names[i]

#     model_lit   = SAFTgammaMie(components[i])
#     model_fitted = SAFTgammaMie(components[i]; userlocations = [
#         param_dir * "singledata_SAFTgammaMie_5th.csv",
#         param_dir * "pairdata_SAFTgammaMie_5th.csv",
#         # param_dir * "assocdata_SAFTgammaMie.csv",
#     ])

#     P_exp    = CSV.read(data_dir * "$(name)_sat_p.csv",DataFrame; header=3)

#     T_min = 300
#     T_max = 900
#     T_range = LinRange(T_min, T_max, 200)

#     all_curves["$(name)_Papaioannou"] = sat_envelope(model_lit, collect(T_range))
#     all_curves["$(name)_fitted"] = sat_envelope(model_fitted, collect(T_range))

#     p_sat_lit = sat_envelope(model_lit, P_exp.T).p
#     p_sat_fitted = sat_envelope(model_fitted, P_exp.T).p
#     p_sat = P_exp.out_p 

#     AAD_p_lit = AAD(p_sat_lit, p_sat)
#     AAD_p_fitted = AAD(p_sat_fitted, p_sat)   

#     for (prop, paad_lit, paad_fit) in [
#         ("Saturated Pressure", AAD_p_lit[2], AAD_p_fitted[2]),
#     ]
#         push!(AAD_results, (name, prop, paad_lit, paad_fit))
#     end
# end

# summary = combine(groupby(AAD_results, :property), 
#     :percent_AAD_lit => mean => :percent_AAD_lit,
#     :percent_AAD_fit => mean => :percent_AAD_fit
# )
# insertcols!(summary, 1, :component => fill("Average", nrow(summary)))

# CSV.write("long_chain_AAD.csv", vcat(AAD_results, summary); append=false)

# apply_theme!()
# colours = line_colour(comp_names)

# fig = Figure(size = (700, 500))
# ax = Axis(fig[1,1];
#     xlabel = "Temperature / K",
#     ylabel = "log₁₀(Pressure / Pa)"
# )

# for i in 1:N
#     name = comp_names[i]
#     marker = MARKERS[mod1(i, length(MARKERS))]

#     P_exp    = CSV.read(data_dir * "$(name)_sat_p.csv",DataFrame; header=3)
#     c_lit = all_curves["$(name)_Papaioannou"]
#     c_fit = all_curves["$(name)_fitted"]

#     lines!(ax, c_lit.T, log10.(c_lit.p);
#     color = PALETTE[1], linewidth = LINEWIDTH,
#     label = "Papaioannou")

#     lines!(ax, c_fit.T, log10.(c_fit.p);
#     color = PALETTE[2], linewidth = LINEWIDTH,
#     label = "fitted")

#     scatter!(ax, P_exp.T, log10.(P_exp.out_p);
#     color = EXP_COLOR,
#     marker = marker, 
#     markersize = MARKERSIZE,
#     strokewidth = 1.0,
#     strokecolor = :black,
#     label = name)
# end

# Legend(fig[1,2], ax; merge = true, unique = true)
# save("long_chain_sat_p.png", fig; px_per_unit = 3)

name = comp_names[1]
model_lit   = SAFTgammaMie(components[1])
model_fitted = SAFTgammaMie(components[1]; userlocations = [
    param_dir * "singledata_SAFTgammaMie_5th.csv",
    param_dir * "pairdata_SAFTgammaMie_5th.csv",
    # param_dir * "assocdata_SAFTgammaMie.csv",
])


u_exp = CSV.read(data_dir * "$(name)_u.csv", DataFrame; header=3)

u_curves = Dict(
    "lit"      => u_plot(model_lit, u_exp),
    "fitted"   => u_plot(model_fitted, u_exp),
)

fig_u = plot_u(u_curves;
    exp_T   = u_exp.T,
    exp_p   = u_exp.p,
    exp_u   = u_exp.out_u,
)

save("$(name)_u.png", fig_u; px_per_unit=3)

# components = groups_from_smiles(["C"^4, "C"^10])
# name = ["butane", "decane"]

# model_lit   = SAFTgammaMie(components)
# model_fitted = SAFTgammaMie(components; userlocations = [
#     param_dir * "singledata_SAFTgammaMie_5th.csv",
#     param_dir * "pairdata_SAFTgammaMie_5th.csv",
#     # param_dir * "assocdata_SAFTgammaMie.csv",
# ])

# Tc_lit,   pc_lit,   vc_lit   = crit_pure(model_lit)
# Tc_fitted, pc_fitted, vc_fitted = crit_pure(model_fitted)

# crits = Dict(
#     "Papaioannou" => (Tc_lit,   pc_lit,   vc_lit),
#     "fitted" => (Tc_fitted, pc_fitted, vc_fitted),
# )

# pxy_exp = CSV.read(data_dir * "butane-decane_Pxy.csv", DataFrame; header=3)

# pxy_curves = Dict(
#     "lit"      => _plot(model_lit, u_exp),
#     "fitted"   => u_plot(model_fitted, u_exp),
# )

# fig_u = plot_u(u_curves;
#     exp_T   = u_exp.T,
#     exp_p   = u_exp.p,
#     exp_u   = u_exp.out_u,
# )

# save("$(name)_u.png", fig_u; px_per_unit=3)