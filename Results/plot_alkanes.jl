using GCIdentifier, ChemicalIdentifiers, Clapeyron, PyCall, CSV, DataFrames
import PyPlot; const plt = PyPlot
plt.matplotlib.use("TkAgg")
println("Set up complete")

ethane = get_groups_from_smiles("CC", SAFTgammaMieGroups)       
propane = get_groups_from_smiles("CCC", SAFTgammaMieGroups) 
pentane = get_groups_from_smiles("CCCCC", SAFTgammaMieGroups) 
octane = get_groups_from_smiles("CCCCCCCC", SAFTgammaMieGroups)  
nonane = get_groups_from_smiles("CCCCCCCCC", SAFTgammaMieGroups) 
decane = get_groups_from_smiles("CCCCCCCCCC", SAFTgammaMieGroups) 
components = [ethane, propane, pentane, octane, nonane, decane]
components_names = ["ethane", "propane" , "pentane", "octane", "nonane", "decane"]
println(components)

ethane_exp = CSV.read("Fitting Data/ethane_sat_p.csv", DataFrame; header=3)
propane_exp = CSV.read("Fitting Data/propane_sat_p.csv", DataFrame; header=3)
pentane_exp = CSV.read("Fitting Data/pentane_sat_p.csv", DataFrame; header=3)
octane_exp = CSV.read("Fitting Data/octane_sat_p.csv", DataFrame; header=3)
nonane_exp = CSV.read("Fitting Data/nonane_sat_p.csv", DataFrame; header=3)
decane_exp = CSV.read("Fitting Data/decane_sat_p.csv", DataFrame; header=3)
exp_data = [ethane_exp, propane_exp, pentane_exp, octane_exp, nonane_exp, decane_exp]
println("experimental data loaded")

N = length(components)

# for i in 1:N
#     component = components[i]
#     name = components_names[i]
#     model = SAFTgammaMie(component)
#     exp = exp_data[i]
#     pressures = unique(exp.p)
    
#     plt.clf()
#     for p_val in pressures
#         subset = filter(row -> row.p == p_val, exp)
        
#         T_rho = subset.T
#         rho_exp = subset.out_rhol
#         N_rho = length(T_rho)
#         rho_calc = zeros(N_rho)

#         println(subset)
#         for j in 1:N_rho
#             T = T_rho[j]
#             rho_calc[j] = molar_density(model, p_val, T; phase=:liquid)  
#         end

#         plt.plot([T_rho], [rho_exp], marker = "o", color = "black")
#         plt.plot(T_rho, rho_calc, label = "p = $(p_val)")
#     end
#     plt.legend(loc="upper right",frameon=false,fontsize=11)
#     plt.title(name)
#     plt.show()
# end


# for i in 1:N
#     component = components[i]
#     name = components_names[i]
#     model = SAFTgammaMie(component)
#     exp = exp_data[i]
    
#     plt.clf()

#     T_rho = exp.T
#     rho_exp = exp.out_rhol
#     N_rho = length(T_rho)
#     rho_calc = zeros(N_rho)

#     for j in 1:N_rho
#         T = T_rho[j]
#         sat = saturation_pressure(model,T)
#         rho_calc[j] = 1/sat[2] 
#     end

#     plt.plot([T_rho], [rho_exp], marker = "o", color = "black")
#     plt.plot(T_rho, rho_calc)
#     plt.title(name)
#     plt.show()
# end

for i in 1:N
    component = components[i]
    name = components_names[i]
    model = SAFTgammaMie(component)
    exp = exp_data[i]
    
    plt.clf()

    T_p = exp.T
    p_exp = exp.out_P
    N_p = length(T_p)
    p_calc = zeros(N_p)

    for j in 1:N_p
        T = T_p[j]
        sat = saturation_pressure(model,T)
        p_calc[j] = sat[1] 
    end

    plt.plot([T_p], [p_exp], marker = "o", color = "black")
    plt.plot(T_p, p_calc)
    plt.title(name)
    plt.show()
end