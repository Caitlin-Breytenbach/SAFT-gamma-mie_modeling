using GCIdentifier, ChemicalIdentifiers, Clapeyron, Metaheuristics
println("Set up complete")
# Alkanes n-alkanes: ethane to n-decane, using saturated vapour pressure, liquid saturated density and liquid density
# w1=w2=w3=1
# Properties from NIST
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
println(components)

method = SHADE(;options=Options(iterations=200, debug=true));
model = SAFTgammaMie(components)
println("Model built successfully")

toestimate = [
    Dict( #epsilon CH3
        :param => :epsilon,
        :indices => (1,1),
        :lower => 200.,
        :upper => 400.,
        :guess => 250.
    ),
    Dict( #epsilon CH2
        :param => :epsilon,
        :indices => (2,2),
        :lower => 200.,
        :upper => 500.,
        :guess => 470.
    ),
    Dict( #sigma CH3
        :param => :sigma,
        :indices => (1,1),
        :recombine => true,
        :factor => 1e-10,
        :lower => 2.,
        :upper => 5.,
        :guess => 4.
    ),
    Dict( #sigma CH2
        :param => :sigma,
        :indices => (2,2),
        :recombine => true,
        :factor => 1e-10,
        :lower => 2.,
        :upper => 5.,
        :guess => 4.
    ),
    Dict( #Sk CH3
        :param => :shapefactor,
        :indices => (1,1),
        :lower => 0.1,
        :upper => 1.,
        :guess => 0.5
    ),
    Dict( #Sk CH2
        :param => :shapefactor,
        :indices => (2,2),
        :lower => 0.1,
        :upper => 1.,
        :guess => 0.2
    ),
    Dict( #lambda_r CH3
        :param => :lambda_r,
        :indices => (1,1),
        :recombine => true,
        :lower => 8.,
        :upper => 30.,
        :guess => 10.
    ),
    Dict( #lambda_r CH2
        :param => :lambda_r,
        :indices => (2,2),
        :recombine => true,
        :lower => 8.,
        :upper => 30.,
        :guess => 20.
    ),
    Dict( #epsilon CH3-CH2
        :param => :epsilon,
        :indices => (1,2),
        :lower => 200.,
        :upper => 500.,
        :guess => 350.
    ),

];

function rhol(model::EoSModel, T, p)    #K, Pa
    try
        rhol = molar_density(model, p, T; phase=:liquid)
        return rhol             #mol/m3
    catch
        return NaN
    end
end

function saturation_p(model::EoSModel,T)    #K
    try
        sat = saturation_pressure(model,T)
        return sat[1]       #Pa
    catch
        return NaN
    end
end

function saturation_rhol(model::EoSModel,T)   #K
    try
        sat = saturation_pressure(model,T)
        return 1/sat[2]     #mol/m3
    catch
        return NaN
    end    
end


# function saturation_rhov(model::EoSModel,T)   #K
#     try
#         sat = saturation_pressure(model,T)
#         return 1/sat[3]     #mol/m3
#     catch
#         return NaN
#     end
# end

# function en_vap(model::EoSModel,T)           #K
#     try
#         sat = enthalpy_vap(model,T)
#         return sat          #J
#     catch
#         return NaN
#     end
# end

# function bubble_p(model::EoSModel,x,T)   #mole frac, K
#     try
#         bub = bubble_pressure(model,T,[x,1-x])
#         return bub[1], bub[4][1]           #Pa, mol frac 
#     catch
#         return NaN, NaN
#     end
# end

# function bubble_t(model::EoSModel,x,T)   #mole frac, K
#     try
#         bub = bubble_temperature(model,T,[x,1-x])
#         return bub[1], bub[4][1]           #Pa, mol frac 
#     catch
#         return NaN, NaN
#     end
# end

# function binary_density(model::EoSModel,x,P,T)    # mole frac, Pa, K
#     try
#         rho = mass_density(model,P,T,[x,1-x])        #kg/m3
#         return rho
#     catch
#         return NaN
#     end
# end

# function binary_he(model::EoSModel,x,P,T)         # mole frac, Pa, K
#     try
#         he = excess(model,P,T,[x,1-x],enthalpy)       #J/mol
#         return he
#     catch
#         return NaN
#     end
# end

estimator,objective,initial,upper,lower = Estimation(model,toestimate,[
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/ethane_rhol.csv", 
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/propane_rhol.csv", 
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/pentane_rhol.csv",
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/octane_rhol.csv", 
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/nonane_rhol.csv",  
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/decane_rhol.csv",
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/ethane_sat_p.csv", 
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/propane_sat_p.csv", 
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/pentane_sat_p.csv",
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/octane_sat_p.csv", 
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/nonane_sat_p.csv",  
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/decane_sat_p.csv",
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/ethane_sat_rhol.csv", 
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/propane_sat_rhol.csv", 
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/pentane_sat_rhol.csv",
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/octane_sat_rhol.csv", 
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/nonane_sat_rhol.csv",  
"C:/Users/jhbre/Downloads/masters/SAFT-gamma-mie_modeling/Parameter fitting/Data/decane_sat_rhol.csv"
],[:vrmodel]);

println("Estimator built. Starting optimization")

params, model = optimize(objective, estimator, method);

export_model(model);