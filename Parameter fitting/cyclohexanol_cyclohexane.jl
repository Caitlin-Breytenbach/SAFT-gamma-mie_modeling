using GCIdentifier, ChemicalIdentifiers, Clapeyron, Metaheuristics
println("Set up complete")
#cyclopentanol = get_groups_from_name("cyclopentanol", SAFTgammaMieGroups; check=true)
cyclohexanol = get_groups_from_name("cyclohexanol", SAFTgammaMieGroups; check=true)
#cycloheptanol = get_groups_from_name("cycloheptanol", SAFTgammaMieGroups; check=true)
cyclohexane = get_groups_from_name("cyclohexane", SAFTgammaMieGroups; check=true)
# two_methylcyclohexanol = get_groups_from_name("2-methylcyclohexanol", SAFTgammaMieGroups; check=true)
# menthol = get_groups_from_name("menthol", SAFTgammaMieGroups; check=true)
components = [cyclohexanol, cyclohexane]
method = SHADE(;options=Options(iterations=10, debug=true));
model = SAFTgammaMie(components)
println("Model built successfully")
# println(components)


toestimate = [
    Dict( #epsilon cCHOH
        :param => :epsilon,
        :indices => (2,2),
        :lower => 200.,
        :upper => 400.,
        :guess => 300.
    ),
    Dict( #sigma cCHOH
        :param => :sigma,
        :indices => (2,2),
        :recombine => true,
        :factor => 1e-10,
        :lower => 2.,
        :upper => 5.,
        :guess => 4.
    ),
    Dict( #Sk cCHOH
        :param => :shapefactor,
        :indices => (2,2),
        :lower => 0.1,
        :upper => 1.,
        :guess => 0.3
    ),
    Dict( #lambda_r cCHOH
        :param => :lambda_r,
        :indices => (2,2),
        :recombine => true,
        :lower => 8.,
        :upper => 30.,
        :guess => 10.
    ),
    Dict( #epsilon cCHOH-cCH2
        :param => :epsilon,
        :indices => (1,2),
        :lower => 200.,
        :upper => 500.,
        :guess => 300.
    ),

    # cCHOH H >=< e1
    Dict(
        :param => :epsilon_assoc,
        :lower => 1000.,
        :upper => 3500.,
        :guess => 1000.
    ),
    Dict(
        :param => :bondvol,
        :lower => 1e-30,
        :upper => 10e-28,
        :guess => 1e-28
    )
];


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

function saturation_rhov(model::EoSModel,T)   #K
    try
        sat = saturation_pressure(model,T)
        return 1/sat[3]     #mol/m3
    catch
        return NaN
    end
end

function en_vap(model::EoSModel,T)           #K
    try
        sat = enthalpy_vap(model,T)
        return sat          #J
    catch
        return NaN
    end
end

function bubble_p(model::EoSModel,x,T)   #mole frac, K
    try
        bub = bubble_pressure(model,T,[x,1-x])
        return bub[1], bub[4][1]           #Pa, mol frac 
    catch
        return NaN, NaN
    end
end

# function bubble_t(model::EoSModel,x,T)   #mole frac, K
#     try
#         bub = bubble_temperature(model,T,[x,1-x])
#         return bub[1], bub[4][1]           #Pa, mol frac 
#     catch
#         return NaN, NaN
#     end
# end

function binary_density(model::EoSModel,x,P,T)    # mole frac, Pa, K
    try
        rho = mass_density(model,P,T,[x,1-x])        #kg/m3
        return rho
    catch
        return NaN
    end
end

function binary_he(model::EoSModel,x,P,T)         # mole frac, Pa, K
    try
        he = excess(model,P,T,[x,1-x],enthalpy)       #J/mol
        return he
    catch
        return NaN
    end
end

# function binary_u(model::EoSModel,x,P,T)        # mole frac, Pa, K
#     try
#         u = speed_of_sound(model,P,T,[x,1-x])       #m/s
#         return u
#     catch
#         return NaN
#     end
# end

estimator,objective,initial,upper,lower = Estimation(model,toestimate,[#"Fitting Data/cyclohexanol-cyclohexane_speed_of_sound.csv",
"Fitting Data/cyclohexanol-cyclohexane_density.csv", 
"Fitting Data/cyclohexanol-cyclohexane_Pxy.csv", 
# "Fitting Data/cyclohexanol-cyclohexane_Txy.csv",
"Fitting Data/cyclohexanol_sat_rhol.csv", 
"Fitting Data/cyclohexanol-cyclohexane_excess_enthalpy.csv",  
"Fitting Data/cyclohexanol_sat_rhov.csv",
"Fitting Data/cyclohexanol_sat_p.csv", 
"Fitting Data/cyclohexanol_enthalpy.csv"
],[:vrmodel]);

println("Estimator built. Starting optimization")

params, model = optimize(objective, estimator, method);

export_model(model);