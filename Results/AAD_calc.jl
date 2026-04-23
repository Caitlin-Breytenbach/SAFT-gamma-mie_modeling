using GCIdentifier, ChemicalIdentifiers, Clapeyron, PyCall, CSV, DataFrames
println("Set up complete")

cyclohexanol = get_groups_from_name("cyclohexanol", SAFTgammaMieGroups; check=true)
cyclohexane = get_groups_from_name("cyclohexane", SAFTgammaMieGroups; check=true)
components = [cyclohexanol, cyclohexane]
println(components)

cyclohexanol_B = SAFTgammaMie(cyclohexanol)
cyclohexanol_ECA100 = SAFTgammaMie(cyclohexanol; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_100ECA.csv", 
"Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_100ECA.csv",
"Parameters/assocdata_SAFTgammaMie_cCHOH_100ECA.csv"])
cyclohexanol_ECA1000 = SAFTgammaMie(cyclohexanol; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_1000_2ECA.csv", 
"Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_1000_2ECA.csv",
"Parameters/assocdata_SAFTgammaMie_cCHOH_1000_2ECA.csv"])
cyclohexanol_SHADE100 = SAFTgammaMie(cyclohexanol; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_100SHADE.csv", 
"Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_100SHADE.csv",
"Parameters/assocdata_SAFTgammaMie_cCHOH_100SHADE.csv"])
# cyclohexanol_SHADE500 = SAFTgammaMie(cyclohexanol; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_500SHADE.csv", 
# "Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_500SHADE.csv",
# "Parameters/assocdata_SAFTgammaMie_cCHOH_500SHADE.csv"])
cyclohexanol_models = [cyclohexanol_B, cyclohexanol_ECA100, cyclohexanol_ECA1000, cyclohexanol_SHADE100] #, cyclohexanol_SHADE500]

# binary_B = SAFTgammaMie(components)
# binary_ECA100 = SAFTgammaMie(components; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_100ECA.csv", 
# "Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_100ECA.csv",
# "Parameters/assocdata_SAFTgammaMie_cCHOH_100ECA.csv"])
# binary_ECA1000 = SAFTgammaMie(components; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_1000_2ECA.csv", 
# "Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_1000_2ECA.csv",
# "Parameters/assocdata_SAFTgammaMie_cCHOH_1000_2ECA.csv"])
# binary_models = [binary_B, binary_ECA100, binary_ECA1000]
println("model set up complete")

P_exp = CSV.read("Fitting Data/cyclohexanol_sat_p.csv", DataFrame; header=3)
rhov_exp = CSV.read("Fitting Data/cyclohexanol_sat_rhov.csv", DataFrame; header=3)
rhol_exp = CSV.read("Fitting Data/cyclohexanol_sat_rhol.csv", DataFrame; header=3)
hvap_exp = CSV.read("Fitting Data/cyclohexanol_enthalpy.csv", DataFrame; header=3)
# Pxy_exp = CSV.read("Fitting Data/cyclohexanol-cyclohexane_Pxy.csv", DataFrame; header=3)
# Pxy_exp423 = filter(row -> row.T == 423.56, Pxy_exp)
# Txy_exp = CSV.read("Fitting Data/cyclohexanol-cyclohexane_Txy.csv", DataFrame; header=3)
# rhobinary_exp = CSV.read("Fitting Data/cyclohexanol-cyclohexane_density.csv", DataFrame; header=3)
# rhobinary_exp293 = filter(row -> row.T ==293.15, rhobinary_exp)
# Hbinary_exp = CSV.read("Fitting Data/cyclohexanol-cyclohexane_excess_enthalpy.csv",  DataFrame; header=3)
# ubinary_exp = CSV.read("Fitting Data/cyclohexanol-cyclohexane_speed_of_sound.csv", DataFrame; header=3)
# ubinary_exp298 = filter(row -> row.T ==298.15, ubinary_exp)
println("experimental data loaded")

function sat(model, T)

    N = length(T)
    psat = zeros(N)
    vl   = zeros(N)
    vv   = zeros(N)
    hvap = zeros(N)

    for i in 1:N
        if i==1
            sat = saturation_pressure(model, T[i])
            psat[i] = sat[1]
            vl[i] = sat[2]
            vv[i] = sat[3]
        else
            v0 = [vl[i-1],vv[i-1]]
            sat = saturation_pressure(model, T[i]; v0=v0)
            psat[i] = sat[1]
            vl[i] = sat[2]
            vv[i] = sat[3]
        end
        hvap[i] = enthalpy_vap(model, T[i])
    end
    rhov = 1 ./vv
    rhol = 1 ./vl
    return psat, rhov, rhol, hvap
end

T_P = P_exp.T
N_P = length(T_P)
psat_exp = P_exp.out_p

println("calculating ADD")
ADD_psat = zeros(length(cyclohexanol_models))
percentADD_psat = zeros(length(cyclohexanol_models))

for j in 1:length(cyclohexanol_models)

    psat_calc = sat(cyclohexanol_models[j],T_P)[1]
    
    ADD_psat_values = zeros(N_P)
    for i in 1:N_P
        ADD_psat_values[i] = abs(psat_exp[i]-psat_calc[i])
    end
    ADD_psat[j] = (1/N_P)*(sum(ADD_psat_values))

    percentADD_psat_values = zeros(N_P)
    for i in 1:N_P
        percentADD_psat_values[i] = abs((psat_exp[i]-psat_calc[i])/psat_exp[i])*100
    end
    percentADD_psat[j] = (1/N_P)*(sum(percentADD_psat_values))
end

model_names = ["Bernet", "ECA100", "ECA1000", "SHADE100"] #"SHADE500"]

# results = DataFrame(
#     Model = model_names,
#     ADD_psat_kPa = ADD_psat ./1e3,
#     Percent_ADD_psat = percentADD_psat
# )

# println(results)

T_rhov = rhov_exp.T
N_rhov = length(T_rhov)
rhov_exp = rhov_exp.out_rhov

println("calculating ADD")
ADD_rhov = zeros(length(cyclohexanol_models))
percentADD_rhov = zeros(length(cyclohexanol_models))

for j in 1:length(cyclohexanol_models)

    rhov_calc = sat(cyclohexanol_models[j],T_rhov)[2]
    
    ADD_rhov_values = zeros(N_rhov)
    for i in 1:N_rhov
        ADD_rhov_values[i] = abs(rhov_exp[i]-rhov_calc[i])
    end
    ADD_rhov[j] = (1/N_rhov)*(sum(ADD_rhov_values))

    percentADD_rhov_values = zeros(N_rhov)
    for i in 1:N_rhov
        percentADD_rhov_values[i] = abs((rhov_exp[i]-rhov_calc[i])/rhov_exp[i])*100
    end
    percentADD_rhov[j] = (1/N_rhov)*(sum(percentADD_rhov_values))
end

# results = DataFrame(
#     Model = model_names,
    # ADD_rhov = ADD_rhov,
    # Percent_ADD_rhov = percentADD_rhov
# )

# println(results)

T_rhol = rhol_exp.T
N_rhol = length(T_rhol)
rhol_exp = rhol_exp.out_rhol

println("calculating ADD")
ADD_rhol = zeros(length(cyclohexanol_models))
percentADD_rhol = zeros(length(cyclohexanol_models))

for j in 1:length(cyclohexanol_models)

    rhol_calc = sat(cyclohexanol_models[j],T_rhol)[3]
    
    ADD_rhol_values = zeros(N_rhol)
    for i in 1:N_rhol
        ADD_rhol_values[i] = abs(rhol_exp[i]-rhol_calc[i])
    end
    ADD_rhol[j] = (1/N_rhol)*(sum(ADD_rhol_values))

    percentADD_rhol_values = zeros(N_rhol)
    for i in 1:N_rhol
        percentADD_rhol_values[i] = abs((rhol_exp[i]-rhol_calc[i])/rhol_exp[i])*100
    end
    percentADD_rhol[j] = (1/N_rhol)*(sum(percentADD_rhol_values))
end

# results = DataFrame(
#     Model = model_names,
    # ADD_rhol = ADD_rhol,
    # Percent_ADD_rhol = percentADD_rhol
# )

# println(results)

T_hvap = hvap_exp.T
N_hvap = length(T_hvap)
hvap_exp = hvap_exp.out_en_vap

println("calculating ADD")
ADD_hvap = zeros(length(cyclohexanol_models))
percentADD_hvap = zeros(length(cyclohexanol_models))

for j in 1:length(cyclohexanol_models)

    hvap_calc = sat(cyclohexanol_models[j],T_hvap)[4]
    
    ADD_hvap_values = zeros(N_hvap)
    for i in 1:N_hvap
        ADD_hvap_values[i] = abs(hvap_exp[i]-hvap_calc[i])
    end
    ADD_hvap[j] = (1/N_hvap)*(sum(ADD_hvap_values))

    percentADD_hvap_values = zeros(N_hvap)
    for i in 1:N_hvap
        percentADD_hvap_values[i] = abs((hvap_exp[i]-hvap_calc[i])/hvap_exp[i])*100
    end
    percentADD_hvap[j] = (1/N_hvap)*(sum(percentADD_hvap_values))
end

results = DataFrame(
    Model = model_names,
    ADD_psat_kPa = ADD_psat ./1e3,
    Percent_ADD_psat = percentADD_psat,
    ADD_rhov_molm3 = ADD_rhov,
    Percent_ADD_rhov = percentADD_rhov,
    ADD_rhol_molm3 = ADD_rhol,
    Percent_ADD_rhol = percentADD_rhol,
    ADD_hvap_kJmol = ADD_hvap,
    Percent_ADD_hvap = percentADD_hvap
)

println(results)
CSV.write("Cyclohexanol_ADD.csv", results)
println("Results successfully exported.")

# T_Pxy = [423.56, 443.6, 463.62]

# function Pxy(model, T, x)
#     N = length(x)
#     p_bub    = zeros(N)
#     y    = zeros(N)
#     p_dew = zeros(N)

#     for i in 1:N
#         bub = bubble_pressure(model, T, [x[i], 1-x[i]])
#         p_bub[i] = bub[1]
#         y[i] = bub[4][1]

#         dew = dew_pressure(model, T, [y[i], 1-y[i]])
#         p_dew[i] = dew[1]
#     end
#     return p_bub, p_dew  
# end

# for k in 1:length(T_Pxy)
#     T = T_Pxy[k]
#     Pxy_expT = filter(row -> row.T == T, Pxy_exp)
#     x = Pxy_expT.x
#     N_Pxy = length(x)
#     P_exp = Pxy_expT.out_p
#     y_exp = Pxy_expT.out_y

#     for j in 1:length(binary_models)
#         ADD_Pdew = zeros(N_Pxy)
#         ADD_Pbub = zeros(N_Pxy)
#         percentADD_Pdew = zeros(N_Pxy)


#     end

# end