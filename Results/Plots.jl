using GCIdentifier, ChemicalIdentifiers, Clapeyron, PyCall, CSV, DataFrames
import PyPlot; const plt = PyPlot
plt.matplotlib.use("TkAgg")
println("Set up complete")

# cyclohexanol = get_groups_from_name("cyclohexanol", SAFTgammaMieGroups; check=true)
# cyclohexane = get_groups_from_name("cyclohexane", SAFTgammaMieGroups; check=true)
ethane = get_groups_from_smiles("CC", SAFTgammaMieGroups)       #T_triple = 90.35K, 0.9*T_c = 274.79K
propane = get_groups_from_smiles("CCC", SAFTgammaMieGroups) 
pentane = get_groups_from_smiles("CCCCC", SAFTgammaMieGroups) 
octane = get_groups_from_smiles("CCCCCCCC", SAFTgammaMieGroups)  #T_triple = 216.418K, 0.9*T_c = 569.57K
nonane = get_groups_from_smiles("CCCCCCCCC", SAFTgammaMieGroups) #T_triple = 219.68K, 0.9*T_c = 535.09K
decane = get_groups_from_smiles("CCCCCCCCCC", SAFTgammaMieGroups) #T_triple = 243.536K, 0.9*T_c = 555.93K
components = [ethane, propane, pentane, octane, nonane, decane]

println(components)

cyclohexanol_B = SAFTgammaMie(cyclohexanol)
cyclohexanol_SHADE100 = SAFTgammaMie(cyclohexanol; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_100SHADE.csv", 
"Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_100SHADE.csv",
"Parameters/assocdata_SAFTgammaMie_cCHOH_100SHADE.csv"])
cyclohexanol_ECA100 = SAFTgammaMie(cyclohexanol; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_100ECA.csv", 
"Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_100ECA.csv",
"Parameters/assocdata_SAFTgammaMie_cCHOH_100ECA.csv"])
cyclohexanol_ECA1000 = SAFTgammaMie(cyclohexanol; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_1000_2ECA.csv", 
"Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_1000_2ECA.csv",
"Parameters/assocdata_SAFTgammaMie_cCHOH_1000_2ECA.csv"])
# cyclohexanol_SHADE500 = SAFTgammaMie(cyclohexanol; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_500SHADE.csv", 
# "Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_500SHADE.csv",
# "Parameters/assocdata_SAFTgammaMie_cCHOH_500SHADE.csv"])
# cyclohexanol_DE500 = SAFTgammaMie(cyclohexanol; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_1000DE.csv", 
# "Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_1000DE.csv",
# "Parameters/assocdata_SAFTgammaMie_cCHOH_1000DE.csv"])

binary_B = SAFTgammaMie(components)
binary_SHADE100 = SAFTgammaMie(components; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_100SHADE.csv", 
"Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_100SHADE.csv",
"Parameters/assocdata_SAFTgammaMie_cCHOH_100SHADE.csv"])
binary_ECA100 = SAFTgammaMie(components; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_100ECA.csv", 
"Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_100ECA.csv",
"Parameters/assocdata_SAFTgammaMie_cCHOH_100ECA.csv"])
binary_ECA1000 = SAFTgammaMie(components; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_1000_2ECA.csv", 
"Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_1000_2ECA.csv",
"Parameters/assocdata_SAFTgammaMie_cCHOH_1000_2ECA.csv"])
# binary_SHADE500 = SAFTgammaMie(components; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_500SHADE.csv", 
# "Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_500SHADE.csv",
# "Parameters/assocdata_SAFTgammaMie_cCHOH_500SHADE.csv"])
# binary_DE500 = SAFTgammaMie(components; userlocations = ["Parameters/singledata_SAFTgammaMie_cCHOH_1000DE.csv", 
# "Parameters/pairdata_SAFTgammaMie_cCHOH-cCH2_1000DE.csv",
# "Parameters/assocdata_SAFTgammaMie_cCHOH_1000DE.csv"])
println("model set up complete")

P_exp = CSV.read("Fitting Data/cyclohexanol_sat_p.csv", DataFrame; header=3)
rhov_exp = CSV.read("Fitting Data/cyclohexanol_sat_rhov.csv", DataFrame; header=3)
rhol_exp = CSV.read("Fitting Data/cyclohexanol_sat_rhol.csv", DataFrame; header=3)
H_exp = CSV.read("Fitting Data/cyclohexanol_enthalpy.csv", DataFrame; header=3)
Pxy_exp = CSV.read("Fitting Data/cyclohexanol-cyclohexane_Pxy.csv", DataFrame; header=3)
Pxy_exp423 = filter(row -> row.T == 423.56, Pxy_exp)
Txy_exp = CSV.read("Fitting Data/cyclohexanol-cyclohexane_Txy.csv", DataFrame; header=3)
rhobinary_exp = CSV.read("Fitting Data/cyclohexanol-cyclohexane_density.csv", DataFrame; header=3)
rhobinary_exp293 = filter(row -> row.T ==293.15, rhobinary_exp)
Hbinary_exp = CSV.read("Fitting Data/cyclohexanol-cyclohexane_excess_enthalpy.csv",  DataFrame; header=3)
ubinary_exp = CSV.read("Fitting Data/cyclohexanol-cyclohexane_speed_of_sound.csv", DataFrame; header=3)
ubinary_exp298 = filter(row -> row.T ==298.15, ubinary_exp)
println("experimental data loaded")

(Tc_B, pc_B, vc_B) = crit_pure(cyclohexanol_B)
(Tc_ECA100, pc_ECA100, vc_ECA100) = crit_pure(cyclohexanol_ECA100)
(Tc_SHADE100, pc_SHADE100, vc_SHADE100) = crit_pure(cyclohexanol_SHADE100)
(Tc_ECA1000, pc_ECA1000, vc_ECA1000) = crit_pure(cyclohexanol_ECA1000)
# (Tc_SHADE500, pc_SHADE500, vc_SHADE500) = crit_pure(cyclohexanol_SHADE500)
(Tc, pc) = (647.1, 4303.912*1000,)
# (Tc_DE500, pc_DE500, vc_DE500) = crit_pure(cyclohexanol_ECA1000)

N    = 100

T_B    = LinRange(300, Tc_B,  N)
psat_B = zeros(N)
vl_B   = zeros(N)
vv_B   = zeros(N)

hL_B   = zeros(N)
hV_B   = zeros(N)
hvap_B = zeros(N)
cpL_B  = zeros(N)
cpV_B  = zeros(N)

for i in 1:N
    if i==1
        sat = saturation_pressure(cyclohexanol_B, T_B[i])
        psat_B[i] = sat[1]
        vl_B[i] = sat[2]
        vv_B[i] = sat[3]
    else
        v0 = [vl_B[i-1],vv_B[i-1]]
        sat = saturation_pressure(cyclohexanol_B, T_B[i]; v0=v0)
        psat_B[i] = sat[1]
        vl_B[i] = sat[2]
        vv_B[i] = sat[3]
    end
    hL_B[i]  = Clapeyron.VT_enthalpy(cyclohexanol_B,vl_B[i],T_B[i],[1.])
    hV_B[i]  = Clapeyron.VT_enthalpy(cyclohexanol_B,vv_B[i],T_B[i],[1.])
    hvap_B[i] = enthalpy_vap(cyclohexanol_B,T_B[i])
    cpL_B[i] = Clapeyron.VT_isobaric_heat_capacity(cyclohexanol_B,vl_B[i],T_B[i],[1.])
    cpV_B[i] = Clapeyron.VT_isobaric_heat_capacity(cyclohexanol_B,vv_B[i],T_B[i],[1.])
end
N    = 100

T_SHADE100    = LinRange(300, Tc_SHADE100,  N)
psat_SHADE100 = zeros(N)
vl_SHADE100   = zeros(N)
vv_SHADE100   = zeros(N)

hL_SHADE100   = zeros(N)
hV_SHADE100   = zeros(N)
cpL_SHADE100  = zeros(N)
cpV_SHADE100  = zeros(N)

for i in 1:N
    if i==1
        sat = saturation_pressure(cyclohexanol_SHADE100, T_SHADE100[i])
        psat_SHADE100[i] = sat[1]
        vl_SHADE100[i] = sat[2]
        vv_SHADE100[i] = sat[3]
    else
        v0 = [vl_SHADE100[i-1],vv_SHADE100[i-1]]
        sat = saturation_pressure(cyclohexanol_SHADE100, T_SHADE100[i]; v0=v0)
        psat_SHADE100[i] = sat[1]
        vl_SHADE100[i] = sat[2]
        vv_SHADE100[i] = sat[3]
    end
    hL_SHADE100[i]  = Clapeyron.VT_enthalpy(cyclohexanol_SHADE100,vl_SHADE100[i],T_SHADE100[i],[1.])
    hV_SHADE100[i]  = Clapeyron.VT_enthalpy(cyclohexanol_SHADE100,vv_SHADE100[i],T_SHADE100[i],[1.])
    cpL_SHADE100[i] = Clapeyron.VT_isobaric_heat_capacity(cyclohexanol_SHADE100,vl_SHADE100[i],T_SHADE100[i],[1.])
    cpV_SHADE100[i] = Clapeyron.VT_isobaric_heat_capacity(cyclohexanol_SHADE100,vv_SHADE100[i],T_SHADE100[i],[1.])
end

T_ECA100    = LinRange(300, Tc_ECA100,  N)
psat_ECA100 = zeros(N)
vl_ECA100   = zeros(N)
vv_ECA100   = zeros(N)

hL_ECA100   = zeros(N)
hV_ECA100   = zeros(N)
cpL_ECA100  = zeros(N)
cpV_ECA100  = zeros(N)

for i in 1:N
    if i==1
        sat = saturation_pressure(cyclohexanol_ECA100, T_ECA100[i])
        psat_ECA100[i] = sat[1]
        vl_ECA100[i] = sat[2]
        vv_ECA100[i] = sat[3]
    else
        v0 = [vl_ECA100[i-1],vv_ECA100[i-1]]
        sat = saturation_pressure(cyclohexanol_ECA100, T_ECA100[i]; v0=v0)
        psat_ECA100[i] = sat[1]
        vl_ECA100[i] = sat[2]
        vv_ECA100[i] = sat[3]
    end
    hL_ECA100[i]  = Clapeyron.VT_enthalpy(cyclohexanol_ECA100,vl_ECA100[i],T_ECA100[i],[1.])
    hV_ECA100[i]  = Clapeyron.VT_enthalpy(cyclohexanol_ECA100,vv_ECA100[i],T_ECA100[i],[1.])
    cpL_ECA100[i] = Clapeyron.VT_isobaric_heat_capacity(cyclohexanol_ECA100,vl_ECA100[i],T_ECA100[i],[1.])
    cpV_ECA100[i] = Clapeyron.VT_isobaric_heat_capacity(cyclohexanol_ECA100,vv_ECA100[i],T_ECA100[i],[1.])
end

T_ECA1000    = LinRange(300, Tc_ECA1000,  N)
psat_ECA1000 = zeros(N)
vl_ECA1000   = zeros(N)
vv_ECA1000   = zeros(N)

hL_ECA1000   = zeros(N)
hV_ECA1000   = zeros(N)
cpL_ECA1000  = zeros(N)
cpV_ECA1000  = zeros(N)

for i in 1:N
    if i==1
        sat = saturation_pressure(cyclohexanol_ECA1000, T_ECA1000[i])
        psat_ECA1000[i] = sat[1]
        vl_ECA1000[i] = sat[2]
        vv_ECA1000[i] = sat[3]
    else
        v0 = [vl_ECA1000[i-1],vv_ECA1000[i-1]]
        sat = saturation_pressure(cyclohexanol_ECA1000, T_ECA1000[i]; v0=v0)
        psat_ECA1000[i] = sat[1]
        vl_ECA1000[i] = sat[2]
        vv_ECA1000[i] = sat[3]
    end
    hL_ECA1000[i]  = Clapeyron.VT_enthalpy(cyclohexanol_ECA1000,vl_ECA1000[i],T_ECA1000[i],[1.])
    hV_ECA1000[i]  = Clapeyron.VT_enthalpy(cyclohexanol_ECA1000,vv_ECA1000[i],T_ECA1000[i],[1.])
    cpL_ECA1000[i] = Clapeyron.VT_isobaric_heat_capacity(cyclohexanol_ECA1000,vl_ECA1000[i],T_ECA1000[i],[1.])
    cpV_ECA1000[i] = Clapeyron.VT_isobaric_heat_capacity(cyclohexanol_ECA1000,vv_ECA1000[i],T_ECA1000[i],[1.])
end

# T_SHADE500    = LinRange(300, Tc_SHADE500,  N)
# psat_SHADE500 = zeros(N)
# vl_SHADE500   = zeros(N)
# vv_SHADE500   = zeros(N)

# hL_SHADE500   = zeros(N)
# hV_SHADE500   = zeros(N)
# cpL_SHADE500  = zeros(N)
# cpV_SHADE500  = zeros(N)

# for i in 1:N
#     if i==1
#         sat = saturation_pressure(cyclohexanol_SHADE500, T_SHADE500[i])
#         psat_SHADE500[i] = sat[1]
#         vl_SHADE500[i] = sat[2]
#         vv_SHADE500[i] = sat[3]
#     else
#         v0 = [vl_SHADE500[i-1],vv_SHADE500[i-1]]
#         sat = saturation_pressure(cyclohexanol_SHADE500, T_SHADE500[i]; v0=v0)
#         psat_SHADE500[i] = sat[1]
#         vl_SHADE500[i] = sat[2]
#         vv_SHADE500[i] = sat[3]
#     end
#     hL_SHADE500[i]  = Clapeyron.VT_enthalpy(cyclohexanol_SHADE500,vl_SHADE500[i],T_SHADE500[i],[1.])
#     hV_SHADE500[i]  = Clapeyron.VT_enthalpy(cyclohexanol_SHADE500,vv_SHADE500[i],T_SHADE500[i],[1.])
#     cpL_SHADE500[i] = Clapeyron.VT_isobaric_heat_capacity(cyclohexanol_SHADE500,vl_SHADE500[i],T_SHADE500[i],[1.])
#     cpV_SHADE500[i] = Clapeyron.VT_isobaric_heat_capacity(cyclohexanol_SHADE500,vv_SHADE500[i],T_SHADE500[i],[1.])
# end
# N    = 100

# T_DE500    = LinRange(300, Tc_DE500,  N)
# psat_DE500 = zeros(N)
# vl_DE500   = zeros(N)
# vv_DE500   = zeros(N)

# hL_DE500   = zeros(N)
# hV_DE500   = zeros(N)
# cpL_DE500  = zeros(N)
# cpV_DE500  = zeros(N)

# for i in 1:N
#     if i==1
#         sat = saturation_pressure(cyclohexanol_DE500, T_DE500[i])
#         psat_DE500[i] = sat[1]
#         vl_DE500[i] = sat[2]
#         vv_DE500[i] = sat[3]
#     else
#         v0 = [vl_DE500[i-1],vv_DE500[i-1]]
#         sat = saturation_pressure(cyclohexanol_DE500, T_DE500[i]; v0=v0)
#         psat_DE500[i] = sat[1]
#         vl_DE500[i] = sat[2]
#         vv_DE500[i] = sat[3]
#     end
#     hL_DE500[i]  = Clapeyron.VT_enthalpy(cyclohexanol_DE500,vl_DE500[i],T_DE500[i],[1.])
#     hV_DE500[i]  = Clapeyron.VT_enthalpy(cyclohexanol_DE500,vv_DE500[i],T_DE500[i],[1.])
#     cpL_DE500[i] = Clapeyron.VT_isobaric_heat_capacity(cyclohexanol_DE500,vl_DE500[i],T_DE500[i],[1.])
#     cpV_DE500[i] = Clapeyron.VT_isobaric_heat_capacity(cyclohexanol_DE500,vv_DE500[i],T_DE500[i],[1.])
# end

plt.clf()
plt.plot(T_B, log10.(psat_B), color="blue", label = "Bernet")
plt.plot([Tc_B], log10.([pc_B]), marker="o",color="blue")
plt.plot(T_SHADE100, log10.(psat_SHADE100), color="red", label = "SHADE100")
plt.plot([Tc_SHADE100],log10.([pc_SHADE100]), marker="o",color="red")
plt.plot(T_ECA1000, log10.(psat_ECA1000), color="green", label = "ECA1000")
plt.plot([Tc_ECA1000],log10.([pc_ECA1000]), marker="o",color="green")
plt.plot(T_ECA100, log10.(psat_ECA100), color="purple", label = "ECA100")
plt.plot([Tc_ECA100],log10.([pc_ECA100]), marker="o",color="purple")
# plt.plot(T_SHADE500, log10.(psat_SHADE500), color="orange", label = "SHADE500")
# plt.plot([Tc_SHADE500],log10.([pc_SHADE500]), marker="o",color="orange")
# plt.plot(T_DE500, log10.(psat_DE500), color="orange", label = "DE500")
# plt.plot([Tc_DE500],log10.([pc_DE500]), marker="o",color="orange")
plt.plot(P_exp.T, log10.(P_exp.out_p), marker="o", linestyle="none", color="black", label="Experimental")
plt.plot([Tc],log10.([pc]), marker="*",color="black")
plt.legend(loc="lower right",frameon=false,fontsize=11)
plt.xlabel("Temperature / K",fontsize=16)
plt.ylabel("log(Pressure / Pa)",fontsize=16)
plt.xlim([300,700])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

plt.clf()
plt.plot(1e-3 ./vl_B, T_B, color="blue", label = "Bernet")
plt.plot(1e-3 ./vv_B, T_B, color="blue")
plt.plot([1e-3/vc_B],[Tc_B], marker="o",color="blue")
plt.plot(1e-3 ./vl_SHADE100, T_SHADE100, color="red", label = "SHADE100")
plt.plot(1e-3 ./vv_SHADE100, T_SHADE100, color="red")
plt.plot([1e-3/vc_SHADE100],[Tc_SHADE100], marker="o",color="red")
plt.plot(1e-3 ./vl_ECA1000, T_ECA1000, color="green", label = "ECA1000")
plt.plot(1e-3 ./vv_ECA1000, T_ECA1000, color="green")
plt.plot([1e-3/vc_ECA1000],[Tc_ECA1000], marker="o",color="green")
plt.plot(1e-3 ./vl_ECA100, T_ECA100, color="purple", label = "ECA100")
plt.plot(1e-3 ./vv_ECA100, T_ECA100, color="purple")
plt.plot([1e-3/vc_ECA100],[Tc_ECA100], marker="o",color="purple")
# plt.plot(1e-3 ./vl_SHADE500, T_SHADE500, color="orange", label = "SHADE500")
# plt.plot(1e-3 ./vv_SHADE500, T_SHADE500, color="orange")
# plt.plot([1e-3/vc_SHADE500],[Tc_SHADE500], marker="o",color="orange")
# plt.plot(1e-3 ./vl_DE500, T_DE500, color="orange", label = "DE500")
# plt.plot(1e-3 ./vv_DE500, T_DE500, color="orange")
# plt.plot([1e-3/vc_DE500],[Tc_DE500], marker="o",color="orange")
plt.plot(rhov_exp.out_rhov*1e-3, rhov_exp.T, marker="o", linestyle="none", color="black", label="Experimental")
plt.plot(rhol_exp.out_rhol*1e-3, rhol_exp.T, marker="o", linestyle="none", color="black")
plt.legend(loc="upper right",frameon=false,fontsize=11)
plt.ylabel("Temperature / K",fontsize=16)
plt.xlabel("Density / (mol/dm³)",fontsize=16)
plt.ylim([300,690])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

plt.clf()
plt.plot(T_B, (hvap_B)./1e3, color="blue", label = "Bernet")
plt.plot([Tc_B],[0], marker="o",color="blue")
plt.plot(T_SHADE100, (hV_SHADE100.-hL_SHADE100)./1e3, color="red", label = "SHADE100")
plt.plot([Tc_SHADE100],[0], marker="o",color="red")
plt.plot(T_ECA1000, (hV_ECA1000.-hL_ECA1000)./1e3, color="green", label = "ECA1000")
plt.plot([Tc_ECA1000],[0], marker="o",color="green")
plt.plot(T_ECA100, (hV_ECA100.-hL_ECA100)./1e3, color="purple", label = "ECA100")
plt.plot([Tc_ECA100],[0], marker="o",color="purple")
# plt.plot(T_SHADE500, (hV_SHADE500.-hL_SHADE500)./1e3, color="orange", label = "SHADE500")
# plt.plot([Tc_SHADE500],[0], marker="o",color="orange")
# plt.plot(T_DE500, (hV_DE500.-hL_DE500)./1e3, color="orange", label = "DE500")
# plt.plot([Tc_DE500],[0], marker="o",color="orange")
plt.plot(H_exp.T, H_exp.out_en_vap/1e3, marker="o", linestyle="none", color="black", label="Experimental")
plt.legend(loc="upper right",frameon=false,fontsize=11)
plt.xlabel("Temperature / K",fontsize=16)
plt.ylabel("Enthalpy of Vapourisation / (kJ/mol)",fontsize=16)
#plt.xlim([300,690])
#plt.ylim([0,60])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

N    = 201
T    = 423.56

x    = LinRange(0., 1.,  N)
p_B    = zeros(N)
y_B    = zeros(N)

v0 = nothing
for i in 1:N
    bub = bubble_pressure(binary_B, T, [x[i], 1-x[i]])
    p_B[i] = bub[1]
    y_B[i] = bub[4][1]
end

p_SHADE100    = zeros(N)
y_SHADE100    = zeros(N)

v0 = nothing
for i in 1:N
    bub = bubble_pressure(binary_SHADE100, T, [x[i], 1-x[i]])
    p_SHADE100[i] = bub[1]
    y_SHADE100[i] = bub[4][1]
end

p_ECA1000    = zeros(N)
y_ECA1000    = zeros(N)

v0 = nothing
for i in 1:N
    bub = bubble_pressure(binary_ECA1000, T, [x[i], 1-x[i]])
    p_ECA1000[i] = bub[1]
    y_ECA1000[i] = bub[4][1]
end

p_ECA100    = zeros(N)
y_ECA100    = zeros(N)

v0 = nothing
for i in 1:N
    bub = bubble_pressure(binary_ECA100, T, [x[i], 1-x[i]])
    p_ECA100[i] = bub[1]
    y_ECA100[i] = bub[4][1]
end

# p_SHADE500    = zeros(N)
# y_SHADE500    = zeros(N)

# v0 = nothing
# for i in 1:N
#     bub = bubble_pressure(binary_SHADE500, T, [x[i], 1-x[i]])
#     p_SHADE500[i] = bub[1]
#     y_SHADE500[i] = bub[4][1]
# end
# p_DE500    = zeros(N)
# y_DE500    = zeros(N)

# v0 = nothing
# for i in 1:N
#     bub = bubble_pressure(binary_DE500, T, [x[i], 1-x[i]])
#     p_DE500[i] = bub[1]
#     y_DE500[i] = bub[4][1]
# end

plt.clf()
plt.plot(x, p_B./1e3, color="blue", label = "Bernet")
plt.plot(y_B, p_B./1e3, color="blue")
plt.plot(x, p_SHADE100./1e3, color="red", label = "SHADE100")
plt.plot(y_SHADE100, p_SHADE100./1e3, color="red")
plt.plot(x, p_ECA1000./1e3, color="green", label = "ECA1000")
plt.plot(y_ECA1000, p_ECA1000./1e3, color="green")
plt.plot(x, p_ECA100./1e3, color="purple", label = "ECA100")
plt.plot(y_ECA100, p_ECA100./1e3, color="purple")
# plt.plot(x, p_SHADE500./1e3, color="purple", label = "SHADE500")
# plt.plot(y_SHADE500, p_SHADE500./1e3, color="purple")
# plt.plot(x, p_DE500./1e3, color="orange", label = "DE500")
# plt.plot(y_DE500, p_DE500./1e3, color="orange")
plt.plot(Pxy_exp423.x, Pxy_exp423.out_p/1e3, marker="o", linestyle="none", color="black", label="Experimental")
plt.plot(Pxy_exp423.out_y, Pxy_exp423.out_p/1e3, marker="o", linestyle="none", color="black")
plt.legend(loc="upper right",frameon=true,fontsize=11)
plt.xlabel("composition / (mol/mol)",fontsize=16)
plt.ylabel("Pressure / kPa",fontsize=16)
plt.xlim([0,1])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

# N    = 201
# p    = 100e3

# x    = LinRange(0., 1.,  N)
# T_B    = zeros(N)
# y_B    = zeros(N)

# v0 = nothing
# for i in 1:N
#     bub = bubble_temperature(binary_B, p, [x[i], 1-x[i]])
#     T_B[i] = bub[1]
#     y_B[i] = bub[4][1]
# end


# plt.clf()
# plt.plot(x, T_B, color="blue", label = "Bernet")
# plt.plot(y_B, T_B, color="blue")
# plt.plot(x, p_SHADE100./1e3, color="red", label = "SHADE100")
# plt.plot(y_SHADE100, p_SHADE100./1e3, color="red")
# plt.plot(x, p_ECA1000./1e3, color="green", label = "ECA1000")
# plt.plot(y_ECA1000, p_ECA1000./1e3, color="green")
# plt.plot(x, p_SHADE500./1e3, color="purple", label = "SHADE500")
# plt.plot(y_SHADE500, p_SHADE500./1e3, color="purple")
# plt.plot(x, p_DE500./1e3, color="orange", label = "DE500")
# plt.plot(y_DE500, p_DE500./1e3, color="orange")
# plt.plot(Txy_exp.x, Txy_exp.out_T, marker="o", linestyle="none", color="black", label="Experimental")
# plt.plot(Txy_exp.out_y, Txy_exp.out_T, marker="o", linestyle="none", color="black")
# plt.legend(loc="upper right",frameon=true,fontsize=11)
# plt.xlabel("composition / (mol/mol)",fontsize=16)
# plt.ylabel("Temperature / K",fontsize=16)
# plt.xlim([0,1])
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.show()

p = 101.325e3     #Pa
T = 293.15        #K
rho_B = zeros(N)

for i in 1:N
    rho_B[i] = mass_density(binary_B, p, T, [x[i], 1-x[i]])
end

rho_SHADE100 = zeros(N)

for i in 1:N
    rho_SHADE100[i] = mass_density(binary_SHADE100, p, T, [x[i], 1-x[i]])
end

rho_ECA1000 = zeros(N)

for i in 1:N
    rho_ECA1000[i] = mass_density(binary_ECA1000, p, T, [x[i], 1-x[i]])
end

rho_ECA100 = zeros(N)

for i in 1:N
    rho_ECA100[i] = mass_density(binary_ECA100, p, T, [x[i], 1-x[i]])
end

# rho_SHADE500 = zeros(N)

# for i in 1:N
#     rho_SHADE500[i] = mass_density(binary_SHADE500, p, T, [x[i], 1-x[i]])
# end

# rho_DE500 = zeros(N)

# for i in 1:N
#     rho_DE500[i] = mass_density(binary_DE500, p, T, [x[i], 1-x[i]])
# end

plt.clf()
plt.plot(x, rho_B, color="blue", label = "Bernet")
plt.plot(x, rho_SHADE100, color="red", label = "SHADE100")
plt.plot(x, rho_ECA1000, color="green", label = "ECA1000")
plt.plot(x, rho_ECA100, color="purple", label = "ECA100")
# plt.plot(x, rho_SHADE500, color="orange", label = "SHADE500")
# plt.plot(x, rho_DE500, color="orange", label = "DE500")
plt.plot(rhobinary_exp293.x, rhobinary_exp293.out_mass_density, marker="o", linestyle="none", color="black", label="Experimental")
plt.legend(loc="upper left",frameon=false,fontsize=11)
plt.xlabel("composition / (mol/mol)",fontsize=16)
plt.ylabel("ρ / (kg/m³)",fontsize=16)
plt.xlim([0,1])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

p = 101.325e3     #Pa
T = 298.15        #K
u_B = zeros(N)

for i in 1:N
    u_B[i] = speed_of_sound(binary_B, p, T, [x[i], 1-x[i]])
end

u_SHADE100 = zeros(N)

for i in 1:N
    u_SHADE100[i] = speed_of_sound(binary_SHADE100, p, T, [x[i], 1-x[i]])
end

u_ECA1000 = zeros(N)

for i in 1:N
    u_ECA1000[i] = speed_of_sound(binary_ECA1000, p, T, [x[i], 1-x[i]])
end

u_ECA100 = zeros(N)

for i in 1:N
    u_ECA100[i] = speed_of_sound(binary_ECA100, p, T, [x[i], 1-x[i]])
end

# u_SHADE500 = zeros(N)

# for i in 1:N
#     u_SHADE500[i] = speed_of_sound(binary_SHADE500, p, T, [x[i], 1-x[i]])
# end
u_DE500 = zeros(N)

# for i in 1:N
#     u_DE500[i] = speed_of_sound(binary_DE500, p, T, [x[i], 1-x[i]])
# end

plt.clf()
plt.plot(x, u_B, color="blue", label = "Bernet")
plt.plot(x, u_SHADE100, color="red", label = "SHADE100")
plt.plot(x, u_ECA1000, color="green", label = "ECA1000")
plt.plot(x, u_ECA100, color="purple", label = "ECA100")
# plt.plot(x, u_SHADE500, color="orange", label = "SHADE500")
# plt.plot(x, u_DE500, color="orange", label = "DE500")
plt.plot(ubinary_exp298.x, ubinary_exp298.out_u, marker="o", linestyle="none", color="black", label="Experimental")
plt.legend(loc="upper left",frameon=false,fontsize=11)
plt.xlabel("composition / (mol/mol)",fontsize=16)
plt.ylabel("u / (m/s)",fontsize=16)
plt.xlim([0,1])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

p = 101.325e3     #Pa
T = 298.15        #K
Hexcess_B = zeros(N)

for i in 1:N
    Hexcess_B[i] = excess(binary_B, p, T, [x[i], 1-x[i]], enthalpy)
end

Hexcess_SHADE100 = zeros(N)

for i in 1:N
    Hexcess_SHADE100[i] = excess(binary_SHADE100, p, T, [x[i], 1-x[i]], enthalpy)
end

Hexcess_ECA1000 = zeros(N)

for i in 1:N
    Hexcess_ECA1000[i] = excess(binary_ECA1000, p, T, [x[i], 1-x[i]], enthalpy)
end

Hexcess_ECA100 = zeros(N)

for i in 1:N
    Hexcess_ECA100[i] = excess(binary_ECA100, p, T, [x[i], 1-x[i]], enthalpy)
end

# Hexcess_SHADE500 = zeros(N)

# for i in 1:N
#     Hexcess_SHADE500[i] = excess(binary_SHADE500, p, T, [x[i], 1-x[i]], enthalpy)
# end

# Hexcess_DE500 = zeros(N)

# for i in 1:N
#     Hexcess_DE500[i] = excess(binary_DE500, p, T, [x[i], 1-x[i]], enthalpy)
# end

plt.clf()
plt.plot(x, Hexcess_B, color="blue", label = "Bernet")
plt.plot(x, Hexcess_SHADE100, color="red", label = "SHADE100")
plt.plot(x, Hexcess_ECA1000, color="green", label = "ECA1000")
plt.plot(x, Hexcess_ECA100, color="purple", label = "ECA100")
# plt.plot(x, Hexcess_SHADE500, color="orange", label = "SHADE500")
# plt.plot(x, Hexcess_DE500, color="orange", label = "DE500")
plt.plot(Hbinary_exp.x, Hbinary_exp.out_he, marker="o", linestyle="none", color="black", label="Experimental")
plt.legend(loc="upper right",frameon=false,fontsize=11)
plt.xlabel("composition / (mol/mol)",fontsize=16)
plt.ylabel("Excess enthalpy / (J/mol)",fontsize=16)
plt.xlim([0,1])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()