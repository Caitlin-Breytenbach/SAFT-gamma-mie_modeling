cd(@__DIR__)
include("../src/ThermoProps.jl")
using GCIdentifier, ChemicalIdentifiers, Clapeyron, PyCall, CSV, DataFrames

export AAD


function AAD(calc, exp)
    N = length(exp)
    AAD_values = zeros(N)
    
    for i in 1:N 
        AAD_values[i] = abs(exp[i] - calc[i])        
    end
    AAD = (1/N)*(sum(AAD_values))

    AAD_percent_values = zeros(N)

    for i in 1:N 
        AAD_percent_values[i] = abs((exp[i] - calc[i])/exp[i])*100      
    end
    AAD_percent = (1/N)*(sum(AAD_percent_values))

    return AAD, AAD_percent
end
