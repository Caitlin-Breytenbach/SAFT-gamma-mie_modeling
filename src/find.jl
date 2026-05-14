using GCIdentifier, ChemicalIdentifiers, Clapeyron

function groups_from_names(inputs::Vector{String})
    comps = []

    for i in eachindex(inputs)
        input = inputs[i]

        try
            comp = get_groups_from_name(input, SAFTgammaMieGroups)
            push!(comps, comp)
        catch
            println("Use SMILES for '$input' ")
        end
    end
    println(comps)
    return comps
end

function groups_from_smiles(inputs::Vector{String})
    comps = []

    for i in eachindex(inputs)
        input = inputs[i]

        try
            comp = get_groups_from_smiles(input, SAFTgammaMieGroups)
            push!(comps, comp)
        catch
            println("Groups missing for '$input' ")
        end
    end
    println(comps)
    return comps
end
