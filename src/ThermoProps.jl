using Clapeyron
 
export saturation_p, saturation_rhol, saturation_rhov, en_vap
export sat_envelope
export rhol, rhov_compressed
export bubble_p, dew_p, bubble_t, dew_t
export binary_density, binary_he, binary_u
export calc_ADD, calc_percentADD
export load_model

function load_model(components; userlocations::Vector{String}=String[])
    if isempty(userlocations)
        return SAFTgammaMie(components)
    else
        return SAFTgammaMie(components; userlocations=userlocations)
    end
end

function saturation_p(model::EoSModel, T::Real)
    try
        return saturation_pressure(model, T)[1]
    catch
        return NaN
    end
end

function sat_p_plot(model, exp_data, crits)
    T_min = minimum(exp_data.T) 
    Tc = crits[1]
    sat_p = zeros(N)
    T = LinRange(T_min, Tc, 100)
    N = length(T)
    for i in 1:N
        sat_p[i] = saturation_p(model, T[i])
    end
    return T, sat_p
end

function saturation_rhol(model::EoSModel, T::Real)
    try
        return 1.0 / saturation_pressure(model, T)[2]
    catch
        return NaN
    end
end

function saturation_rhov(model::EoSModel, T::Real)
    try
        return 1.0 / saturation_pressure(model, T)[3]
    catch
        return NaN
    end
end

function sat_rho_plot(model, exp_data, crits)
    T_min = minimum(exp_data.T) 
    Tc = crits[1]
    sat_rho = zeros(N)
    T = LinRange(T_min, Tc, 100)
    N = length(T)
    for i in 1:N
        sat_rhol[i] = saturation_rhol(model, T[i])
        sat_rhov[i] = saturation_rhov(model, T[i])
    end
    return T, sat_rhol, sat_rhov
end

function en_vap(model::EoSModel, T::Real)
    try
        return enthalpy_vap(model, T)
    catch
        return NaN
    end
end

function sat_envelope(model::EoSModel, T_range::AbstractVector)
    N    = length(T_range)
    psat = zeros(N)
    vl   = zeros(N)
    vv   = zeros(N)
    hvap = zeros(N)
 
    for i in 1:N
        try
            sat = (i == 1) ?
                saturation_pressure(model, T_range[i]) :
                saturation_pressure(model, T_range[i]; v0=[vl[i-1], vv[i-1]])
            psat[i] = sat[1]
            vl[i]   = sat[2]
            vv[i]   = sat[3]
        catch
            psat[i] = NaN
            vl[i]   = (i > 1 && !isnan(vl[i-1])) ? vl[i-1] : NaN
            vv[i]   = (i > 1 && !isnan(vv[i-1])) ? vv[i-1] : NaN
        end
        hvap[i] = en_vap(model, T_range[i])
    end
 
    rhol = 1.0 ./ vl
    rhov = 1.0 ./ vv
 
    return (T=T_range, p=psat, vl=vl, vv=vv, Δh=hvap)
end

function rhol(model::EoSModel, T::Real, p::Real)
    try
        return molar_density(model, p, T; phase=:liquid)
    catch
        return NaN
    end
end

function rhol_curve(model, exp_data)
    T_vals   = Float64[]
    p_vals   = Float64[]
    rho_vals = Float64[]
    for row in eachrow(exp_data)
        push!(T_vals,   row.T)
        push!(p_vals,   row.p)
        push!(rho_vals, rhol(model, row.T, row.p))
    end
    return (T_vals=T_vals, p_vals=p_vals, rho_vals=rho_vals)
end

function bubble_p(model::EoSModel, x::Real, T::Real)
    try
        bub = bubble_pressure(model, T, [x, 1 - x])
        return bub[1], bub[4][1]
    catch
        return NaN, NaN
    end
end

function dew_p(model::EoSModel, y::Real, T::Real)
    try
        dew = dew_pressure(model, T, [y, 1 - y])
        return dew[1], dew[4][1]
    catch
        return NaN, NaN
    end
end

function bubble_t(model::EoSModel, x::Real, p::Real)
    try
        bub = bubble_temperature(model, p, [x, 1 - x])
        return bub[1], bub[4][1]
    catch
        return NaN, NaN
    end
end

function dew_t(model::EoSModel, y::Real, p::Real)
    try
        dew = dew_temperature(model, p, [y, 1 - y])
        return dew[1], dew[4][1]
    catch
        return NaN, NaN
    end
end

function binary_density(model::EoSModel, x::Real, p::Real, T::Real)
    try
        return mass_density(model, p, T, [x, 1 - x])
    catch
        return NaN
    end
end

function binary_he(model::EoSModel, x::Real, p::Real, T::Real)
    try
        return excess(model, p, T, [x, 1 - x], enthalpy)
    catch
        return NaN
    end
end

function binary_u(model::EoSModel, x::Real, p::Real, T::Real)
    try
        return speed_of_sound(model, p, T, [x, 1 - x])
    catch
        return NaN
    end
end

