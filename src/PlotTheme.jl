module PlotTheme
 
using CairoMakie, Colors
 
export apply_theme!, PALETTE, line_colour, EXP_COLOR, LINEWIDTH, MARKERSIZE, FONTSIZE, MARKERS
export plot_saturation_pressure, plot_VLE_envelope, plot_enthalpy_vap, plot_rhol, plot_Cp, plot_Cv_isobaric, plot_Cv_isothermal, plot_u
export plot_Pxy, plot_binary_property
export exp_scatter!, model_line!, crit_point!

const PALETTE = [
    colorant"#4878CF",   # blue
    colorant"#D65F5F",   # red
    colorant"#6ACC65",   # green
    colorant"#B47CC7",   # purple
    colorant"#F0A30A",   # amber  
    colorant"#3BBDBD",   # teal   
]
 
const EXP_COLOR  = colorant"#222222"
const LINEWIDTH  = 2.0
const MARKERSIZE = 10.0
const FONTSIZE   = 14
const MARKERS = [:circle, :rect, :diamond, :utriangle, :dtriangle, :pentagon, :hexagon, :star5]

function line_colour(names::Vector{String})
    Dict(name => PALETTE[mod1(i, length(PALETTE))] for (i, name) in enumerate(names))    
end

function apply_theme!()
    update_theme!(
        fonts    = Attributes(
            regular = "TeX Gyre Heros",
            bold    = "TeX Gyre Heros Bold",
        ),
        fontsize = FONTSIZE,
        Axis = (
            spinewidth         = 1.2,
            xgridvisible       = false,
            ygridvisible       = false,
            xminorticksvisible = true,
            yminorticksvisible = true,
            xminorticks        = IntervalsBetween(5),
            yminorticks        = IntervalsBetween(5),
            xminorticksize     = 3,
            yminorticksize     = 3,
            xticksize          = 5,
            yticksize          = 5,
            xtickwidth         = 1.2,
            ytickwidth         = 1.2,
            rightspinevisible  = false,
            topspinevisible    = false,
            xlabelsize         = FONTSIZE,
            ylabelsize         = FONTSIZE,
            xticklabelsize     = FONTSIZE - 2,
            yticklabelsize     = FONTSIZE - 2,
        ),
        Legend = (
            framevisible = false,
            labelsize    = FONTSIZE - 1,
            patchsize    = (20, 2),
            rowgap       = 2,
        ),
    )
end

function model_line!(ax, x, y; name::String, colors::Dict, kwargs...)
    lines!(ax, x, y;
        color     = get(colors, name, PALETTE[1]),
        linewidth = LINEWIDTH,
        label     = name,
        kwargs...)
end

function exp_scatter!(ax, x, y; label="Experiment", kwargs...)
    scatter!(ax, x, y;
        color       = EXP_COLOR,
        marker      = :circle,
        markersize  = MARKERSIZE,
        strokewidth = 1.,
        strokecolor = :black,
        label       = label,
        kwargs...)
end

function crit_point!(ax, x, y; name::String, colors::Dict, marker=:circle)
    scatter!(ax, [x], [y];
        color       = get(colors, name, PALETTE[1]),
        marker      = marker,
        markersize  = MARKERSIZE + 3,
        strokewidth = 1.0,
        strokecolor = :white)
end

function plot_saturation_pressure(
    model_curves::Dict; model_crits = nothing,
    exp_T, exp_p,
    exp_crit = nothing,
    Tlims    = (nothing, nothing),
    size     = (500, 420),
)
    names  = collect(keys(model_curves))
    colors = line_colour(names)
 
    fig = Figure(; size)
    ax  = Axis(fig[1,1];
        xlabel = "Temperature / K",
        ylabel = "log₁₀(Pressure / Pa)",
        limits = (Tlims, nothing),
    )
 
    for name in names
        c      = model_curves[name]
        if !isnothing(model_crits)
        Tc, pc = model_crits[name]
        crit_point!(ax, Tc, log10(pc); name, colors)            
        end
        model_line!(ax, c.T, log10.(c.p); name, colors)
    end
 
    exp_scatter!(ax, exp_T, log10.(exp_p))
    if !isnothing(exp_crit)
        scatter!(ax, [exp_crit[1]], [log10(exp_crit[2])];
            color=EXP_COLOR, marker=:star5, markersize=MARKERSIZE+4,
            strokewidth=0.5, strokecolor=:white)
    end
 
    Legend(fig[1,2], ax)
    return fig
end

function plot_VLE_envelope(
    model_curves::Dict, model_crits::Dict;
    exp_rhol_T, exp_rhol,
    exp_rhov_T = Float64[],
    exp_rhov = Float64[],
    Tlims = (nothing, nothing),
    size  = (500, 420),
)
    names  = collect(keys(model_curves))
    colors = line_colour(names)
 
    fig = Figure(; size)
    ax  = Axis(fig[1,1];
        xlabel = "Density / (mol dm⁻³)",
        ylabel = "Temperature / K",
        limits = (nothing, Tlims),
    )
 
    for name in names
        c          = model_curves[name]
        Tc, pc, vc = model_crits[name]
        col        = get(colors, name, PALETTE[1])
        lines!(ax, c.rhol .*1e-3, c.T; color=col, linewidth=LINEWIDTH, label=name)
        lines!(ax, c.rhov .*1e-3, c.T; color=col, linewidth=LINEWIDTH)
        crit_point!(ax, 1e-3/vc, Tc; name, colors)
    end
 
    exp_scatter!(ax, exp_rhol .* 1e-3, exp_rhol_T)
    exp_scatter!(ax, exp_rhov .* 1e-3, exp_rhov_T; label="")   # no duplicate legend entry
 
    Legend(fig[1,2], ax)
    return fig
end

function plot_enthalpy_vap(
    model_curves::Dict, model_crits::Dict;
    exp_T, exp_hvap,
    size = (500, 420),
)
    names  = collect(keys(model_curves))
    colors = line_colour(names)
 
    fig = Figure(; size)
    ax  = Axis(fig[1,1];
        xlabel = "Temperature / K",
        ylabel = "Enthalpy of vaporisation / (kJ mol⁻¹)",
    )
 
    for name in names
        c      = model_curves[name]
        Tc     = first(model_crits[name])
        model_line!(ax, c.T, c.Δh ./ 1e3; name, colors)
        crit_point!(ax, Tc, 0.0; name, colors)
    end
 
    exp_scatter!(ax, exp_T, exp_hvap ./ 1e3)
 
    Legend(fig[1,2], ax)
    return fig
end

function plot_Pxy(
    model_curves::Dict;
    exp_x, exp_p,
    exp_y    = nothing,
    T_label  = nothing,
    size     = (500, 420),
)
    names  = collect(keys(model_curves))
    colors = line_colour(names)
 
    fig = Figure(; size)
    ax  = Axis(fig[1,1];
        xlabel = "Composition x, y (mol mol⁻¹)",
        ylabel = "Pressure / kPa",
        limits = ((0, 1), nothing),
        xticks = 0:0.2:1,
    )
 
    for name in names
        c   = model_curves[name]
        col = get(colors, name, PALETTE[1])
        lines!(ax, c.x, c.p ./ 1e3; color=col, linewidth=LINEWIDTH, label=name)
        if hasproperty(c, :y)
            lines!(ax, c.y, c.p ./ 1e3; color=col, linewidth=LINEWIDTH, linestyle=:dash)
        end
    end
 
    exp_scatter!(ax, exp_x, exp_p ./ 1e3)
    if !isnothing(exp_y)
        exp_scatter!(ax, exp_y, exp_p ./ 1e3; label="")
    end
 
    if !isnothing(T_label)
        text!(ax, 0.02, 0.97; text=T_label, align=(:left,:top),
              space=:relative, fontsize=FONTSIZE-2)
    end
 
    Legend(fig[1,2], ax)
    return fig
end

function plot_binary_property(
    model_curves::Dict;
    exp_x, exp_y,
    xlabel::String,
    ylabel::String,
    annotation       = nothing,
    legend_position  = :rt,
    size             = (500, 420),
)
    names  = collect(keys(model_curves))
    colors = line_colour(names)
 
    fig = Figure(; size)
    ax  = Axis(fig[1,1];
        xlabel = xlabel,
        ylabel = ylabel,
        limits = ((0, 1), nothing),
        xticks = 0:0.2:1,
    )
 
    for name in names
        c = model_curves[name]
        model_line!(ax, c.x, c.y; name, colors)
    end
 
    exp_scatter!(ax, exp_x, exp_y)
 
    if !isnothing(annotation)
        halign, valign, tx, ty = if legend_position == :rt
            (:left, :top, 0.02, 0.97)
        elseif legend_position == :lt
            (:left, :top, 0.02, 0.97)
        elseif legend_position == :rb
            (:left, :bottom, 0.02, 0.03)
        else   # :lb
            (:left, :bottom, 0.02, 0.03)
        end
        text!(ax, tx, ty; text=annotation, align=(halign, valign),
              space=:relative, fontsize=FONTSIZE-2)
    end
 
    Legend(fig[1,2], ax)
    return fig
end


function plot_rhol(
    model_curves::Dict;
    exp_T, exp_p, exp_rho,
    size  = (500, 420),
)
    names  = collect(keys(model_curves))
    colors = line_colour(names)
 
    fig = Figure(; size)
    ax  = Axis(fig[1,1];
        xlabel = "Temperature / K",
        ylabel = "Density / (mol m⁻³)",
    )

    all_p = sort(unique(exp_p))
    n_p = length(all_p)
    temp_marker = Dict(p => MARKERS[mod1(k, length(MARKERS))] for (k, p) in enumerate(all_p))

    for p_val in all_p
        mask = exp_p .== p_val
        scatter!(ax, exp_T[mask], exp_rho[mask];
            color = EXP_COLOR,
            marker = temp_marker[p_val],
            markersize = MARKERSIZE,
            strokewidth = 1.0,
            strokecolor = :black,
            label = "$(round(Int, p_val)/1e6) MPa",
            )
        
    end
 
    for name in names
        c          = model_curves[name]
        col        = get(colors, name, PALETTE[1])
        for p_val in all_p
            mask = c.p_vals .== p_val
            T_sorted = c.T_vals[mask]
            rho_sorted = c.rho_vals[mask]
            idx = sortperm(T_sorted)
            lines!(ax, T_sorted[idx], rho_sorted[idx];
                color = col,
                linewidth = LINEWIDTH,
                label = name,
                ) 
        end
    end
 
    Legend(fig[1,2], ax; merge=true, unique=true)
    return fig
end

function plot_Cp(
    model_curves::Dict;
    exp_T, exp_p, exp_Cp,
    size  = (500, 420),
)
    names  = collect(keys(model_curves))
    colors = line_colour(names)
 
    fig = Figure(; size)
    ax  = Axis(fig[1,1];
        xlabel = "Temperature / K",
        ylabel = "Cp / (J mol⁻¹ K⁻¹)",
    )

    all_p = sort(unique(exp_p))
    n_p = length(all_p)
    press_marker = Dict(P => MARKERS[mod1(k, length(MARKERS))] for (k, P) in enumerate(all_p))

    for p_val in all_p
        mask = exp_p .== p_val
        scatter!(ax, exp_T[mask], exp_Cp[mask];
            color = EXP_COLOR,
            marker = press_marker[p_val],
            markersize = MARKERSIZE,
            strokewidth = 1.0,
            strokecolor = :black,
            label = "$(round(Int, p_val/1e3)) kPa",
            )
        
    end
 
    for name in names
        c          = model_curves[name]
        col        = get(colors, name, PALETTE[1])
        for p_val in all_p
            mask = c.p_vals .== p_val
            T_sorted = c.T_vals[mask]
            Cp_sorted = c.Cp_vals[mask]
            idx = sortperm(T_sorted)
            lines!(ax, T_sorted[idx], Cp_sorted[idx];
                color = col,
                linewidth = LINEWIDTH,
                label = name,
                ) 
        end
    end
 
    Legend(fig[1,2], ax; merge=true, unique=true)
    return fig
end

function plot_Cv_isobaric(
    model_curves::Dict;
    exp_T, exp_p, exp_Cv, 
    size  = (500, 420),
)
    names  = collect(keys(model_curves))
    colors = line_colour(names)
 
    fig = Figure(; size)
    ax  = Axis(fig[1,1];
        xlabel = "Temperature / K",
        ylabel = "Cv / (J mol⁻¹ K⁻¹)",
    )

    all_p = sort(unique(exp_p))
    n_p = length(all_p)
    press_marker = Dict(P => MARKERS[mod1(k, length(MARKERS))] for (k, P) in enumerate(all_p))

    for p_val in all_p
        mask = exp_p .== p_val
        scatter!(ax, exp_T[mask], exp_Cv[mask];
            color = EXP_COLOR,
            marker = press_marker[p_val],
            markersize = MARKERSIZE,
            strokewidth = 1.0,
            strokecolor = :black,
            label = "$(round(Int, p_val/1e3)) kPa",
            )
        
    end
 
    for name in names
        c          = model_curves[name]
        col        = get(colors, name, PALETTE[1])
        for p_val in all_p
            mask = c.p_vals .== p_val
            T_sorted = c.T_vals[mask]
            Cv_sorted = c.Cv_vals[mask]
            idx = sortperm(T_sorted)
            lines!(ax, T_sorted[idx], Cv_sorted[idx];
                color = col,
                linewidth = LINEWIDTH,
                label = name,
                ) 
        end
    end
 
    Legend(fig[1,2], ax; merge=true, unique=true)
    return fig
end

function plot_Cv_isothermal(
    model_curves::Dict;
    exp_T, exp_p, exp_Cv,
    size  = (500, 420),
)
    names  = collect(keys(model_curves))
    colors = line_colour(names)
 
    fig = Figure(; size)
    ax  = Axis(fig[1,1];
        xlabel = "Pressure / MPa",
        ylabel = "Cv / (J mol⁻¹ K⁻¹)",
    )

    all_T = sort(unique(exp_T))
    n_T = length(all_T)
    press_marker = Dict(T => MARKERS[mod1(k, length(MARKERS))] for (k, T) in enumerate(all_T))

    for T_val in all_T
        mask = exp_T .== T_val
        scatter!(ax, exp_p[mask], exp_Cv[mask];
            color = EXP_COLOR,
            marker = press_marker[T_val],
            markersize = MARKERSIZE,
            strokewidth = 1.0,
            strokecolor = :black,
            label = "$(round(Int, T_val)) K",
            )
        
    end
 
    for name in names
        c          = model_curves[name]
        col        = get(colors, name, PALETTE[1])
        for T_val in all_T
            mask = c.T_vals .== T_val
            p_sorted = c.p_vals[mask]
            Cv_sorted = c.Cv_vals[mask]
            idx = sortperm(p_sorted)
            lines!(ax, p_sorted[idx], Cv_sorted[idx];
                color = col,
                linewidth = LINEWIDTH,
                label = name,
                ) 
        end
    end
 
    Legend(fig[1,2], ax; merge=true, unique=true)
    return fig
end

function plot_u(
    model_curves::Dict;
    exp_T, exp_p, exp_u,
    size  = (500, 420),
)
    names  = collect(keys(model_curves))
    colors = line_colour(names)
 
    fig = Figure(; size)
    ax  = Axis(fig[1,1];
        xlabel = "Temperature / K",
        ylabel = "u / (m s⁻¹)",
    )

    all_p = sort(unique(exp_p))
    n_p = length(all_p)
    press_marker = Dict(P => MARKERS[mod1(k, length(MARKERS))] for (k, P) in enumerate(all_p))

    for p_val in all_p
        mask = exp_p .== p_val
        scatter!(ax, exp_T[mask], exp_u[mask];
            color = EXP_COLOR,
            marker = press_marker[p_val],
            markersize = MARKERSIZE,
            strokewidth = 1.0,
            strokecolor = :black,
            label = "$(round(Int, p_val/1e3)) kPa",
            )
        
    end
 
    for name in names
        c          = model_curves[name]
        col        = get(colors, name, PALETTE[1])
        for p_val in all_p
            mask = c.p_vals .== p_val
            T_sorted = c.T_vals[mask]
            u_sorted = c.u_vals[mask]
            idx = sortperm(T_sorted)
            lines!(ax, T_sorted[idx], u_sorted[idx];
                color = col,
                linewidth = LINEWIDTH,
                label = name,
                ) 
        end
    end
 
    Legend(fig[1,2], ax; merge=true, unique=true)
    return fig
end

function plot_pxy(
    model_curves::Dict;
    exp_T, exp_p, exp_x, exp_y,
    size  = (500, 420),
)
    names  = collect(keys(model_curves))
    colors = line_colour(names)
 
    fig = Figure(; size)
    ax  = Axis(fig[1,1];
        xlabel = "x,y",
        ylabel = "Pressure / (MPa)",
    )

    all_T = sort(unique(exp_T))
    n_T = length(all_T)
    press_marker = Dict(T => MARKERS[mod1(k, length(MARKERS))] for (k, T) in enumerate(all_T))

    for T_val in all_T
        mask = exp_T .== T_val
        scatter!(ax, exp_p[mask] ./1e6, exp_x[mask];
            color = EXP_COLOR,
            marker = press_marker[T_val],
            markersize = MARKERSIZE,
            strokewidth = 1.0,
            strokecolor = :black,
            label = "$(round(Int, T_val)) K",
            )
        
            scatter!(ax, exp_p[mask] ./1e6, exp_y[mask];
            color = EXP_COLOR,
            marker = press_marker[T_val],
            markersize = MARKERSIZE,
            strokewidth = 1.0,
            strokecolor = :black,
            label = "$(round(Int, T_val)) K",
            )

    end
 
    for name in names
        c          = model_curves[name]
        col        = get(colors, name, PALETTE[1])
        for T_val in all_T
            mask = c.T_vals .== T_val
            p_sorted = c.p_vals[mask]
            x_sorted = c.x_vals[mask]
            y_sorted = c.y_vals[mask]
            idx = sortperm(p_sorted)
            lines!(ax, p_sorted[idx], x_sorted[idx];
                color = col,
                linewidth = LINEWIDTH,
                label = name,
            ) 
            lines!(ax, p_sorted[idx], y_sorted[idx];
                color = col,
                linewidth = LINEWIDTH,
                label = name,
            ) 
        end
    end
 
    Legend(fig[1,2], ax; merge=true, unique=true)
    return fig
end
end