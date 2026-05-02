cd(@__DIR__)
include("../src/ThermoProps.jl")

using NLopt, Clapeyron
export optimiser

function optimiser(estimator, objective, initial, upper, lower )
    x0_norm = @. (initial - lower)/(upper - lower)
    iter = Ref(0)

    function objective_fn(x_norm::Vector, grad::Vector)
        x = @. x_norm*(upper - lower) + lower
        val = objective(x)

        iter[] += 1 
        if iter[] % 10 == 0
            println("Iteration: $(iter[])")
            println("objective: $(round(val; digits = 6))")
            println("parameters: $(round.(x; sigdigits = 5))")
        end
        if length(grad)> 0
            del = 1e-6
            for i in eachindex(x_norm)
                xp = copy(x_norm); xp[i] += del 
                xm = copy(x_norm); xm[i] -= del 
                grad[i] = (objective(@. xp*(upper - lower) + lower) - objective(@. xm*(upper - lower) + lower))/(2*del*(upper[i] - lower[i]))    
            end
        end
        return val
    end
    n = length(initial)
    opt = NLopt.Opt(:LD_SLSQP, n)

    NLopt.lower_bounds!(opt, zeros(n))
    NLopt.upper_bounds!(opt, ones(n))
    NLopt.xtol_rel!(opt, 1e-8)
    NLopt.min_objective!(opt, objective_fn)
    min_f, min_xnorm, ret = NLopt.optimize(opt, x0_norm)
    println(
        """
        objective value       : $min_f
        solution status       : $ret
        """
    )

    min_x = @. min_xnorm*(upper - lower) + lower
    model_opt = Clapeyron.return_model(estimator, model, min_x)
    return model_opt
end
