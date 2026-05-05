cd(@__DIR__)
include("../src/ThermoProps.jl")

using NLopt, Clapeyron, FiniteDiff
export optimiser

function optimiser(estimator, objective, initial, upper, lower )
    iter = Ref(0)
    best_f = Ref(Inf)
    best_x = copy(initial)
    val = objective(initial)
    x = initial

    function objective_fn(x_norm::Vector, grad::Vector)

        if val < best_f[]
            best_f[] = val 
            best_x .= x 
        end 

        x = @. x_norm*(upper - lower) + lower
        val = objective(x)

        iter[] += 1 
        if iter[] % 10 == 0
            println("Iteration: $(iter[])")
            println("objective: $(round(val; digits = 6))")
            println("parameters: $(round.(x; sigdigits = 5))")
        end

        if length(grad)> 0
            f = x_norm -> objective(@. x_norm*(upper-lower) + lower)
            FiniteDiff.finite_difference_gradient!(grad, f, x_norm)
            println("  Gradient      : $(round.(grad; sigdigits=3))")
        end
        return val
    end

    n = length(initial)

    function run_opt(x0)
        x0_norm = @. (x0 - lower)/(upper - lower)
        opt = NLopt.Opt(:LD_SLSQP, n)
        NLopt.lower_bounds!(opt, zeros(n))
        NLopt.upper_bounds!(opt, ones(n))
        NLopt.xtol_rel!(opt, 1e-8)
        NLopt.min_objective!(opt, objective_fn)
        return NLopt.optimize(opt, x0_norm)
    end

    min_f, min_xnorm, ret = run_opt(initial)
    println(
        """
        Iterations            : $(iter[])
        objective value       : $min_f
        solution status       : $ret
        """
    )

    restart = 1
    improvement = Inf
    while min_f > 1e-2 && improvement > 1e-4
        prev_f = min_f 
        restart += 1
    
        println("Restart from best point")
        min_f, min_xnorm, ret = run_opt(best_x)
        println(
            """
            Run $restart
            Iterations            : $(iter[])
            objective value       : $min_f
            solution status       : $ret
            """
        )
        improvement = prev_f - min_f
         
    end

    if min_f > 0.1 && improvement <= 1e-4
        println("No improvement. Best objective $(best_f[])")
    end

    min_x = @. min_xnorm*(upper - lower) + lower
    model_opt = Clapeyron.return_model(estimator, model, min_x)
    return model_opt
end
