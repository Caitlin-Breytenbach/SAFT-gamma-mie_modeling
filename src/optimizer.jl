cd(@__DIR__)
include("../src/ThermoProps.jl")

using NLopt, Clapeyron, FiniteDiff, SqpSolver, JuMP, Ipopt
export optimizer

function optimizer(estimator, objective, initial, upper, lower, starts)
    iter = Ref(0)
    function objective_fn(x_norm::Vector, grad::Vector)

        x = @. x_norm*(upper - lower) + lower
        val = objective(x)

        iter[] += 1 
        if iter[] % 10 == 0 || iter[] == 1
            println("Iteration: $(iter[])")
            println("objective: $(round(val; digits = 6))")
            println("parameters: $(round.(x; sigdigits = 5))")
        end

        if length(grad)> 0
            f = x_norm -> objective(@. x_norm*(upper-lower) + lower)
            FiniteDiff.finite_difference_gradient!(grad, f, x_norm)
            # println("  Gradient      : $(round.(grad; sigdigits=3))")
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

    min_x = @. min_xnorm*(upper - lower) + lower

    for i in 1:starts
        x0 = @. rand()*(upper - lower) + lower
        println("Restarting optimization with $(x0)")

        min_f_random, min_xnorm_random, ret_random = run_opt(x0)

        println(
            """
            Iterations            : $(iter[])
            objective value       : $min_f_random
            solution status       : $ret_random
            """
        )

        if min_f_random < min_f
            min_f = min_f_random
            min_xnorm = min_xnorm_random
            min_x = @. min_xnorm*(upper - lower) + lower
            println("New best solution found with objective value: $min_f")
        end
    end

    println("Optimized parameters: $(round.(min_x; sigdigits=5))")
    model_opt = Clapeyron.return_model(estimator, model, min_x)
    return model_opt
end
