function plot_glc_limiting_study(simid, Ds, ϵs, cg, is_glclim_stst)


    ## -----------------------------------------------------
    # simulation
    Dyn_ugs, Dyn_Xs = [], []
    net_Uugs, net_UXs = [], []
    ϵs_ = []
    
    @showprogress for (ϵi, ϵ) in enumerate(ϵs)
            
        for D in Ds

            ## ------------------------------------------------------
            # Dynamic
            simparams = (;D, ϵ, cg)
            status = get_status(simid, simparams)
            if (status != STST_SIM_STATUS) 
                @warn("Wrong status ", (D, ϵ, cg), status)
                continue
            end
            datfile = Dyn.simdat_file(simid, simparams)
            if !isfile(datfile) 
                @warn("File missing ", (D, ϵ, cg), datfile)
                continue
            end
            S = ldat(datfile)
            dyn_ug = S.ug_av
            dyn_X = S.X
            
            ## ------------------------------------------------------
            # LP 
            net = S.V.net
            Dyn.fix!(net, :z, D)
            Dyn.ub!(net, :ug, cg * D / S.X)
            Dyn.lb!(net, :ug, 0.0)
            Dyn.reducebox!(net)
            Lug, Uug = Dyn.fva(net, :ug)
            UX = cg * D / Lug
            if !is_glclim_stst[simparams]
                @warn("Not glc-lim stst ", (D, ϵ, cg))
                continue
            end

            push!(Dyn_ugs, dyn_ug)
            push!(Dyn_Xs, dyn_X)
            push!(net_Uugs, Uug)
            push!(net_UXs, UX)

            push!(ϵs_, ϵ)
            
        end

    end

    isempty(ϵs_) && error("Not data recorded!!!")
    max_ϵs_ = maximum(ϵs)
    ms0 = 5

    ps = Plots.Plot[]
    for (p, dyn_vs, net_vs) in [
            (plot(;xlabel = _textbf("Dynamic (u_g)"), ylabel = _textbf("c_{g}D/X")), Dyn_ugs, net_Uugs),
            (plot(;xlabel = _textbf("Dynamic (X)"), ylabel = _textbf("Theo. maximum (X)")), Dyn_Xs, net_UXs)
        ]

        all_ = [dyn_vs; net_vs] |> sort!
        lims = (last(all_) - first(all_))
        lims = (first(all_) - (lims * 0.1), last(all_) + (lims * 0.1))
        plot!(p, all_, all_; 
            ls = :dash, lw = 3, c = :black, alpha = 0.7, label = "", xlim = lims, ylim = lims
        )

        scatter!(p, dyn_vs, net_vs; 
            label = "",
            markerstrokecolor = :white,
            color = :black,
            ms = ms0 .+ (2 * ms0) .* (ϵs_ ./ max_ϵs_),
        )
        push!(ps, p)
    end
    
    return ps

end