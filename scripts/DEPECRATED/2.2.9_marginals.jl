## -------------------------------------------------------------------
# marginal distributions
let 
    objider = iJR.BIOMASS_IDER
    size = [300, 250]
    method2 = ME_MAX_POL

    # Iders
    model_iders, Fd_iders = [objider], ["D"]
    for Fd_met in EXCH_FLX_IDERS
        model_met = Fd_mets_map[Fd_met]
        model_exch = Fd_rxns_map[Fd_met]
        push!(model_iders, model_exch)
        push!(Fd_iders, string("u", Fd_met))
    end
    
    for (model_ider, Fd_ider) in zip(model_iders, Fd_iders)
        ps = Plots.Plot[]
        ps_bs = Plots.Plot[]
        for exp in EXPS
            p = plot(title = string(Fd_ider, " exp: ", exp))
            p_bs = plot(title = string(Fd_ider, " exp: ", exp))
            margin, m, M = -Inf, Inf, -Inf
            Fd_av = Fd.val(Fd_ider, exp)
            
            # EP
            for method in ALL_METHODS
                color = method_colors[method]    

                datfile = INDEX[method, :DFILE, exp]
                dat = deserialize(datfile)
                model = dat[:model]
                objidx = ChU.rxnindex(model, objider)
                epouts = dat[:epouts]
                exp_beta = maximum(keys(epouts))
                epout = epouts[exp_beta]
                ep_av = ChU.av(model, epout, model_ider)
                ep_va = sqrt(ChU.va(model, epout, model_ider))
                        
                ChP.plot_marginal!(p, model, [epout], model_ider; 
                    legend = false, color, alpha = 0.6, lw = 5)
                
                m = minimum([m, ep_av, Fd_av])
                M = maximum([M, ep_av, Fd_av])
                margin = maximum([margin, 3 * ep_va])

                if method == method2
                    for (beta, epout) in sort(epouts; by = first)
                        ep_av = ChU.av(model, epout, model_ider)
                        ep_va = sqrt(ChU.va(model, epout, model_ider))
                        
                        alpha = 0.15
                        color = method_colors[method]
                        ChP.plot_marginal!(p_bs, model, epout, model_ider; 
                            legend = false, color, alpha, lw = 1)

                        if beta == exp_beta
                            ChP.plot_marginal!(p_bs, model, epout, model_ider; 
                                legend = false, color, 
                                alpha = 1.0, lw = 3
                            )
                            break
                        end
                    end
                    push!(ps_bs, p_bs)
                end

            end
            # Experimental
            vline!(p, [Fd_av]; label = "", lw = 6, color = :black, alpha = 0.3)
            vline!(p_bs, [Fd_av]; label = "", lw = 6, color = :black, alpha = 0.3)
            
            plot!(p; xlim = [m - margin, M + margin], size)
            plot!(p_bs; xlim = [m - margin, M + margin], size)
            push!(ps, p)
        end

        for k in [:xi, :D, :uGLC]
            p = plot(;title = Fd_ider, size)
            xticks =  (EXPS, string.(EXPS))
            vals = [Fd.val(k, exp) for exp in EXPS]
            p = bar!(p, EXPS, vals; title = k, label = "", xticks)
            push!(ps, p)
            push!(ps_bs, p)
        end

        # legend
        p = plot(;title = "Legend", size)
        for (method, color) in method_colors
            color = method_colors[method]
            bar!(p, [string(method)], [1]; yaxis = nothing, 
                color, xrotation = 35, label = ""
            )
        end
        push!(ps, p)
        push!(ps_bs, p)

        pname = string(Fd_ider, "_marginals")
        mysavefig(ps, pname)

        method = method2
        pname = string(Fd_ider, "_marginals_vs_beta")
        mysavefig(ps_bs, pname; method)
    end

end 

## -------------------------------------------------------------------
# marginals v2
let 
    objider = iJR.BIOMASS_IDER
    size = [300, 250]
    Fd_mets_map = iJR.load_mets_map()
    exch_met_map = iJR.load_exch_met_map()

    # Iders
    model_iders, Fd_iders = [objider], ["D"]
    for Fd_met in EXCH_FLX_IDERS
        model_met = Fd_mets_map[Fd_met]
        model_exch = Fd_rxns_map[Fd_met]
        push!(model_iders, model_exch)
        push!(Fd_iders, string("u", Fd_met))
    end
    
    for (model_ider, Fd_ider) in zip(model_iders, Fd_iders)
        marg_params = (;xlabel = string(Fd_ider), yaxis = nothing, ylabel = "prob")

        epps = Plots.Plot[]
        exps = Plots.Plot[]
        for method in ALL_METHODS
            expp = plot(;title = string("Experimental"), marg_params...)
            epp = plot(;title = string(" MaxEnt: ", method), marg_params...)
            margin, m, M = -Inf, Inf, -Inf
            
            # EP
            for exp in EXPS
                Fd_av = Fd.val(Fd_ider, exp)
                color = exp_colors[exp]    

                datfile = INDEX[method, :DFILE, exp]
                dat = deserialize(datfile)
                model = dat[:model]
                objidx = ChU.rxnindex(model, objider)
                epouts = dat[:epouts]
                exp_beta = maximum(keys(epouts))
                epout = epouts[exp_beta]
                ep_av = ChU.av(model, epout, model_ider)
                ep_va = sqrt(ChU.va(model, epout, model_ider))
                        
                ChP.plot_marginal!(epp, model, [epout], model_ider; 
                    legend = false, color, alpha = 0.8, lw = 3)
                
                m = minimum([m, ep_av, Fd_av])
                M = maximum([M, ep_av, Fd_av])
                margin = maximum([margin, 3 * ep_va])

                # Experimental
                vline!(expp, [Fd_av]; label = "", lw = 3, color, alpha = 0.8)
                
            end
            
            map([expp, epp]) do p
                plot!(p; xlim = [m - margin, M + margin], size)
            end

            push!(epps, epp)
            push!(exps, expp)
        end

        extras = Plots.Plot[]
        for k in [:xi, :X, Fd_ider] |> unique
            p = plot(;title = "Experimental", size, 
                xlabel = "rep", ylabel = string(k))
            xticks =  (EXPS, string.(EXPS))
            vals = [Fd.val(k, exp) for exp in EXPS]
            color = [exp_colors[exp] for exp in EXPS]
            p = bar!(p, EXPS, vals; label = "", xticks, color)
            push!(extras, p)
        end

        # legend
        p = plot(;title = "Legend", size)
        for (method, color) in method_colors
            color = method_colors[method]
            bar!(p, [string(method)], [1]; yaxis = nothing, 
                color, xrotation = 35, label = ""
            )
        end
        push!(extras, p)

        ps = Plots.Plot[exps; epps; extras]
        layout = (3, 4)
        pname = string(Fd_ider, "_marginals_v2")
        mysavefig(ps, pname; layout)

    end # for (model_ider, Fd_ider)

end 
