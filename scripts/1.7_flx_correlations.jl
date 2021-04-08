## -------------------------------------------------------------------
# EP biomass corr
let
    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData), 
        (ChF.iJR904, ChF.FolsomData), 
        (ChH.iJR904, ChH.HeerdenData), 
    ]
        src = nameof(Data)

        datfile = iJR.procdir("dat.bson")
        DAT = UJL.load_data(datfile; verbose = false)

        ME_METHODS = [:ME_MAX_POL]
        EXPS = DAT[:EXPS]

        ps = Plots.Plot[]
        for method in ME_METHODS
            p = plot(title = string(src, " method: ", method), 
                xlabel = "exp biom", ylabel = "model biom"
            )
            model_vals = DAT[method, :flx, "D", EXPS]
            eperr_vals = DAT[method, :err, "D", EXPS]
            exp_vals = DAT[:exp, :flx, "D", EXPS]
            color = [exp_colors[exp] for exp in EXPS]
            m, M = myminmax([exp_vals; model_vals])
            margin = abs(M - m) * 0.1
            scatter!(p, exp_vals, model_vals; 
                yerr = eperr_vals,
                label = "", color,
                alpha = 0.7, ms = 7,
                xlim = [m - margin, M + margin],
                ylim = [m - margin, M + margin],
            )
            push!(ps, p)
        end
        layout = (1, length(ps))
        mysavefig(ps, "obj_val_ep_corr"; src, layout)
    end
end

## -------------------------------------------------------------------
# flx correlations
let

    p0 = plot(;
        ylabel = "model abs flx",
        xlabel = "exp abs flx", 
    )
    ps_pool = Dict()

    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData),
        (ChF.iJR904, ChF.FolsomData),
        (ChH.iJR904, ChH.HeerdenData),
    ]
        src = nameof(Data)
        marker = (8, source_markers[Data])

        datfile = iJR.procdir("dat.bson")
        DAT = UJL.load_data(datfile; verbose = false)

        METHODS = [:FBA_Z_FIX_MIN_COST, :ME_MAX_POL]
        FLX_IDERS = DAT[:FLX_IDERS]
        EXPS = DAT[:EXPS]

        tot_ps = Plots.Plot[]
        for method in METHODS    
            mp0 = plot!(deepcopy(p0); title = string(method))    
            # total corr
            let            
                model_vals = DAT[method, :flx, FLX_IDERS, EXPS] .|> abs
                model_errs = haskey(DAT[method], :err) ? 
                    DAT[method, :err, FLX_IDERS, EXPS] .|> abs : 0.0
                exp_vals = DAT[:exp, :flx, FLX_IDERS, EXPS] .|> abs
                
                color = [ider_colors[ider] for ider in FLX_IDERS, exp in EXPS]
                # ep corr
                local_p = plot!(deepcopy(p0); title = "$(src): $(method)")
                top_p = get!(ps_pool, (;method), deepcopy(mp0))

                for p in [local_p, top_p]
                    scatter!(p, exp_vals, model_vals; yerr = model_errs, 
                        marker, label = "", color, alpha = 0.7
                    )
                    all_vals = [model_vals; exp_vals] |> sort!
                    plot!(p, all_vals, all_vals; 
                        ls = :solid, color = :black, label = ""
                    )
                end
                mysavefig(local_p, "corr"; src, method)
            end

            # per ider
            let       
                for ider in FLX_IDERS
                    model_vals = DAT[method, :flx, ider, EXPS] .|> abs
                    model_errs = haskey(DAT[method], :err) ? 
                        DAT[method, :err, ider, EXPS] .|> abs : 0.0
                    exp_vals = DAT[:exp, :flx, ider, EXPS]  .|> abs
                    
                    color = ider_colors[ider]
                    # ep corr
                    local_p = plot!(deepcopy(p0); title = "$(src): $(method)")
                    top_p = get!(ps_pool, (;ider, method), deepcopy(mp0))

                    for p in [local_p, top_p]
                        scatter!(p, exp_vals, model_vals; yerr = model_errs, 
                        marker, label = "", color, alpha = 0.7
                        )
                        all_vals = [model_vals; exp_vals] |> sort!
                        plot!(p, all_vals, all_vals; 
                            ls = :solid, color = :black, label = ""
                        )
                    end
                    mysavefig(local_p, "corr"; src, ider, method)
                end
            end
        end
    end

    for (kws, p) in ps_pool
        mysavefig(p, "corr"; kws...)
    end
end 