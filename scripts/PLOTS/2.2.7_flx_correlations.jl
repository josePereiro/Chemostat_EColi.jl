## -------------------------------------------------------------------
# EP biomass corr
let
    objider = iJR.BIOMASS_IDER
    ps = Plots.Plot[]
    for method in ALL_METHODS
        p = plot(title = string(iJR.PROJ_IDER, " method: ", method), 
            xlabel = "exp biom", ylabel = "model biom")
        ep_vals = DAT[method, :ep, :flx, objider, EXPS]
        eperr_vals = DAT[method, :eperr, :flx, objider, EXPS]
        Fd_vals = DAT[method, :Fd, :flx, objider, EXPS]
        Fderr_vals = DAT[method, :Fderr, :flx, objider, EXPS]
        color = [exp_colors[exp] for exp in EXPS]
        m, M = myminmax([Fd_vals; ep_vals])
        margin = abs(M - m) * 0.1
        scatter!(p, Fd_vals, ep_vals; 
            xerr = Fderr_vals,
            yerr = eperr_vals,
            label = "", color,
            alpha = 0.7, ms = 7,
            xlim = [m - margin, M + margin],
            ylim = [m - margin, M + margin],
        )
        push!(ps, p)
    end
    layout = (1, length(ps))
    mysavefig(ps, "obj_val_ep_corr"; layout)
end

## -------------------------------------------------------------------
# correlations
let

    tot_ps = Plots.Plot[]
    for method in ALL_METHODS        
        # total corr
        let            
            ep_vals = DAT[method, :ep, :flx, FLX_IDERS, EXPS] .|> abs
            ep_errs = DAT[method, :eperr, :flx, FLX_IDERS, EXPS] .|> abs
            Fd_vals = DAT[method, :Fd, :flx, FLX_IDERS, EXPS] .|> abs
            
            color = [ider_colors[ider] for ider in FLX_IDERS, exp in EXPS]
            scatter_params = (;label = "", color, ms = 7, alpha = 0.7)
            # ep corr
            p = plot(title = "$(iJR.PROJ_IDER) (EP) $method", 
                ylabel = "model abs flx",
                xlabel = "exp abs flx", 
            )
            scatter!(p, Fd_vals, ep_vals; yerr = ep_errs, scatter_params...)
            all_vals = [ep_vals; Fd_vals] |> sort!
            plot!(p, all_vals, all_vals; ls = :dash, color = :black, label = "")
            push!(tot_ps, deepcopy(p))
        end

        # per ider
        let       
            for ider in FLX_IDERS
                ep_vals = DAT[method, :ep, :flx, ider, EXPS] .|> abs
                ep_errs = DAT[method, :eperr, :flx, ider, EXPS] .|> abs
                Fd_vals = DAT[method, :Fd, :flx, ider, EXPS] .|> abs
                
                color = ider_colors[ider]
                # ep corr
                p = plot(title = "$(iJR.PROJ_IDER) (EP) $method", 
                    ylabel = "model abs flx",
                    xlabel = "exp abs flx", 
                )
                scatter!(p, Fd_vals, ep_vals; yerr = ep_errs, 
                    label = "", color, ms = 7, alpha = 0.7
                )
                bounds = DAT[method, :bounds, :flx, ider, EXPS]
                all_vals = [ep_vals; Fd_vals] |> sort!
                plot!(p, all_vals, all_vals; 
                    ls = :dash, color = :black, label = "", 
                    # xlim = [lb, ub], ylim = [lb, ub]
                )
                mysavefig(p, "corr"; ider, method)
            end
        end
    end

    layout = (1, length(tot_ps))
    mysavefig(tot_ps, "flx_tot_corr"; layout)

end 