let
    method = ME_MAX_POL
    biom_ider = iJR.BIOMASS_IDER

    ps_pool = Dict()
    for exp in EXPS

        for Fd_ider in FLX_IDERS

            # 2D Projection
            p = plot(;title = string("Folsom2014 exp:", exp), 
                xlabel = string(biom_ider), ylabel = string(Fd_ider),
                legend = :left
            )
            proj = DAT[method, :ep, :proj, Fd_ider, exp]
            ChP.plot_projection2D!(p, proj; l = 50)

            # bounds
            lb, ub = DAT[method, :bounds, :flx, Fd_ider, exp]
            hline!(p, [lb]; lw = 3, 
                label = "fva lb",
                color = :blue, ls = :solid
            )
            hline!(p, [ub]; lw = 3,
                label = "fva ub", 
                color = :red, ls = :solid
            )

            # EXPERIMENTAL FLXS
            exp_biom = DAT[method, :Fd, :flx, biom_ider, exp]
            exp_exch = DAT[method, :Fd, :flx, Fd_ider, exp]
            scatter!(p, [exp_biom], [exp_exch]; 
                m = 8, color = :red, label = "exp"
            )
            
            # MAXENT FLXS
            ep_biom = DAT[method, :ep, :flx, biom_ider, exp]
            ep_biom_err = DAT[method, :eperr, :flx, biom_ider, exp]
            ep_exch = DAT[method, :ep, :flx, Fd_ider, exp]
            ep_exch_err = DAT[method, :eperr, :flx, Fd_ider, exp]
            scatter!(p, [ep_biom], [ep_exch]; 
                xerr = [ep_biom_err], yerr = [ep_exch_err],
                m = 8, color = :blue, label = "maxent"
            )

            # mysavefig(p, "polytope"; Fd_ider, exp, method)
            ps_pool[(exp, Fd_ider)] = deepcopy(p)
        end
    end

    # collect 
    for exp in EXPS
        ps = Plots.Plot[ps_pool[(exp, Fd_ider)] for Fd_ider in FLX_IDERS]
        mysavefig(ps, "polytope"; exp, method)
    end

    for Fd_ider in FLX_IDERS
        ps = Plots.Plot[ps_pool[(exp, Fd_ider)] for exp in EXPS]
        mysavefig(ps, "polytope"; Fd_ider, method)
    end
end