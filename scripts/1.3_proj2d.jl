
## -------------------------------------------------------------------
# proj 2D
let
    me_method = :ME_MAX_POL
    fba_method = :FBA_Z_FIX_MIN_COST


    for (iJR, Data) in [
            (ChK.iJR904, ChK.KayserData), 
            (ChF.iJR904, ChF.FolsomData), 
            (ChH.iJR904, ChH.HeerdenData), 
        ]

        src = nameof(Data)

        datfile = iJR.procdir("dat.bson")
        DAT = UJL.load_data(datfile; verbose = false)
        FLX_IDERS = DAT[:FLX_IDERS]
        EXPS = DAT[:EXPS]

        ps_pool = Dict()
        for exp in EXPS

            for ider in FLX_IDERS

                # 2D Projection
                p = plot(;title = string(src, " exp:", exp), 
                    xlabel = "D", ylabel = ider,
                    legend = :right
                )
                proj = DAT[me_method, :proj, ider, exp]
                ChP.plot_projection2D!(p, proj; l = 50)

                # bounds
                lb, ub = DAT[:bounds, ider, exp]
                hline!(p, [lb]; lw = 3, 
                    label = "fva lb",
                    color = :blue, ls = :solid
                )
                hline!(p, [ub]; lw = 3,
                    label = "fva ub", 
                    color = :red, ls = :solid
                )

                # EXPERIMENTAL FLXS
                exp_biom = DAT[:exp, :flx, "D", exp]
                exp_exch = DAT[:exp, :flx, ider, exp]
                scatter!(p, [exp_biom], [exp_exch]; 
                    m = 8, color = :black, label = "exp"
                )
                
                # MAXENT FLXS
                me_biom = DAT[me_method, :flx, "D", exp]
                me_biom_err = DAT[me_method, :err, "D", exp]
                me_exch = DAT[me_method, :flx, ider, exp]
                me_exch_err = DAT[me_method, :err, ider, exp]
                scatter!(p, [me_biom], [me_exch]; 
                    xerr = [me_biom_err], yerr = [me_exch_err],
                    m = 8, color = :red, label = "maxent"
                )
                
                # FBA FLXS
                fba_biom = DAT[fba_method, :flx, "D", exp]
                fba_exch = DAT[fba_method, :flx, ider, exp]
                scatter!(p, [fba_biom], [fba_exch]; 
                    m = 8, color = :blue, label = "fba"
                )

                ps_pool[(exp, ider)] = deepcopy(p)
            end
        end

        # collect 
        for exp in EXPS
            ps = Plots.Plot[ps_pool[(exp, ider)] for ider in FLX_IDERS]
            mysavefig(ps, "polytope"; exp, src, me_method, fba_method)
        end

        for ider in FLX_IDERS
            ps = Plots.Plot[ps_pool[(exp, ider)] for exp in EXPS]
            mysavefig(ps, "polytope"; ider, src, me_method, fba_method)
        end
    end # for (iJR, Data)
end
