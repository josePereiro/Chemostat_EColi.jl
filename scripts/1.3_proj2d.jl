
## -------------------------------------------------------------------
# proj 2D
let
    me_method = :ME_MAX_POL
    fba_method = :FBA_Z_FIX_MIN_COST


    for (iJR, Data) in [
            (ChK.iJR904, ChK.KayserData), 
            (ChN.iJR904, ChN.NanchenData), 
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

## -------------------------------------------------------------------
# TODO: update ChemostatPlots
function plot_projection2D!(p, proj; l = 10, spkwargs...)

    def_spkwargs = (;label = "", alpha = 0.5)

    flx1s = sort(collect(keys(proj)))
    rflx1s = reverse(flx1s)
    flx2lbs = [first(proj[flx1]) for flx1 in flx1s]
    flx2ubs = [last(proj[flx1]) for flx1 in rflx1s]

    s = Shape([flx1s; rflx1s], [flx2lbs; flx2ubs])
    plot!(p, [s]; label = "", def_spkwargs..., spkwargs...)

    return p
end

## -------------------------------------------------------------------
let
    me_method = :ME_MAX_POL
    # fba_method = :FBA_Z_FIX_MIN_COST
    
    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData), 
        (ChN.iJR904, ChN.NanchenData), 
        (ChF.iJR904, ChF.FolsomData), 
        # (ChH.iJR904, ChH.HeerdenData),  
    ]
        src = nameof(Data)

        datfile = iJR.procdir("dat.bson")
        DAT = UJL.load_data(datfile; verbose = false)
        FLX_IDERS = DAT[:FLX_IDERS]
        EXPS = DAT[:EXPS]

        
        for ider in FLX_IDERS
            # 2D Projection
            exp = EXPS[1]
            fontsize = 13
            p = plot(;
                # title = string(src, " exp:", exp), 
                xlabel = "D", ylabel = ider,
                legend = :right, 
                dpi = 1000,
                thickness_scaling = 1.3, 
                xguidefontsize = fontsize, yguidefontsize = fontsize
            )
            proj = DAT[me_method, :proj, ider, exp]
            # ChP.plot_projection2D!(p, proj)
            plot_projection2D!(p, proj)

            # EXPERIMENTAL FLXS
            for exp in EXPS
                exp_biom = DAT[:exp, :flx, "D", exp]
                exp_exch = DAT[:exp, :flx, ider, exp]
                marker = (8, source_markers[Data])
                scatter!(p, [exp_biom], [exp_exch]; 
                    marker, color = :white, label = "", 
                )
            end

            # plot!(p; ylim = [-10, 0.0])
            mysavefig(p, "proj2D"; src, ider)
        end
    end
end

