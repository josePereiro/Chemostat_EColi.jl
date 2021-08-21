## -------------------------------------------------------------------
# flx correlations
let

    MEname(c) = "\$\\textbf{ME}^{\\langle \\textbf{$c} \\rangle} \\textbf{ (abs. flux)} \$"

    me_method = :ME_MAX_POL
    fba_method = :FBA_Z_FIX_MIN_COST
    me_label = MEname(1)
    fba_label = "\$ \\textbf{FBA} \\textbf{ (abs. flux)} \$"

    for (get_flxs, flxs) in [
        ((DAT) -> DAT[:EXCH_FLX_IDERS], "all"),
        ((DAT) -> ["GLC", "AC"], "common"),
    ]
        
        fontsize = 13
        p = plot(;
            title = "exchanges",
            xlabel = fba_label, ylabel = me_label, 
            thickness_scaling = 1.3, dpi = 1000,
            xguidefontsize = fontsize, yguidefontsize = fontsize
        )
    
        all_vals = []
        for (iJR, Data) in [
            (ChK.iJR904, ChK.KayserData),
            (ChN.iJR904, ChN.NanchenData), 
            (ChF.iJR904, ChF.FolsomData),
            # (ChH.iJR904, ChH.HeerdenData),
        ]
            src = nameof(Data)
            marker = (6, source_markers[Data])

            datfile = iJR.procdir("dat.bson")
            DAT = UJL.load_data(datfile; verbose = false)

            # EXCH_FLX_IDERS = DAT[:EXCH_FLX_IDERS]
            # EXCH_FLX_IDERS = ["GLC", "AC"]
            EXCH_FLX_IDERS = get_flxs(DAT)
            EXPS = DAT[:EXPS]

            tot_ps = Plots.Plot[]

            # total corr
            let            
                ep_vals = DAT[me_method, :flx, EXCH_FLX_IDERS, EXPS] .|> abs
                ep_errs = DAT[me_method, :err, EXCH_FLX_IDERS, EXPS] .|> abs

                fba_vals = DAT[fba_method, :flx, EXCH_FLX_IDERS, EXPS] .|> abs
                
                color = [ider_colors[ider] for ider in EXCH_FLX_IDERS, exp in EXPS]
                
                # ep corr
                scatter!(p, fba_vals, ep_vals; yerr = ep_errs, 
                    marker, label = "", color, alpha = 0.7
                )
                push!(all_vals, [ep_vals; fba_vals])
            end
        end

        sort!(all_vals)
        plot!(p, all_vals, all_vals; alpha = 0.7,
            ls = :solid, color = :black, label = ""
        )
        mysavefig(p, "fba_vs_me_exch_corr"; flxs)
    end
end 

