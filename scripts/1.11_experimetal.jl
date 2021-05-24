function do_exp_plot(id1, id2)

    tot_p = nothing
    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData), 
        (ChF.iJR904, ChF.FolsomData), 
        (ChN.iJR904, ChN.NanchenData), 
        # (ChH.iJR904, ChH.HeerdenData), 
    ]

        src = nameof(Data)
        # unit1 = Data.unit(id1)
        # unit2 = Data.unit(id2)

        datfile = iJR.procdir("dat.bson")
        DAT = UJL.load_data(datfile; verbose = false)
        EXPS = DAT[:EXPS]
        
        local_p = plot(;
            xlabel = string(id1),
            ylabel = string(id2)
        )
        tot_p = isnothing(tot_p) ? deepcopy(local_p) : tot_p

        xs = abs.([Data.val(id1, exp, 0.0) for exp in EXPS])
        sidxs = sortperm(xs)
        xs = xs[sidxs]
        ys = abs.([Data.val(id2, exp, 0.0) for exp in EXPS][sidxs])

        marker = (8, source_markers[Data])
        for p in [local_p, tot_p]
            scatter!(p, xs, ys; 
                label = "", marker, 
                alpha = 0.8, color = :black
            )
            # plot!(p, xs, ys; 
            #     label = "", lw = 2, ls = :dash, 
            #     alpha = 0.8, color = :black
            # )
        end

        plot!(local_p; title = string(src))
        mysavefig(local_p, "experimental_corr"; id1, id2, src)
    end

    mysavefig(tot_p, "experimental_corr"; id1, id2)
end

## ---------------------------------------------------------------------------
let
    # ids = [:sGLC, :sAC, :D, :X, :uGLC, :uAC]
    ids = [:D, :sGLC]
    for i in 1:length(ids)
        for j in (i + 1):length(ids)
            do_exp_plot(ids[i], ids[j])
        end
    end
end

## ---------------------------------------------------------------------------
let
    fontsize = 13
    commom_params = (; 
        dpi = 1000,
        thickness_scaling = 1.3, 
        xguidefontsize = fontsize, yguidefontsize = fontsize, 
        ylim = [0.0, Inf]
    )

    tot_p = nothing
    for (iJR, Data) in [
        # (ChK.iJR904, ChK.KayserData), 
        # (ChF.iJR904, ChF.FolsomData), 
        # (ChN.iJR904, ChN.NanchenData), 
        (ChH.iJR904, ChH.HeerdenData), 
    ]

        src = nameof(Data)

        
        DAT = ChE.load_DAT(src)
        EXPS = DAT[:EXPS]
        # EXPS = Data.EXPS
        
        local_p = plot(;
            xlabel = "D (1/ h)",
            ylabel = "sGLC / cGLC",
            commom_params...
        )
        tot_p = isnothing(tot_p) ? deepcopy(local_p) : tot_p

        Ds = abs.([Data.val(:D, exp, 0.0) for exp in EXPS])
        sidxs = sortperm(Ds)
        Ds = Ds[sidxs]
        sGLC = abs.([Data.val(:sGLC, exp, 0.0) for exp in EXPS][sidxs])
        cGLC = abs.([Data.val(:cGLC, exp, 0.0) for exp in EXPS][sidxs])

        marker = (8, source_markers[Data])
        for p in [local_p, tot_p]

            # sGLC
            scatter!(p, Ds, sGLC ./ cGLC; 
                label = "", marker, 
                alpha = 0.8, color = :black
            )
            plot!(p, Ds, sGLC ./ cGLC; 
                label = "", lw = 2, ls = :dash, 
                alpha = 0.8, color = :black
            )
        end

        plot!(local_p; title = string(src))
        mysavefig(local_p, "exp_GLC_vs_D"; src)
    end

    mysavefig(tot_p, "exp_GLC_vs_D")
end