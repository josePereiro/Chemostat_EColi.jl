function do_exp_plot(id1, id2)

    tot_p = nothing
    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData), 
        (ChF.iJR904, ChF.FolsomData), 
        (ChH.iJR904, ChH.HeerdenData), 
    ]

        src = nameof(Data)
        # unit1 = Data.unit(id1)
        # unit2 = Data.unit(id2)

        local_p = plot(;
            xlabel = string(id1),
            ylabel = string(id2)
        )
        tot_p = isnothing(tot_p) ? deepcopy(local_p) : tot_p

        xs = abs.([Data.val(id1, exp, 0.0) for exp in Data.EXPS])
        sidxs = sortperm(xs)
        xs = xs[sidxs]
        ys = abs.([Data.val(id2, exp, 0.0) for exp in Data.EXPS][sidxs])

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
    ids = [:sGLC, :sAC, :D, :X, :uGLC, :uAC]
    for i in 1:length(ids)
        for j in (i + 1):length(ids)
            do_exp_plot(ids[i], ids[j])
        end
    end
end