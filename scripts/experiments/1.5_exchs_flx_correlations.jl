## -------------------------------------------------------------------
# flx correlations
function plot_exch_corrs(METHODS, METHODS_LABELS = string.(METHODS))

    ps = Plots.Plot[]

    for (method_label, method) in zip(METHODS_LABELS, METHODS)

        p = plot(;
            xlabel = _textbf("exp abs. flux"), 
            ylabel = _textbf("model abs. flux")
        )

        all_vals = []

        for (iJR, Data) in [
            (ChN.iJR904, ChN.NanchenData), 
            (ChK.iJR904, ChK.KayserData),
            (ChF.iJR904, ChF.FolsomData),
        ]
            src = nameof(Data)
            marker = (10, source_markers[Data])

            @show src
            DAT = ChE.load_DAT(src)

            EXCH_FLX_IDERS = DAT[:EXCH_FLX_IDERS]
            EXPS = DAT[:EXPS]
            
            for ider in EXCH_FLX_IDERS, exp in EXPS

                model_val = DAT[method, :flx, ider, exp] |> abs
                model_err = DAT[method, :err, ider, exp] |> abs 
                exp_val = DAT[:exp, :flx, ider, exp] |> abs
                
                # corr
                color = ider_colors[ider]
                scatter!(p, [exp_val], [model_val]; yerr = [model_err], 
                    marker, label = "", color, alpha = 0.7
                )

                push!(all_vals, exp_val, model_val)
            end

        end # for (iJR, Data)
        
        # ref
        sort!(all_vals) 
        plot!(p, all_vals, all_vals; label = "",
            ls = :dash, alpha = 0.5, lw = 3, 
            color = :black,
        )   

        # title
        m, M = minimum(all_vals), maximum(all_vals)
        x = m + (M - m) * 0.01
        y = M - (M - m) * 0.05
        annotate!(p, [(x, y, (_textbf(method_label), 26, :left, :top, :black))])

        push!(ps, p)

    end # for (method_label, method)

    sfig(ChE, ps, 
        "exchs_corr", ".png", 
        layout = (1, length(ps)), 
    )

    return ps
end 