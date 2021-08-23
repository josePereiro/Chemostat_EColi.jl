## -------------------------------------------------------------------
# Xexp Xmax
function plot_Xexp_Xmax_study()
    methodL = :FBA_Z_FIX_MIN_VG_MIN_COST
    methodU = :FBA_Z_FIX_MAX_VG_MIN_COST
    
    p = plot(; 
        xlabel = _textbf("exp. X (log scale)"), 
        ylabel = _textbf("theo. max. X (log scale)"), 
    )

    all_vals = []
    for (src_name, iJR, Data) in [
        ("Kayser", ChK.iJR904, ChK.KayserData), 
        ("Nanchen", ChN.iJR904, ChN.NanchenData), 
        ("Folsom", ChF.iJR904, ChF.FolsomData), 
    ]
        src = nameof(Data)
        marker = (9, source_markers[Data])

        @info("Doing", src)
        
        DAT = ChE.load_DAT(src)
        EXPS = DAT[:EXPS]
        ug_mins = DAT[methodL, :flx, "GLC", EXPS]

        Xexps = [Data.val(:X, exp) for exp in EXPS]
        cgs = [Data.cval(:GLC, exp) for exp in EXPS]
        Ds = [Data.val(:D, exp) for exp in EXPS]
        Xmaxs = ((cgs .* Ds) ./ ug_mins) .|> abs

        Xexps = log10.(Xexps)
        Xmaxs = log10.(Xmaxs)

        scatter!(p, Xexps, Xmaxs; 
            label = _textbf(src_name), 
            color = :black,
            marker, alpha = 0.8,
        )
        # # fake label for spacing
        # scatter!(p, [Xexps[end]], [Xmaxs[end]]; label = " ", alpha = 0.0)

        plot!(p; legend = :topleft)
        
        push!(all_vals, Xexps...)
    end

    sort!(all_vals)
    plot!(p, all_vals, all_vals; 
        label = "", lw = 3, ls = :dash, 
        alpha = 0.7, color = :black
    )
    plot!(;legend = :bottomright)
    # sfig(ChE, p, 
    #     "Xexp_Xmax_study", ".png"
    # )

    return p
end 