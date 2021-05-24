## -------------------------------------------------------------------
let
    p = plot()
    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData), 
        (ChN.iJR904, ChN.NanchenData), 
        (ChF.iJR904, ChF.FolsomData), 
        # (ChH.iJR904, ChH.HeerdenData), 
    ]
        src = nameof(Data)
        marker = (6, source_markers[Data])
        scatter!(p, Data.val(:D), abs.(Data.val(:D) ./ Data.val(:uGLC));
            label = "", marker
        )
    end
    mysavefig(p, "exp_yield_vs_D_coor")
end

## -------------------------------------------------------------------
# EP biomass corr
let
    method = FBA_Z_FIX_MIN_VG_MAX_COST

    fontsize = 13
    common_params = (;
        thickness_scaling = 1.3, 
        xguidefontsize = fontsize, yguidefontsize = fontsize,
        dpi = 1000
    )
    
    p = plot(; 
        xlabel = "\$ \\textbf{exp. yield} \$", 
        ylabel = "\$ \\textbf{model max. yield} \$", 
        common_params...
    )

    all_vals = []
    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData), 
        (ChN.iJR904, ChN.NanchenData), 
        (ChF.iJR904, ChF.FolsomData), 
        (ChH.iJR904, ChH.HeerdenData), 
    ]
        src = nameof(Data)
        marker = (6, source_markers[Data])

        DAT = ChE.load_DAT(src)

        EXPS = DAT[:EXPS]
        # EXPS = Data.EXPS
    
        @info("Doing", src, method)
        
        model_vals = DAT[method, :flx, "D", EXPS] ./ DAT[method, :flx, "GLC", EXPS]

        exp_vals = DAT[:exp, :flx, "D", EXPS] ./ DAT[:exp, :flx, "GLC", EXPS]
        # color = [exp_colors[exp] for exp in EXPS]
        scatter!(p, abs.(exp_vals), abs.(model_vals); 
            label = "",
            color = :black,
            marker, alpha = 0.8,
        )
        
        push!(all_vals, exp_vals..., model_vals...)
    end

    sort!(all_vals)
    plot!(p, abs.(all_vals), abs.(all_vals); 
        label = "", lw = 3, ls = :dash, 
        alpha = 0.7, color = :black
    )
    plot!(;legend = :bottomright)
    mysavefig(p, "max_yield_coor")
end 