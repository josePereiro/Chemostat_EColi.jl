## -------------------------------------------------------------------
# Glc uptake
function plot_Glc_uptake_study()
    
    methodL = :FBA_Z_FIX_MIN_VG_MIN_COST
    methodU = :FBA_Z_FIX_MAX_VG_MIN_COST
    
    p = plot(; 
        xlabel = "\$ \\textbf{exp. vg} \$", 
        ylabel = "\$ \\textbf{model. vg} \$", 
    )

    all_vals = []
    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData), 
        (ChN.iJR904, ChN.NanchenData), 
        (ChF.iJR904, ChF.FolsomData), 
        # (ChH.iJR904, ChH.HeerdenData), 
    ]
        src = nameof(Data)
        marker = (9, source_markers[Data])

        DAT = ChE.load_DAT(src)

        EXPS = DAT[:EXPS]
    
        @info("Doing", src)
        
        exp_vals = DAT[:exp, :flx, "GLC", EXPS]
        model_Uvals = DAT[methodU, :flx, "GLC", EXPS]
        model_Lvals = DAT[methodL, :flx, "GLC", EXPS]

        # color = [exp_colors[exp] for exp in EXPS]
        scatter!(p, abs.(exp_vals), abs.(exp_vals); 
            label = "",
            color = :red,
            marker, alpha = 0.8,
        )
        scatter!(p, abs.(exp_vals), abs.(model_Lvals); 
            label = "",
            color = :blue,
            marker, alpha = 0.8,
        )
        
        push!(all_vals, exp_vals...)
    end

    sort!(all_vals)
    plot!(p, abs.(all_vals), abs.(all_vals); 
        label = "", lw = 3, ls = :dash, 
        alpha = 0.7, color = :black
    )
    plot!(;legend = :bottomright)
    sfig(ChE, p, 
        "glc_uptake_study", ".png"
    )

    p
end 