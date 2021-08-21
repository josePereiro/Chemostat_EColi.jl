## -------------------------------------------------------------------
# Yield corr
let
    method = FBA_Z_FIX_ADJ_MIN_VG_MIN_COST

    fontsize = 13
    common_params = (;
        thickness_scaling = 1.3, 
        xguidefontsize = fontsize, yguidefontsize = fontsize,
        dpi = 1000
    )
    
    p = plot(; 
        xlabel = "D", 
        ylabel = "biom. adj. factor", 
        common_params...
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

        
        DAT = ChE.load_DAT(src)
        EXPS = DAT[:EXPS]
        LP_DAT = ChE.load_LP_DAT(src)
    
        @info("Doing", src, method)
        
        D = DAT[:exp, :flx, "D", EXPS] 
        exp_vals = LP_DAT[method, :Sfac, EXPS]

        scatter!(p, abs.(D), abs.(exp_vals); 
            label = "",
            color = :black,
            marker, alpha = 0.8,
        )
        
    end

    plot!(;legend = :bottomright)
    mysavefig(p, "biom_adj_study")
end 