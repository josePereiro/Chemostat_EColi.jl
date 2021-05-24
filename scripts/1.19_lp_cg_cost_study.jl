## -------------------------------------------------------------------
let
    METHODS = [

        FBA_Z_FIX_MAX_VG_MIN_COST, 
        FBA_Z_FIX_MIN_VG_MIN_COST, 
        FBA_Z_FIX_MAX_VG_MAX_COST, 
        FBA_Z_FIX_MIN_VG_MAX_COST,
        FBA_Z_VG_FIX_MAX_COST, 
        FBA_Z_VG_FIX_MIN_COST, 
        FBA_Z_FIX_MAX_COST, 
        FBA_Z_FIX_MIN_COST, 
        FBA_MAX_Z_MIN_COST, 
        FBA_MAX_Z_MAX_COST,

        ME_MAX_POL,
        ME_Z_EXPECTED_G_BOUNDED, 
        # ME_Z_OPEN_G_OPEN
    ]

    fontsize = 13
    common_params = (;
        thickness_scaling = 1.3, 
        xguidefontsize = fontsize, yguidefontsize = fontsize,
        # dpi = 1000
    )
    
    for method in METHODS
        p = plot(; 
            xlabel = "vg", 
            ylabel = "cost", 
            common_params...
        )

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
            
            @info("Doing", src, method)
            
            vg_vals = DAT[method, :flx, "GLC", EXPS]
            cost_vals = DAT[method, :flx, "cost", EXPS]

            scatter!(p, abs.(vg_vals), abs.(cost_vals); 
                label = "",
                color = :black,
                marker, alpha = 0.8,
            )
        end

        plot!(;legend = :bottomright, ylim = [0.0, Inf])
        mysavefig(p, "vg_cost_study"; method)
    end

end