## -------------------------------------------------------------------
# Glc uptake
function plot_residual_Glc_study()
    
    p = plot(; 
        xlabel = _textbf("experiment "), 
        ylabel = _textbf("s_g", "/", "c_g"), 
    )

    for Ed in [ChK.KayserData, ChN.NanchenData, ChF.FolsomData]
        
        src = nameof(Ed)
        marker = (9, source_markers[Ed])

        DAT = ChE.load_DAT(src)
        EXPS = DAT[:EXPS]
    
        @info("Doing", src)
        
        sGLC = Ed.val.(:sGLC, EXPS, 0.0)
        cGLC = Ed.val.(:cGLC, EXPS, 0.0)

        scatter!(p, EXPS, abs.(sGLC ./ cGLC); 
            label = "",
            color = :black,
            marker,
        )
        # plot!(p, EXPS, abs.(sGLC ./ cGLC); 
        #     label = "",
        #     color = :black,
        #     ls = :dash, 
        #     alpha = 0.8,
        # )
        
    end

    # plot!(;legend = :bottomright)
    # sfig(ChE, p, 
    #     "glc_residual_study", ".png"
    # )

    p
end