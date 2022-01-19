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
        
    end

    p
end