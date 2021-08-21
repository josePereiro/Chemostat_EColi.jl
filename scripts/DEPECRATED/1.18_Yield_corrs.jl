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

