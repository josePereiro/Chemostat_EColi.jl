## ---------------------------------------------------------------------------------
function _box_vol(model, idxs)
    try
        L, U = ChLP.fva(model, idxs)
        vol = prod(abs.(L .- U))
        max(0.0, log10(vol + 1e-50))
    catch; NaN end
end

## ---------------------------------------------------------------------------------
function plot_pol_box_size(iJR, Data; 
        Dlim = nothing, 
        cgD_Xlim = nothing,
        bins = 50
    )    
    src = nameof(Data)

    model0 = iJR.load_model("max_model")
    rxns_map = iJR.load_rxns_map()
    offsetf = 1.1
    
    exglcidx = ChU.rxnindex(model0, iJR.EX_GLC_IDER)
    biomidx = ChU.rxnindex(model0, iJR.BIOMASS_IDER)
    exglcL, exglcU = ChU.bounds(model0, exglcidx)
    
    maxD = maximum(Data.val(:D)) * offsetf 
    max_cgD_X = maximum(Data.ciD_X(:GLC)) * offsetf
    
    Dlim = isnothing(Dlim) ? [0.01, maxD] : Dlim
    Ds = range(Dlim...; length = bins)
    @show Ds
    
    # ug negative means uptake !!!
    cgD_Xlim = isnothing(cgD_Xlim) ? [abs(exglcU), max_cgD_X] : cgD_Xlim
    cgD_Xs = range(cgD_Xlim...; length = bins) # make possitive
    @show cgD_Xs

    cid = (:BOX_VOL, src, Dlim, cgD_Xlim, bins)
    box_vols = lcache(ChE, cid) do

        box_vols_ = zeros(bins, bins) 
    
        ids = ["GLC"] 

        model_ids = [rxns_map[id] for id in ids]
        model_idxs = [ChU.rxnindex(model0, model_id) for model_id in model_ids]

        # feeding task
        Ch = Channel(nthreads()) do Ch_
            @showprogress for (Di, D) in enumerate(Ds)
                for (cgD_Xi, cgD_X) in enumerate(cgD_Xs)
                    put!(Ch_, (Di, D, cgD_Xi, cgD_X))
                end
            end 
        end

        # compute volume map
        @threads for _ in 1:nthreads()
            model = deepcopy(model0)
            for (Di, D, cgD_Xi, cgD_X) in Ch
                
                # Reduce Pol
                # ug negative means uptake !!!
                ChU.lb!(model, exglcidx, -cgD_X)
                ChU.bounds!(model, biomidx, D, D)
                box_vols_[Di, cgD_Xi] = _box_vol(model, model_idxs)
            end
        end

        return box_vols_
    end

    # vol map
    p = heatmap(Ds, cgD_Xs, box_vols'; 
        label = "", 
        xlabel = _textbf("D"), 
        ylabel = _textbf("c_g", "D/X"), 
        color = :greys,
    )

    fontsize = 20
    text_ = _textbf("polytope volume")
    plot!(p; rightmargin = 19mm)
    annotate!(p, 0.66, 20.0, text(text_, fontsize, :black, rotation = 90))

    # exp vals
    marker = (10, source_markers[Data])
    scatter!(p, Data.val(:D), abs.(Data.ciD_X(:GLC));
        label = "", color = :white, marker
    )

    # sfig(ChE, p, 
    #     @fileid, "pol_box_volume_with_dat", (;src, bins), ".png"
    # )

    return p
end
