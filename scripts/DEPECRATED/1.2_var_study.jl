let
    p0 = plot(;xlabel = "politope box volume (log)", ylabel = "std (mmol/ gDW hr)")
    all_p = deepcopy(p0)
    zeroth = 1e-2
    for (iJR, Data) in [
                (ChK.iJR904, ChK.KayserData), 
                (ChN.iJR904, ChN.NanchenData), 
                (ChF.iJR904, ChF.FolsomData), 
                (ChH.iJR904, ChH.HeerdenData), 
            ]
        
        model = iJR.load_model("max_model")
        datfile = iJR.procdir("dat.bson")
        DAT = UJL.load_data(datfile; verbose = false)
        
        method = :ME_MAX_POL
        # EXCH_FLX_IDERS = DAT[:EXCH_FLX_IDERS]
        EXCH_FLX_IDERS = ["GLC", "AC"]
        EXPS = DAT[:EXPS]
        rxns_map = iJR.load_rxns_map()
        
        exglcidx = ChU.rxnindex(model, iJR.EX_GLC_IDER)
        biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
        
        model_ids = [rxns_map[id] for id in EXCH_FLX_IDERS]
        model_idxs = [ChU.rxnindex(model, model_id) for model_id in model_ids]
        
        
        # Reduce Pol
        p = deepcopy(p0)
        vols = []

        # --------------------------------------------------- 
        # ME_MAX_POL
        for exp in EXPS
            # box vol
            cgD_X = -Data.ciD_X(:GLC, exp)
            D = Data.val(:D, exp)

            ChU.lb!(model, exglcidx, cgD_X)
            ChU.bounds!(model, biomidx, D, D)
            vol = box_vol(model, model_idxs) .+ zeroth
            isnan(vol) && continue

            push!(vols, vol)
        end
        sidxs = sortperm(vols)

        for ider in EXCH_FLX_IDERS
            stds = DAT[method, :err, ider, EXPS]
            color = ider_colors[ider]
            marker = source_markers[Data]

            for p_ in [p, all_p]
                scatter!(p_, vols[sidxs], stds[sidxs]; 
                    label = "", m = (8, marker), color
                )
            end
        end

        # ---------------------------------------------------
        # ME_MAX_POL_B0
        method = :ME_MAX_POL_B0
        model = iJR.load_model("max_model")
        vol = box_vol(model, model_idxs) + zeroth

        for ider in EXCH_FLX_IDERS

            std = DAT[method, :err, ider]
            marker = source_markers[Data]

            for (m, color) in [
                    [(10, marker), :black], 
                    [(8, marker), ider_colors[ider]], 
                    [(3, :circle), :black], 
                ]
                scatter!(all_p, [vol], [std]; label = "", m, color)
            end
        end

        src = nameof(Data)
        mysavefig(p, "pol_box_vs_var"; src)
    end
    plot!(all_p; xscale = :log10)

    mysavefig(all_p, "pol_box_vs_var")
end  