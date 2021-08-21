let
    model0 = iJR.load_model("max_model")

    bins = 50
    Fd_rxns_map = iJR.load_rxns_map()
    offsetf = 1.1
    
    exglcidx = ChU.rxnindex(model0, iJR.EX_GLC_IDER)
    biomidx = ChU.rxnindex(model0, iJR.BIOMASS_IDER)
    exglcL, exglcU = ChU.bounds(model0, exglcidx)

    
    maxD = maximum(Fd.val(:D)) * offsetf 
    max_cgD_X = -maximum(Fd.ciD_X(:GLC)) * offsetf
    
    Ds = range(0.01, maxD; length = bins)
    cgD_Xs = range(max_cgD_X, exglcU; length = bins)

    box_vols = zeros(bins, bins) 
    
    Fd_ids = ["GLC", "PYR", "SUCC", "LAC", "FORM", "AC"]
    model_ids = [Fd_rxns_map[Fd_id] for Fd_id in Fd_ids]
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
        thid = threadid()
        thid == 1 && continue
        model = deepcopy(model0)
        for (Di, D, cgD_Xi, cgD_X) in Ch
            
            # Reduce Pol
            ChU.lb!(model, exglcidx, cgD_X)
            ChU.bounds!(model, biomidx, D, D)
            try
                L, U = ChLP.fva(model, model_idxs)
                vol = prod(abs.(L .- U))
                box_vols[Di, cgD_Xi] = max(0.0, log10(vol + 1e-50))
            catch
                box_vols[Di, cgD_Xi] = NaN
            end
        end
    end

    # vol map
    p = heatmap(Ds, cgD_Xs, box_vols'; 
        title = "Polytope Box volume", label = "", 
        xlabel = "D", ylabel = "cgD_X"
    )

    # exp vals
    scatter!(p, Fd.val(:D), -Fd.ciD_X(:GLC);
        label = "exp", color = :white, 
        m = 8
    )

    mysavefig(p, "pol_box_volume"; bins)
end
