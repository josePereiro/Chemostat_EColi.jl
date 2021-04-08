function plot_pol_box_size(iJR, Data; bins = 50)    
    
    model0 = iJR.load_model("max_model")
    rxns_map = iJR.load_rxns_map()
    offsetf = 1.1
    
    exglcidx = ChU.rxnindex(model0, iJR.EX_GLC_IDER)
    biomidx = ChU.rxnindex(model0, iJR.BIOMASS_IDER)
    exglcL, exglcU = ChU.bounds(model0, exglcidx)

    
    maxD = maximum(Data.val(:D)) * offsetf 
    max_cgD_X = -maximum(Data.ciD_X(:GLC)) * offsetf
    
    Ds = range(0.01, maxD; length = bins)
    cgD_Xs = range(max_cgD_X, exglcU; length = bins)

    box_vols = zeros(bins, bins) 
    
    ids = ["GLC", "AC", "NH4"]
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
        # thid = threadid()
        # thid == 1 && continue
        model = deepcopy(model0)
        for (Di, D, cgD_Xi, cgD_X) in Ch
            
            # Reduce Pol
            ChU.lb!(model, exglcidx, cgD_X)
            ChU.bounds!(model, biomidx, D, D)
            box_vols[Di, cgD_Xi] = box_vol(model, model_idxs)
        end
    end

    # vol map
    p = heatmap(Ds, cgD_Xs, box_vols'; 
        title = "Polytope Box volume", label = "", 
        xlabel = "D", ylabel = "cgD_X"
    )

    # exp vals
    scatter!(p, Data.val(:D), -Data.ciD_X(:GLC);
        label = "exp", color = :white, 
        m = 8
    )

    src = nameof(Data)
    mysavefig(p, "pol_box_volume"; src, bins)
end

## ------------------------------------------------------------------------
let
    # setup
    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData), 
        (ChF.iJR904, ChF.FolsomData), 
        (ChH.iJR904, ChH.HeerdenData), 
        ]
        plot_pol_box_size(iJR, Data; bins = 10)
    end
end

## -------------------------------------------------------------------
# method_diffs vs vol
let
    me_method = :ME_MAX_POL
    fba_method = :FBA_Z_FIX_MIN_COST

    p0 = plot(;xlabel = "politope box volume (log)", ylabel = "abs(ME - FBA)")
    all_p = deepcopy(p0)
    zeroth = 1e-2
    for (iJR, Data) in [
                (ChK.iJR904, ChK.KayserData), 
                (ChF.iJR904, ChF.FolsomData), 
                (ChH.iJR904, ChH.HeerdenData), 
            ]
        
        src = nameof(Data)

        model = iJR.load_model("max_model")
        datfile = iJR.procdir("dat.bson")
        DAT = UJL.load_data(datfile; verbose = false)
        
        FLX_IDERS = DAT[:FLX_IDERS] #["GLC", "AC"]
        EXPS = DAT[:EXPS]
        rxns_map = iJR.load_rxns_map()
        
        exglcidx = ChU.rxnindex(model, iJR.EX_GLC_IDER)
        biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
        
        model_ids = [rxns_map[id] for id in FLX_IDERS]
        model_idxs = [ChU.rxnindex(model, model_id) for model_id in model_ids]
        
        
        # Reduce Pol
        p = deepcopy(p0)
        vols = []
        method_diffs_pool = Dict()
        # fba_vals_pool = Dict()

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

            # MSE
            exp_exglc = DAT[:exp, :flx, "GLC", exp] 
            for ider in FLX_IDERS

                # MDiff
                diffs = get!(method_diffs_pool, ider, [])
                me_flx = DAT[me_method, :flx, ider, exp]
                fba_flx = DAT[fba_method, :flx, ider, exp]
                push!(diffs, abs(me_flx - fba_flx) / abs(exp_exglc))

                # # ME
                # me_vals = get!(me_vals_pool, ider, [])
                # me_flx = DAT[me_method, :flx, ider, exp]
                # push!(me_vals, me_flx)
                
                # # FBA
                # fba_vals = get!(fba_vals_pool, ider, [])
                # fba_flx = DAT[fba_method, :flx, ider, exp]
                # push!(fba_vals, fba_flx)
            end
            
        end
        sidxs = sortperm(vols)

        marker = source_markers[Data]
        for _p in [p, all_p]

            for ider in FLX_IDERS
                diffs = method_diffs_pool[ider]
                scatter!(_p, vols[sidxs], diffs[sidxs]; 
                    label = "", m = (8, marker), 
                    color = ider_colors[ider]
                )
                
                # me_vals = me_vals_pool[ider]
                # scatter!(_p, vols[sidxs], me_vals[sidxs]; 
                #     label = "", m = (8, marker), 
                #     color = method_colors[me_method]
                # )
                    
                # fba_vals = fba_vals_pool[ider]
                # scatter!(_p, vols[sidxs], fba_vals[sidxs]; 
                #     label = "", m = (8, marker), 
                #     color = method_colors[fba_method]
                # )
            end
        end
    
        mysavefig(p, "pol_box_vs_method_diff"; src)
    end
    plot!(all_p; xscale = :log10)
    mysavefig(all_p, "pol_box_vs_method_diff")
end  

## ------------------------------------------------------------------------
# vol vs exp
function vol_vs_id(id)
    me_method = :ME_MAX_POL
    fba_method = :FBA_Z_FIX_MIN_COST

    p0 = plot(;xlabel = "politope box volume (log)", ylabel = "$id")
    all_p = deepcopy(p0)
    zeroth = 1e-2
    for (iJR, Data) in [
                (ChK.iJR904, ChK.KayserData), 
                (ChF.iJR904, ChF.FolsomData), 
                (ChH.iJR904, ChH.HeerdenData), 
            ]
        
        src = nameof(Data)

        model = iJR.load_model("max_model")
        datfile = iJR.procdir("dat.bson")
        DAT = UJL.load_data(datfile; verbose = false)
        
        FLX_IDERS = ["GLC", "AC"]
        EXPS = DAT[:EXPS]
        rxns_map = iJR.load_rxns_map()
        
        exglcidx = ChU.rxnindex(model, iJR.EX_GLC_IDER)
        biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
        
        model_ids = [rxns_map[id] for id in FLX_IDERS]
        model_idxs = [ChU.rxnindex(model, model_id) for model_id in model_ids]
        
        
        # Reduce Pol
        p = deepcopy(p0)
        vols = []
        sGLCs = []

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


            # MSE
            sGLC = Data.val(id, exp, 0.0)
            push!(sGLCs, sGLC)
            
        end
        sidxs = sortperm(vols)

        marker = source_markers[Data]
        for _p in [p, all_p]

            scatter!(_p, vols[sidxs], sGLCs[sidxs]; 
                label = "", m = (8, marker), 
                color = method_colors[me_method]
            )
        end
    
        mysavefig(p, "pol_box_vs_$id"; src)
    end

    plot!(all_p; xscale = :log10)
    mysavefig(all_p, "pol_box_vs_$id")
end  

let
    for id in [:sGLC, :D, :sAC, :uGLC, :xi]
        vol_vs_id(id)
    end
end