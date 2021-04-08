## -------------------------------------------------------------------
# Collect
DAT = ChU.DictTree()
let 
    
    DATfile = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, "2.1_DAT.jls")
    # CACHE
    if isfile(DATfile) 
        global DAT = deserialize(DATfile) 
        @info("DAT CACHE LOADED")
        return
    end

    objider = iJR.BIOMASS_IDER
    DAT[:FLX_IDERS] = FLX_IDERS
    DAT[:EXPS] = []

    # Find exps
    for exp in 1:4
        ok = false
        for method in ALL_METHODS
            ok = haskey(INDEX, method, :DFILE, exp) &&
                INDEX[method, :DFILE, exp] != :unfeasible
            !ok && break
        end
        !ok && continue
        push!(DAT[:EXPS], exp)
    end
    max_model = iJR.load_model("max_model"; uncompress = false)

    for exp in DAT[:EXPS], method in ALL_METHODS
            
        # !haskey(INDEX, method, :DFILE, exp) && continue
        datfile = INDEX[method, :DFILE, exp]
        dat = deserialize(datfile)
        
        model = dat[:model]
        objidx = ChU.rxnindex(model, objider)
        epouts = dat[:epouts]
        exp_beta = maximum(keys(epouts))
        epout = epouts[exp_beta]
        exp_xi = Fd.val(:xi, exp)
        fva_model = iJR.load_model("fva_models", exp; uncompress = false)

        println()
        @info("Doing", exp, method, length(dat[:epouts]), epout.iter);

        # Biomass
        ep_biom = ChU.av(model, epout, objidx)
        ep_std = sqrt(ChU.va(model, epout, objidx))
        Fd_biom = Fd.val("D", exp)
        Fd_biom_err = Fd.err("D", exp)
        max_lb, max_ub = ChU.bounds(max_model, objidx)
        fva_lb, fva_ub = ChU.bounds(fva_model, objidx)
        lb = max(max_lb, fva_lb)
        ub = min(max_ub, fva_ub)
        
        # store
        DAT[method, :ep   , :flx, objider, exp] = ep_biom
        DAT[method, :eperr, :flx, objider, exp] = ep_std
        DAT[method, :Fd   , :flx, objider, exp] = Fd_biom
        DAT[method, :Fderr, :flx, objider, exp] = Fd_biom_err
        DAT[:Fd   , :flx, objider, exp] = Fd_biom
        DAT[:Fderr, :flx, objider, exp] = Fd_biom_err
        DAT[method, :bounds, :flx, objider, exp] = (lb, ub)
        
        # fluxes
        for Fd_met in FLX_IDERS

                model_met = Fd_mets_map[Fd_met]
                model_exch = Fd_rxns_map[Fd_met]
                model_exchi = ChU.rxnindex(model, model_exch)

                proj = ChLP.projection2D(model, objider, model_exchi; l = 50)
                ep_av = ChU.av(model, epout, model_exchi)
                ep_std = sqrt(ChU.va(model, epout, model_exchi))
                Fd_flx = Fd.val("u$Fd_met", exp)
                Fd_err = Fd.err("u$Fd_met", exp)

                max_lb, max_ub = ChU.bounds(max_model, Fd_rxns_map[Fd_met])
                fva_lb, fva_ub = ChU.bounds(fva_model, Fd_rxns_map[Fd_met])
                lb = max(max_lb, fva_lb)
                ub = min(max_ub, fva_ub)
                
                DAT[method, :ep, :proj, Fd_met, exp] = proj
                DAT[method, :Fd, :flx, Fd_met, exp] = Fd_flx
                DAT[method, :Fderr, :flx, Fd_met, exp] = Fd_err
                DAT[:Fd, :flx, Fd_met, exp] = Fd_flx
                DAT[:Fderr, :flx, Fd_met, exp] = Fd_err
                DAT[method, :ep, :flx, Fd_met, exp] = ep_av
                DAT[method, :eperr, :flx, Fd_met, exp] = ep_std
                
                DAT[method, :bounds, :flx, Fd_met, exp] = (lb, ub)

        end # for Fd_met
    end # for exp in EXPS, for method

    DAT[:EXPS] |> unique! |> sort!
    serialize(DATfile, DAT)
end;

# -------------------------------------------------------------------
# Inter project comunication
let
    CORR_DAT = isfile(iJR.CORR_DAT_FILE) ? ChU.load_data(iJR.CORR_DAT_FILE) : Dict()
    CORR_DAT[:MAXENT_EP] = DAT
    ChU.save_data(iJR.CORR_DAT_FILE, CORR_DAT)
end
