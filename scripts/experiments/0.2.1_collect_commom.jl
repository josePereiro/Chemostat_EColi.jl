function collect_commom_DAT!(DAT, Ed, iJR)

    # maps
    # ----------------------------------------------------------------------------------
    
    # ----------------------------------------------------------------------------------
    max_model = iJR.load_model("max_model"; uncompress = false)
    objider = iJR.BIOMASS_IDER
    
    # ----------------------------------------------------------------------------------
    for exp in DAT[:EXPS]
        # exp dat
        Ed_biom = Ed.val("D", exp)
        DAT[:exp, :flx, "D", exp] = Ed_biom
        DAT[:exp, :err, "D", exp] = 0.0

        # bounds
        fva_model = iJR.load_model("fva_models", exp; uncompress = false)
        max_lb, max_ub = ChU.bounds(max_model, objider)
        fva_lb, fva_ub = ChU.bounds(fva_model, objider)
        lb = max(max_lb, fva_lb)
        ub = min(max_ub, fva_ub)
        DAT[:bounds, "D", exp] = (lb, ub)

        # exchanges
        for Ed_ider in DAT[:EXCH_FLX_IDERS]

            # exp dat
            Ed_flx = Ed.uval(Ed_ider, exp)
            Ed_err = try; Ed.uerr(Ed_ider, exp) catch; 0.0 end
            DAT[:exp, :flx, Ed_ider, exp] = Ed_flx
            DAT[:exp, :err, Ed_ider, exp] = Ed_err

        end

        # inner
        ider_map = iJR.load_rxns_map()
        for model_ider in DAT[:INNER_FLX_IDERS]

            Ed_ider = ider_map[model_ider]

            # exp dat
            Ed_flx = Ed.val(Ed_ider, exp)
            Ed_err = try; Ed.err(Ed_ider, exp) catch; 0.0 end
            DAT[:exp, :flx, Ed_ider, exp] = Ed_flx
            DAT[:exp, :err, Ed_ider, exp] = Ed_err

        end
    end

    return DAT
end