## ----------------------------------------------------------------------------------
function collect_lp_DAT!(DAT, Ed, iJR)

    proj = string(parentmodule(iJR))
    src = string(nameof(Ed))

    rxns_map = iJR.load_rxns_map()

    objider = iJR.BIOMASS_IDER
    costider = iJR.COST_IDER

    LP_DAT = ChE.load_LP_DAT(src)

    for method in DAT[:LP_METHODS]
            
        for exp in DAT[:EXPS]

            @info("Doing", 
                proj, exp, method
            ); println()

            model = LP_DAT[method, :model, exp]
            fbaout = LP_DAT[method, :fbaout, exp]

            # Biomass
            fba_biom = ChU.av(model, fbaout, objider)
            DAT[method, :flx, "D", exp] = fba_biom
            DAT[method, :err, "D", exp] = 0.0

            fba_cost = ChU.av(model, fbaout, costider)
            DAT[method, :flx, "cost", exp] = fba_cost
            DAT[method, :err, "cost", exp] = 0.0

            for Ed_ider in DAT[:EXCH_FLX_IDERS]
                model_ider = rxns_map[Ed_ider]

                fba_flx = ChU.av(model, fbaout, model_ider)
                DAT[method, :flx, Ed_ider, exp] = fba_flx
                DAT[method, :err, Ed_ider, exp] = 0.0
            end

            idermap = merge(iJR.load_inners_idermap(), iJR.load_inners_idermap())
            for (exglcider, model_iders) in idermap
                # flxs
                fba_flx = ChU.av(model, fbaout, model_iders[1])
                if length(model_iders) == 2 # reversible
                    # r = r+ - r-
                    fba_flx -= ChU.av(model, fbaout, model_iders[2])
                end
                        
                DAT[method, :flx, exglcider, exp] = fba_flx
                DAT[method, :err, exglcider, exp] = 0.0
            end

        end

    end # for method

    return DAT
end