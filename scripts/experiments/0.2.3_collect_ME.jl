## ----------------------------------------------------------------------------------
function collect_ME_DAT!(DAT, iJR; nths = nthreads())

    project = string(parentmodule(iJR))

    costider = iJR.COST_IDER
    objider = iJR.BIOMASS_IDER

    rxns_map = iJR.load_rxns_map()

    # util stuff
    WLOCK = ReentrantLock()
    isexpdep(method) = !(method in DAT[:EXP_INDEPENDENT])
    depks(method, typ, Ed_met, exp) = 
        isexpdep(method) ? (method, typ, Ed_met, exp) : (method, typ, Ed_met)

    function dat_file(;method, exp)
        kwargs = isexpdep(method) ? (;method, exp) : (;method)
        procdir(iJR, "maxent_ep_dat", kwargs, ".jls")
    end

    # Feed jobs
    Ch = Channel(nths) do ch
        for method in DAT[:ME_METHODS]
            for exp in DAT[:EXPS]
                put!(ch, (exp, method))
                !isexpdep(method) && break # Do only once
            end
        end
    end

    @threads for thid in 1:nths
        thid = threadid()
        for (exp, method) in Ch
            
            # ME data
            datfile = dat_file(;method, exp)
            # !isfile(datfile) && continue
            dat = ldat(datfile)
            
            model = dat[:model]
            epouts = dat[:epouts]
            exp_beta = dat[:exp_beta] # maximum(keys(epouts)) 
            epout = epouts[exp_beta]
            
            lock(WLOCK) do
                @info("Doing",
                    project,
                    exp, method, 
                    length(epouts), 
                    epout.iter, thid
                ); println()
            end

            # Biomass and cost
            ep_biom = ChU.av(model, epout, objider)
            ep_biomstd = sqrt(ChU.va(model, epout, objider))
            ep_cost = ChU.av(model, epout, costider)
            ep_coststd = sqrt(ChU.va(model, epout, costider))
            
            # store
            lock(WLOCK) do
                DAT[depks(method, :flx, "D", exp)...] = ep_biom
                DAT[depks(method, :err, "D", exp)...] = ep_biomstd
                DAT[depks(method, :flx, "cost", exp)...] = ep_cost
                DAT[depks(method, :err, "cost", exp)...] = ep_coststd
            end

            # exchanges
            for Ed_met in DAT[:EXCH_FLX_IDERS]
                model_exch = rxns_map[Ed_met]

                # flxs
                ep_av = ChU.av(model, epout, model_exch)
                ep_std = sqrt(ChU.va(model, epout, model_exch))
                        
                lock(WLOCK) do
                    DAT[depks(method, :flx, Ed_met, exp)...] = ep_av
                    DAT[depks(method, :err, Ed_met, exp)...] = ep_std
                end
            end

            # inner flxs
            idermap = iJR.load_inners_idermap()
            for exider in DAT[:INNER_FLX_IDERS]
                
                model_iders = idermap[exider]

                # flxs
                ep_av = ChU.av(model, epout, model_iders[1])
                ep_std = sqrt(ChU.va(model, epout, model_iders[1]))
                if length(model_iders) == 2 # reversible
                    # r = r+ - r-
                    ep_av -= ChU.av(model, epout, model_iders[2])
                    # ep_std += sqrt(ChU.va(model, epout, model_iders[2]))
                    ep_std = NaN
                end
                
                lock(WLOCK) do
                    DAT[depks(method, :flx, exider, exp)...] = ep_av
                    DAT[depks(method, :err, exider, exp)...] = ep_std
                end
            end
        
        end # for (exp, method)
    end # for thid

    return DAT
end