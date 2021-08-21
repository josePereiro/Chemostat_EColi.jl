## -------------------------------------------------------------------
# TODO: work better inter project comunication
# Full corr FBA vs ME
let
    
    DAT_FILE_PREFFIX =  "maxent_ep_dat"
    me_method = :ME_MAX_POL
    fba_method = :FBA_Z_FIX_MIN_COST

    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData),
        (ChN.iJR904, ChN.NanchenData), 
        (ChF.iJR904, ChF.FolsomData),
        (ChH.iJR904, ChH.HeerdenData),
    ]
        src = nameof(Data)
        marker = (8, source_markers[Data])

        datfile = iJR.procdir("dat.bson")
        DAT = UJL.load_data(datfile; verbose = false)

        # Preload FBA Data
        LP_DAT = ChE.load_LP_DAT(src)

        EXPS = DAT[:EXPS]
        ME_STOI_ERRS = []
        FBA_STOI_ERRS = []
        for exp in EXPS

            # norm
            exp_glc = abs(DAT[:exp, :flx, "GLC", exp])

            # ME
            ME_DAT = UJL.mysavename(DAT_FILE_PREFFIX, "jls"; 
                method = me_method, exp) |> iJR.procdir |> deserialize
            model = ME_DAT[:model]
            epouts = ME_DAT[:epouts]
            exp_beta = maximum(keys(epouts))
            epout = epouts[exp_beta]
            
            ERRS = ChU.stoi_err(model, epout) ./ exp_glc .* 100.0
            push!(ME_STOI_ERRS, ERRS)

            # FBA
            model = LP_DAT[fba_method, :model, exp]
            fbaout = LP_DAT[fba_method, :fbaout, exp]

            ERRS = ChU.stoi_err(model, fbaout) ./ exp_glc .* 100.0
            push!(FBA_STOI_ERRS, ERRS)
        end

        p = plot(;title = string(src), 
            xlabel = "experiments", ylabel = "Stoi Err (% exp_glc)"
        )
        
        for (method, ERRS) in [
            (me_method, ME_STOI_ERRS), 
            (fba_method, FBA_STOI_ERRS)
        ]
            color = method_colors[method]
            scatter!(p, EXPS, mean.(ERRS); label = string(method),
                fillrange = [minimum.(ERRS), maximum.(ERRS)],
                fillcolor = color,
                fillalpha = 0.3,
                alpha = 0.5, color, m = 6, 
            )
        end

        for (method, ERRS) in [
            (me_method, ME_STOI_ERRS), 
            (fba_method, FBA_STOI_ERRS)
        ]
            color = method_colors[method]
            plot!(p, EXPS, mean.(ERRS); label = "",
                lw = 2, color, ls = :dash, alpha = 0.5, 
                xticks = (EXPS, string.(EXPS))
            )
        end

        mysavefig(p, "stoi_err"; src, me_method, fba_method)
    end
end