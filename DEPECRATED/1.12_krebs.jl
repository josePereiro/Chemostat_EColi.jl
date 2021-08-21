## -------------------------------------------------------------------
# krebs per ider per exp per src
let

    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData), 
        (ChN.iJR904, ChN.NanchenData), 
        (ChF.iJR904, ChF.FolsomData), 
        (ChH.iJR904, ChH.HeerdenData), 
    ]

        src = nameof(Data)

        datfile = iJR.procdir("dat.bson")
        DAT = UJL.load_data(datfile; verbose = false)

        METHODS = [:FBA_Z_FIX_MIN_COST, :ME_MAX_POL]
        EXPS = DAT[:EXPS]
        EXCH_FLX_IDERS = iJR.load_krebs_iders()

        for exp in EXPS
            p = plot(;xlabel = "flx ider", ylabel = "flx val / exp glcex")
            glc_exp_val = abs(DAT[:exp, :flx, "GLC", exp])

            for method in METHODS
                model_vals = abs.(DAT[method, :flx, EXCH_FLX_IDERS, exp]) ./ glc_exp_val

                color = method_colors[method]
                scatter!(p, EXCH_FLX_IDERS, model_vals; 
                    color, label = string(method), 
                    m = 8
                )
                plot!(p, EXCH_FLX_IDERS, model_vals; label = "",
                    color, ls = :dash, alpha = 0.5, lw = 3
                )
            end

            plot!(;xrotation = 35, ylim = [0.0, Inf])
            mysavefig(p, "krebs"; src, exp)
        end
    end
end


## -------------------------------------------------------------------
# krebs avs
let

    fba_method = :FBA_Z_FIX_MIN_COST
    me_method = :ME_MAX_POL
    EXCH_FLX_IDERS = KiJR.krebs_iders

    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData), 
        (ChN.iJR904, ChN.NanchenData), 
        (ChF.iJR904, ChF.FolsomData), 
        (ChH.iJR904, ChH.HeerdenData), 
    ]
        me_avs, me_stds = [], []
        fba_avs, fba_stds = [], []
        src = nameof(Data)
        for ider in EXCH_FLX_IDERS

            datfile = iJR.procdir("dat.bson")
            DAT = UJL.load_data(datfile; verbose = false)

            EXPS = DAT[:EXPS]

            glc_exp_vals = abs.(DAT[:exp, :flx, "GLC", EXPS])
            fba_vals = abs.(DAT[fba_method, :flx, ider, EXPS]) ./ glc_exp_vals
            me_vals = abs.(DAT[me_method, :flx, ider, EXPS]) ./ glc_exp_vals
            push!(fba_avs, mean(fba_vals)); push!(fba_stds, std(fba_vals))
            push!(me_avs, mean(me_vals)); push!(me_stds, std(me_vals))
        end

        p = plot(;title = string(src), xlabel = "flx ider", ylabel = "flx val / exp glcex")
       
        for (avs, stds, method) in [
            (me_avs, me_stds, me_method),
            (fba_avs, fba_stds, fba_method),
        ]
            color = method_colors[method]
            lb, ub = avs .- stds, avs .+ stds
            lb, ub = avs .- stds, avs .+ stds
            scatter!(p, EXCH_FLX_IDERS, [avs avs]; label = "", 
                fillrange = [lb ub], fillalpha = 0.2, color,
                markeralpha = 0.0
            )
        end

        marker = (8, source_markers[Data])
        for (avs, method) in [
            (me_avs, me_method),
            (fba_avs, fba_method),
        ]
            color = method_colors[method]
            scatter!(p, EXCH_FLX_IDERS, avs; label = string(method), 
                color, marker
            )
            plot!(p, EXCH_FLX_IDERS, avs; label = "", 
                color, lw = 1, alpha = 0.3, ls = :dash
            )
        end
       
        plot!(p; xrotation = 35, ylim = [0.0, Inf])
        mysavefig(p, "krebs_ave_comparison"; src)
    end
end

  