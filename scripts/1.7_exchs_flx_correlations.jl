## -------------------------------------------------------------------
# X vs D
let
    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData), 
        (ChN.iJR904, ChN.NanchenData), 
        (ChF.iJR904, ChF.FolsomData), 
        # (ChH.iJR904, ChH.HeerdenData), 
    ]
        src = nameof(Data)
        xid = :D
        yid = :X
        p = plot(;xlabel = string(xid), ylabel = string(yid))

        xs = Data.val(xid)
        sidxs = sortperm(xs)
        ys = try; Data.val(yid) .|> abs catch; zero(xs) end
        plot!(p, xs[sidxs], ys[sidxs]; 
            label = "", ls = :dash, color = :black, lw = 3
        )
        scatter!(p, xs[sidxs], ys[sidxs]; 
            label = "", m = 7, color = :black, 
            ylim = [0.0, Inf]
        )
        mysavefig(p, "exp_corrs"; src)
    end
end

## -------------------------------------------------------------------
# EP biomass corr
let
    METHODS = [
        # FBA_Z_FIX_MIN_COST,
        FBA_MAX_Z_MIN_COST, 
        # FBA_Z_FIX_MAX_COST, 
        # FBA_Z_FIX_MIN_VG_COST, 
        # FBA_Z_VG_FIX_MIN_COST,
        # FBA_Z_FIX_MAX_VG_MIN_COST,
        
        # ME_MAX_POL,
        # ME_Z_EXPECTED_G_BOUNDED
    ]

    common_params = (;thickness_scaling = 1.3, dpi = 1000)
    
    for method in METHODS
        p = plot(; 
            xlabel = "\$ \\textrm{exp. D (1/ h)} \$", 
            ylabel = "\$ \\textrm{model max. } {\\langle z \\rangle}_{P_X} \\textrm{ (1/ h)} \$", 
            common_params...
        )

        all_vals = []
        for (iJR, Data) in [
            (ChK.iJR904, ChK.KayserData), 
            (ChN.iJR904, ChN.NanchenData), 
            (ChF.iJR904, ChF.FolsomData), 
            # (ChH.iJR904, ChH.HeerdenData), 
        ]
            src = nameof(Data)
            marker = (6, source_markers[Data])

            datfile = iJR.procdir("dat.bson")
            DAT = UJL.load_data(datfile; verbose = false)

            EXPS = DAT[:EXPS]
        
            @info("Doing", src, method)
            
            model_vals = DAT[method, :flx, "D", EXPS]
            model_errs = haskey(DAT[method], :err) ? 
                DAT[method, :err, "D", EXPS] : fill(0.0, length(EXPS))

            exp_vals = DAT[:exp, :flx, "D", EXPS]
            # color = [exp_colors[exp] for exp in EXPS]
            m, M = myminmax([exp_vals; model_vals])
            margin = abs(M - m) * 0.1
            scatter!(p, exp_vals, model_vals; 
                # yerr = model_errs,
                label = string(source_labels[Data]), color = :black,
                marker, alpha = 0.7,
                xlim = [m - margin, M + margin],
                ylim = [m - margin, M + margin]
            )
            
            push!(all_vals, exp_vals..., model_vals...)
        end

        sort!(all_vals)
        plot!(p, all_vals, all_vals; label = "", lw = 3, ls = :dash, alpha = 0.7, color = :black)
        plot!(;legend = :bottomright)
        mysavefig(p, "biom_corr"; method)
    end
end

## -------------------------------------------------------------------
# flx correlations
let

    # TODO: make this global
    MEname(c) = "\$\\textbf{ME}^{\\langle \\textbf{$c} \\rangle} \\textbf{ (abs. flux)} \$"

    METHODS = [

        FBA_Z_FIX_MAX_VG_MIN_COST, 
        FBA_Z_FIX_MIN_VG_MIN_COST, 
        FBA_Z_FIX_MAX_VG_MAX_COST, 
        FBA_Z_FIX_MIN_VG_MAX_COST,
        FBA_Z_VG_FIX_MAX_COST, 
        FBA_Z_VG_FIX_MIN_COST, 
        FBA_Z_FIX_MAX_COST, 
        FBA_Z_FIX_MIN_COST, 
        FBA_MAX_Z_MIN_COST, 
        FBA_MAX_Z_MAX_COST,
        
        ME_MAX_POL,
        # ME_Z_EXPECTED_G_BOUNDED
    ]

    # MEname(c) = "\$\\textbf{ME}^{\\langle \\textbf{$c} \\rangle} \\textbf{ (abs. flux)} \$"

    METHODS_LABELS = [
        # "\$ \\textbf{FBA} \$",
        "FBA_Z_FIX_MAX_VG_MIN_COST", 
        FBA_Z_FIX_MIN_VG_MIN_COST, 
        "FBA_Z_FIX_MAX_VG_MAX_COST", 
        "\$ \\textbf{FBA} \\textbf{ (abs. flux)} \$",
        "FBA_Z_VG_FIX_MAX_COST", 
        "FBA_Z_VG_FIX_MIN_COST", 
        "FBA_Z_FIX_MAX_COST", 
        "FBA_Z_FIX_MIN_COST", 
        "FBA_MAX_Z_MIN_COST", 
        "FBA_MAX_Z_MAX_COST",

        MEname(1), # "ME unbiased",
        # MEname(0)
    ] 

    for (get_flxs, flxs) in [
        ((DAT) -> DAT[:FLX_IDERS], "all"),
        # ((DAT) -> ["GLC", "AC"], "common"),
    ]
        ps_pool = Dict()
        fontsize = 13
        p0 = plot(;
            title = "Exchanges", 
            dpi = 1000,
            thickness_scaling = 1.3, 
            xguidefontsize = fontsize, yguidefontsize = fontsize
        )

        for (iJR, Data) in [
            (ChK.iJR904, ChK.KayserData),
            # (ChN.iJR904, ChN.NanchenData), 
            # (ChF.iJR904, ChF.FolsomData),
            # (ChH.iJR904, ChH.HeerdenData),
        ]
            src = nameof(Data)
            marker = (6, source_markers[Data])

            @show src
            DAT = ChE.load_DAT(src)

            # FLX_IDERS = DAT[:FLX_IDERS]
            # FLX_IDERS = ["GLC", "AC"]
            FLX_IDERS = get_flxs(DAT)
            EXPS = DAT[:EXPS]

            tot_ps = Plots.Plot[]
            for (method_label, method) in zip(METHODS_LABELS, METHODS)
                mp0 = plot!(deepcopy(p0); 
                    xlabel = "\$ \\textbf{exp. (abs. flux)}\$", ylabel = method_label
                )    
                # total corr
                let            
                    model_vals = DAT[method, :flx, FLX_IDERS, EXPS] .|> abs
                    model_errs = haskey(DAT[method], :err) ? 
                        DAT[method, :err, FLX_IDERS, EXPS] .|> abs : 
                        zeros(prod(length.([FLX_IDERS, EXPS])))
                    exp_vals = DAT[:exp, :flx, FLX_IDERS, EXPS] .|> abs
                    
                    color = [ider_colors[ider] for ider in FLX_IDERS, exp in EXPS]
                    # ep corr
                    local_p = plot!(deepcopy(p0); xlabel = "\$ \\textbf{exp.}\$", ylabel = method_label)
                    top_p = get!(ps_pool, (;method), deepcopy(mp0))

                    for p in [local_p, top_p]
                        scatter!(p, exp_vals, model_vals; yerr = model_errs, 
                            marker, label = "", color, alpha = 0.7
                        )
                        all_vals = [model_vals; exp_vals] |> sort!
                        plot!(p, all_vals, all_vals; 
                            ls = :solid, color = :black, label = ""
                        )
                    end
                    # mysavefig(local_p, "exchs_corr"; src, method, flxs)
                end

                # per ider
                let       
                    for ider in FLX_IDERS
                        model_vals = DAT[method, :flx, ider, EXPS] .|> abs
                        model_errs = haskey(DAT[method], :err) ? 
                            DAT[method, :err, ider, EXPS] .|> abs : 0.0
                        exp_vals = DAT[:exp, :flx, ider, EXPS]  .|> abs
                        
                        color = ider_colors[ider]
                        # ep corr
                        local_p = plot!(deepcopy(p0); title = "$(src): $(method_label)")
                        top_p = get!(ps_pool, (;ider, method), deepcopy(mp0))

                        for p in [local_p, top_p]
                            scatter!(p, exp_vals, model_vals; yerr = model_errs, 
                                marker, label = "", color, alpha = 0.7
                            )
                            all_vals = [model_vals; exp_vals] |> sort!
                            plot!(p, all_vals, all_vals; 
                                ls = :solid, color = :black, label = ""
                            )
                        end
                        # mysavefig(local_p, "exchs_corr"; src, ider, method, flxs)
                    end
                end
            end
        end

        for (kws, p) in ps_pool
            length(kws) > 1 && continue # store just globals
            mysavefig(p, "exchs_corr"; kws..., flxs)
        end
    end
end 

## -------------------------------------------------------------------
# Full corr FBA vs ME
let
    
    DAT_FILE_PREFFIX =  "maxent_ep_dat"
    me_method = :ME_MAX_POL
    fba_method = :FBA_Z_FIX_MIN_COST
    base = 1e-5
    function f(x) 
        x += base
        s = sign(x)
        return s * log(abs(x))
    end

    p0 = plot(;xlabel = "ME flxs (log)", ylabel = "FBA flxs (log)", 
        thickness_scaling = 1.3
    )
    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData),
        (ChN.iJR904, ChN.NanchenData), 
        (ChF.iJR904, ChF.FolsomData),
        (ChH.iJR904, ChH.HeerdenData),
    ]
        src = nameof(Data)
        marker = (8, source_markers[Data])

        DAT = ChE.load_DAT(src)

        # Preload FBA Data
        LP_DAT = ChE.load_LP_DAT(src)

        EXPS = DAT[:EXPS]
        src_p = deepcopy(p0)
        for exp in EXPS

            # load ME dat
            ME_DAT = UJL.mysavename(DAT_FILE_PREFFIX, "jls"; 
                method = me_method, exp
            ) |> iJR.procdir |> deserialize

            # common rxns
            me_model = ME_DAT[:model]
            fba_model = LP_DAT[fba_method, :model, exp]
            RXNS = intersect(me_model.rxns, fba_model.rxns)

            # ME flxs
            epouts = ME_DAT[:epouts]
            exp_beta = maximum(keys(epouts))
            epout = epouts[exp_beta]
            me_flxs = f.(ChU.av.([me_model], [epout], RXNS))

            # FBA
            fbaout = LP_DAT[fba_method, :fbaout, exp]
            fba_flxs = f.(ChU.av.([fba_model], [fbaout], RXNS))

            local_p = deepcopy(p0)
            for p in [local_p, src_p]
                scatter!(p, me_flxs, fba_flxs; label = "",
                    alpha = 0.5, color = :black, m = 6
                )
                all_vals = [me_flxs; fba_flxs] |> sort!
                plot!(p, all_vals, all_vals; label = "",
                    lw = 2, color = :black
                )
            end
            # mysavefig(local_p, "fba_vs_me_corr"; src, exp)
        end

        mysavefig(src_p, "fba_vs_me_corr"; src)
    end
end