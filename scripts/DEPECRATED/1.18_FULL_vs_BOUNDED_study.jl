MEname(c) = "\$\\textbf{ME}^{\\langle \\textbf{$c} \\rangle}\$"

## ----------------------------------------------------------------------------
# legend
let
    me_name1 = "\$MEP^2\$" #MEname(1)
    me_name2 = "\$MEP^1\$" #MEname(0)

    p = plot()
    params = (;lw = 15, legendfont = 15, legendlinewidth = 15, dpi=1000)
    plot!(p, fill(2, 10); label = me_name2, color = :blue, params...)
    plot!(p, fill(4, 10); label = me_name1, color = :red, params...)
    plot!(p, fill(6, 10); label = "", alpha = 0.0)
    mysavefig(p, "ME_legend")
end

## ----------------------------------------------------------------------------
# exchanges ave var correlations
let
    DAT_FILE_PREFFIX =  "maxent_ep_dat"
    me_method1 = ME_MAX_POL
    me_method2 = ME_Z_EXPECTED_G_BOUNDED
    me_name1 = "on_PX"
    me_name2 = "on_V"
    me_name1 = "\$MEP^2\$" #MEname(1)
    me_name2 = "\$MEP^1\$" #MEname(0)

    common_params = (;thickness_scaling = 1.3, 
        # dpi=1000
    )
    
    diff_th = 0.2
    diff_va_p0 = plot(; title = "Flux variance", 
        xlabel = "exp", ylabel = "biased flux fraction", common_params...
    )
    diff_va_p = deepcopy(diff_va_p0)
    
    fontsize = 13
    corr_av_p0 = plot(;title = "Flux average", common_params...,
        xlabel = me_name1, ylabel = me_name2, 
        xguidefontsize = fontsize, yguidefontsize = fontsize
    )
    corr_va_p0 = plot(;title = "Flux variance (log scale)", common_params...,
        xlabel = me_name1, ylabel = me_name2,
        xguidefontsize = fontsize, yguidefontsize = fontsize
    )
    corr_av_p, corr_va_p = deepcopy.([corr_av_p0, corr_va_p0])
    corr_ref_params = (;label = "", lw = 3, alpha = 0.8, ls = :dash)
    
    vas, avs = [], []
    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData),
        (ChN.iJR904, ChN.NanchenData), 
        (ChF.iJR904, ChF.FolsomData),
        # (ChH.iJR904, ChH.HeerdenData),
    ]
        src = string(nameof(Data))
        marker = (5, source_markers[Data])

        corr_lav_p, corr_lva_p = deepcopy.([corr_av_p0, corr_va_p0])
        diff_lav_p =  deepcopy(diff_va_p0)

        # load resume
        DAT = ChE.load_DAT(src)
        EXPS = DAT[:EXPS]

        lavs, lvas = [], []
        counts_ = []
        for exp in EXPS

            # norm
            exp_glc = abs(DAT[:exp, :flx, "GLC", exp])

            ## -------------------------------------------------------------
            # Common iders
            ME1_DAT = UJL.mysavename(DAT_FILE_PREFFIX, "jls"; 
                method = me_method1, exp) |> iJR.procdir |> deserialize
            model1 = ME1_DAT[:model]
            ME2_DAT = UJL.mysavename(DAT_FILE_PREFFIX, "jls"; 
                method = me_method2, exp) |> iJR.procdir |> deserialize
            model2 = ME2_DAT[:model]

            iders = intersect(model1.rxns, model2.rxns)

            ## -------------------------------------------------------------
            # method1
            epouts = ME1_DAT[:epouts]
            exp_beta = maximum(keys(epouts))
            epout = epouts[exp_beta]
            
            me1_av = ChU.av(model1, epout, iders)
            me1_va = ChU.va(model1, epout, iders)

            ## -------------------------------------------------------------
            # method2
            epouts = ME2_DAT[:epouts]
            exp_beta = maximum(keys(epouts))
            epout = epouts[exp_beta]

            me2_av = ChU.av(model2, epout, iders)
            me2_va = ChU.va(model2, epout, iders)

            ## -------------------------------------------------------------
            for p in [corr_lav_p, corr_av_p]
                scatter!(p, me1_av, me2_av; label = "", marker, 
                    color = :black, alpha = 0.5
                )
            end

            for p in [corr_lva_p, corr_va_p]
                scatter!(p, log.(me1_va), log.(me2_va); label = "", marker, 
                    color = :black, alpha = 0.5
                )
            end

            push!(lavs, me1_av..., me2_av...)
            push!(lvas, me1_va..., me2_va...)
            sort!(lavs); sort!(lvas)

            ## -------------------------------------------------------------
            count_ = count(me2_va ./ me1_va) do r
                r < (1.0 - diff_th)
            end
            push!(counts_, count_/ length(iders))
        end

        ## -------------------------------------------------------------
        for p in [diff_va_p, diff_lav_p]
            # diffs = 1.0 .- me2_va ./ me1_va
            scatter!(p, EXPS, counts_; 
                label = "", marker, 
                color = :black, alpha = 0.8
            )
            plot!(p, EXPS, counts_; 
                label = "", lw = 2, ls = :dash, 
                color = :black, alpha = 0.7
            )
        end
        
        push!(avs, lavs...)
        push!(vas, lvas...)
        sort!(avs); sort!(vas)

        # local plots
        plot!(corr_lav_p, lavs, lavs; corr_ref_params...)
        plot!(corr_lva_p, log.(lvas), log.(lvas); corr_ref_params...)

        # mysavefig([corr_lav_p, corr_lva_p], "corrs_me_max_vs_bounded"; src)
        # mysavefig(diff_lav_p, "diff_me_max_vs_bounded"; src)
    end
    
    plot!(corr_av_p, avs, avs; corr_ref_params...)
    plot!(corr_va_p, log.(vas), log.(vas); corr_ref_params...)
    
    mysavefig([corr_av_p, corr_va_p], "corrs_me_max_vs_bounded")
    mysavefig(diff_va_p, "diff_me_max_vs_bounded")
end 

## ----------------------------------------------------------------------------
# glc, biom marginals
let
    DAT_FILE_PREFFIX =  "maxent_ep_dat"

    me_method1 = ME_MAX_POL
    me_method2 = ME_Z_EXPECTED_G_BOUNDED
    me_method3 = ME_MAX_POL_B0
    
    fontsize = 16
    common_params = (;
        thickness_scaling = 1.6, 
        xguidefontsize = fontsize, yguidefontsize = fontsize,
        # dpi = 1000
    )

    _method_pair(iJR, method; kwargs...) = method => UJL.mysavename(
        DAT_FILE_PREFFIX, "jls"; method, kwargs...
    ) 

    datnames(iJR, exp) = Dict(
        _method_pair(iJR, me_method1; exp),
        _method_pair(iJR, me_method2; exp),
        # _method_pair(me_method3),
    )
    allfound = false

    colors = Dict(
        me_method1 => :red,
        me_method2 => :blue,
        me_method3 => :black,
    )

    for (iJR, Data, EXPS) in [
        (ChK.iJR904, ChK.KayserData, 5:5),
        (ChN.iJR904, ChN.NanchenData, 4:4), 
        (ChF.iJR904, ChF.FolsomData, 1:1),
        (ChH.iJR904, ChH.HeerdenData, 4:4),
    ]

        src = nameof(Data)
        rxns_map = iJR.load_rxns_map()
        rxns_map["cost"] = iJR.COST_IDER

        for exp in EXPS
            @info("Doing", src, exp)

            ps = Plots.Plot[]
            for exp_ider in [Data.msd_mets; "D"; "cost"]
                model_ider = rxns_map[exp_ider]

                units = exp_ider == "D" ? "1/ h" : "mmol/ gCDW h"
                title = exp_ider == "D" ? "z" : exp_ider
                p = plot(;title, common_params...,
                    xlabel = "Flux", 
                    ylabel = "norm. pdf"
                )
                avs, vas = [], []
                for (method, datname) in datnames(iJR, exp)
                    datfile = iJR.procdir(datname)
                    allfound = isfile(datfile)
                    !allfound && break

                    ME_DAT = datfile |> deserialize
                    model = ME_DAT[:model]
                    epouts = ME_DAT[:epouts]
                    exp_beta = maximum(keys(epouts))
                    epout = epouts[exp_beta]

                    ChP.plot_marginal!(p, model, epout, model_ider;
                        color = colors[method], label = "", lw = 4,
                        normalize = true
                    )

                    av = ChU.av(model, epout, model_ider)
                    va = ChU.va(model, epout, model_ider)
                    push!(avs, av); push!(vas, va) 
                end
                !allfound && break

                xlim = [minimum(avs) - 5 * sqrt(maximum(vas)), 
                    maximum(avs) + 5 * sqrt(maximum(vas))]
                xticks = range(xlim...; length = 4)
                xticks = round.(xticks; sigdigits = 2)

                plot!(p; xlim, xticks)
                push!(ps, p)
                # mysavefig(p, "exchange_marginals"; src, exp, exp_ider)
            end
            (isempty(ps) || !allfound) && continue
            mysavefig(ps, "exchange_marginals"; src, exp)
        end
    end

end