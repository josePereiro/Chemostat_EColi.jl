
## ----------------------------------------------------------------------------
_load_ME_dat(iJR, method, exp) = ldat(iJR, "maxent_ep_dat", (;method, exp), ".jls")

## ----------------------------------------------------------------------------
# exchanges ave var correlations
function plot_ME_methods_comparizon()

    me_method1 = :ME_MAX_POL
    me_method2 = :ME_Z_EXPECTED_G_BOUNDED
    me_name1 = "ME^2" 
    me_name2 = "ME^1"
    
    corr_av_p0 = plot(;
        # title = _textbf("Flux average"),
        xlabel = _textbf(me_name1, " flux ave."), 
        ylabel = _textbf(me_name2, " flux ave."), 
    )
    corr_va_p0 = plot(;
        # title = _textbf("Flux variance"), 
        xlabel = _textbf(me_name1, " flux var. (log scale)"), 
        ylabel = _textbf(me_name2, " flux var. (log scale)"), 
    )
    corr_av_p, corr_va_p = deepcopy.([corr_av_p0, corr_va_p0])
    corr_ref_params = (;label = "", lw = 3, alpha = 0.8, ls = :dash)
    
    vas, avs = [], []
    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData),
        (ChN.iJR904, ChN.NanchenData), 
        (ChF.iJR904, ChF.FolsomData),
    ]
        src = string(nameof(Data))
        marker = (8, source_markers[Data])

        # load resume
        DAT = ChE.load_DAT(src)
        EXPS = DAT[:EXPS]

        lavs, lvas = [], []
        for exp in EXPS

            ## -------------------------------------------------------------
            # Common iders
            
            ME1_DAT = _load_ME_dat(iJR, me_method1, exp)
            model1 = ME1_DAT[:model]
            ME2_DAT = _load_ME_dat(iJR, me_method2, exp)
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
            scatter!(corr_av_p, me1_av, me2_av; label = "", marker, 
                color = :black, alpha = 0.5
            )

            scatter!(corr_va_p, log.(me1_va), log.(me2_va); label = "", marker, 
                color = :black, alpha = 0.5
            )

            push!(lavs, me1_av..., me2_av...)
            push!(lvas, me1_va..., me2_va...)
            sort!(lavs); sort!(lvas)
        end

        ## -------------------------------------------------------------
        push!(avs, lavs...)
        push!(vas, lvas...)
        sort!(avs); sort!(vas)

    end
    
    plot!(corr_av_p, avs, avs; corr_ref_params...)
    plot!(corr_va_p, log.(vas), log.(vas); corr_ref_params...)
    
    sfig(ChE, [corr_av_p, corr_va_p], 
        @fileid, "corrs_me_max_vs_bounded", ".png", 
        layout = (1, 2)
    )
    
    return corr_av_p, corr_va_p
end 

## ----------------------------------------------------------------------------
# glc, biom marginals
function plot_me_marginals()

    iJR = ChN.iJR904
    Data = ChN.NanchenData
    exp = 4
    src = nameof(Data)

    me_method1 = :ME_MAX_POL
    me_method2 = :ME_Z_EXPECTED_G_BOUNDED

    colors = Dict(
        me_method1 => :red,
        me_method2 => :blue,
    )
    
    labels = Dict(
        me_method1 => "ME^1",
        me_method2 => "ME^2",
    )
    
    rxns_map = iJR.load_rxns_map()

    @info("Doing", src, exp)

    ps = Plots.Plot[]
    for (exp_ider, xlim) in [
            ("GLC", (-2.6, -0.8)), 
            ("AC", (-0.015, 0.075)), 
            ("D", (0.0, 0.2)),
        ]

        model_ider = rxns_map[exp_ider]

        p = plot(; 
            xlabel = _textbf(exp_ider == "D" ? "z" : exp_ider),
            ylabel = _textbf("norm. pdf")
        )
        avs, vas = [], []

        for method in [me_method1, me_method2]
            
            ME_DAT = _load_ME_dat(iJR, method, exp)
            model = ME_DAT[:model]
            epouts = ME_DAT[:epouts]
            exp_beta = maximum(keys(epouts))
            epout = epouts[exp_beta]

            ChP.plot_marginal!(p, model, epout, model_ider;
                color = colors[method], 
                label = _textbf(labels[method]), 
                xlim,
                lw = 4,
                normalize = true
            )

            av = ChU.av(model, epout, model_ider)
            va = ChU.va(model, epout, model_ider)
            push!(avs, av); push!(vas, va) 
        end
        push!(ps, p)

    end # for exp_ider
    sfig(ChE, ps, 
        @fileid, "exchange_marginals", (;src, exp), ".png", 
        layout = (3, 1)
    )

    return ps

end