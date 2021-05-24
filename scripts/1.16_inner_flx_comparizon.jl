## -------------------------------------------------------------------
let
    iJR = NiJR
    Ed = Nd

    objider = iJR.BIOMASS_IDER
    costider = iJR.COST_IDER
    exglcider = iJR.EX_GLC_IDER
    ids = [objider, costider, exglcider, "GLCpts", "PGI_fwd"]

    exp = 4
    # methods = [:FBA_Z_FIX_MIN_VG_MIN_COST, :FBA_Z_FIX_MIN_VG_MAX_COST]
    methods = [:FBA_Z_FIX_MAX_VG_MIN_COST, :FBA_Z_FIX_MAX_VG_MAX_COST]
    LP_DAT = ChE.load_LP_DAT(nameof(Ed));

    # for method in methods
    #     LP_DAT |> keys
    # end
    method1 = methods[1]
    model1 = LP_DAT[method1][:model][exp]
    fbaout1 = LP_DAT[method1][:fbaout][exp]

    method2 = methods[2]
    model2 = LP_DAT[method2][:model][exp]
    fbaout2 = LP_DAT[method2][:fbaout][exp]

    lim = [-15, 15]
    p = scatter(fbaout1.v, fbaout2.v; 
        label = "", m = 8, alpha = 0.7, 
        xlim = lim, ylim = lim,
        xlabel = string(methods[1]),
        ylabel = string(methods[2])
    )
    mysavefig(p, "test")
end

## -------------------------------------------------------------------
# model vs exp comparison
let
    METHODS = [:FBA_Z_FIX_MIN_VG_MIN_COST, :FBA_Z_FIX_MIN_VG_MAX_COST, :ME_MAX_POL]
    # METHODS = [LP_METHODS; :ME_MAX_POL]

    iJR = ChN.iJR904
    Data = ChN.NanchenData
    
    Nd_rxns_map = iJR.load_rxns_map2() # inner reacts
    iJR_ider_subs = iJR.load_inner_rxns_subs()

    sort_iders(iders) = sort!(deepcopy(iders), 
        by = (ider) -> ider == "GLC" ? "AAA" : ider
    )

    # Collect
    src = nameof(Data)

    DAT = ChE.load_DAT(src)

    FLX_IDERS = sort_iders(DAT[:FLX_IDERS])
    # EXPS = DAT[:EXPS]
    
    # experiments with approx D
    Dintervals = Dict(
        0.05 => 1:3,
        0.1 => 4:9,
        0.15 => 4:11,
        0.2 => 10:11,
        0.25 => 10:15,
        0.3 => 12:15,
        0.35 => 12:17,
        0.4 => 16:17,
        :tot => 1:17
    )

    @info("Doing", src)

    for (avD, exps) in Dintervals
        dat_pool = Dict()
        # Experiments
        edat = get!(dat_pool, :exp, Dict())
        for exp in exps
            !haskey(DAT, :exp, :flx, "GLC", exp) && continue
            glc_exp_val = DAT[:exp, :flx, "GLC", exp]
            
            for (sub, iJR_iders) in iJR_ider_subs
                for iJR_ider in iJR_iders
                    Nd_ider = Nd_rxns_map[iJR_ider]
                    exp_val = abs(DAT[:exp, :flx, Nd_ider, exp] / glc_exp_val)

                    @info(:exp, exp, iJR_ider, exp_val); println()

                    eidat = get!(edat, iJR_ider, [])
                    push!(eidat, exp_val)
                end
            end
        end
                
        # methods
        for method in METHODS
            mdat = get!(dat_pool, method, Dict())

            for exp in exps
                !haskey(DAT, :exp, :flx, "GLC", exp) && continue
                glc_exp_val = DAT[:exp, :flx, "GLC", exp]

                for (sub, iJR_iders) in iJR_ider_subs
                    for iJR_ider in iJR_iders
                        !haskey(DAT, method, :flx, iJR_ider, exp) && continue
                        model_val = abs(DAT[method, :flx, iJR_ider, exp] / glc_exp_val)
                        
                        # @info(method, exp, ider, model_val); println()
                        
                        midat = get!(mdat, iJR_ider, [])
                        push!(midat, model_val)
                    end
                end
            end
        end

        # Plots
        # Plot
        xticks_pos = []
        xticks_iders = []
        fontsize = 13
        p = plot(;xlabel = "Reactions", 
            ylabel = "flux / glc exch.", xgrid = false, 
            dpi = 1000,
            thickness_scaling = 1.3, 
            xguidefontsize = fontsize, yguidefontsize = fontsize
        )
        sep1 = 5
        sep2 = 3 * sep1

        # backgound
        gc = 1
        last_gc = gc
        ci = 1
        for sub in NiJR.subs()
            iJR_iders = iJR_ider_subs[sub]

            for iJR_ider in iJR_iders
                vline!(p, [gc - sep1]; 
                    label = "", color = :darkgray, alpha = 0.8, ls = :dash
                )
                gc += sep1
                for method in METHODS
                    gc += sep1
                end
                gc += sep2
            end
            fillalpha = isodd(ci) ? 0.6 : 0.3
            vspan!(p, [last_gc - sep1, gc - sep1]; label = "",
                linecolor = :darkgray, fillcolor = :darkgray, 
                fillalpha, linealpha = 0.6
            )
            last_gc = gc
            ci += 1
        end

        # fluxes
        gc = 1
        # f(v) = iszero(v) ? -6.0 : log(v)
        f(v) = v
        for sub in NiJR.subs()
            iJR_iders = iJR_ider_subs[sub]
            for iJR_ider in iJR_iders

                # Experimental
                edat = dat_pool[:exp]
                !haskey(edat, iJR_ider) && continue
                m = mean(abs.(edat[iJR_ider]))
                e = std(abs.(edat[iJR_ider]))

                @info(iJR_ider, m, e); println()

                scatter!(p, [gc + sep1], [f(m)]; 
                    yerr = [e], 
                    color = :black, label = "", ms = 6, m = :star,
                    alpha = 0.8
                )
                push!(xticks_pos, gc + 2 * sep1)
                push!(xticks_iders, iJR_ider)
                gc += sep1

                # modeled 
                for method in METHODS
                    mdat = dat_pool[method]
                    color = method_colors[method]
                    μ = mean(abs.(mdat[iJR_ider]))
                    e = std(abs.(mdat[iJR_ider]))
                    m = (6, get(method_markers, method, :square))
                    scatter!(p, [gc + sep1], [f(μ)]; 
                        yerr = [e], 
                        color, label = "", m, alpha = 0.8
                    )
                    gc += sep1
                end
                gc += sep2
            end
        end    

        xticks = (xticks_pos, xticks_iders)
        xrotation = 35
        size = (1100, 500)
        thickness_scaling = 1.3
        plot!(p; xticks, xrotation, size, thickness_scaling)
        mysavefig(p, "inner_flx_comparizon"; avD)
    end
end