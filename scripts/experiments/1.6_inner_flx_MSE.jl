## -------------------------------------------------------------------
# model vs exp comparison
function _plot_inner_fluxs_MSE(METHODS, METHODS_LABELS, METHODS_COLORS, METHODS_SHAPES)

    iJR = ChN.iJR904
    Data = ChN.NanchenData
    
    Nd_rxns_map = iJR.load_rxns_map()
    iJR_ider_subs = iJR.load_inner_rxns_subs()
    
    sort_iders(iders) = sort!(deepcopy(iders), 
        by = (ider) -> ider == "GLC" ? "AAA" : ider
    )
    
    # Collect
    src = nameof(Data)

    DAT = ChE.load_DAT(src)
    
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
    
    ps = Dict()

    for (avD, exps) in Dintervals

        # ----------------------------------------------------------------
        # Collect data pool
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

                @info(method, exp); println()

                for (sub, iJR_iders) in iJR_ider_subs
                    (sub in ["others"]) && continue
                    for iJR_ider in iJR_iders
                        !haskey(DAT, method, :flx, iJR_ider, exp) && continue
                        model_val = abs(DAT[method, :flx, iJR_ider, exp] / glc_exp_val)
                        
                        @info(method, exp, iJR_ider, model_val); println()
                        
                        midat = get!(mdat, iJR_ider, [])
                        push!(midat, model_val)
                    end
                end
            end
        end

        # ----------------------------------------------------------------
        # Plots
        xticks_pos = []
        xticks_iders = []
        ms = 8
        p = plot(;xlabel = _textbf("Reactions"), 
            ylabel = _textbf("log MSE"), xgrid = false
        )
        sep1 = 5
        sep2 = 3 * sep1

        # backgound
        gc = 1
        last_gc = gc
        ci = 1
        for sub in NiJR.subs()
            (sub in ["others"]) && continue
            iJR_iders = iJR_ider_subs[sub]

            for iJR_ider in iJR_iders
                vline!(p, [gc - sep1]; 
                    label = "", color = :darkgray, 
                    alpha = 0.8, ls = :dash
                )

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
        minval, maxval = Inf, -Inf
        f(x) = x <= 0 ? -5 : log10(x)
        for sub in NiJR.subs()
            (sub in ["others"]) && continue
            
            iJR_iders = iJR_ider_subs[sub]
            for iJR_ider in iJR_iders

                # Experimental
                edat = dat_pool[:exp]
                !haskey(edat, iJR_ider) && continue

                push!(xticks_pos, gc + 2 * sep1)
                push!(xticks_iders, iJR_ider)

                # modeled 
                for method in METHODS
                    mdat = dat_pool[method]
                    color = METHODS_COLORS[method]

                    SE = (mdat[iJR_ider] .- edat[iJR_ider]).^2
                    MSE = mean(SE)
                    SSE = std(SE)

                    m = (ms, get(METHODS_SHAPES, method, :square))
                    scatter!(p, [gc + sep1], [f.(MSE)]; 
                        color, label = "", m, alpha = 0.8
                    )
                    gc += sep1
                    maxval = max(maxval, MSE)
                    minval = min(minval, MSE)
                end
                gc += sep2
            end
        end    

        # labels
        for method in METHODS
            label =  METHODS_LABELS[method]
            color =  METHODS_COLORS[method]
            shape =  METHODS_SHAPES[method]
            scatter!(p, [1], [-10];
                # markerstrokewidth = 0.3,
                label, color, m = (ms, shape),
                legend = :outertopright
            )
            # fake label for spacing
            # !islast && scatter!(p, [Ds[end ÷ 2]], [ϵs[end ÷ 2]]; label = _textbf(" "), alpha = 0.0)
        end

        # ticks and margin
        xticks = (xticks_pos, xticks_iders)
        xrotation = 45
        minval, maxval = f(minval), f(maxval)
        margin = (maxval - minval) * 0.1
        ylim = [minval - margin, maxval + margin]
        plot!(p; xticks, xrotation, ylim)

        ps[avD] = p

        # return ps

    end # for (avD, exps) in Dintervals

    return ps
end

# -------------------------------------------------------------------
let

    # [:none, :auto, :circle, :rect, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, 
    # :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star6, 
    # :star7, :star8, :vline, :hline, :+, :x].
    METHODS_SHAPES = Dict(
        :FBA_MAX_Z_MIN_COST => :utriangle,
        :FBA_Z_FIX_MAX_VG_MIN_COST => :dtriangle,
        :FBA_Z_FIX_MIN_VG_MIN_COST => :rtriangle,
        :FBA_Z_FIX_MAX_ATP_MAX_COST => :ltriangle,
        :ME_MAX_POL => :circle,
    )

    METHODS = [
        :FBA_MAX_Z_MIN_COST,
        :FBA_Z_FIX_MAX_VG_MIN_COST, 
        :FBA_Z_FIX_MIN_VG_MIN_COST, 
        :FBA_Z_FIX_MAX_ATP_MAX_COST,
        :ME_MAX_POL
    ]

    METHODS_COLORS = Dict(
        :FBA_MAX_Z_MIN_COST => :red,
        :FBA_Z_FIX_MAX_VG_MIN_COST => :pink, 
        :FBA_Z_FIX_MIN_VG_MIN_COST => :yellow, 
        :FBA_Z_FIX_MAX_ATP_MAX_COST => :brown,
        :ME_MAX_POL => :blue
    )

    METHODS_LABELS = Dict(
        :FBA_MAX_Z_MIN_COST => "FBA max. z",
        :FBA_Z_FIX_MAX_VG_MIN_COST => "FBA min. ug",
        :FBA_Z_FIX_MIN_VG_MIN_COST => "FBA max. ug",
        :FBA_Z_FIX_MAX_ATP_MAX_COST => "FBA max. ATP",
        :ME_MAX_POL => "ME",
    )

    params = (;
       guidefont = font(16),
       xtickfont = font(13),
       ytickfont = font(13),
       legendfont = font(12),
       thickness_scaling = 1.7,
       size = (1700, 900),
       bottom_margin = 7mm
    )

    ps = _plot_inner_fluxs_MSE(METHODS, METHODS_LABELS, METHODS_COLORS, METHODS_SHAPES)

    for (avD, p) in ps
        plot!(p; params...)
        sfig(ChE, p,
            "inner_flx_MSE", (;avD), ".png"
        )
    end
end