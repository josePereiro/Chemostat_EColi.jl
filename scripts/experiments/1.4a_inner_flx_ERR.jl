## -------------------------------------------------------------------
# model vs exp comparison
function _plot_inner_fluxs_ERR(METHODS, METHODS_LABELS, METHODS_COLORS, METHODS_SHAPES, REF_METHOD)

    iJR = ChN.iJR904
    Data = ChN.NanchenData

    SUBS = ["glc tranport", "biomass", "glycolysis", "pentose phosphate", "krebs", "glyoxylate shunt"]
    
    Nd_rxns_map = iJR.load_rxns_map()
    iJR_ider_subs = iJR.load_inner_rxns_subs()
    iJR_ider_subs["biomass"] = ["D"]
    
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

        # avD != :tot && continue

        # ----------------------------------------------------------------
        # Collect data pool
        dat_pool = Dict()

        # Experiments
        edat = get!(dat_pool, :exp, Dict())
        for exp in exps
            !haskey(DAT, :exp, :flx, "GLC", exp) && continue
            # glc_exp_val = DAT[:exp, :flx, "GLC", exp]
            
            for sub in SUBS
                for iJR_ider in iJR_ider_subs[sub]
                    Nd_ider = (iJR_ider == "D") ? "D" : Nd_rxns_map[iJR_ider]
                    # exp_val = abs(DAT[:exp, :flx, Nd_ider, exp] / glc_exp_val)
                    exp_val = abs(DAT[:exp, :flx, Nd_ider, exp])

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
                # glc_exp_val = DAT[:exp, :flx, "GLC", exp]

                @info(method, exp); println()

                for sub in SUBS
                    (sub in ["others"]) && continue
                    iJR_iders = iJR_ider_subs[sub]
                    for iJR_ider in iJR_iders
                        !haskey(DAT, method, :flx, iJR_ider, exp) && continue
                        # model_val = abs(DAT[method, :flx, iJR_ider, exp] / glc_exp_val)
                        model_val = abs(DAT[method, :flx, iJR_ider, exp])
                        
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
        p = plot(;xlabel = "", 
            ylabel = _textbf("ARE"), xgrid = false
        )
        sep1 = 5
        sep2 = 3 * sep1

        # backgound
        gc = 1
        last_gc = gc
        ci = 1
        for sub in SUBS
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
        # f(x) = x <= 0 ? -5 : log10(x)
        f(x) = x
        for sub in SUBS
            (sub in ["others"]) && continue
            
            iJR_iders = iJR_ider_subs[sub]
            for iJR_ider in iJR_iders

                # Experimental
                edat = dat_pool[:exp]
                !haskey(edat, iJR_ider) && continue

                push!(xticks_pos, gc + 2 * sep1)
                push!(xticks_iders, iJR_ider)

                # reference
                ref_mdat = dat_pool[REF_METHOD]
                # ref_SE = (ref_mdat[iJR_ider] .- edat[iJR_ider]).^2
                # ref_SE = (ref_mdat[iJR_ider] ./ edat[iJR_ider] .- 1.0).^2
                ref_SE = (ref_mdat[iJR_ider] ./ edat[iJR_ider] .- 1.0) .|> abs
                ref_MERR = mean(ref_SE)
                ref_MERR = clamp(ref_MERR, 0, 100) # Inf = 10
                ref_SSE = std(ref_SE)
                
                # modeled 
                for method in METHODS
                    mdat = dat_pool[method]
                    color = METHODS_COLORS[method]
                    
                    # SE = (mdat[iJR_ider] .- edat[iJR_ider]).^2
                    # SE = (mdat[iJR_ider] ./ edat[iJR_ider] .- 1.0).^2
                    SE = (mdat[iJR_ider] ./ edat[iJR_ider] .- 1.0) .|> abs
                    MERR = mean(SE)
                    MERR = clamp(MERR, 0, 100) # Inf = 10
                    SSE = std(SE)

                    m = (ms, get(METHODS_SHAPES, method, :square))
                    # val = MERR - ref_MERR
                    
                    val = MERR
                    val = clamp(val, -1, 1)

                    scatter!(p, [gc + sep1], [f.(val)]; 
                        color, label = "", m, alpha = 0.8
                    )
                    gc += sep1
                    maxval = max(maxval, val)
                    minval = min(minval, val)
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
        xticks_iders[xticks_iders .== "D"] .= "BiomassEcoli"
        xticks = (xticks_pos, xticks_iders)
        xrotation = 45
        minval, maxval = f(minval), f(maxval)
        margin = (maxval - minval) * 0.1
        ylim = [minval - margin, maxval + margin]
        plot!(p; xticks, xrotation, ylim)
        # plot!(p; xticks, xrotation)

        ps[avD] = p

    end # for (avD, exps) in Dintervals

    return ps
end