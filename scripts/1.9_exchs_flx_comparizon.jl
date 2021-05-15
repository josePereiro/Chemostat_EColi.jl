## -------------------------------------------------------------------
# model vs exp comparison
let
    METHODS = [:FBA_Z_FIX_MIN_VG_MAX_COST, :ME_MAX_POL]

    sort_iders(iders) = sort!(deepcopy(iders), 
        by = (ider) -> ider == "GLC" ? "AAA" : ider
    )

    # Plot
    gc = 1
    xticks_pos = []
    xticks_iders = []
    p = plot(;xlabel = "exchange", ylabel = "flx/exp_glc_flx")

    for srcs in [
        [
            (ChK.iJR904, ChK.KayserData)], 
            [(ChN.iJR904, ChN.NanchenData)], 
            [(ChF.iJR904, ChF.FolsomData)
        ],
        [
            (ChH.iJR904, ChH.HeerdenData)
        ]
    ]

        # Collect
        dat_pool = Dict()
        for (iJR, Data) in srcs

            src = nameof(Data)

            datfile = iJR.procdir("dat.bson")
            DAT = UJL.load_data(datfile; verbose = false)

            FLX_IDERS = sort_iders(DAT[:FLX_IDERS])
            EXPS = DAT[:EXPS]
            
            @info("Doing", src)

            # Experiments
            edat = get!(dat_pool, :exp, Dict())
            for exp in EXPS
                !haskey(DAT, :exp, :flx, "GLC", exp) && continue
                glc_exp_val = DAT[:exp, :flx, "GLC", exp]
                
                for ider in FLX_IDERS
                    exp_val = abs(DAT[:exp, :flx, ider, exp] / glc_exp_val)

                    @info(:exp, exp, ider, exp_val); println()

                    eidat = get!(edat, ider, [])
                    push!(eidat, exp_val)
                end
            end
                
            # methods
            for method in METHODS
                mdat = get!(dat_pool, method, Dict())

                for exp in EXPS
                    !haskey(DAT, :exp, :flx, "GLC", exp) && continue
                    glc_exp_val = DAT[:exp, :flx, "GLC", exp]

                    for ider in FLX_IDERS
                        !haskey(DAT, method, :flx, ider, exp) && continue
                        model_val = abs(DAT[method, :flx, ider, exp] / glc_exp_val)
                        
                        # @info(method, exp, ider, model_val); println()
                        
                        midat = get!(mdat, ider, [])
                        push!(midat, model_val)
                    
                    end
                end
            end
        end # for dsource

        # Plos
        vline!(p, [gc]; label = "", color = :black, alpha = 0.7, ls = :dash)
        for (ii, ider) in sort_iders(ALL_IDERS) |> enumerate

            # Experimental
            edat = dat_pool[:exp]
            !haskey(edat, ider) && continue
            m = mean(abs.(edat[ider]))
            e = std(abs.(edat[ider]))

            @info(ider, m, e); println()

            scatter!(p, [gc], [m]; yerr = [e], 
                color = :black, label = "", ms = 6, m = :star,
                alpha = 0.8
            )
            push!(xticks_pos, gc)
            push!(xticks_iders, ider)
            gc += 5

            # modeled 
            for method in METHODS
                mdat = dat_pool[method]
                color = method_colors[method]
                μ = mean(abs.(mdat[ider]))
                e = std(abs.(mdat[ider]))
                m = (6, method_markers[method])
                scatter!(p, [gc], [μ]; yerr = [e], 
                    color, label = "", m, alpha = 0.8
                )
                gc += 5
            end
            gc += 15
        end    
    end # for srcs

    xticks = (xticks_pos, xticks_iders)
    xrotation = 35
    size = (1100, 500)
    thickness_scaling = 1.3
    plot!(p; xticks, xrotation, size, thickness_scaling)
    mysavefig(p, "exchs_flx_comparizon")
end
