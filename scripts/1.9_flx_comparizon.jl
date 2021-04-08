## -------------------------------------------------------------------
# model vs exp comparison
let

    # Collect
    dat_pool = Dict()
    for (iJR, Data) in [
        (ChK.iJR904, ChK.KayserData), 
        (ChF.iJR904, ChF.FolsomData), 
        (ChH.iJR904, ChH.HeerdenData), 
    ]

        src = nameof(Data)

        datfile = iJR.procdir("dat.bson")
        DAT = UJL.load_data(datfile; verbose = false)

        METHODS = [:FBA_Z_FIX_MIN_COST, :ME_MAX_POL]
        FLX_IDERS = DAT[:FLX_IDERS]
        EXPS = DAT[:EXPS]
            
        for method in METHODS
            mdat = get!(dat_pool, method, Dict())
            for exp in EXPS
                !haskey(DAT, method, :flx, "GLC", exp) && continue
                glc_exp_val = abs(DAT[method, :flx, "GLC", exp])
                for ider in FLX_IDERS
                    idat = get!(mdat, ider, Dict())
                    !haskey(DAT, method, :flx, ider, exp) && continue
                    exp_val = DAT[:exp, :flx, ider, exp] / glc_exp_val
                    model_val = DAT[method, :flx, ider, exp] / glc_exp_val

                    model_vals = get!(idat, :model_vals, [])
                    push!(model_vals, model_val)
                    exp_vals = get!(idat, :exp_vals, [])
                    push!(exp_vals, exp_val)
                end
            end
        end
    end # for dsource

    # Plot
    p = plot(;xlabel = "exchange", ylabel = "flx/exp_glc_flx")
    gc = 1
    xticks_pos = []
    for (ii, ider) in ALL_IDERS |> enumerate

        # Experimental
        _, edat = first(dat_pool)
        l = length(edat)
        dat = edat[ider]
        m = mean(abs.(dat[:exp_vals]))
        e = std(abs.(dat[:exp_vals]))
        scatter!(p, [gc], [m]; yerr = [e], 
            color = :black, label = "", ms = 8, m = :star,
            alpha = 0.8
        )
        push!(xticks_pos, gc)
        gc += 1

        # modeled 
        for (mi, (method, mdat)) in dat_pool |> enumerate
            color = method_colors[method]
            l = length(mdat)
            dat = mdat[ider]
            μ = mean(abs.(dat[:model_vals]))
            e = std(abs.(dat[:model_vals]))
            m = (8, method_markers[method])
            scatter!(p, [gc], [μ]; yerr = [e], 
                color, label = "", m, alpha = 0.8
            )
            gc += 1
        end
        gc += 5
         
    end
    xticks = (xticks_pos, ALL_IDERS)
    plot!(p; xticks)

    mysavefig(p, "rel_flx_per_ider")
end
