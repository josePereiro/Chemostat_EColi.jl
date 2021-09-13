## ----------------------------------------------------------------------------
function plot_inner_plots(METHODS, METHODS_LABELS = string.(METHODS))

    # This is just for Nanchen
    Data = Nd
    iJR = NiJR
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

    int_ps = Dict()
    for (avD, EXPS) in Dintervals

        Nd_rxns_map = iJR.load_rxns_map() 
        iJR_ider_subs = iJR.load_inner_rxns_subs()

        internals = iJR.load_inners_idermap() |> keys
        color_pool = Plots.distinguishable_colors(length(internals))
        ider_colors = Dict(id => c for (id, c) in zip(internals, color_pool))    

        ps = Plots.Plot[]
        for (method_label, method) in zip(METHODS_LABELS, METHODS)

            p = plot(;
                xlabel = _textbf("exp abs. flux"), 
                ylabel = _textbf("model abs. flux")
            )
            
            vals = []
            
            for sub in NiJR.subs()
                iJR_iders = iJR_ider_subs[sub]

                (sub in ["others"]) && continue
                marker = (35, subs_markers[sub])

                
                for iJR_ider in iJR_iders

                    Nd_ider = Nd_rxns_map[iJR_ider]

                    @info "Doing" Nd_ider iJR_ider

                    exp_exglcs = abs.(DAT[:exp, :flx, "GLC", EXPS])
                    exp_vals = abs.(DAT[:exp, :flx, Nd_ider, EXPS]) ./ exp_exglcs
                    model_vals = abs.(DAT[method, :flx, iJR_ider, EXPS]) ./ exp_exglcs
                    model_errs = abs.(get.([DAT], 0.0, method, :err, iJR_ider, EXPS)) ./ exp_exglcs
                    model_errs = map(model_errs) do err
                        isnan(err) && return 0.0
                        err > 1.0 && return 0.0
                        err
                    end
                    
                    push!(vals, exp_vals..., model_vals...)

                    scatter!(p, exp_vals , model_vals;
                        label = "", 
                        yerr = model_errs,
                        color = ider_colors[iJR_ider], 
                        alpha = 0.8,
                        marker,
                    )

                end # for iJR_ider
            
            end # for (sub, iJR_iders)

            sort!(vals)
            plot!(p, [vals], [vals]; label = "", 
                ls = :dash, alpha = 0.8, lw = 3, 
                color = :black,
            )
 
            # title
            m, M = minimum(vals), maximum(vals)
            x = m + (M - m) * 0.01
            y = M - (M - m) * 0.05
            fontsize = 50
            Dstr = (avD == :tot) ? "all" : round(avD; sigdigits=2)
            text = _textbf(method_label)
            annotate!(p, [(x, y, (text, fontsize, :left, :top, :black))])
            text = _textbf("D=", Dstr)
            annotate!(p, [(x, y - y * 0.13, (text, fontsize, :left, :top, :black))])

            push!(ps, p)
        
        end # for (method_label

        # sfig(ChE, ps, 
        #     @fileid, "inner_flx_corr", (;avD), ".png", 
        #     layout = (1, length(ps))
        # )

        int_ps[avD] = ps
    end

    return int_ps
end