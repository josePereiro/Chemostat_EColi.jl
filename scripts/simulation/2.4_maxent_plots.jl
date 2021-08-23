function _fix_outlayer!(v, th)
    last = first(v)
    for i in eachindex(v)
        err = abs(v[i] - last)
        ref = abs((v[i] + last)/ 2) * th
        if (err > ref)
            v[i] = last
        end
        last = v[i]
    end
    v
end

_normalize!(v) = v .= v ./ maximum(v)

function _extract_hist_vecs(h::Dict)
    xs = collect(keys(h))
    sort!(xs)
    ys = [h[x] for x in xs]
    _fix_outlayer!(ys, 1.1)
    n_ave_conv!(ys, 25)
    _normalize!(ys)
    return xs, ys
end

## ------------------------------------------------------
function _plot_marginal!(plt, P, Vv; pltargs...)
    v_av = sum(P .* Vv)
    vline!(plt, [v_av]; pltargs..., label = "", ls = :dash)

    Pv = hist(Vv, P)
    xs, ys = _extract_hist_vecs(Pv)
    m, M = minimum(ys), maximum(ys)
    ylim = (m, M) .+ ((M - m) .* 0.1 .* [-1, 1])
    plot!(plt, xs, ys; label = "", ylim,
        pltargs...,
    )

    plt
end

## ------------------------------------------------------
function plot_flux_marginals(simid, D, ϵs, cg)
    MEmode = "ME_AVZ_EXPECTED_AVUG_BOUNDED"
    Dstr = string(round(D; sigdigits = 3))

    ϵcolors = palette(:thermal, length(ϵs))

    # Plots
    dyn_pz, dyn_pug = plot(), plot()
    me_pz, me_pug = plot(), plot()
    ps = [dyn_pz, dyn_pug, me_pz, me_pug]
    Vz, Vug = nothing, nothing
    for (ϵi, ϵ) in enumerate(ϵs)
        
        simparams = (;D, ϵ, cg)

        # Check and load
        status = Dyn.get_status(simid, simparams)
        if (status != STST_SIM_STATUS) 
            @warn("Wrong status ", status)
            continue
        end

        datfile = Dyn.maxent_file(simid, MEmode, simparams)
        if !isfile(datfile) 
            @warn("File missing ", datfile)
            continue
        end

        dat = Dyn.load_maxent(simid, MEmode, simparams)
        PME = dat.PME

        S = Dyn.load_simdat(simid, simparams)
        PX = S.P
        V = S.V
        Vz = Vi(V, :z)
        Vug = Vi(V, :ug)

        for (P, pz, pug) in [
                (PX, dyn_pz, dyn_pug), 
                (PME, me_pz, me_pug)
            ]
            ylim = [0.0, 1.2]
            ylabel = _textbf("norm. pdf")
            _plot_marginal!(pz, P, Vz; 
                xlabel = _textbf("z"), 
                ylabel, 
                ylim,
                lw = 3, 
                color = ϵcolors[ϵi], 
                zcolor = ϵ,
                colorbar_title = _textbf("\\epsilon")
            )

            _plot_marginal!(pug, P, Vug; 
                xlabel = _textbf("u_g"), 
                ylabel, 
                ylim,
                xlim = (-1.0, Inf),
                lw = 3, color = ϵcolors[ϵi], 
                zcolor = ϵ,
                colorbar_title = _textbf("\\epsilon")
            )

        end
    end

    # anotations
    x1s = maximum.([Vz, Vug])
    x0s = (0.0, -0.5)
    fontsize = 16
    for (title, ps) in [
            ("Dynamic", [dyn_pz, dyn_pug]), 
            ("ME", [me_pz, me_pug]), 
        ]
        for (p, x0, x1) in zip(ps, x0s, x1s)
            text = _textbf(title)
            annotate!(p, [(x0 + x1 * 0.03, 1.0, (text, fontsize, :left, :top, :black))])
            text = _textbf("D=", Dstr)
            annotate!(p, [(x0 + x1 * 0.03, 0.92, (text, fontsize, :left, :top, :black))])
        end
    end

    return ps

end