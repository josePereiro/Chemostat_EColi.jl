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

function _normalize!(v)
    v .= v ./ maximum(v)
end

function _extract_hist_vecs(h::Dict)
    xs = collect(keys(h))
    sort!(xs)
    ys = [h[x] for x in xs]
    _fix_outlayer!(ys, 1.0)
    n_ave_conv!(ys, 8)
    _normalize!(ys)
    return xs, ys
end

## ------------------------------------------------------
function _plot_marginal!(plt, P, Vv; pltargs...)
    Pv = hist(Vv, P)
    xs, ys = _extract_hist_vecs(Pv)
    plot!(plt, xs, ys; label = "", 
        pltargs...,
    )
    v_av = sum(P .* Vv)
    vline!(plt, [v_av]; pltargs..., label = "", ls = :dash)
    plt
end

## ------------------------------------------------------
function plot_flux_marginals(simid, D, ϵs, cg)
    MEmode = "ME_AVZ_EXPECTED_AVUG_BOUNDED"

    ϵcolors = palette(:thermal, length(ϵs))

    # Plots
    dyn_pz, dyn_pug = plot(), plot()
    me_pz, me_pug = plot(), plot()
    ps = [dyn_pz, dyn_pug, me_pz, me_pug]
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

        Dstr = string(round(D; sigdigits = 3))
        for (title, P, pz, pug) in [
                (string("Dynamic D=", Dstr), PX, dyn_pz, dyn_pug), 
                (string("ME D=", Dstr), PME, me_pz, me_pug)
            ]
            _plot_marginal!(pz, P, Vz; 
                xlabel = "z", ylabel = "proj", 
                title, lw = 3, 
                color = ϵcolors[ϵi], 
                zcolor = ϵ,
                colorbar_title = "ϵ"
            )

            _plot_marginal!(pug, P, Vug; 
                xlabel = "ug", ylabel = "proj", 
                title, lw = 3, color = ϵcolors[ϵi], 
                zcolor = ϵ,
                colorbar_title = "ϵ"
            )

        end
    end

    return ps

end