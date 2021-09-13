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
    # _fix_outlayer!(ys, 1.0)
    n_ave_conv!(ys, 25)
    _normalize!(ys)
    return xs, ys
end

## ------------------------------------------------------
function _plot_marginal!(plt, P, Vv; pltargs...)
    v_av = sum(P .* Vv)
    vline!(plt, [v_av]; pltargs..., label = "", ls = :dot)

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
function plot_flux_marginals3D(simid, D, ϵs, cg, is_glclim_stst)
    MEmode = "ME_AVZ_EXPECTED_AVUG_BOUNDED"

    # Plots
    dyn_pz, dyn_pug, dyn_puo = plot(), plot(), plot()
    me_pz, me_pug, me_puo = plot(), plot(), plot()
    ps = [dyn_pz, dyn_pug, dyn_puo, me_pz, me_pug, me_puo]
    Vz, Vug, Vuo = nothing, nothing, nothing
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

        if !is_glclim_stst[simparams]
            @warn("Not glc-lim STST ", status)
            continue
        end

        dat = Dyn.load_maxent(simid, MEmode, simparams)
        PME = dat.PME

        S = Dyn.load_simdat(simid, simparams)
        PX = S.P
        V = S.V
        Vz = Vi(V, :z)
        Vug = Vi(V, :ug)
        Vuo = Vi(V, :uo)

        z_avPME = sum(PME .* Vz)
        if (abs(z_avPME - D) / D) > 0.05
            @warn("D and z not converged", status, z_avPME, D)
            continue
        end

        max_ϵs_ = maximum(ϵs)
        ls0 = 3

        for (P, pz, pug, puo) in [
                (PX, dyn_pz, dyn_pug, dyn_puo), 
                (PME, me_pz, me_pug, me_puo)
            ]
            comparams = (;
                ylim = [0.0, 1.2],
                ylabel = _textbf("non-normalized pdf"),
                lw =  ls0 + (2 * ls0) * (ϵ / max_ϵs_),
                color = :black, 
                alpha = 0.7
            )
            _plot_marginal!(pz, P, Vz; 
                xlabel = _textbf("z"), 
                comparams...
            )

            _plot_marginal!(pug, P, Vug; 
                xlabel = _textbf("u_g"), 
                comparams...
            )

            _plot_marginal!(puo, P, Vuo; 
                xlabel = _textbf("u_o"), 
                comparams...
            )
        end
    end

    return ps

end