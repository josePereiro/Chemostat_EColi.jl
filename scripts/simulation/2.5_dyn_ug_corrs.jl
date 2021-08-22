function _plot_dyn_ug_corrs(simid, Ds, ϵs, cg)

    MEmode = "ME_AVZ_EXPECTED_AVUG_BOUNDED"
    dyn_zs, dyn_ugs = Float64[], Float64[]
    me_zs, me_ugs = Float64[], Float64[]
    lp_mug_ugs, lp_Mug_ugs = Float64[], Float64[]
    coll_ϵis, coll_ϵs = Int[], Float64[]

    for D in Ds, (ϵi, ϵ) in enumerate(ϵs)
        simparams = (;D, ϵ, cg)

        ## ------------------------------------------------------
        # Dynamic
        simparams = (;D, ϵ, cg)
        status = get_status(simid, simparams)
        if (status != STST_SIM_STATUS) 
            @warn("Wrong status ", status)
            continue
        end

        datfile = Dyn.simdat_file(simid, simparams)
        if !isfile(datfile) 
            @warn("File missing ", datfile)
            continue
        end
        Dyn_Dat = ldat(datfile)

        ## ------------------------------------------------------
        # MaxEnt
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
        ME_Dat = ldat(datfile)

        ## ------------------------------------------------------
        # LP max ug
        net = Dyn_Dat.V.net
        Dyn.fix!(net, :z, Dyn_Dat.D)
        Dyn.ub!(net, :ug, Dyn_Dat.cgD_X)
        Lug, Uug = Dyn.fva(net, :ug)

        ## ------------------------------------------------------
        # Collect
        push!(dyn_zs, Dyn_Dat.z_av)
        push!(dyn_ugs, Dyn_Dat.ug_av)
        push!(me_zs, ME_Dat.z_avPME)
        push!(me_ugs, ME_Dat.ug_avPME)
        push!(lp_mug_ugs, Lug)
        push!(lp_Mug_ugs, Uug)
        push!(coll_ϵis, ϵi)
        push!(coll_ϵs, ϵ)

    end

    ## ------------------------------------------------------
    # Plots

    ϵcolors = palette(:thermal, length(ϵs))

    me_pz = plot(;xlabel = _textbf("Dynamic z"), ylabel = _textbf("ME z"))
    me_pug = plot(;xlabel = _textbf("Dynamic u_g"), ylabel = _textbf("ME u_g"))
    lp_mug_pug = plot(;xlabel = _textbf("Dynamic u_g"), ylabel = _textbf("FBA min. u_g"))
    lp_Mug_pug = plot(;xlabel = _textbf("Dynamic u_g"), ylabel = _textbf("FBA max. u_g"))
    
    for (p, v1, v2) in [
            (me_pz, dyn_zs, me_zs), 
            (me_pug, dyn_ugs, me_ugs), 
            (lp_mug_pug, dyn_ugs, lp_mug_ugs), 
            (lp_Mug_pug, dyn_ugs, lp_Mug_ugs), 
        ]
        scatter!(p, v1, v2; 
            label = "", 
            color = ϵcolors[coll_ϵis], 
            zcolor = coll_ϵs,
            # colorbar_title = _textbf("\\epsilon"), 
            alpha = 0.8, 
            ms = 8
        )
        all_ = sort!([v1; v2])
        plot!(p, all_, all_; label = "", ls = :dash, c = :black, lw = 3)
    end

    return [lp_Mug_pug, lp_mug_pug, me_pug]
    
end