function plot_flx_corrs3D(simid, Ds, ϵs, cg)

    MEmode = "ME_AVZ_EXPECTED_AVUG_BOUNDED"
    dyn_zs, dyn_ugs, dyn_uos = Float64[], Float64[], Float64[]
    me_zs, me_ugs, me_uos = Float64[], Float64[], Float64[]
    lp_min_zs, lp_min_ugs, lp_min_uos = Float64[], Float64[], Float64[]
    lp_max_zs, lp_max_ugs, lp_max_uos = Float64[], Float64[], Float64[]
    coll_ϵis, coll_ϵs = Int[], Float64[]

    for D in Ds, (ϵi, ϵ) in enumerate(ϵs)
        simparams = (;D, ϵ, cg)

        ## ------------------------------------------------------
        # Dynamic
        simparams = (;D, ϵ, cg)
        status = get_status(simid, simparams)
        if (status != STST_SIM_STATUS) 
            @warn("status ", status)
            continue
        end

        datfile = Dyn.simdat_file(simid, simparams)
        if !isfile(datfile) 
            @warn("File missing ", datfile)
            continue
        end
        Dyn_Dat = ldat(datfile)
        PX = Dyn_Dat.P

        ## ------------------------------------------------------
        # MaxEnt
        status = Dyn.get_status(simid, simparams)
        if (status != STST_SIM_STATUS) 
            @warn("status ", status)
            continue
        end

        datfile = Dyn.maxent_file(simid, MEmode, simparams)
        if !isfile(datfile) 
            @warn("File missing ", datfile)
            continue
        end
        ME_Dat = ldat(datfile)
        PME = ME_Dat.PME

        ## ------------------------------------------------------
        # LP 
        net = Dyn_Dat.V.net
        Dyn.fix!(net, :z, Dyn_Dat.D)
        Dyn.ub!(net, :ug, Dyn_Dat.cgD_X)
        Dyn.lb!(net, :ug, 0.0)
        Lug, Uug = Dyn.fva(net, :ug)
        
        Dyn.fix!(net, :ug, Lug)
        Lz, _ = Dyn.fva(net, :z)
        Luo, _ = Dyn.fva(net, :uo)
        
        Dyn.fix!(net, :ug, Uug)
        _, Uz = Dyn.fva(net, :z)
        _, Uuo = Dyn.fva(net, :uo)

        ## ------------------------------------------------------
        # Collect
        V = Dyn_Dat.V
        Vz = Vi(V, :z)
        Vug = Vi(V, :ug)
        Vuo = Vi(V, :uo)

        push!(dyn_zs, sum(PX .* Vz))
        push!(dyn_ugs, sum(PX .* Vug))
        push!(dyn_uos, sum(PX .* Vuo))
        push!(me_zs, sum(PME .* Vz))
        push!(me_ugs, sum(PME .* Vug))
        push!(me_uos, sum(PME .* Vuo))
        push!(lp_min_zs, Lz)
        push!(lp_min_ugs, Lug)
        push!(lp_min_uos, Luo)
        push!(lp_max_zs, Uz)
        push!(lp_max_ugs, Uug)
        push!(lp_max_uos, Uuo)
        push!(coll_ϵis, ϵi)
        push!(coll_ϵs, ϵ)

    end

    ## ------------------------------------------------------
    # Plots

    ϵcolors = palette(:thermal, length(ϵs))

    me_zp = plot(;title =  _textbf("z"), xlabel = _textbf("Dynamic"), ylabel = _textbf("ME"))
    me_ugp = plot(;title =  _textbf("u_g"), xlabel = _textbf("Dynamic"), ylabel = _textbf("ME"))
    me_uop = plot(;title =  _textbf("u_o"), xlabel = _textbf("Dynamic"), ylabel = _textbf("ME"))
    lp_min_zp = plot(;title =  _textbf("z"), xlabel = _textbf("Dynamic"), ylabel = _textbf("FBA min. u_g"))
    lp_min_ugp = plot(;title =  _textbf("u_g"), xlabel = _textbf("Dynamic"), ylabel = _textbf("FBA min. u_g"))
    lp_min_uop = plot(;title =  _textbf("u_o"), xlabel = _textbf("Dynamic"), ylabel = _textbf("FBA min. u_g"))
    lp_max_zp = plot(;title =  _textbf("z"), xlabel = _textbf("Dynamic"), ylabel = _textbf("FBA max. u_g"))
    lp_max_ugp = plot(;title =  _textbf("u_g"), xlabel = _textbf("Dynamic"), ylabel = _textbf("FBA max. u_g"))
    lp_max_uop = plot(;title =  _textbf("u_o"), xlabel = _textbf("Dynamic"), ylabel = _textbf("FBA max. u_g"))
    
    for (p, v1, v2) in [
            (me_zp, dyn_zs, me_zs), 
            (me_ugp, dyn_ugs, me_ugs), 
            (me_uop, dyn_uos, me_uos), 
            (lp_min_zp, dyn_zs, lp_min_zs), 
            (lp_min_ugp, dyn_ugs, lp_min_ugs), 
            (lp_min_uop, dyn_uos, lp_min_uos), 
            (lp_max_zp, dyn_zs, lp_min_zs), 
            (lp_max_ugp, dyn_ugs, lp_max_ugs), 
            (lp_max_uop, dyn_uos, lp_max_uos), 
        ]
        scatter!(p, v1, v2; 
            label = "", 
            color = ϵcolors[coll_ϵis], 
            zcolor = coll_ϵs,
            colorbar_title = _textbf("\\epsilon"), 
            alpha = 0.8, 
            # ms = 8, Test
            ms = 12
        )
        all_ = sort!([v1; v2])
        plot!(p, all_, all_; label = "", ls = :dash, c = :black, lw = 3)
    end

    return [
        # lp_min_zp, lp_max_zp, me_zp, 
        lp_min_ugp, lp_max_ugp, me_ugp, 
        lp_min_uop, lp_max_uop, me_uop, 
    ]
    
end