## -----------------------------------------------------------------------------------------------
function plot_pol_vol_study(simid, Ds, ϵs, cg , ; bins = 150)
    ## -----------------------------------------------------
    # heatmap
    net = lglob(Dyn, :net3D)

    Lz, Uz = Dyn.bounds(net, :z)
    Ds_ = range(Lz, Uz; length = bins)
    Lug, Uug = Dyn.bounds(net, :ug)
    cgD_Xs = range(Lug, Uug; length = bins)

    mat = fill(NaN, bins, bins)

    for (Di, D) in enumerate(Ds_)
        Dyn.fix!(net, :z, D)
        ugl, ugu = Dyn.fva(net, :ug)
        
        for (cgD_Xi, cgD_X) in enumerate(cgD_Xs)
            !(ugl <= cgD_X <= ugu) && continue
            mat[Di, cgD_Xi] = (cgD_X - ugl)
        end
    end

    p = heatmap(Ds_, cgD_Xs, mat';
        label = "",
        color = :greys,
        title = "", 
        xlabel = _textbf("D"), 
        ylabel = _textbf("c_g", "D/X"), 
        # colorbar_title = _textbf("polytope volume")
    )

    
    fontsize = 19
    text_ = _textbf("polytope volume")
    plot!(p; rightmargin = 15mm)
    annotate!(p, 0.65, 13.0 / 2, text(text_, fontsize, :black, rotation = 90))

    ## -----------------------------------------------------
    # simulation
    max_ϵs_ = maximum(ϵs)
    ms0 = 5

    @showprogress for (ϵi, ϵ) in enumerate(ϵs)
        zs = []
        ugs = []
            
        for D in Ds
            simparams = (;D, ϵ, cg)
            status = get_status(simid, simparams)
            (status != STST_SIM_STATUS) && continue
            fn = Dyn.simdat_file(simid, simparams)
            !isfile(fn) && continue
            S = ldat(fn)
            push!(zs, S.z_av)
            push!(ugs, S.ug_av)
        end

        isempty(zs) && continue
        scatter!(p, zs, ugs; 
            # label = _textbf("\\epsilon", "=", round(ϵ; sigdigits=3)),
            label = "",
            # m = (9, shape[ϵi]), 
            markerstrokecolor = :white,
            color = :black,
            ms = ms0 + (2 * ms0) * (ϵ / max_ϵs_),
            # color = :red
        )
    end
    plot!(p; tite = _textbf("A"), legend = :topleft)
    
    return p
    
end

     