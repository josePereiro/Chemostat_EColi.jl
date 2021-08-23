## -----------------------------------------------------------------------------------------------
function plot_pol_vol_study(simid, Ds, ϵs, cg , ; bins = 150)
    ## -----------------------------------------------------
    # heatmap
    net = lglob(Dyn, :net2D)

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
        xlabel = _textbf("D"), 
        ylabel = _textbf("c_g", "D/X"), 
        colorbar_title = _textbf("polytope volume")
    )

    ## -----------------------------------------------------
    # simulation
    shape = [:circle, :rect, :star5, :diamond, :hexagon, :cross, :xcross, 
        :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, 
        :heptagon, :octagon, :star4, :star6, :star7, :star8, :vline, :hline, :+, :x]

    @showprogress for (ϵi, ϵ) in enumerate(ϵs)
            zs = []
            ugs = []
            ϵs_ = []
            
        for D in Ds
            simparams = (;D, ϵ, cg)
            status = get_status(simid, simparams)
            (status != STST_SIM_STATUS) && continue
            fn = Dyn.simdat_file(simid, simparams)
            !isfile(fn) && continue
            S = ldat(fn)
            push!(zs, S.z_av)
            push!(ugs, S.ug_av)
            push!(ϵs_, ϵ)
        end

        scatter!(p, zs, ugs; 
            label = _textbf("\\epsilon", "=", round(ϵ; sigdigits=3)),
            m = (9, shape[ϵi]), 
            color = :red
        )
        # fake label for spacing
        # scatter!(p, [zs[end]], [ugs[end]]; label = _textbf(" "), alpha = 0.0)
    end
    plot!(p; legend = :topleft)
    
    return p
    
end

     