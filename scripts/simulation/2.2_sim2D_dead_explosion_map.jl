## ---------------------------------------------------------------
function plot_dead_explosion_map()
    
    params = lglob(Dyn, :dyn, :params, :infinite_cg)
    @extract params: Ds ϵs cg simid
    
    Ds_idxs = Dict(D => i for (i, D) in enumerate(Ds))
    ϵs_idxs = Dict(ϵ => i for (i, ϵ) in enumerate(ϵs))

    colGRAD = cgrad([colorant"black", colorant"blue", colorant"yellow"]);
    mat = fill(3, length(Ds_idxs), length(ϵs_idxs))
    
    for (D, ϵ) in Iterators.product(Ds, ϵs)
        simparams = (;D, ϵ, cg)
        
        status = get_status(simid, simparams)

        if status == DEAD_SIM_STATUS
            mat[Ds_idxs[D], ϵs_idxs[ϵ]] = 1
        elseif status == STST_SIM_STATUS
            mat[Ds_idxs[D], ϵs_idxs[ϵ]] = 2
        elseif status == EXPLODED_SIM_STATUS
            mat[Ds_idxs[D], ϵs_idxs[ϵ]] = 3
        else
            @warn("Unexpected status: ", simparams, status)
        end
    end

    @info("At", cg)
    p = plot(;
        title = _textbf("c_{g}=Inf"),
        xlabel = _textbf("D"), 
        ylabel = _textbf("\\epsilon")
    )
    # labels
    for (label, color, islast) = [
            (_textbf("dead"), :black, false),
            (_textbf("steady state"), :blue, false),
            (_textbf("explosion"), :yellow, true),
        ]
        scatter!(p, [Ds[end ÷ 2]], [ϵs[end ÷ 2]];
            markerstrokewidth=0.3,
            label, color, m = :square
        )
        # fake label for spacing
        # !islast && scatter!(p, [Ds[end ÷ 2]], [ϵs[end ÷ 2]]; label = _textbf(" "), alpha = 0.0)
    end

    # map
    heatmap!(p, Ds, ϵs, mat'; 
        color = colGRAD,
        colorbar = false,
        label = "", 
    )
    
    # sfig(ChE, p, 
    #     simid, "D_ϵ_heatmap", (;cg), ".png"
    # )

    return p
end    