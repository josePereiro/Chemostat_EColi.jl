let

    METHODS = [
        :FBA_MAX_Z_MIN_COST,
        :FBA_Z_FIX_MAX_VG_MIN_COST, 
        # :FBA_Z_FIX_MAX_VG_MAX_COST, 
        :FBA_Z_FIX_MIN_VG_MIN_COST, 
        # :FBA_Z_FIX_MIN_VG_MAX_COST, 
        # :FBA_Z_FIX_ADJ_MAX_VG_MIN_COST, 
        # :FBA_Z_FIX_ADJ_MAX_VG_MAX_COST, 
        # :FBA_Z_FIX_ADJ_MIN_VG_MIN_COST, 
        # :FBA_Z_FIX_ADJ_MIN_VG_MAX_COST, 
        # :FBA_Z_FIX_MAX_ATP_MIN_COST, 
        :FBA_Z_FIX_MAX_ATP_MAX_COST,
        :ME_MAX_POL, 
        # :ME_MAX_POL_B0
    ]

    fontsize = 10
    commom_params = (; 
        # dpi = 1000,
        # thickness_scaling = 1.0, 
        # xlabel = "exp", 
        # ylabel = "model",
        titlefont = 15,
        # aspectatio = 1.0,
        # xguidefontsize = fontsize, 
        # yguidefontsize = fontsize
        # xaxis = false, 
        # yaxis = false
    )

    avDs = [:exch, 0.05, 0.1, 0.2, 0.3, 0.4, :tot]
    ps = Plots.Plot[]
    for avD in avDs, method in METHODS
        title = string(method, "\n ------- \nD = ", avD)
        if avD == :exch
            k = (method, :exchs_flx_cor, :all_iders)
            p = plot!(deepcopy(PS[k...]); title, commom_params...)
        else
            # corrs
            k = (method, :inner_flx_cor, :all_subs, avD)
            p = plot!(deepcopy(PS[k...]); title, commom_params...)
        end
        push!(ps, p)
    end
    layout = (length(avDs), length(METHODS))
    mysavefig(ps, "inner_corr_callage"; layout)
    return

    # comparizon
    # push!(ps, PS[:inner_flx_comparizon, :all_subs, avD])

    # MSE
    push!(ps, PS[:inner_flx_MSE, :all_subs, avD])

    layout = (2, 4)
    margin = -3
    mysavefig(ps, "inner_callage"; avD, layout, margin)
end 