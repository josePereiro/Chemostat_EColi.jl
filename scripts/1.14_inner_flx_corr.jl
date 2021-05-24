## ----------------------------------------------------------------------------
let

    # FBA_Z_FIX_MIN_VG_MAX_COST
    METHODS = [

        # FBA_Z_FIX_MAX_VG_MIN_COST, 
        # FBA_Z_FIX_MIN_VG_MIN_COST, 
        # FBA_Z_FIX_MAX_VG_MAX_COST, 
        FBA_Z_FIX_MIN_VG_MAX_COST,
        # FBA_Z_VG_FIX_MAX_COST, 
        # FBA_Z_VG_FIX_MIN_COST, 
        # FBA_Z_FIX_MAX_COST, 
        # FBA_Z_FIX_MIN_COST, 
        # FBA_MAX_Z_MIN_COST, 
        # FBA_MAX_Z_MAX_COST,
        
        ME_MAX_POL,
        # ME_Z_EXPECTED_G_BOUNDED
    ]

    MEname(c) = "\$\\textbf{ME}^{\\langle \\textbf{$c} \\rangle} \\textbf{ (abs. flux)} \$"

    METHODS_LABELS = [
        # "FBA_Z_FIX_MAX_VG_MIN_COST", 
        # "FBA_Z_FIX_MIN_VG_MIN_COST", 
        # "FBA_Z_FIX_MAX_VG_MAX_COST", 
        "\$ \\textbf{FBA} \\textbf{ (abs. flux)} \$",
        # "FBA_Z_VG_FIX_MAX_COST", 
        # "FBA_Z_VG_FIX_MIN_COST", 
        # "FBA_Z_FIX_MAX_COST", 
        # "FBA_Z_FIX_MIN_COST", 
        # "FBA_MAX_Z_MIN_COST", 
        # "FBA_MAX_Z_MAX_COST",

        MEname(1), # "ME unbiased",
        # MEname(0)
    ] 
    
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

    # EXPS = DAT[:EXPS]
    # EXPS = [Dintervals[0.1]; Dintervals[0.2]] |> unique |> sort

    for (Dav, EXPS) in Dintervals
    
        fontsize = 13
        commom_params = (; 
            title = "Internals", 
            dpi = 1000,
            thickness_scaling = 1.3, 
            xguidefontsize = fontsize, yguidefontsize = fontsize
        )

        Nd_rxns_map = iJR.load_rxns_map2() # inner reacts
        iJR_ider_subs = iJR.load_inner_rxns_subs()

        internals = iJR.load_inners_idermap() |> keys
        color_pool = Plots.distinguishable_colors(length(internals))
        ider_colors = Dict(id => c for (id, c) in zip(internals, color_pool))    

        ps = Plots.Plot[]
        for (method_label, method) in zip(METHODS_LABELS, METHODS)

            global_p = plot(;
                commom_params..., 
                xlabel = "\$ \\textbf{exp. (abs. flux)}\$", ylabel = method_label,
            )
            
            vals = []
            
            for sub in NiJR.subs()
                iJR_iders = iJR_ider_subs[sub]

                (sub in ["others"]) && continue
                marker = (6, subs_markers[sub])

                local_p = plot(;
                    commom_params...
                )
                
                for iJR_ider in iJR_iders

                    # try
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

                        for p in [local_p, global_p]
                            scatter!(p, exp_vals , model_vals;
                                label = "", 
                                yerr = model_errs,
                                color = ider_colors[iJR_ider], 
                                alpha = 0.8,
                                marker, 
                                xlim = [0.0, 2.3],
                                ylim = [0.0, 2.3],
                            )
                        end
                    # catch err; @warn(err) end
                end # for iJR_ider

                sort!(vals)
                plot!(local_p, [vals], [vals]; label = "",
                    ls = :dash, alpha = 0.5, lw = 3, color = :black,
                )
                # mysavefig(local_p, "inner_flx_corr"; method, sub, Dav)
            
            end # for (sub, iJR_iders)

            sort!(vals)
            plot!(global_p, [vals], [vals]; label = "",
                ls = :dash, alpha = 0.5, lw = 3, color = :black,
            )
            # mysavefig(global_p, "inner_flx_corr"; method, Dav)
            push!(ps, deepcopy(global_p))
        end

        mysavefig(ps, "inner_flx_corr"; Dav)
    end
end

# ## ----------------------------------------------------------------------------
#  let
#     METHODS = [
#         FBA_Z_FIX_MIN_COST,
#         # FBA_MAX_Z_MIN_COST, 
#         # FBA_Z_FIX_MAX_COST, 
#         # FBA_Z_FIX_MIN_VG_COST, 
#         # FBA_Z_VG_FIX_MIN_COST,
        
#         ME_MAX_POL,
#         ME_Z_EXPECTED_G_BOUNDED
#     ]

#     METHODS_LABELS = [
#         "FBA",
#         "FBA_MAX_Z_MIN_COST",
#         "FBA_Z_FIX_MAX_COST",
#         "FBA_Z_FIX_MIN_VG_COST",
#         "FBA_Z_VG_FIX_MIN_COST",

#         "ME unbiased",
#         "ME biased"
#     ]
    
#     # This is just for Nanchen
#     Data = Nd
#     iJR = NiJR
#     src = nameof(Data)

#     datfile = iJR.procdir("dat.bson")
#     DAT = UJL.load_data(datfile; verbose = false)
#     EXPS = DAT[:EXPS]

#     global_p = plot(;
#         xlabel = "square error", ylabel = "prob"
#     )

    
#     Nd_rxns_map = iJR.load_rxns_map2() # inner reacts
#     iJR_iders = iJR.load_inner_iders()
#     # iJR_iders = iJR.load_krebs_iders()
#     iJR_ider_subs = iJR.load_inner_rxns_subs()
    
#     # Nd_iders = getindex.([Nd_rxns_map], iJR_iders)
    
#     for (method_label, method) in zip(METHODS_LABELS, METHODS)
        
#         local_p = plot(;
#             title = method_label,
#             xlabel = "square error", ylabel = "prob"
#         )
        
#         diffs = []

#         exp_exglcs = abs.(DAT[:exp, :flx, "GLC", EXPS])

#         for (sub, iJR_iders) in iJR_ider_subs
#             for iJR_ider in iJR_iders
#                 Nd_ider = Nd_rxns_map[iJR_ider]
#                 exp_vals = abs.(DAT[:exp, :flx, Nd_ider, EXPS]) ./ exp_exglcs
#                 model_vals = abs.(DAT[method, :flx, iJR_ider, EXPS]) ./ exp_exglcs
#                 push!(diffs, abs.(exp_vals .- model_vals)...)
#             end
#         end

#         for p in [local_p, global_p]
#             color = method_colors[method]
#             histogram!(p, diffs; label = "", 
#                 color, bins = 80
#             )
#         end
    
#         mysavefig(local_p, "me_fba_diff_histogram"; method)
#     end
#     mysavefig(global_p, "me_fba_diff_histogram")
# end