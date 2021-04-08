# ## -------------------------------------------------------------------
# # MSE per method
# let

#     for (iJR, Data) in [
#             (ChK.iJR904, ChK.KayserData), 
#             (ChF.iJR904, ChF.FolsomData), 
#             (ChH.iJR904, ChH.HeerdenData), 
#         ]

#         model = iJR.load_model("max_model")
#         datfile = iJR.procdir("dat.bson")
#         DAT = UJL.load_data(datfile; verbose = false)

#         FLX_IDERS = DAT[:FLX_IDERS]
#         EXPS = DAT[:EXPS]

#         p = plot(;xlabel = "experiment", ylabel = "MSE")
#         for method in [:ME_MAX_POL]
#             MSEs = []
#             for exp in EXPS

#                 sum = 0.0
#                 N = 0
#                 glc_flx = DAT[:exp, :flx, "GLC", exp]
#                 for ider in FLX_IDERS
#                     model_val = DAT[method, :flx, ider, exp]
#                     exp_val = DAT[:exp, :flx, ider, exp]
#                     sum += (model_val/glc_flx - exp_val/glc_flx)^2
#                     N += 1
#                 end
#                 push!(MSEs, sum / N)
#             end
#             scatter!(p, EXPS, MSEs; color = method_colors[method],
#                 label = string(method), m = 8, alpha = 0.8, 
#                 legend = :topleft
#             )
#             plot!(p, EXPS, MSEs; color = method_colors[method],
#                 label = "", ls = :dash, alpha = 0.8
#             )
#         end
#         src = nameof(Data)
#         mysavefig(p, "MSE_per_method"; src)
#     end
# end

## -------------------------------------------------------------------
let
    me_method = :ME_MAX_POL
    fba_method = :FBA_Z_FIX_MIN_COST

    p0 = plot(;xlabel = "politope box volume (log)", ylabel = "MSE")
    all_p = deepcopy(p0)
    zeroth = 1e-2
    for (iJR, Data) in [
                (ChK.iJR904, ChK.KayserData), 
                (ChF.iJR904, ChF.FolsomData), 
                (ChH.iJR904, ChH.HeerdenData), 
            ]
        
        src = nameof(Data)

        model = iJR.load_model("max_model")
        datfile = iJR.procdir("dat.bson")
        DAT = UJL.load_data(datfile; verbose = false)
        
        FLX_IDERS = ["GLC", "AC"]
        EXPS = DAT[:EXPS]
        rxns_map = iJR.load_rxns_map()
        
        exglcidx = ChU.rxnindex(model, iJR.EX_GLC_IDER)
        biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
        
        model_ids = [rxns_map[id] for id in FLX_IDERS]
        model_idxs = [ChU.rxnindex(model, model_id) for model_id in model_ids]
        
        
        # Reduce Pol
        p = deepcopy(p0)
        vols = []
        me_MSEs = []
        fba_MSEs = []

        # --------------------------------------------------- 
        # ME_MAX_POL
        for exp in EXPS
            # box vol
            cgD_X = -Data.ciD_X(:GLC, exp)
            D = Data.val(:D, exp)

            ChU.lb!(model, exglcidx, cgD_X)
            ChU.bounds!(model, biomidx, D, D)
            vol = box_vol(model, model_idxs) .+ zeroth
            isnan(vol) && continue

            push!(vols, vol)

            # MSE
            exp_glc = DAT[:exp, :flx, "GLC", exp]
            exp_flxs = DAT[:exp, :flx, FLX_IDERS, exp] 

            # ME
            me_flxs = DAT[me_method, :flx, FLX_IDERS, exp]
            MSEsum = sum(((exp_flxs ./ exp_glc) .- (me_flxs ./ exp_glc)).^2)
            MSE = MSEsum / length(FLX_IDERS)
            push!(me_MSEs, MSE)

            any(MSE .> 0.3) && let
                println(); @info("Outlayer", exp, src); println()
            end
            
            # FBA
            fba_flxs = DAT[fba_method, :flx, FLX_IDERS, exp]
            MSEsum = sum(((exp_flxs ./ exp_glc) .- (fba_flxs ./ exp_glc)).^2)
            MSE = MSEsum / length(FLX_IDERS)
            push!(fba_MSEs, MSE)
            
        end
        sidxs = sortperm(vols)

        marker = source_markers[Data]
        for _p in [p, all_p]

            scatter!(_p, vols[sidxs], me_MSEs[sidxs]; 
                label = "", m = (8, marker), 
                color = method_colors[me_method]
            )
                
            scatter!(_p, vols[sidxs], fba_MSEs[sidxs]; 
                label = "", m = (8, marker), 
                color = method_colors[fba_method]
            )
        end
    
        mysavefig(p, "pol_box_vs_MSE"; src)
    end
    plot!(all_p; xscale = :log10)
    mysavefig(all_p, "pol_box_vs_MSE")
end  

## -------------------------------------------------------------------
# # MSE per ider
# let
#     p = plot(;xlabel = "experiment", ylabel = "MSE")
#     for method in ALL_METHODS
#         MSEs = []

#         for ider in FLX_IDERS
#             sum = 0.0
#             N = 0
#             for exp in EXPS
#                 glc_flx = DAT[method, :Fd, :flx, "GLC", exp]
#                 model_val = DAT[method, :ep, :flx, ider, exp]
#                 exp_val = DAT[method, :Fd, :flx, ider, exp]
#                 sum += (model_val/glc_flx - exp_val/glc_flx)^2
#                 N += 1
#             end
#             push!(MSEs, sum / N)
#         end

#         scatter!(p, FLX_IDERS, MSEs; color = method_colors[method],
#             label = string(method), m = 8, alpha = 0.8, 
#             legend = :topleft
#         )
#         plot!(p, FLX_IDERS, MSEs; color = method_colors[method],
#             label = "", ls = :dash, alpha = 0.8
#         )
#     end
#     mysavefig(p, "MSE_per_ider")
# end

# ## -------------------------------------------------------------------
# # MSE per beta
# let
#     method = ME_MAX_POL

#     ps = Plots.Plot[]
#     for exp in EXPS
#         p = plot(;title = string("exp: ", exp), 
#             xlabel = "beta", ylabel = "MSE"
#         )

#         datfile = INDEX[method, :DFILE, exp]
#         dat = deserialize(datfile)
#         epouts = dat[:epouts]
#         betas = epouts |> keys |> collect |> sort
#         exp_beta = maximum(betas) # dat[:exp_beta]
#         model = dat[:model]
        
#         MSEs = []
#         for beta in betas
#             epout = epouts[beta]

#             sum = 0.0
#             N = 0

#             glc_flx = Fd.uval(:GLC, exp)
#             for ider in FLX_IDERS

#                 model_met = Fd_mets_map[ider]
#                 model_exch = Fd_rxns_map[ider]
#                 model_exchi = ChU.rxnindex(model, model_exch)

#                 model_flx = ChU.av(model, epout, model_exchi)
#                 exp_flx = Fd.uval(ider, exp)

#                 sum += (model_flx/glc_flx - exp_flx/glc_flx)^2
#                 N += 1
#             end
            
#             push!(MSEs, sum / N)
#         end

#         scatter!(p, first.(betas), MSEs; color = :black,
#             label = "", m = 8, alpha = 0.8
#         )
#         plot!(p, first.(betas), MSEs; color = :black,
#             label = "", ls = :dash, alpha = 0.8
#         )
#         vline!(p, [first(exp_beta)]; color = :black, 
#             label = "", ls = :dot, lw = 3, alpha = 0.9
#         )
#         push!(ps, p)
#     end
#     mysavefig(ps, "MSE_vs_beta")
# end
