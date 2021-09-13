function plot_flx_corrs(simid, Ds, ϵs, cg, is_glclim_stst)

    MEmode = "ME_AVZ_EXPECTED_AVUG_BOUNDED"

    # dat pool
    datpool = Dict()

    _push!(d, k, v) = push!(get!(d, k, typeof(v)[]), v)

    # Collect
    for D in Ds, (ϵi, ϵ) in enumerate(ϵs)
        simparams = (;D, ϵ, cg)

        ## ------------------------------------------------------
        # Dynamic
        simparams = (;D, ϵ, cg)
        status = get_status(simid, simparams)
        if (status != STST_SIM_STATUS) 
            # @warn("Wrong status ", (D, ϵ, cg), status)
            continue
        end

        datfile = Dyn.simdat_file(simid, simparams)
        if !isfile(datfile) 
            # @warn("File missing ", (D, ϵ, cg), datfile)
            continue
        end

        if !is_glclim_stst[simparams]
            @warn("Not glc-lim STST ", status)
            continue
        end

        S = ldat(datfile)
        PX = S.P
        V = S.V
        Vz = Vi(V, :z)
        Vug = Vi(V, :ug)
        Vuo = Vi(V, :uo)

        # Dyn dat
        z_avPX = sum(PX .* Vz)
        ug_avPX = sum(PX .* Vug)
        uo_avPX = sum(PX .* Vuo)
        cgD_X = S.cgD_X

        ## ------------------------------------------------------
        # LP 
        net0 = S.V.net
        
        # fmm & fmM
        net = deepcopy(net0)
        Dyn.fix!(net, :z, D)
        Dyn.ub!(net, :ug, cgD_X)
        Lug, _ = Dyn.fva(net, :ug)
        Dyn.fix!(net, :ug, Lug)
        Dyn.reducebox!(net)
        Luo, Uuo = Dyn.fva(net, :uo)

        z_fba_fmm, ug_fba_fmm, uo_fba_fmm = D, Lug, Luo
        z_fba_fmM, ug_fba_fmM, uo_fba_fmM = D, Lug, Uuo

        # fMm & fMM
        net = deepcopy(net0)
        Dyn.fix!(net, :z, D)
        Dyn.ub!(net, :ug, cgD_X)
        _, Uug = Dyn.fva(net, :ug)
        Dyn.fix!(net, :ug, Uug)
        Dyn.reducebox!(net)
        Luo, Uuo = Dyn.fva(net, :uo)

        z_fba_fMm, ug_fba_fMm, uo_fba_fMm = D, Uug, Luo
        z_fba_fMM, ug_fba_fMM, uo_fba_fMM = D, Uug, Uuo

        # M00
        net = deepcopy(net0)
        Dyn.ub!(net, :ug, cgD_X)
        Dyn.reducebox!(net)
        _, Uz = Dyn.fva(net, :z)
        Dyn.fix!(net, :z, Uz)
        Dyn.reducebox!(net)
        Lug, Uug = Dyn.fva(net, :ug)
        if !isapprox(Lug, Uug) 
            @warn("Lug != Uug", (D, ϵ, cg), Lug, Uug)
            continue
        end
        Dyn.fix!(net, :ug, Lug)
        Dyn.reducebox!(net)
        Luo, Uuo = Dyn.fva(net, :uo)
        !isapprox(Luo, Uuo) && @warn("Luo != Uuo", (D, ϵ, cg), Luo, Uuo)

        z_fba_M00, ug_fba_M00, uo_fba_M00 = Uz, Lug, Luo
        
        ## ------------------------------------------------------
        # MaxEnt
        status = Dyn.get_status(simid, simparams)
        if (status != STST_SIM_STATUS) 
            @warn("Wrong status ", (D, ϵ, cg), status)
            continue
        end

        datfile = Dyn.maxent_file(simid, MEmode, simparams)
        if !isfile(datfile) 
            @warn("File missing ", (D, ϵ, cg), datfile)
            continue
        end
        ME_Dat = ldat(datfile)
        PME = ME_Dat.PME

        # me dat
        z_avPME = sum(PME .* Vz)
        ug_avPME = sum(PME .* Vug)
        uo_avPME = sum(PME .* Vuo)

        if (abs(z_avPME - D) / D) > 0.05
            @warn("D and z not converged", status, z_avPME, D)
            continue
        end

        
        # collect
        for (plot_id, dyn_v, model_v) in [
            (:fba_fmm_z, z_avPX, z_fba_fmm),
            (:fba_fmm_ug, ug_avPX, ug_fba_fmm),
            (:fba_fmm_uo, uo_avPX, uo_fba_fmm),

            (:fba_fmM_z, z_avPX, z_fba_fmM),
            (:fba_fmM_ug, ug_avPX, ug_fba_fmM),
            (:fba_fmM_uo, uo_avPX, uo_fba_fmM),

            (:fba_fMM_z, z_avPX, z_fba_fMM),
            (:fba_fMM_ug, ug_avPX, ug_fba_fMM),
            (:fba_fMM_uo, uo_avPX, uo_fba_fMM),

            (:fba_fMm_z, z_avPX, z_fba_fMm),
            (:fba_fMm_ug, ug_avPX, ug_fba_fMm),
            (:fba_fMm_uo, uo_avPX, uo_fba_fMm),

            (:fba_M00_z, z_avPX, z_fba_M00),
            (:fba_M00_ug, ug_avPX, ug_fba_M00),
            (:fba_M00_uo, uo_avPX, uo_fba_M00),

        ]
            plt_dat = get!(datpool, plot_id, Dict())
            _push!(plt_dat, :dyn, dyn_v)
            _push!(plt_dat, :model, model_v)
            _push!(plt_dat, :ϵi, ϵi)
            _push!(plt_dat, :ϵ, ϵ)
        end

        for (plot_id, dyn_v, model_v) in [
            (:me_z, z_avPX, z_avPME)
            (:me_ug, ug_avPX, ug_avPME)
            (:me_uo, uo_avPX, uo_avPME)
        ]
            plt_dat = get!(datpool, plot_id, Dict())
            _push!(plt_dat, :dyn, dyn_v)
            _push!(plt_dat, :model, model_v)
            _push!(plt_dat, :ϵi, ϵi)
            _push!(plt_dat, :ϵ, ϵ)
        end
        
    end

    # # Test
    # datpool = lcache(ChE, :TEMP_CACHE) do
    #     datpool
    # end

    # sort
    _sort(f::Function, arr) = sort(arr; by = f)
    sdatpool = _sort(datpool) do elm
        model_label = split(string(elm), "_")
        replace!(model_label, "z" => "1z", "ug" => "2ug", "uo" => "3uo")
        join(model_label)
    end
    
    # plot
    ps = Plots.Plot[]
    for (plot_id, dat) in sdatpool
        @info("At", plot_id)

        model_label = split(string(plot_id), "_")
        replace!(model_label, "ug" => "(u_g)", "uo" => "(u_o)", "z" => "(z)", "fba" => "FBA", "me" => "ME")
        flux_label = last(model_label)
        model_label = join(model_label, " ")
        model_label = _textbf(model_label)
        dyn_label = _textbf("Dynamic", " ", flux_label)
        
        p = plot(;
            xlabel = dyn_label, 
            ylabel = model_label
        )
        push!(ps, p)

        dyn_vs = dat[:dyn]
        model_vs = dat[:model]
        ϵs_ = dat[:ϵ]
        # ϵis_ = dat[:ϵi]

        max_ϵs_ = maximum(ϵs)
        ms0 = 6

        scatter!(p, dyn_vs, model_vs; 
            label = "", 
            color = :black,
            alpha = 0.65, 
            m = :circle,
            ms = ms0 .+ (2 * ms0) .* (ϵs_ ./ max_ϵs_)
        )

        all_ = sort!([dyn_vs; model_vs])
        lims = (last(all_) - first(all_))
        lims = (first(all_) - (lims * 0.1), last(all_) + (lims * 0.1))
        plot!(p, all_, all_; 
            label = "", ls = :dash, c = :black, lw = 3, xlim = lims, ylim = lims
        )

        # annotation
        # if plot_id in [:fba_M00_ug, :fba_fmm_ug, :fba_fmM_ug, :fba_fMm_ug, :fba_fMM_ug, :me_ug]
            # plot!(p; leftmargin = 18mm)
            # m_, M_ = sort(all_)[[1, end]]
            # text_ = split(string(plot_id), "_")
            # replace!(text_, "ug" => "(u_g)", "uo" => "(u_o)", "z" => "(z)", "fba" => "FBA", "me" => "ME")
            # text_ = _textbf(text_...)
            # annotate!(p, m_ - (M_ - m_) * 0.25, m_ + (M_ - m_) * 0.5, text(text_, 25, :black, rotation = 90))
        # end

        # # title
        # if plot_id == :fba_M00_ug
        #     plot!(p; title = _textbf("u_g"))
        # elseif plot_id == :fba_M00_uo
        #     plot!(p; title = _textbf("u_o"))
        # elseif plot_id == :fba_M00_z
        #     plot!(p; title = _textbf("z"))
        # else
        #     plot!(p; title = _textbf(""))
        # end
        # plot!(p; title = "$plot_id") # Test

    end

    return ps

    # return [lp_Mug_pug, lp_mug_pug, me_pug]
    
end