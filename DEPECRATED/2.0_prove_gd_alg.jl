using Base: Float64
import ProjAssistant
const PjAss = ProjAssistant
PjAss.@quickactivate

@time begin
    import Chemostat_EColi
    const ChE = Chemostat_EColi

    import Chemostat_Kayser2005
    const ChK = Chemostat_Kayser2005
    const KiJR = ChK.iJR904
    const Kd = ChK.KayserData

    import ChemostatPlots
    const ChP = ChemostatPlots
    
    import Chemostat
    import Chemostat.LP.MathProgBase
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChLP = Ch.LP
    const ChEP = Ch.MaxEntEP
    const ChT = Ch.Test
    using ProgressMeter
    using Base.Threads
    using SparseArrays

    using Serialization
    PjAss.set_verbose(false)

    using Base.Threads
    using Statistics
    using Plots
    import GR
    !isinteractive() && GR.inline("png")

    import SimTools
    import SimTools: gd_value, grad_desc!
    const ST = SimTools
    
end

## --------------------------------------------------------------------
# z_beta vs z_av
function beta_vs_av(model, mov_ider, mov_betas, fix_ider, fix_betas; 
        # MaxEnt
        alpha = Inf,
        maxiter = 5000,
        verbose = false,
        epsconv = 1e-5,
    )
    
    M, N = size(model)
    mov_idx = ChU.rxnindex(model, mov_ider)
    fix_idx = ChU.rxnindex(model, fix_ider)
    
    p = plot(; xlabel = string(mov_ider, " beta"), ylabel = string(mov_ider, " av"))
    
    lk = ReentrantLock()
    @threads for fix_beta in collect(fix_betas)

        @info("At", fix_beta)
        
        epout = nothing
        avs = Float64[]
        betas = Float64[]
        beta_vec = zeros(N)
        for mov_beta in mov_betas

            try
                beta_vec[mov_idx] = mov_beta
                beta_vec[fix_idx] = fix_beta
                epout = ChEP.maxent_ep(model; 
                    alpha, beta_vec, maxiter, verbose, epsconv, 
                    solution = epout
                )
                av = ChU.av(model, epout, mov_idx)
                (isnan(av) || isinf(av)) && break

                push!(avs, av)
                push!(betas, mov_beta)
            catch err; break; end

        end

        isempty(avs) && continue
        lock(lk) do
            plot!(p, betas, avs; label = "", lw = 3)
        end
    end

    ChE.sfig(p, 
        "2.0", "beta_vs_av", (;mov_ider, fix_ider), "png"
    )
end

## --------------------------------------------------------------------
# ToyModel
let
    model = ChT.toy_model();
    mov_betas = range(0.0, 1e3; length = 500)
    fix_betas = range(0.0, 1e3; length = 50)

    mov_ider = "biom"
    fix_ider = "gt"

    beta_vs_av(model, mov_ider, mov_betas, fix_ider, fix_betas)

    mov_ider = "gt"
    fix_ider = "biom"

    beta_vs_av(model, mov_ider, mov_betas, fix_ider, fix_betas)
end

## --------------------------------------------------------------------
# beta_z vs bate_ug
function make_heatmap(f::Function, model, 
        ider1, ider1_betas, ider2, ider2_betas;
        alpha = Inf,
        maxiter = 5000,
        verbose = false,
        epsconv = 1e-5,
    )

    M, N = size(model)

    idx1 = ChU.rxnindex(model, ider1)
    idx2 = ChU.rxnindex(model, ider2)

    mat = zeros(length(ider1_betas), length(ider2_betas))
    @threads for (ider1_betai, ider1_beta) in collect(enumerate(ider1_betas))
        
        beta_vec = zeros(N)
        @info("At", ider1_beta)
        epout = nothing
        
        beta_vec[idx1] = ider1_beta

        for (ider2_betai, ider2_beta) in enumerate(ider2_betas)
            beta_vec[idx2] = ider2_beta
            epout = ChEP.maxent_ep(model; 
                alpha, beta_vec, maxiter, verbose, epsconv, 
                solution = epout
            )
            mat[ider1_betai, ider2_betai] =  f(epout, ider1_beta, ider2_beta)
        end
    end

    p = plot(;xlabel = string(ider1, " beta"), ylabel = string(ider2, " beta"))
    heatmap!(p, ider1_betas, ider2_betas, mat'; label = "")

end

## --------------------------------------------------------------------
# "z_valid", "betas_heatmap"
let
    model = ChT.toy_model();
    
    Ds = range(0.0, 3.0; length = 12)
    
    bins = 50
    ider1 = "biom"
    ider1_betas = range(-1e2, 1e3; length = bins)
    ider2 = "gt"
    ider2_betas = range(-1e3, 1e3; length = bins)

    ps = map(Ds) do D
        
        z_valid_th = max(D * 0.1, 0.01)

        function get_valid(epout, args...)
            ider1_av = ChU.av(model, epout, ider1)
            ider2_av = ChU.av(model, epout, ider2)
            (isnan(ider1_av) || isinf(ider1_av)) && return NaN
            (isnan(ider2_av) || isinf(ider2_av)) && return NaN

            is_valid = (abs(ider1_av - D) <= z_valid_th)
            # is_valid = (ider2_av <= cgDX)
            is_valid ? ider1_av : 0.0
        end

        p =  make_heatmap(get_valid, model, ider1, ider1_betas, ider2, ider2_betas)
        plot!(p; title = string("D: ", round(D; sigdigits = 3)))
        return p
    end

    ChE.sfig(ps, 
        "2.0", "z_valid", "betas_heatmap", (;ider1, ider2), "png"
    )

end

## --------------------------------------------------------------------
# "ug_valid", "betas_heatmap"
let
    model = ChT.toy_model();

    cgDXs = range(0.0, 10.0; length = 12)
    
    bins = 50
    ider1 = "biom"
    ider1_betas = range(-1e2, 1e3; length = bins)
    ider2 = "gt"
    ider2_betas = range(-1e3, 1e3; length = bins)

    ps = map(cgDXs) do cgDX

        function get_valid(epout, args...)
            ider1_av = ChU.av(model, epout, ider1)
            ider2_av = ChU.av(model, epout, ider2)
            (isnan(ider1_av) || isinf(ider1_av)) && return NaN
            (isnan(ider2_av) || isinf(ider2_av)) && return NaN

            is_valid = (ider2_av <= cgDX)
            is_valid ? ider2_av : 0.0
        end

        p =  make_heatmap(get_valid, model, ider1, ider1_betas, ider2, ider2_betas)
        plot!(p; title = string("cgDX: ", round(cgDX; sigdigits = 3)))
        return p
    end

    ChE.sfig(ps, 
        "2.0", "ug_valid", "betas_heatmap", (;ider1, ider2), "png"
    )

end

## --------------------------------------------------------------------
# ## --------------------------------------------------------------------
# # EColi core
# let
#     model = ChT.ecoli_core_model();
#     mov_betas = range(0.0, 1e4; length = 500)
#     fix_betas = range(0.0, 1e2; length = 15)

#     mov_ider = "Biomass_Ecoli_core_w_GAM"
#     fix_ider = "EX_glc(e)"

#     beta_vs_av(model, mov_ider, mov_betas, fix_ider, fix_betas)
# end