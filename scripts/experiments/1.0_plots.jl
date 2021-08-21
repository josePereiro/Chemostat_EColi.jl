using ProjAssistant: ImgTools
using ProjAssistant
@quickactivate

# ---------------------------------------------------------------------------------
@time begin
    import Chemostat_EColi
    const ChE = Chemostat_EColi
    
    import Chemostat_Kayser2005
    const ChK = Chemostat_Kayser2005
    const KiJR = ChK.iJR904
    const Kd = ChK.KayserData

    import Chemostat_Nanchen2006
    const ChN = Chemostat_Nanchen2006
    const NiJR = ChN.iJR904
    const Nd = ChN.NanchenData

    import Chemostat_Folsom2014
    const ChF = Chemostat_Folsom2014
    const FiJR = ChF.iJR904
    const Fd = ChF.FolsomData

    import ChemostatPlots
    const ChP = ChemostatPlots
    
    import Chemostat
    import Chemostat.LP.MathProgBase
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChLP = Ch.LP

    import ImgTools
    import ImgTools: make_grid, WHITE_PIX

    using ProgressMeter
    using Base.Threads
    using SparseArrays
    using MagicDicts
    using Serialization
    using Statistics
    using Plots

    import GR
    !isinteractive() && GR.inline("png")
    
end

## ---------------------------------------------------------------------------------
ider_colors = Dict(
    "GLC" => :black, 
    "CO2" => :blue, 
    "O2" => :red, 
    "AC" => :pink, 
    "NH4" => :orange, 
    "FORM" => :purple, 
    "SUCC" => :yellow,
)

exp_colors = let
    EXPS = 1:20
    len = length(EXPS)
    colors = Plots.distinguishable_colors(len)
    Dict(ider => color for (ider, color) in zip(EXPS, colors))
end

subs_markers = Dict(
    "glycolysis" => :circle,
    "krebs" => :diamond,
    "pentose phosphate" => :utriangle,
    "others" => :star,
    "glyoxylate shunt" => :square,
    "glc tranport" => :hex
)

source_markers = Dict(
    Kd => :square, 
    Fd => :circle, 
    Nd => :utriangle, 
)

## ------------------------------------------------------
function pltparams()
    return (;
       titlefont = 16,
       axisfont = 16,
       guidefont = 16,
       colorbar_titlefont = 16,
       xtickfont = 12,
       ytickfont = 12,
       legendfont = 12,
       thickness_scaling = 1.6,
       size = (1220, 940)
    )
end

## ---------------------------------------------------------------------------------
# figure 8
include("1.1_pol_box_size.jl")
include("1.2_Xexp_Xmax_study.jl")
include("1.3_Glc_uptake_study.jl")

let
    # panel 1
    Dlim = [0.01, 0.5]
    cgD_Xlim = [-40.0, 0.0]
    bins = 60
    panel1 = plot_pol_box_size(
        ChN.iJR904, ChN.NanchenData; 
        Dlim, cgD_Xlim, bins
    )

    # panel 2
    panel2 = plot_Xexp_Xmax_study()
    
    # panel 3
    panel3 = plot_Glc_uptake_study()

    for p in [panel1, panel2, panel3]
        plot!(p; pltparams()...)
    end
    
    sfig(ChE, [panel1, panel2, panel3], 
        @fileid, "exp_space_location", ".png";
        layout = (1, 3)
    )

end

## ---------------------------------------------------------------------------------
# figure 9
include("1.4_inner_flx_corr.jl")
include("1.5_exchs_flx_correlations.jl")

let

    METHODS = [
        :FBA_MAX_Z_MIN_COST,
        :FBA_Z_FIX_MAX_VG_MIN_COST, 
        :FBA_Z_FIX_MIN_VG_MIN_COST, 
        :FBA_Z_FIX_MAX_ATP_MAX_COST,
        :ME_MAX_POL
    ]

    METHODS_LABELS = [
        "FBA max. z",
        "FBA min. ug",
        "FBA max. ug",
        "FBA max. ATP",
        "ME",
    ]

    exch_ps = plot_exch_corrs(METHODS, METHODS_LABELS)
    inner_ps = plot_inner_plots(METHODS, METHODS_LABELS)

    ps = []
    push!(ps, exch_ps...)
    push!(ps, inner_ps[0.05]...)
    push!(ps, inner_ps[0.1]...)
    push!(ps, inner_ps[0.2]...)
    push!(ps, inner_ps[0.3]...)
    push!(ps, inner_ps[0.4]...)
    push!(ps, inner_ps[:tot]...)

    params = (;
       titlefont = 26,
       axisfont = 26,
       guidefont = 26,
       xtickfont = 22,
       ytickfont = 22,
       legendfont = 12,
       thickness_scaling = 1.6,
       size = (1220, 940)
    )

    for p in ps
        plot!(p; params...)
    end

    sfig(ChE, ps,
        "flux_corrs_collage", ".png";
        layout = (7, 5)
    )

end;

## ---------------------------------------------------------------------------------
# figure 10
include("1.6_var_study_and_marginals.jl")

let    
    corr_av_p, corr_va_p = plot_ME_methods_comparizon()
    pGLC, pAC, pz = plot_me_marginals()

    size0 = (1220, 940)
    for p in [pGLC, pAC, pz]
        plot!(p; 
            titlefont = 26,
            axisfont = 26,
            guidefont = 26,
            xtickfont = 22,
            ytickfont = 22,
            legendfont = 22,
            thickness_scaling = 1.6,
            size = size0
        )
    end

    margs_file = sfig(ChE, [pGLC, pAC, pz],
        "margs.png";
        layout = (3, 1)
    )

    for p in [corr_av_p, corr_va_p]
        plot!(p; 
            titlefont = 26,
            axisfont = 26,
            guidefont = 26,
            xtickfont = 22,
            ytickfont = 22,
            legendfont = 22,
            thickness_scaling = 1.6,
            size = 1.5 .* size0
        )
    end

    corrs_file = sfig(ChE, [corr_av_p, corr_va_p],
        "corrs.png";
        layout = (2, 1)
    )

    grid_ = make_grid(
        lfig.([margs_file, corrs_file]); 
        layout = (1, 2), 
        tofill = WHITE_PIX
    )

    sfig(ChE, grid_,
        "bias_study.png";
    )

end;
