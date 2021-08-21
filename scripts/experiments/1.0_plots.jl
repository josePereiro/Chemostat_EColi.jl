using Plots: length
using Serialization: include
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
myminmax(a) = isempty(a) ? (0.0, 0.0) : (minimum(a), maximum(a))

function box_vol(model, idxs)
    try
        L, U = ChLP.fva(model, idxs)
        vol = prod(abs.(L .- U))
        max(0.0, log10(vol + 1e-50))
    catch; NaN end
end

## ---------------------------------------------------------------------------------
LP_METHODS = [
    :FBA_Z_FIX_MAX_VG_MIN_COST, :FBA_Z_FIX_MAX_VG_MAX_COST, 
    :FBA_Z_FIX_MIN_VG_MIN_COST, :FBA_Z_FIX_MIN_VG_MAX_COST,
    :FBA_Z_VG_FIX_MAX_COST, :FBA_Z_VG_FIX_MIN_COST, :FBA_Z_FIX_MAX_COST, 
    :FBA_Z_FIX_MIN_COST, :FBA_MAX_Z_MIN_COST, :FBA_MAX_Z_MAX_COST, 
    :FBA_Z_FIX_ADJ_MIN_VG_MIN_COST, :FBA_Z_FIX_ADJ_MIN_VG_MAX_COST,
    :FBA_Z_FIX_MAX_ATP_MIN_COST, :FBA_Z_FIX_MAX_ATP_MAX_COST
]

ME_METHODS = [
    :ME_MAX_POL,
    :ME_MAX_POL_B0,
    :ME_Z_EXPECTED_G_BOUNDED, 
]

## ---------------------------------------------------------------------------------
ALL_IDERS = ["AC", "CO2", "FORM", "GLC", "NH4", "O2", "SUCC"]

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

# [:none, :auto, :circle, :rect, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, 
# :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star6, 
# :star7, :star8, :vline, :hline, :+, :x].
subs_markers = Dict(
    "glycolysis" => :circle,
    "krebs" => :diamond,
    "pentose phosphate" => :utriangle,
    "others" => :star,
    "glyoxylate shunt" => :square,
    "glc tranport" => :hex
)

subs_colors = Dict(
    "glycolysis" => :red,
    "krebs" => :blue,
    "pentose phosphate" => :yellow,
    "others" => :purple,
    "glyoxylate shunt" => :brown,
    "glc tranport" => :black
)

source_markers = Dict(
    Kd => :square, 
    Fd => :circle, 
    Nd => :utriangle, 
)

source_labels = Dict(
    Kd => "Kayser", 
    Nd => "Nanchen", 
    Fd => "Folsom", 
)

# method_colors = Dict(
#     :ME_MAX_POL => :blue,
#     :ME_Z_EXPECTED_G_BOUNDED => :purple,
#     :FBA_Z_FIX_MIN_COST => :red,
#     :FBA_Z_FIX_MAX_COST => :red,
#     :FBA_Z_FIX_MIN_VG_MAX_COST => :red,
#     :FBA_Z_FIX_MIN_VG_MIN_COST => :red,
# )

method_colors = let
    dict_ = Dict{Symbol, Any}(
        :ME_MAX_POL => :blue,
        :ME_Z_EXPECTED_G_BOUNDED => :purple,
    )
    colors = Plots.distinguishable_colors(length(LP_METHODS))
    for (method, color) in zip(LP_METHODS, colors)
        dict_[method] = color
    end
    dict_[:FBA_Z_FIX_MIN_VG_MAX_COST] = :red
    dict_[:FBA_Z_FIX_MIN_VG_MIN_COST] = :green
    dict_
end

method_markers = Dict(
    :ME_MAX_POL => :circle,
    :ME_Z_EXPECTED_G_BOUNDED => :circle,
    :FBA_Z_FIX_MIN_COST => :square,
    :FBA_Z_FIX_MAX_COST => :square,
    :FBA_Z_FIX_MIN_VG_MAX_COST => :square,
    :FBA_Z_FIX_MIN_VG_MIN_COST => :square,
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
let
    p = plot(1:10)
    annotate!(p, [(7,3,"(7,3)"),(3,7,text("hey", 14, :left, :top, :green))])
    annotate!(p, [(4, 4, ("More text", 8, 45.0, :bottom, :red))])
    p
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
