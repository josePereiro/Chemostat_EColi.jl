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
    using Plots.Measures

    using LaTeXStrings

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
function _textbf(strs...)
   strs = replace.(string.(strs), " " => "~")
   latexstring(string("\\mathbf{", strs..., "}"))
end

## ---------------------------------------------------------------------------------
# figure 8
include("1.1_pol_box_size.jl")
include("1.2_Xexp_Xmax_study.jl")
include("1.3_residual_Glc_study.jl")

let
    # panel 1
    # TODO: Add colorbar_title
    Dlim = [0.01, 0.5]
    cgD_Xlim = [0.0, 40.0]
    bins = 60
    panel1 = plot_pol_box_size(
        ChN.iJR904, ChN.NanchenData; 
        Dlim, cgD_Xlim, bins
    )

    # panel 2
    panel2 = plot_Xexp_Xmax_study()
    
    # panel 3
    panel3 = plot_residual_Glc_study()

    for (title, p) in [
            ("A", panel1), 
            ("B", panel2), 
            ("C", panel3)
        ]
        plot!(p;
            title,
            titlefont = font(22),
            legendfont = font(18),
            guidefont = font(22),
            xtickfont = font(21),
            ytickfont = font(21),
            thickness_scaling = 1.6,
            size = (1320, 940)
        )
    end
    
    sfig(ChE, [panel1, panel2, panel3], 
        "figure_8", ".png";
        layout = (1, 3)
    )

end

## ---------------------------------------------------------------------------------
# figure 9 and 10
include("1.4_inner_flx_corr.jl")
include("1.4a_inner_flx_ERR.jl")
include("1.5_exchs_flx_correlations.jl")

let

    METHODS = [
        :FBA_MAX_Z_MIN_COST,
        :FBA_Z_FIX_MAX_VG_MIN_COST, 
        :FBA_Z_FIX_MIN_VG_MIN_COST, 
        :FBA_Z_FIX_MAX_ATP_MAX_COST,
        :ME_MAX_POL
    ]

    METHODS_LABELS = Dict(
        :FBA_MAX_Z_MIN_COST => "FBA max. z",
        :FBA_Z_FIX_MAX_VG_MIN_COST => "FBA min. ug",
        :FBA_Z_FIX_MIN_VG_MIN_COST => "FBA max. ug",
        :FBA_Z_FIX_MAX_ATP_MAX_COST => "FBA max. ATP",
        :ME_MAX_POL => "ME",
    )

    METHODS_SHAPES = Dict(
        :FBA_MAX_Z_MIN_COST => :utriangle,
        :FBA_Z_FIX_MAX_VG_MIN_COST => :dtriangle,
        :FBA_Z_FIX_MIN_VG_MIN_COST => :rtriangle,
        :FBA_Z_FIX_MAX_ATP_MAX_COST => :ltriangle,
        :ME_MAX_POL => :circle,
    )

    METHODS_COLORS = Dict(
        :FBA_MAX_Z_MIN_COST => :red,
        :FBA_Z_FIX_MAX_VG_MIN_COST => :pink, 
        :FBA_Z_FIX_MIN_VG_MIN_COST => :yellow, 
        :FBA_Z_FIX_MAX_ATP_MAX_COST => :brown,
        :ME_MAX_POL => :blue
    )

    ## ------------------------------------------------------
    # figure_9

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

    for p in ps
        plot!(p; 
            titlefont = font(26),
            axisfont = font(26),
            guidefont = font(26),
            xtickfont = font(22),
            ytickfont = font(22),
            legendfont = font(12),
            thickness_scaling = 1.6,
            size = (1220, 940)
        )
    end

    # add margin
    for ps in [exch_ps, values(inner_ps)...]
        plot!(first(ps); leftmargin = 5mm)
        plot!(last(ps); rightmargin = 7mm)
    end

    sfig(ChE, ps,
        "figure_9", ".png";
        layout = (7, 5)
    )

    ## ------------------------------------------------------
    # figure_10

    ps = _plot_inner_fluxs_ERR(METHODS, METHODS_LABELS, METHODS_COLORS, METHODS_SHAPES, :ME_MAX_POL)

    tot_p = ps[:tot]
    plot!(tot_p; 
        guidefont = font(16),
        xtickfont = font(13),
        ytickfont = font(13),
        legendfont = font(12),
        thickness_scaling = 1.7,
        size = (1700, 900),
        bottom_margin = 7mm
    )

    sfig(ChE, tot_p,
        "figure_10", ".png"
    )

end

## ---------------------------------------------------------------------------------
# figure 12
include("1.6_var_study_and_marginals.jl")

let    
    corr_av_p, corr_va_p = plot_ME_methods_comparizon()
    pGLC, pAC, pz = plot_me_marginals()

    size0 = (1220, 940)
    for p in [pGLC, pAC, pz]
        plot!(p; 
            titlefont = font(26),
            axisfont = font(26),
            xtickfont = font(22),
            ytickfont = font(22),
            legendfont = font(22),
            guidefont = font(26),
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
            titlefont = font(26),
            axisfont = font(26),
            guidefont = font(26),
            xtickfont = font(22),
            ytickfont = font(22),
            legendfont = font(22),
            thickness_scaling = 1.6,
            size = 1.5 .* size0
        )
    end

    corrs_file = sfig(ChE, [corr_av_p, corr_va_p],
        "corrs.png";
        layout = (2, 1)
    )

    sfig(ChE, lfig.([margs_file, corrs_file]),
        "figure_12", ".png";
        layout = (1, 2)
    )

    rm.([margs_file, corrs_file]; force = true)

    nothing
end
