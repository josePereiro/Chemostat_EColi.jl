import DrWatson
const DW = DrWatson
DW.quickactivate(@__DIR__, "Chemostat_EColi")

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

    import Chemostat_Heerden2013
    const ChH = Chemostat_Heerden2013
    const HiJR = ChH.iJR904
    const Hd = ChH.HeerdenData

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

    using Serialization
    import UtilsJL
    const UJL = UtilsJL
    UJL.set_verbose(false)

    using Statistics
    using Plots
    import GR
    GR.inline("png")
    
end

# ---------------------------------------------------------------------------------
fileid = "1.0"
function mysavefig(p, pname; params...) 
    name = string(fileid, "_", pname)
    fname = UJL.mysavefig(p, name, ChE.plotsdir(); params...)
    @info "Plotting" fname
end
myminmax(a) = isempty(a) ? (0.0, 0.0) : (minimum(a), maximum(a))

# ---------------------------------------------------------------------------------
 function box_vol(model, idxs)
    try
        L, U = ChLP.fva(model, idxs)
        vol = prod(abs.(L .- U))
        max(0.0, log10(vol + 1e-50))
    catch; NaN end
end

# ---------------------------------------------------------------------------------
# ME methods
const ME_Z_OPEN_G_OPEN          = :ME_Z_OPEN_G_OPEN
const ME_MAX_POL                = :ME_MAX_POL
const ME_MAX_POL_B0             = :ME_MAX_POL_B0
const ME_Z_EXPECTED_G_MOVING    = :ME_Z_EXPECTED_G_MOVING
const ME_Z_EXPECTED_G_BOUNDED   = :ME_Z_EXPECTED_G_BOUNDED
const ME_Z_FIXXED_G_BOUNDED     = :ME_Z_FIXXED_G_BOUNDED

# LP methods
const FBA_Z_FIX_MAX_VG_MIN_COST = :FBA_Z_FIX_MAX_VG_MIN_COST
const FBA_Z_FIX_MIN_VG_MIN_COST = :FBA_Z_FIX_MIN_VG_MIN_COST
const FBA_Z_FIX_MAX_VG_MAX_COST = :FBA_Z_FIX_MAX_VG_MAX_COST
const FBA_Z_FIX_MIN_VG_MAX_COST = :FBA_Z_FIX_MIN_VG_MAX_COST
const FBA_Z_VG_FIX_MAX_COST = :FBA_Z_VG_FIX_MAX_COST
const FBA_Z_VG_FIX_MIN_COST = :FBA_Z_VG_FIX_MIN_COST
const FBA_Z_FIX_MAX_COST = :FBA_Z_FIX_MAX_COST
const FBA_Z_FIX_MIN_COST = :FBA_Z_FIX_MIN_COST
const FBA_MAX_Z_MIN_COST = :FBA_MAX_Z_MIN_COST
const FBA_MAX_Z_MAX_COST = :FBA_MAX_Z_MAX_COST

LP_METHODS = [
    FBA_Z_FIX_MAX_VG_MIN_COST, FBA_Z_FIX_MAX_VG_MAX_COST, 
    FBA_Z_FIX_MIN_VG_MIN_COST, FBA_Z_FIX_MIN_VG_MAX_COST,
    FBA_Z_VG_FIX_MAX_COST, FBA_Z_VG_FIX_MIN_COST, FBA_Z_FIX_MAX_COST, 
    FBA_Z_FIX_MIN_COST, FBA_MAX_Z_MIN_COST, FBA_MAX_Z_MAX_COST
]

ME_METHODS = [
    # ME_Z_OPEN_G_OPEN, 
    ME_MAX_POL,
    ME_MAX_POL_B0,
    # ME_Z_FIXXED_G_BOUNDED, 
    ME_Z_EXPECTED_G_BOUNDED, 
    # ME_Z_EXPECTED_G_MOVING
]

# ---------------------------------------------------------------------------------
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
    # Hd => :circle, 
    Kd => :square, 
    Fd => :circle, 
    Nd => :utriangle, 
)

source_labels = Dict(
    Kd => "Kayser", 
    Nd => "Nanchen", 
    # Hd => "Heerden", 
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
    dict_[FBA_Z_FIX_MIN_VG_MAX_COST] = :red
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

## ---------------------------------------------------------------------------------
# pol box size
# include("1.1_pol_box_size.jl")

## ---------------------------------------------------------------------------------
# var study
# include("1.2_var_study.jl")

## ---------------------------------------------------------------------------------
