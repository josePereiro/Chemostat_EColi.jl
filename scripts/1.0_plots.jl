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
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChLP = Ch.LP
    using ProgressMeter
    using Base.Threads
    using SparseArrays

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

# ## ---------------------------------------------------------------------------------
# let
#     n = 100
#     vols = rand(n)
#     vars = rand(n)
#     xticks = UJL.get_ticks(log.(vols); l = 4) do vol
#         round(exp(vol))
#     end
#     p = scatter(log.(vols), vars; label = "", xticks)
#     mysavefig(p, "test") 
# end

# ---------------------------------------------------------------------------------
 function box_vol(model, idxs)
    try
        L, U = ChLP.fva(model, idxs)
        vol = prod(abs.(L .- U))
        max(0.0, log10(vol + 1e-50))
    catch; NaN end
end

# ---------------------------------------------------------------------------------
ALL_IDERS = ["AC", "CO2", "FORM", "GLC", "NH4", "O2", "SUCC"]

ider_colors = Dict(
    "AC" => :pink, 
    "CO2" => :blue, 
    "FORM" => :orange, 
    "GLC" => :black, 
    "NH4" => :orange, 
    "O2" => :red, 
    "SUCC" => :green
)

exp_colors = let
    EXPS = 1:20
    len = length(EXPS)
    colors = Plots.distinguishable_colors(len)
    Dict(ider => color for (ider, color) in zip(EXPS, colors))
end

source_markers = Dict(
    Hd => :circle, 
    Kd => :square, 
    Fd => :star, 
)

method_colors = Dict(
    :ME_MAX_POL => :blue,
    :FBA_Z_FIX_MIN_COST => :red
)

method_markers = Dict(
    :ME_MAX_POL => :circle,
    :FBA_Z_FIX_MIN_COST => :square
)


# ---------------------------------------------------------------------------------
# pol box size
# include("1.1_pol_box_size.jl")

## ---------------------------------------------------------------------------------
# var study
include("1.2_var_study.jl")

## ---------------------------------------------------------------------------------
