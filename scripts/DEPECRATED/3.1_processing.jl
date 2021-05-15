import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Folsom2014")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    # -------------------------------------------------------------------
    import Chemostat_Folsom2014
    const ChF = Chemostat_Folsom2014

    const iJR = ChF.iJR904
    const Fd = ChF.FolsomData # experimental data
    const Bd = ChF.BegData    # cost data

    #  ----------------------------------------------------------------------------
    # run add "https://github.com/josePereiro/Chemostat" in the 
    # julia Pkg REPL for installing the package
    import Chemostat
    import Chemostat.LP.MathProgBase
    const Ch = Chemostat
    const ChU = Ch.Utils

    using Serialization
    import UtilsJL
    const UJL = UtilsJL
    
    import FileIO
    using Plots
    import GR
    GR.inline("png")
end

## -----------------------------------------------------------------------------------------------
LPDAT = ChU.load_data(iJR.LP_DAT_FILE)

const FBA_Z_FIX_MIN_COST = :FBA_Z_FIX_MIN_COST
const FBA_MAX_Z_MIN_COST = :FBA_MAX_Z_MIN_COST
const FBA_Z_FIX_MIN_VG_COST = :FBA_Z_FIX_MIN_VG_COST
const FBA_Z_VG_FIX_MIN_COST = :FBA_Z_VG_FIX_MIN_COST

# -----------------------------------------------------------------------------------------------
fileid = "3.1"
function mysavefig(p, pname; params...) 
    fname = UJL.mysavefig(p, string(fileid, "_", pname), iJR.MODEL_FIGURES_DIR; params...)
    @info "Plotting" fname
end
myminmax(a::Vector) = (minimum(a), maximum(a))  
FLX_IDERS = ["GLC", "PYR", "SUCC", "LAC", "FORM", "AC", "O2", "CO2"]

EXPS = 1:4 # experiments that have both concentration and flx data

exp_colors = let
    colors = Plots.distinguishable_colors(length(EXPS))
    Dict(exp => color for (exp, color) in zip(EXPS, colors))
end

ider_colors = let
    iders = [FLX_IDERS; "D"]
    colors = Plots.distinguishable_colors(length(iders))
    Dict(ider => color for (ider, color) in zip(iders, colors))
end

method_colors = Dict(
    FBA_Z_FIX_MIN_COST => :red,
    FBA_MAX_Z_MIN_COST => :orange,
    FBA_Z_FIX_MIN_VG_COST => :blue,
    FBA_Z_VG_FIX_MIN_COST => :purple,
)

ALL_METHODS = [
    FBA_Z_FIX_MIN_COST,
    FBA_MAX_Z_MIN_COST, 
    FBA_Z_FIX_MIN_VG_COST,
    FBA_Z_VG_FIX_MIN_COST
]

Fd_mets_map = iJR.load_mets_map()
Fd_rxns_map = iJR.load_rxns_map()

## -----------------------------------------------------------------------------------------------
# correlations
DAT = UJL.DictTree()
DAT[:FLX_IDERS] = FLX_IDERS
DAT[:EXPS] = EXPS

## -----------------------------------------------------------------------------------------------
# Flx correlations
let
    ps = Plots.Plot[]
    for method in ALL_METHODS
        p = plot(;title = string(method), xlabel = "exp flx", ylabel = "model flx")
        for Fd_ider in FLX_IDERS
            model_ider = Fd_rxns_map[Fd_ider]
            for exp in EXPS
                color = ider_colors[Fd_ider]
                # every body is possitive here
                Fd_flx = abs(Fd.uval(Fd_ider, exp))

                model = LPDAT[method, :model, exp]
                fbaout = LPDAT[method, :fbaout, exp]
                
                fba_flx = abs(ChU.av(model, fbaout, model_ider))
                DAT[method, :Fd, :flx, Fd_ider, exp] = Fd_flx
                DAT[method, :lp, :flx, Fd_ider, exp] = fba_flx

                scatter!(p, [Fd_flx], [fba_flx]; 
                    ms = 8, color, label = ""
                )
            end
        end
        xs = DAT[method, [:Fd, :lp], :flx, FLX_IDERS, EXPS] |> sort
        plot!(p, xs, xs; label = "", ls = :dash, 
            alpha = 0.9, color = :black, lw = 3
        )
        push!(ps, p)
    end
    
    pname = "flx_tot_corr"
    mysavefig(ps, pname)

end

## -------------------------------------------------------------------
# Inter project comunication
let
    CORR_DAT = isfile(iJR.CORR_DAT_FILE) ? ChU.load_data(iJR.CORR_DAT_FILE) : Dict()
    CORR_DAT[:LP] = DAT
    ChU.save_data(iJR.CORR_DAT_FILE, CORR_DAT)
end

## -------------------------------------------------------------------
# # # leyends
# # # TODO fix this...
# # let
# #     for (title, colors) in [
# #             ("exp", exp_colors), 
# #             ("iders", ider_colors),
# #             ("method", method_colors)
# #         ]
# #     p = plot(; framestyle = :none)
# #         scolors = sort(collect(colors); by = (p) -> string(first(p)))
# #         for (id, color) in scolors
# #             scatter!(p, [0], [0];
# #                 thickness_scaling = 1,
# #                 color, ms = 8, label = string(id),
# #                 legendfontsize=10, 
# #                 # size = [300, 900],
# #                 # legend = :left
# #             )
# #         end
# #         mysavefig(p, "$(title)_color_legend")
# #     end
# # end