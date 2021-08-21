using ProjAssistant
@quickactivate

# ----------------------------------------------------------------------------
## ARGS
using ArgParse

set = ArgParseSettings()
@add_arg_table! set begin
    "--new-dat"
        help = "ignore disk stored DAT"   
        action = :store_true
    "--skip-me"
        help = "do not recompute MaxEnt part"   
        action = :store_true
    "--skip-lp"
        help = "do not recompute LP part"   
        action = :store_true
    "--dry-run"
        help = "do all but change nothing on disk"   
        action = :store_true
end

if isinteractive()
    # Dev values
    new_dat = false
    skip_lp = false
    skip_me = false
    dry_run = false
else
    parsed_args = parse_args(set)
    new_dat = parsed_args["new-dat"]
    skip_lp = parsed_args["skip-lp"]
    skip_me = parsed_args["skip-me"]
    dry_run = parsed_args["dry-run"]
end

# ----------------------------------------------------------------------------
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

    using Serialization
    using MagicDicts

    using Statistics
end

## ----------------------------------------------------------------------------------
include("0.2.1_collect_commom.jl")
include("0.2.2_collect_LP.jl")
include("0.2.3_collect_ME.jl")

# ----------------------------------------------------------------------------------
let
    nths = nthreads()
    for (iJR, Ed, EXPS, EXCH_FLX_IDERS, INNER_FLX_IDERS) in [
        (NiJR, Nd, 1:17, ["GLC", "AC"], NiJR.load_inner_iders()),
        (KiJR, Kd, 5:12, ["GLC", "CO2", "O2", "AC", "NH4"], []),
        (FiJR, Fd, 1:4, ["GLC", "SUCC", "AC", "FORM"], []),
    ]
        src = string(nameof(Ed))
        
        # ---------------------------------------------------------------------------------
        DAT_FILE = ChE.DATfile(src)
        DAT = (!isfile(DAT_FILE) || new_dat) ? MagicDict() : ChE.load_DAT(src)

        DAT[:EXPS] = EXPS
        DAT[:EXCH_FLX_IDERS] = EXCH_FLX_IDERS
        DAT[:INNER_FLX_IDERS] = INNER_FLX_IDERS

        # LP methods
        DAT[:LP_METHODS] = [
            :FBA_Z_FIX_ADJ_MIN_VG_MIN_COST,
            :FBA_Z_FIX_ADJ_MIN_VG_MAX_COST,
            :FBA_Z_FIX_ADJ_MAX_VG_MIN_COST,
            :FBA_Z_FIX_ADJ_MAX_VG_MAX_COST,
            :FBA_Z_FIX_MIN_VG_MIN_COST,
            :FBA_Z_FIX_MIN_VG_MAX_COST,
            :FBA_Z_FIX_MAX_VG_MIN_COST,
            :FBA_Z_FIX_MAX_VG_MAX_COST,
            :FBA_MAX_Z_MIN_COST,
            :FBA_MAX_Z_MAX_COST,
            :FBA_Z_FIX_MIN_COST,
            :FBA_Z_FIX_MAX_COST,
            :FBA_Z_VG_FIX_MIN_COST,
            :FBA_Z_VG_FIX_MAX_COST,
            :FBA_Z_FIX_MAX_ATP_MIN_COST,
            :FBA_Z_FIX_MAX_ATP_MAX_COST
        ]

        # ME methods
        DAT[:ME_METHODS] = [
            :ME_MAX_POL,
            :ME_MAX_POL_B0,
            :ME_Z_EXPECTED_G_BOUNDED,
        ]

        DAT[:EXP_INDEPENDENT] = [:ME_MAX_POL_B0]
        
        # ---------------------------------------------------------------------------------
        collect_commom_DAT!(DAT, Ed, iJR)
        skip_lp || collect_lp_DAT!(DAT, Ed, iJR)
        skip_me || collect_ME_DAT!(DAT, iJR; nths)
        
        # ---------------------------------------------------------------------------------
        # saving
        dry_run || sdat(DAT, DAT_FILE)
    end

end;