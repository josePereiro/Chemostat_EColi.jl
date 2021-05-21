import DrWatson
const DW = DrWatson
DW.quickactivate(@__DIR__, "Chemostat_EColi")

# ----------------------------------------------------------------------------
## ARGS
using ArgParse

set = ArgParseSettings()
@add_arg_table! set begin
    "--dry-run"
        help = "do all but change nothing on disk"   
        action = :store_true
end

if isinteractive()
    # Dev values
    dry_run = false
else
    parsed_args = parse_args(set)
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
end

# ----------------------------------------------------------------------------
function do_LP(iJR, Ed, EXPS)

    src = string(nameof(Ed))
    LPDAT = UJL.DictTree()

    objider = iJR.BIOMASS_IDER
    costider = iJR.COST_IDER
    exglcider = iJR.EX_GLC_IDER
    max_sense = -1.0
    min_sense = 1.0
 
    for exp in EXPS

        @info("Doing ", exp, src); println()
        
        # FBA_Z_FIX_MIN_VG_MIN_COST
        let
            model = iJR.load_model("fva_models", exp)

            # fix Z
            exp_growth = Ed.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            
            # min vg (vg is negative)
            fbaout1 = ChLP.fba(model, exglcider; sense = max_sense)
            exglc = ChU.av(model, fbaout1, exglcider)
            ChU.bounds!(model, exglcider, exglc, exglc)

            # min cost
            fbaout = ChLP.fba(model, costider; sense = min_sense)

            LPDAT[:FBA_Z_FIX_MIN_VG_MIN_COST, :model, exp] = ChU.compressed_model(model)
            LPDAT[:FBA_Z_FIX_MIN_VG_MIN_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_FIX_MIN_VG_MAX_COST
        let
            model = iJR.load_model("fva_models", exp)

            # fix Z
            exp_growth = Ed.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            
            # min vg (vg is negative)
            fbaout1 = ChLP.fba(model, exglcider; sense = max_sense)
            exglc = ChU.av(model, fbaout1, exglcider)
            ChU.bounds!(model, exglcider, exglc, exglc)

            # max cost
            fbaout = ChLP.fba(model, costider; sense = max_sense)

            LPDAT[:FBA_Z_FIX_MIN_VG_MAX_COST, :model, exp] = ChU.compressed_model(model)
            LPDAT[:FBA_Z_FIX_MIN_VG_MAX_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_FIX_MAX_VG_MIN_COST
        let
            model = iJR.load_model("fva_models", exp)

            # fix Z
            exp_growth = Ed.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            
            # max vg (vg is negative)
            fbaout1 = ChLP.fba(model, exglcider; sense = min_sense)
            exglc = ChU.av(model, fbaout1, exglcider)
            ChU.bounds!(model, exglcider, exglc, exglc)

            # min cost
            fbaout = ChLP.fba(model, costider; sense = min_sense)

            LPDAT[:FBA_Z_FIX_MAX_VG_MIN_COST, :model, exp] = ChU.compressed_model(model)
            LPDAT[:FBA_Z_FIX_MAX_VG_MIN_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_FIX_MAX_VG_MAX_COST
        let
            model = iJR.load_model("fva_models", exp)
            
            # fix Z
            exp_growth = Ed.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            
            # max vg (vg is negative)
            fbaout1 = ChLP.fba(model, exglcider; sense = min_sense)
            exglc = ChU.av(model, fbaout1, exglcider)
            ChU.bounds!(model, exglcider, exglc, exglc)

            # max cost
            fbaout = ChLP.fba(model, costider; sense = max_sense)

            LPDAT[:FBA_Z_FIX_MAX_VG_MAX_COST, :model, exp] = ChU.compressed_model(model)
            LPDAT[:FBA_Z_FIX_MAX_VG_MAX_COST, :fbaout, exp] = fbaout
        end

        # FBA_MAX_Z_MIN_COST
        let
            model = iJR.load_model("fva_models", exp)
            
            # max z
            fbaout1 = ChLP.fba(model, objider; sense = max_sense)
            objval = ChU.av(model, fbaout1, objider)
            ChU.bounds!(model, objider, objval, objval)

            # min cost
            fbaout = ChLP.fba(model, costider; sense = min_sense)
     
            LPDAT[:FBA_MAX_Z_MIN_COST, :model, exp] = ChU.compressed_model(model)
            LPDAT[:FBA_MAX_Z_MIN_COST, :fbaout, exp] = fbaout
        end

        # FBA_MAX_Z_MAX_COST
        let
            model = iJR.load_model("fva_models", exp)
            
            # max z
            fbaout1 = ChLP.fba(model, objider; sense = max_sense)
            objval = ChU.av(model, fbaout1, objider)
            ChU.bounds!(model, objider, objval, objval)

            # max cost
            fbaout = ChLP.fba(model, costider; sense = max_sense)
     
            LPDAT[:FBA_MAX_Z_MAX_COST, :model, exp] = ChU.compressed_model(model)
            LPDAT[:FBA_MAX_Z_MAX_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_FIX_MIN_COST
        let
            model = iJR.load_model("fva_models", exp)
            exp_growth = Ed.val("D", exp)
            
            # fix Z
            exp_growth = Ed.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)

            # min cost
            fbaout = ChLP.fba(model, costider; sense = min_sense)

            LPDAT[:FBA_Z_FIX_MIN_COST, :model, exp] = ChU.compressed_model(model)
            LPDAT[:FBA_Z_FIX_MIN_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_FIX_MAX_COST
        let
            model = iJR.load_model("fva_models", exp)
            exp_growth = Ed.val("D", exp)
            
            # fix Z
            exp_growth = Ed.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)

            # max cost
            fbaout = ChLP.fba(model, costider; sense = max_sense)

            LPDAT[:FBA_Z_FIX_MAX_COST, :model, exp] = ChU.compressed_model(model)
            LPDAT[:FBA_Z_FIX_MAX_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_VG_FIX_MIN_COST
        let
            model = iJR.load_model("fva_models", exp)

            # fix Z
            exp_growth = Ed.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)

            # fix vg
            exp_exglc = -abs(Ed.uval("GLC", exp))
            ChU.bounds!(model, exglcider, exp_exglc, exp_exglc)

            # min cost
            fbaout = ChLP.fba(model, costider; sense = min_sense)

            LPDAT[:FBA_Z_VG_FIX_MIN_COST, :model, exp] = ChU.compressed_model(model)
            LPDAT[:FBA_Z_VG_FIX_MIN_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_VG_FIX_MAX_COST
        let
            model = iJR.load_model("fva_models", exp)
            
            # fix Z
            exp_growth = Ed.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)

            # fix vg
            exp_exglc = -abs(Ed.uval("GLC", exp))
            ChU.bounds!(model, exglcider, exp_exglc, exp_exglc)

            # max cost
            fbaout = ChLP.fba(model, costider; sense = max_sense)

            LPDAT[:FBA_Z_VG_FIX_MAX_COST, :model, exp] = ChU.compressed_model(model)
            LPDAT[:FBA_Z_VG_FIX_MAX_COST, :fbaout, exp] = fbaout
        end

    end

    return LPDAT
end

## ----------------------------------------------------------------------------------
let
    for (iJR, Ed, EXPS) in [
        (NiJR, Nd, Nd.EXPS),
        (KiJR, Kd, 5:13),
        (FiJR, Fd, 1:4),
        (HiJR, Hd, 1:13),
    ]
        src = string(nameof(Ed))

        # ---------------------------------------------------------------------------------
        LP_DAT = do_LP(iJR, Ed, EXPS)
        
        # ---------------------------------------------------------------------------------
        # saving
        dfile = ChE.LP_DATfile(src)
        dry_run || serialize(dfile, LP_DAT)
    end

end;