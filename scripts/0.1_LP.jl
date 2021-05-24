import DrWatson
const DW = DrWatson
DW.quickactivate(@__DIR__, "Chemostat_EColi")

# ----------------------------------------------------------------------------
# ARGS
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

    using Base.Threads
    using Serialization
    import UtilsJL
    const UJL = UtilsJL
    UJL.set_verbose(false)

    using Statistics
end

# ----------------------------------------------------------------------------
function fix_Z_and_readj_model(iJR, Ed, exp; 
        verbose = false, adjth = 0.001,
        proc_fun = (model) -> nothing
    )

    model = iJR.load_model("fva_models", exp) |> deepcopy
    biomider = iJR.BIOMASS_IDER
    biomidx = ChU.rxnindex(model, biomider)

    costider = iJR.COST_IDER
    exglcider = iJR.EX_GLC_IDER
    
    # ----------------------------------------------------
    # fix z
    exp_growth = Ed.val("D", exp)
    
    verbose && @info("Before re-adjustement $("-"^50)")
    ChU.bounds!(model, biomider, exp_growth, exp_growth)
    
    # ----------------------------------------------------
    # find max Yzvg
    vgL, vgU = ChLP.fva(model, exglcider) .|> first
    min_vg = min(abs(vgL), abs(vgU)) # vg is negative
    verbose && @show min_vg
    
    exp_Yzvg = abs(Ed.val("D", exp) / Ed.val("uGLC", exp))
    model_Yzvg = abs(exp_growth/ min_vg)
    verbose && @show exp_Yzvg, model_Yzvg
    verbose && println()

    # ----------------------------------------------------
    # readjust biom reaction
    Sfac = 0.98
    last_Sfac = Sfac # in case of error
    biomS = collect(model.S[:, biomidx]) # backup
    x0 = Sfac
    x1 = 1.0
    gdth = adjth
    minΔx = 0.01
    maxΔx = 0.05
    target = exp_Yzvg
    gdmaxiter = 1000

    function readjust_S!(gdmodel)
        last_Sfac = Sfac
        Sfac = UJL.gd_value(gdmodel)

        verbose && println()
        verbose && @info("After re-adjustement $("-"^50)", Sfac)
        try
            # adjust biom
            model.S[:, biomidx] .= biomS .* Sfac

            # find max Yzvg
            vgL, vgU = ChLP.fva(model, exglcider) .|> first
            verbose && @show vgL, vgU

            min_vg = min(abs(vgL), abs(vgU)) # vg is negative
            verbose && @show min_vg

            model_Yzvg = abs(exp_growth/ min_vg)
            verbose &&  @show exp_Yzvg, model_Yzvg

            ret = proc_fun(model)
            ret === true && 

            # break
            last_Sfac == Sfac && return target

            return model_Yzvg
        catch err
            verbose && println()
            verbose && @warn("Error, using previous last_Sfac", 
                last_Sfac, Sfac, 
                (vgL, vgU),
                (exp_Yzvg, model_Yzvg),
                err
            )

            # use last one
            Sfac = last_Sfac
            model.S[:, biomidx] .= biomS .* Sfac

            return target
        end
    end

    gdmodel = UJL.grad_desc(readjust_S!; 
        x0, x1, gdth, minΔx, maxΔx,
        target, maxiter = gdmaxiter, 
        verbose = false
    )

    proc_fun(model)

    verbose && @show Sfac
    return Sfac
end;

## ----------------------------------------------------------------------------------
#  dev
let
    return
    # F N H
    # 2 1 3

    # iJR = FiJR 
    # Ed = Fd
    # exp = 2 

    # iJR = NiJR 
    # Ed = Nd
    # exp = 1 

    iJR = HiJR 
    Ed = Hd
    exp = 3 

    src = string(nameof(Ed))
    LPDAT = UJL.DictTree()

    biomider = iJR.BIOMASS_IDER
    costider = iJR.COST_IDER
    exglcider = iJR.EX_GLC_IDER
    max_sense = -1.0
    min_sense = 1.0
   
    model = nothing
    fbaout = nothing
    # FBA_Z_FIX_ADJ_MIN_VG_MIN_COST
    function proc_fun(model_)
        model = model_
        # fix Z and readjust biom so exp and model yield (z/vg) match

        # min vg (vg is negative)
        fbaout1 = ChLP.fba(model, exglcider; sense = max_sense)
        exglc = ChU.av(model, fbaout1, exglcider)
        ChU.fixxing(model, exglcider, exglc) do
            # min cost
            fbaout = ChLP.fba(model, costider; sense = min_sense)
        end
    end
    Sfac = fix_Z_and_readj_model(iJR, Ed, exp; proc_fun, verbose = true);
    
    # store
    LPDAT[:FBA_Z_FIX_ADJ_MIN_VG_MIN_COST, :Sfac, exp] = Sfac
    LPDAT[:FBA_Z_FIX_ADJ_MIN_VG_MIN_COST, :model, exp] = ChU.compressed_model(model)
    LPDAT[:FBA_Z_FIX_ADJ_MIN_VG_MIN_COST, :fbaout, exp] = fbaout

end;

## ----------------------------------------------------------------------------
function do_LP(iJR, Ed, EXPS)

    src = string(nameof(Ed))
    LPDAT = UJL.DictTree()

    biomider = iJR.BIOMASS_IDER
    costider = iJR.COST_IDER
    exglcider = iJR.EX_GLC_IDER
    max_sense = -1.0
    min_sense = 1.0
 
    for exp in EXPS

        thid = threadid()
        @info("Doing ", exp, src, thid); println()
        
        # FBA_Z_FIX_ADJ_MIN_VG_MIN_COST
        let
            # fix Z and readjust biom so exp and model yield (z/vg) match
            model = nothing
            fbaout = nothing
            # FBA_Z_FIX_ADJ_MIN_VG_MIN_COST
            function proc_fun(model_)
                model = model_
                
                # min vg (vg is negative)
                fbaout1 = ChLP.fba(model, exglcider; sense = max_sense)
                exglc = ChU.av(model, fbaout1, exglcider)
                ChU.fixxing(model, exglcider, exglc) do
                    # min cost
                    fbaout = ChLP.fba(model, costider; sense = min_sense)
                end
            end
            # fix Z and readjust biom so exp and model yield (z/vg) match
            Sfac = fix_Z_and_readj_model(iJR, Ed, exp; proc_fun, verbose = false);
            
            # store
            LPDAT[:FBA_Z_FIX_ADJ_MIN_VG_MIN_COST, :Sfac, exp] = Sfac
            LPDAT[:FBA_Z_FIX_ADJ_MIN_VG_MIN_COST, :model, exp] = ChU.compressed_model(model)
            LPDAT[:FBA_Z_FIX_ADJ_MIN_VG_MIN_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_FIX_ADJ_MIN_VG_MAX_COST
        let
            # fix Z and readjust biom so exp and model yield (z/vg) match
            model = nothing
            fbaout = nothing
            # FBA_Z_FIX_ADJ_MIN_VG_MIN_COST
            function proc_fun(model_)
                model = model_
                
                # min vg (vg is negative)
                fbaout1 = ChLP.fba(model, exglcider; sense = max_sense)
                exglc = ChU.av(model, fbaout1, exglcider)
                ChU.fixxing(model, exglcider, exglc) do
                    # max cost
                    fbaout = ChLP.fba(model, costider; sense = max_sense)
                end
            end
            # fix Z and readjust biom so exp and model yield (z/vg) match
            Sfac = fix_Z_and_readj_model(iJR, Ed, exp; proc_fun, verbose = false);

            LPDAT[:FBA_Z_FIX_ADJ_MIN_VG_MAX_COST, :model, exp] = ChU.compressed_model(model)
            LPDAT[:FBA_Z_FIX_ADJ_MIN_VG_MAX_COST, :fbaout, exp] = fbaout
            LPDAT[:FBA_Z_FIX_ADJ_MIN_VG_MAX_COST, :Sfac, exp] = Sfac
        end

        # FBA_Z_FIX_MIN_VG_MIN_COST
        let
            model = iJR.load_model("fva_models", exp)

            # fix Z
            exp_growth = Ed.val("D", exp)
            ChU.bounds!(model, biomider, exp_growth, exp_growth)
            
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
            ChU.bounds!(model, biomider, exp_growth, exp_growth)
            
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
            ChU.bounds!(model, biomider, exp_growth, exp_growth)
            
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
            ChU.bounds!(model, biomider, exp_growth, exp_growth)
            
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
            fbaout1 = ChLP.fba(model, biomider; sense = max_sense)
            objval = ChU.av(model, fbaout1, biomider)
            ChU.bounds!(model, biomider, objval, objval)

            # min cost
            fbaout = ChLP.fba(model, costider; sense = min_sense)
     
            LPDAT[:FBA_MAX_Z_MIN_COST, :model, exp] = ChU.compressed_model(model)
            LPDAT[:FBA_MAX_Z_MIN_COST, :fbaout, exp] = fbaout
        end

        # FBA_MAX_Z_MAX_COST
        let
            model = iJR.load_model("fva_models", exp)
            
            # max z
            fbaout1 = ChLP.fba(model, biomider; sense = max_sense)
            objval = ChU.av(model, fbaout1, biomider)
            ChU.bounds!(model, biomider, objval, objval)

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
            ChU.bounds!(model, biomider, exp_growth, exp_growth)

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
            ChU.bounds!(model, biomider, exp_growth, exp_growth)

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
            ChU.bounds!(model, biomider, exp_growth, exp_growth)

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
            ChU.bounds!(model, biomider, exp_growth, exp_growth)

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
    @threads for (iJR, Ed, EXPS) in [
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