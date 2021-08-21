using ProjAssistant
@quickactivate

# ----------------------------------------------------------------------------
# ARGS
using ArgParse

set = ArgParseSettings()
@add_arg_table! set begin
    "--dry-run"
        help = "do all but change nothing on disk"   
        action = :store_true
    "--redo"
        help = "ignore disk version (if exist)"   
        action = :store_true
end

if isinteractive()
    # Dev values
    dry_run = false
    redo = true
else
    parsed_args = parse_args(set)
    dry_run = parsed_args["dry-run"]
    redo = parsed_args["redo"]
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
    const ChSS = Ch.SteadyState
    using ProgressMeter
    using Base.Threads
    using SparseArrays

    import SimTools
    const SimT = SimTools

    using Base.Threads
    using Serialization

    using MagicDicts

    using Statistics
end

# ----------------------------------------------------------------------------
function fix_Z_and_readj_model(iJR, Ed, exp; 
        btol = 0.0,
        verbose = false, adjth = 0.001,
        proc_fun = (model) -> nothing
    )

    model = iJR.load_model("fva_models", exp) |> deepcopy
    biomider = iJR.BIOMASS_IDER
    biomidx = ChU.rxnindex(model, biomider)
    exglcider = iJR.EX_GLC_IDER
    
    # ----------------------------------------------------
    # fix z
    exp_growth = Ed.val("D", exp)
    
    verbose && @info("Before re-adjustement $("-"^50)")
    dflx = abs(exp_growth * btol)
    ChU.bounds!(model, biomider, exp_growth - dflx, exp_growth + dflx)
    
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
        Sfac = SimT.gd_value(gdmodel)

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

    gdmodel = SimT.grad_desc(readjust_S!; 
        x0, x1, gdth, minΔx, maxΔx,
        target, maxiter = gdmaxiter, 
        verbose = false
    )

    proc_fun(model)

    verbose && @show Sfac
    return Sfac
end;

## ----------------------------------------------------------------------------
const wlk = ReentrantLock()
thsv_set!(LPDAT, ks, v) = lock(() -> (LPDAT[ks...] = v), wlk)

## ----------------------------------------------------------------------------
function do_LP!(LPDAT, iJR, Ed, exp)

    biomider = iJR.BIOMASS_IDER
    costider = iJR.COST_IDER
    exglcider = iJR.EX_GLC_IDER
    atpmider = iJR.ATPM_IDER
    max_sense = ChLP.MAX_SENSE
    min_sense = ChLP.MIN_SENSE
    btol = 0.0


    # FBA_Z_FIX_MAX_ATP_MINMAX_COST
    let
        for (FBAmethod, atpm_sense, cost_sense) in [
            (:FBA_Z_FIX_MAX_ATP_MIN_COST, max_sense, min_sense),
            (:FBA_Z_FIX_MAX_ATP_MAX_COST, max_sense, max_sense),
        ]

            !redo && haskey(LPDAT, [FBAmethod, :model, exp]...) && continue

            model = iJR.load_model("base_model")
            exp_xi = Ed.val(:xi, exp)
            exp_intake_info = iJR.intake_info(exp)
            ChSS.apply_bound!(model, exp_xi, exp_intake_info; emptyfirst = true)

            # fix Z
            exp_growth = Ed.val("D", exp)
            dflx = abs(exp_growth * btol)
            ChU.bounds!(model, biomider, exp_growth - dflx, exp_growth + dflx)

            # max atpm
            # rxn[139]: ATPM (ATP maintenance requirement)
            # lb: 7.6, ub: 1000.0
            # (-1.0) atp_c + (-1.0) h2o_c ==> (1.0) adp_c + (1.0) h_c + (1.0) pi_c
            ChU.ub!(model, atpmider, iJR.ABS_MAX_BOUND)
            fbaout1 = ChLP.fba(model, atpmider; sense = atpm_sense)
            atpm = ChU.av(model, fbaout1, atpmider)
            dflx = abs(atpm * btol)
            ChU.bounds!(model, atpmider, atpm - dflx, atpm + dflx)

            # min cost
            fbaout = ChLP.fba(model, costider; sense = cost_sense)

            thsv_set!(LPDAT, [FBAmethod, :model, exp], ChU.compressed_model(model))
            thsv_set!(LPDAT, [FBAmethod, :fbaout, exp], fbaout)
        end
    end

    # FBA_Z_FIX_ADJ_MINMAX_VG_MINMAX_COST
    let
        for (FBAmethod, glc_sense, cost_sense) in [
            (:FBA_Z_FIX_ADJ_MIN_VG_MIN_COST, max_sense, min_sense),
            (:FBA_Z_FIX_ADJ_MIN_VG_MAX_COST, max_sense, max_sense),
            (:FBA_Z_FIX_ADJ_MAX_VG_MIN_COST, min_sense, min_sense),
            (:FBA_Z_FIX_ADJ_MAX_VG_MAX_COST, min_sense, max_sense),
        ]

            !redo && haskey(LPDAT, [FBAmethod, :model, exp]...) && continue

            # fix Z and readjust biom so exp and model yield (z/vg) match
            model = nothing
            fbaout = nothing
            function proc_fun(model_)
                model = model_
                
                # min vg (vg is negative)
                fbaout1 = ChLP.fba(model, exglcider; sense = glc_sense)
                exglc = ChU.av(model, fbaout1, exglcider)
                ChU.fixxing(model, exglcider, exglc; btol) do
                    # min cost
                    fbaout = ChLP.fba(model, costider; sense = cost_sense)
                end
            end
            # fix Z and readjust biom so exp and model yield (z/vg) match
            Sfac = fix_Z_and_readj_model(iJR, Ed, exp; proc_fun, verbose = false);
            
            # store
            thsv_set!(LPDAT, [FBAmethod, :Sfac, exp], Sfac)
            thsv_set!(LPDAT, [FBAmethod, :model, exp], ChU.compressed_model(model))
            thsv_set!(LPDAT, [FBAmethod, :fbaout, exp], fbaout)
        end
    end

    # FBA_Z_FIX_MINMAX_VG_MINMAX_COST
    let
        for (FBAmethod, glc_sense, cost_sense) in [
            (:FBA_Z_FIX_MIN_VG_MIN_COST, max_sense, min_sense),
            (:FBA_Z_FIX_MIN_VG_MAX_COST, max_sense, max_sense),
            (:FBA_Z_FIX_MAX_VG_MIN_COST, min_sense, min_sense),
            (:FBA_Z_FIX_MAX_VG_MAX_COST, min_sense, max_sense),
        ]

            !redo && haskey(LPDAT, [FBAmethod, :model, exp]...) && continue

            model = iJR.load_model("fva_models", exp)

            # fix Z
            exp_growth = Ed.val("D", exp)
            dflx = abs(exp_growth * btol)
            ChU.bounds!(model, biomider, exp_growth - dflx, exp_growth + dflx)
            
            # min vg (vg is negative)
            fbaout1 = ChLP.fba(model, exglcider; sense = glc_sense)
            exglc = ChU.av(model, fbaout1, exglcider)
            dflx = abs(exglc * btol)
            ChU.bounds!(model, exglcider, exglc - dflx, exglc + dflx)

            # min cost
            fbaout = ChLP.fba(model, costider; sense = cost_sense)

            thsv_set!(LPDAT, [FBAmethod, :model, exp], ChU.compressed_model(model))
            thsv_set!(LPDAT, [FBAmethod, :fbaout, exp], fbaout)
        end
    end

    # FBA_MAX_Z_MINMAX_COST
    let
        for (FBAmethod, biom_sense, cost_sense) in [
            (:FBA_MAX_Z_MIN_COST, max_sense, min_sense),
            (:FBA_MAX_Z_MAX_COST, max_sense, max_sense),
        ]
            !redo && haskey(LPDAT, [FBAmethod, :model, exp]...) && continue

            model = iJR.load_model("fva_models", exp)
            
            # max z
            fbaout1 = ChLP.fba(model, biomider; sense = biom_sense)
            objval = ChU.av(model, fbaout1, biomider)
            dflx = abs(objval * btol)
            ChU.bounds!(model, biomider, objval - dflx, objval + dflx)

            # min cost
            fbaout = ChLP.fba(model, costider; sense = cost_sense)
    
            thsv_set!(LPDAT, [FBAmethod, :model, exp], ChU.compressed_model(model))
            thsv_set!(LPDAT, [FBAmethod, :fbaout, exp], fbaout)
        end
    end

    # FBA_Z_FIX_MINMAX_COST
    let
        for (FBAmethod, cost_sense) in [
            (:FBA_Z_FIX_MIN_COST, min_sense),
            (:FBA_Z_FIX_MAX_COST, max_sense),
        ]
            !redo && haskey(LPDAT, [FBAmethod, :model, exp]...) && continue

            model = iJR.load_model("fva_models", exp)
            exp_growth = Ed.val("D", exp)
            
            # fix Z
            exp_growth = Ed.val("D", exp)
            dflx = abs(exp_growth * btol)
            ChU.bounds!(model, biomider, exp_growth - dflx, exp_growth + dflx)

            # min cost
            fbaout = ChLP.fba(model, costider; sense = cost_sense)

            thsv_set!(LPDAT, [FBAmethod, :model, exp], ChU.compressed_model(model))
            thsv_set!(LPDAT, [FBAmethod, :fbaout, exp], fbaout)
        end
    end

    # FBA_Z_VG_FIX_MINMAX_COST
    let
        for (FBAmethod, cost_sense) in [
            (:FBA_Z_VG_FIX_MIN_COST, min_sense),
            (:FBA_Z_VG_FIX_MAX_COST, max_sense),
        ]
            !redo && haskey(LPDAT, [FBAmethod, :model, exp]...) && continue

            model = iJR.load_model("fva_models", exp)

            # fix Z
            exp_growth = Ed.val("D", exp)
            dflx = abs(exp_growth * btol)
            ChU.bounds!(model, biomider, exp_growth - dflx, exp_growth + dflx)

            # fix vg
            exp_exglc = -abs(Ed.uval("GLC", exp))
            dflx = abs(exp_exglc * btol)
            ChU.bounds!(model, exglcider, exp_exglc - dflx, exp_exglc + dflx)

            # min cost
            fbaout = ChLP.fba(model, costider; sense = cost_sense)

            thsv_set!(LPDAT, [FBAmethod, :model, exp], ChU.compressed_model(model))
            thsv_set!(LPDAT, [FBAmethod, :fbaout, exp], fbaout)
        end
    end
    return LPDAT
end


## ----------------------------------------------------------------------------------
# create (iJR, Ed, exp) collection 
function collect_experiments() 
   collprod((iJR, Ed, EXPS)) = vec(collect(Iterators.product([iJR], [Ed], EXPS)))
   mod_dat = [
        (NiJR, Nd, 1:17),
        (KiJR, Kd, 5:13),
        (FiJR, Fd, 1:4),
    ]
    return vcat(map(collprod, mod_dat)...)
end

## ----------------------------------------------------------------------------------
let

    EXP_COLL = collect_experiments() 
    PROGRESS = Dict() # to now when finished
    isfinished(iJR) = count((dat) -> isequal(first(dat), iJR), EXP_COLL) == PROGRESS[iJR]

    Ch = Channel(nthreads()) do Ch_
        for (iJR, Ed, exp) in EXP_COLL
            put!(Ch_, (iJR, Ed, exp))
        end
    end
    
    LPDATs = Dict()
    # @threads 
    for _ in 1:nthreads()
        thid = threadid()
        
        for (iJR, Ed, exp) in Ch

            src = string(nameof(Ed))
            dfile = ChE.LP_DATfile(src)
            
            lock(wlk) do
                # init
                if !haskey(LPDATs, src)
                    LPDATs[src] = isfile(dfile) ? ldat(dfile) : MagicDict()
                end
                get!(PROGRESS, iJR, 0)

                # info
                @info(string("Doing ", "-"^30), 
                    exp, src, thid
                ); println()
            end

            # --------------------------------------------
            do_LP!(LPDATs[src], iJR, Ed, exp)
       
            # --------------------------------------------
            lock(wlk) do
                # progress
                PROGRESS[iJR] += 1

                # finish
                if isfinished(iJR)
                    @info(string("Finished ", "-"^30), 
                        src, thid
                    ); println()

                    # saving
                    dry_run || sdat(LPDATs[src], dfile)
                    
                    # save memory
                    LPDATs[src] = nothing
                    GC.gc()
                end
            end
        end

    end # for thid
end;