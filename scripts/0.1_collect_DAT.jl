import DrWatson
const DW = DrWatson
DW.quickactivate(@__DIR__, "Chemostat_EColi")

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

## ----------------------------------------------------------------------------------
function collect_commom_DAT!(DAT, Ed, iJR)

    # maps
    # ----------------------------------------------------------------------------------
    rxns_map = iJR.load_rxns_map()
    
    # ----------------------------------------------------------------------------------
    max_model = iJR.load_model("max_model"; uncompress = false)
    objider = iJR.BIOMASS_IDER
    
    # ----------------------------------------------------------------------------------
    for exp in DAT[:EXPS]
        # exp dat
        Ed_biom = Ed.val("D", exp)
        DAT[:exp, :flx, "D", exp] = Ed_biom
        DAT[:exp, :err, "D", exp] = 0.0

        # bounds
        fva_model = iJR.load_model("fva_models", exp; uncompress = false)
        max_lb, max_ub = ChU.bounds(max_model, objider)
        fva_lb, fva_ub = ChU.bounds(fva_model, objider)
        lb = max(max_lb, fva_lb)
        ub = min(max_ub, fva_ub)
        DAT[:bounds, "D", exp] = (lb, ub)

        for Ed_ider in DAT[:FLX_IDERS]
            model_exch = rxns_map[Ed_ider]

            # exp dat
            Ed_flx = Ed.uval(Ed_ider, exp)
            Ed_err = try; Ed.uerr(Ed_ider, exp) catch; 0.0 end
            DAT[:exp, :flx, Ed_ider, exp] = Ed_flx
            DAT[:exp, :err, Ed_ider, exp] = Ed_err

            # bounds
            max_lb, max_ub = ChU.bounds(max_model, model_exch)
            fva_lb, fva_ub = ChU.bounds(fva_model, model_exch)
            lb = max(max_lb, fva_lb)
            ub = min(max_ub, fva_ub)
            DAT[:bounds, Ed_ider, exp] = (lb, ub)

        end

        # # all data
        # for Ed_ider in Ed.DAT_IDERS

        #     val = Ed.val(Ed_ider, exp)
        #     err = Ed.err(Ed_ider, exp)

        #     DAT[:exp, :flx, Ed_ider, exp] = val
        #     DAT[:exp, :err, Ed_ider, exp] = err
        # end
    end

    return DAT
end

## ----------------------------------------------------------------------------------
function collect_ME_DAT!(DAT, iJR)

    costider = iJR.COST_IDER
    objider = iJR.BIOMASS_IDER

    rxns_map = iJR.load_rxns_map()

    # util stuff
    WLOCK = ReentrantLock()
    isexpdep(method) = !(method in DAT[:EXP_INDEPENDENT])
    depks(method, typ, Ed_met, exp) = 
        isexpdep(method) ? (method, typ, Ed_met, exp) : (method, typ, Ed_met)
    DAT_FILE_PREFFIX =  "maxent_ep_dat"
    function dat_file(;method, exp)
        kwargs = isexpdep(method) ? (;method, exp) : (;method)
        fname = UJL.mysavename(DAT_FILE_PREFFIX, "jls"; kwargs...)
        iJR.procdir(fname)
    end

    # Feed jobs
    nths = nthreads()
    Ch = Channel(nths) do ch
        for method in DAT[:ME_METHODS]
            for exp in DAT[:EXPS]
                put!(ch, (exp, method))
                !isexpdep(method) && break # Do only once
            end
        end
    end

    @threads for thid in 1:nths
        thid = threadid()
        thid == 1 && nths > 1 && continue
        for (exp, method) in Ch
            
            # ME data
            datfile = dat_file(;method, exp)
            !isfile(datfile) && continue
            dat = deserialize(datfile)
            
            model = dat[:model]
            epouts = dat[:epouts]
            exp_beta = maximum(keys(epouts)) # dat[:exp_beta]
            epout = epouts[exp_beta]
            
            lock(WLOCK) do
                @info("Doing", 
                    exp, method, 
                    length(dat[:epouts]), 
                    epout.iter, thid
                ); println()
            end

            # Biomass and cost
            ep_biom = ChU.av(model, epout, objider)
            ep_biomstd = sqrt(ChU.va(model, epout, objider))
            ep_cost = ChU.av(model, epout, costider)
            ep_coststd = sqrt(ChU.va(model, epout, costider))
            
            # store
            lock(WLOCK) do
                DAT[depks(method, :flx, "D", exp)...] = ep_biom
                DAT[depks(method, :err, "D", exp)...] = ep_biomstd
                DAT[depks(method, :flx, "cost", exp)...] = ep_cost
                DAT[depks(method, :err, "cost", exp)...] = ep_coststd
            end

            for Ed_met in DAT[:FLX_IDERS]
                model_exch = rxns_map[Ed_met]

                # flxs
                ep_av = ChU.av(model, epout, model_exch)
                ep_std = sqrt(ChU.va(model, epout, model_exch))

                # proj 2d
                proj = ChLP.projection2D(model, objider, model_exch; l = 50)
                        
                lock(WLOCK) do
                    DAT[depks(method, :proj, Ed_met, exp)...] = proj
                    DAT[depks(method, :flx, Ed_met, exp)...] = ep_av
                    DAT[depks(method, :err, Ed_met, exp)...] = ep_std
                end
            end

            # inner flxs
            for (exider, model_iders) in iJR.load_inners_idermap()
                # flxs
                ep_av = ChU.av(model, epout, model_iders[1])
                ep_std = sqrt(ChU.va(model, epout, model_iders[1]))
                if length(model_iders) == 2 # reversible
                    # r = r+ - r-
                    ep_av -= ChU.av(model, epout, model_iders[2])
                    # ep_std += sqrt(ChU.va(model, epout, model_iders[2]))
                    ep_std = NaN
                end

                # proj 2d (fwd only)
                proj = ChLP.projection2D(model, objider, model_iders[1]; l = 50)
                
                lock(WLOCK) do
                    DAT[depks(method, :proj, exider, exp)...] = proj
                    DAT[depks(method, :flx, exider, exp)...] = ep_av
                    DAT[depks(method, :err, exider, exp)...] = ep_std
                end
            end
        
        end # for (exp, method)
    end # for thid

    return DAT
end

## ----------------------------------------------------------------------------------
function collect_lp_DAT!(DAT, iJR)

    rxns_map = iJR.load_rxns_map()

    objider = iJR.BIOMASS_IDER
    costider = iJR.COST_IDER

    LP_DAT_FILE = iJR.procdir("lp_dat_file.bson")
    LP_DAT = ChU.load_data(LP_DAT_FILE; verbose = false);

    for method in DAT[:LP_METHODS]
            
        for exp in DAT[:EXPS]

            model = LP_DAT[method, :model, exp]
            fbaout = LP_DAT[method, :fbaout, exp]

            # Biomass
            fba_biom = ChU.av(model, fbaout, objider)
            DAT[method, :flx, "D", exp] = fba_biom
            DAT[method, :err, "D", exp] = 0.0

            fba_cost = ChU.av(model, fbaout, costider)
            DAT[method, :flx, "cost", exp] = fba_cost
            DAT[method, :err, "cost", exp] = 0.0

            for Ed_ider in DAT[:FLX_IDERS]
                model_ider = rxns_map[Ed_ider]

                fba_flx = ChU.av(model, fbaout, model_ider)
                DAT[method, :flx, Ed_ider, exp] = fba_flx
            end

            idermap = merge(iJR.load_inners_idermap(), iJR.load_inners_idermap())
            for (exglcider, model_iders) in idermap
                # flxs
                fba_flx = ChU.av(model, fbaout, model_iders[1])
                if length(model_iders) == 2 # reversible
                    # r = r+ - r-
                    fba_flx -= ChU.av(model, fbaout, model_iders[2])
                end
                        
                DAT[method, :flx, exglcider, exp] = fba_flx
                DAT[method, :err, exglcider, exp] = 0.0
            end
        end

    end # for method

    return DAT
end

## ----------------------------------------------------------------------------------
let
    for (iJR, Ed, EXPS, FLX_IDERS) in [
        (NiJR, Nd, Nd.EXPS, ["GLC", "AC"]),
        (KiJR, Kd, 5:13, ["GLC", "CO2", "O2", "AC", "NH4"]),
        (FiJR, Fd, 1:4, ["GLC", "SUCC", "AC", "FORM"]),
    ]
        src = string(nameof(Ed))
        
        # ---------------------------------------------------------------------------------
        DAT_FILE = ChE.DATfile(src)
        DAT = (!isfile(DAT_FILE) || new_dat) ? ChU.DictTree() : ChE.load_DAT(src)

        DAT[:EXPS] = EXPS
        DAT[:FLX_IDERS] = FLX_IDERS

        # ME methods
        ME_Z_OPEN_G_OPEN           = :ME_Z_OPEN_G_OPEN
        ME_MAX_POL                 = :ME_MAX_POL
        ME_MAX_POL_B0              = :ME_MAX_POL_B0
        ME_Z_EXPECTED_G_MOVING     = :ME_Z_EXPECTED_G_MOVING
        ME_Z_EXPECTED_G_BOUNDED    = :ME_Z_EXPECTED_G_BOUNDED
        ME_Z_FIXXED_G_BOUNDED      = :ME_Z_FIXXED_G_BOUNDED

        # LP methods
        FBA_Z_FIX_MAX_VG_MIN_COST  = :FBA_Z_FIX_MAX_VG_MIN_COST
        FBA_Z_FIX_MIN_VG_MIN_COST  = :FBA_Z_FIX_MIN_VG_MIN_COST
        FBA_Z_FIX_MAX_VG_MAX_COST  = :FBA_Z_FIX_MAX_VG_MAX_COST
        FBA_Z_FIX_MIN_VG_MAX_COST  = :FBA_Z_FIX_MIN_VG_MAX_COST
        FBA_Z_VG_FIX_MAX_COST      = :FBA_Z_VG_FIX_MAX_COST
        FBA_Z_VG_FIX_MIN_COST      = :FBA_Z_VG_FIX_MIN_COST
        FBA_Z_FIX_MAX_COST         = :FBA_Z_FIX_MAX_COST
        FBA_Z_FIX_MIN_COST         = :FBA_Z_FIX_MIN_COST
        FBA_MAX_Z_MIN_COST         = :FBA_MAX_Z_MIN_COST
        FBA_MAX_Z_MAX_COST         = :FBA_MAX_Z_MAX_COST

        DAT[:LP_METHODS] = [
            FBA_Z_FIX_MAX_VG_MIN_COST, FBA_Z_FIX_MAX_VG_MAX_COST, 
            FBA_Z_FIX_MIN_VG_MIN_COST, FBA_Z_FIX_MIN_VG_MAX_COST,
            FBA_Z_VG_FIX_MAX_COST, FBA_Z_VG_FIX_MIN_COST, FBA_Z_FIX_MAX_COST, 
            FBA_Z_FIX_MIN_COST, FBA_MAX_Z_MIN_COST, FBA_MAX_Z_MAX_COST
        ]

        DAT[:ME_METHODS] = [
            # ME_Z_OPEN_G_OPEN, 
            ME_MAX_POL,
            ME_MAX_POL_B0,
            # ME_Z_FIXXED_G_BOUNDED, 
            ME_Z_EXPECTED_G_BOUNDED, 
            # ME_Z_EXPECTED_G_MOVING
        ]

        DAT[:EXP_INDEPENDENT] = [ME_MAX_POL_B0]
        
        # ---------------------------------------------------------------------------------
        collect_commom_DAT!(DAT, Ed, iJR)
        skip_me || collect_ME_DAT!(DAT, iJR)
        skip_lp || collect_lp_DAT!(DAT, iJR)
        
        # ---------------------------------------------------------------------------------
        # saving
        dry_run || serialize(DAT_FILE, DAT)
    end

end;