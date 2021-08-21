module Chemostat_EColi

    try
        import Chemostat_Kayser2005
        import Chemostat_Nanchen2006
        import Chemostat_Folsom2014
    catch err
        @warn("Error loading source package, see scripts/0.0_install.jl")
        rethrow(err)
    end
    const ChK = Chemostat_Kayser2005
    const ChN = Chemostat_Nanchen2006
    const ChF = Chemostat_Folsom2014

    import Plots
    using ProjAssistant
    @gen_top_proj

    using Serialization
    
    function __init__()
        @create_proj_dirs
    end

    # ---------------------------------------------------------------------------------
    # Dat
    DATfile(src) = procdir(Chemostat_EColi, "dat", (;src), ".jls")
    load_DAT(src) = ldat(DATfile(src))
    
    LP_DATfile(src) = procdir(Chemostat_EColi, "lp_dat", (;src), ".jls")
    load_LP_DAT(src) = ldat(LP_DATfile(src))

    ## ------------------------------------------------------
    # ustils
    function lims(d; marginf = 0.1)
        m = minimum(d)
        M = maximum(d)
        margin = (M - m) * marginf
        (m - margin, M + margin)
    end
    
    function lims(d, oldl; marginf = 0.1)
        isempty(oldl) && lims(d; marginf)
        newl = lims(d; marginf)
        [min(first(newl), first(oldl)), max(last(newl), last(oldl))]
    end

end
