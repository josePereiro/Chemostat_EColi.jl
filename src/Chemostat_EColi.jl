module Chemostat_EColi

    import Plots
    using ProjAssistant
    @gen_top_proj

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
