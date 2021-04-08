module Chemostat_EColi

    try
        import Chemostat_Kayser2005
        import Chemostat_Heerden2013
        import Chemostat_Folsom2014
    catch err
        @warn("Error loading source package, see scripts/0.0_install.jl")
        rethrow(err)
    end
    const ChK = Chemostat_Kayser2005
    const ChH = Chemostat_Heerden2013
    const ChF = Chemostat_Folsom2014

    import UtilsJL
    const UJL = UtilsJL
    UJL.gen_top_proj(@__MODULE__)

    function __init__()
        UJL.create_proj_dirs(@__MODULE__)
    end

end
