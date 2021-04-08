import DrWatson
const DW = DrWatson
DW.quickactivate(@__DIR__, "Chemostat_EColi")
import Pkg

## ------------------------------------------------------------------
# Install per source packages 
let
    # Here the source packages where clonned at the same
    # dir than Chemostat_EColi
    proj_dir = DW.projectdir()
    source_pkgs = [
        joinpath(proj_dir, "../Chemostat_Kayser2005"), 
        joinpath(proj_dir, "../Chemostat_Heerden2013"), 
        joinpath(proj_dir, "../Chemostat_Folsom2014")
    ]
    for pkg in source_pkgs
        !isdir(pkg) && 
            error("pkg: ", basename(pkg), " not found in the given path: ", pkg)
        Pkg.develop(path=pkg)
    end
end