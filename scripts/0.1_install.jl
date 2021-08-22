using ProjAssistant
@quickactivate

@time begin
    import Chemostat_EColi
    const ChE = Chemostat_EColi

    import Pkg
    import Pkg: RegistrySpec, PackageSpec
end

## ------------------------------------------------------------------
# globals
PROJ_ROOT = dirname(projdir(ChE))
cd(PROJ_ROOT)
@info("Globals", PROJ_ROOT, pwd())

## ------------------------------------------------------------------
# utils
function _checkout_tag(tag)
    cmd_strs = ["git", "checkout", tag]
    run(Cmd(cmd_strs), wait=true)
end
function _clone_tag(outname, url, tag)
    cmd_strs = ["git", "clone", "--depth", "1", "--branch", tag, "--single-branch", url, outname]
    run(Cmd(cmd_strs), wait=true)
end

## ------------------------------------------------------------------
# Add registers
Pkg.Registry.add(RegistrySpec(url = "https://github.com/josePereiro/CSC_Registry.jl"))

## ------------------------------------------------------------------
# instanciate ChE
Pkg.instantiate()

## ------------------------------------------------------------------
# Clone and dev dependencies
let
    # desable autoprecompilation
    ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0

    for (pkgname, url, tag) in [
            (
                "Chemostat_InSilico", 
                "https://github.com/josePereiro/Chemostat_InSilico.jl", 
                "v0.2.0"
            ),
            (
                "Chemostat_Kayser2005", 
                "https://github.com/josePereiro/Chemostat_Kayser2005.jl",
                "v0.2.0"
            ), 
            (
                "Chemostat_Nanchen2006", 
                "https://github.com/josePereiro/Chemostat_Nanchen2006.jl", 
                "v0.2.0"
            ), 
            (
                "Chemostat_Folsom2014", 
                "https://github.com/josePereiro/Chemostat_Folsom2014.jl",
                "v0.2.0"
            ), 
        ]

        pkg_path = joinpath(PROJ_ROOT, pkgname)

        # clone and checkout
        if !isdir(pkg_path)
            cd(PROJ_ROOT)
            _clone_tag(pkgname, url, tag)
        end

        cd(pkg_path)
        _checkout_tag(tag)
        
        # develop
        cd(PROJ_ROOT)
        Pkg.develop(;path=pkg_path)

        cd(PROJ_ROOT)
    end
end

## ------------------------------------------------------------------
# precompile
Pkg.precompile()

## ------------------------------------------------------------------
# Test import
try
    import Chemostat_Kayser2005
    import Chemostat_Nanchen2006
    import Chemostat_Folsom2014
    import Chemostat_InSilico
catch err
    @warn("Error loading source package, see scripts/0.1_install.jl")
    rethrow(err)
end

# ## ------------------------------------------------------------------
# # clone paper repo
# # TODO: make repo public
# let
#     pkgname = "MaxEnt_EColi_paper"
#     pkg_path = joipath(PROJ_ROOT, pkgname)
#     tag = ""
#     url="https://github.com/josePereiro/MaxEnt_EColi_paper"

#     if !isdir(pkg_path)
#         cd(PROJ_ROOT)
#         _clone_tag(pkgname, url, tag)
#     end
#     cd(pkg_path)
#     _checkout_tag(tag)
# end