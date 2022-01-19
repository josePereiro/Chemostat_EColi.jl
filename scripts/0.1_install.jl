import Pkg
import Pkg: RegistrySpec, PackageSpec, pkgdir, PRESERVE_ALL

## ------------------------------------------------------------------
function _info(msg; kwargs...)
    println()
    msg = rpad(string(msg, " "), 60, "-")
    @info(msg, kwargs...)
    println()
    println()
end

## ------------------------------------------------------------------
# desable autoprecompilation
ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0

## ------------------------------------------------------------------
_info("Add registries")
Pkg.Registry.add(RegistrySpec(url = "https://github.com/josePereiro/CSC_Registry.jl"))

# ----------------------------------------------------------------------------
_info("Activate Chemostat_EColi")
CH_ECOLI_DIR = Base.current_project(@__DIR__) |> dirname
basename(CH_ECOLI_DIR) != "Chemostat_EColi" && error("Not at 'Chemostat_EColi' folder!!!")
Pkg.activate(CH_ECOLI_DIR)

# ----------------------------------------------------------------------------
# ARGS
using ArgParse

set = ArgParseSettings()
@add_arg_table! set begin
    "--git-cmd"
        help = "the git command to use"   
        default = "git"
end

if isinteractive()
    # Dev values
    GIT_CMD = "git"
else
    parsed_args = parse_args(set)
    GIT_CMD = parsed_args["git-cmd"]
end

## ------------------------------------------------------------------
# globals
PROJ_ROOT = dirname(CH_ECOLI_DIR)
cd(PROJ_ROOT)

_info("Globals"; PROJ_ROOT, GIT_CMD)

## ------------------------------------------------------------------
# git utils
function _checkout_tag(tag)
    cmd_strs = [GIT_CMD, "checkout", tag]
    run(Cmd(cmd_strs), wait=true)
end
function _clone_tag(outname, url, tag)
    cmd_strs = [ GIT_CMD, "clone", "--depth", "1", "--branch", tag, "--single-branch", url, outname]
    run(Cmd(cmd_strs), wait=true)
end

## ------------------------------------------------------------------
# Clone and dev dependencies
# TODO: (Chemostat_InSilico) pkg> add BenchmarkTools
let
    _info("Adding .jl devs")
    for (pkgname, url, tag) in [
            (
                "Chemostat_InSilico", 
                "https://github.com/josePereiro/Chemostat_InSilico.jl", 
                "v0.2.3"
            ),
            (
                "Chemostat_Kayser2005", 
                "https://github.com/josePereiro/Chemostat_Kayser2005.jl",
                "v0.2.2"
            ), 
            (
                "Chemostat_Nanchen2006", 
                "https://github.com/josePereiro/Chemostat_Nanchen2006.jl", 
                "v0.2.2"
            ), 
            (
                "Chemostat_Folsom2014", 
                "https://github.com/josePereiro/Chemostat_Folsom2014.jl",
                "v0.2.2"
            ), 
        ]
        
        pkg_path = joinpath(PROJ_ROOT, pkgname)

        _info("clone and checkout"; pkgname, tag)
        if !isdir(pkg_path)
            cd(PROJ_ROOT)
            _clone_tag(pkgname, url, tag)
        end

        cd(pkg_path)
        _checkout_tag(tag)
        
        # add
        _info("develop"; pkgname, tag)
        Pkg.activate(CH_ECOLI_DIR)
        Pkg.develop(;url=relpath(pkg_path), preserve = PRESERVE_ALL)
        
        _info("instanciate"; pkgname, tag)
        Pkg.activate(pkg_path)
        Pkg.instantiate()

        cd(PROJ_ROOT)
    end
end

## ------------------------------------------------------------------
# precompile
_info("precompile Chemostat_EColi")
Pkg.activate(CH_ECOLI_DIR)
Pkg.precompile()

## ------------------------------------------------------------------
# Test import
_info("Testing import")
try
    import Chemostat_Kayser2005
    import Chemostat_Nanchen2006
    import Chemostat_Folsom2014
    import Chemostat_InSilico
catch err
    @warn("Error loading source package, see scripts/0.1_install.jl")
    rethrow(err)
end

## ------------------------------------------------------------------
_info("instantiate Chemostat_EColi")
Pkg.instantiate()

# ## ------------------------------------------------------------------
# # TODO: make repo public
# # clone paper repo
# _info("cloning tex project")
# let
#     pkgname = "MaxEnt_EColi_paper"
#     pkg_path = joipath(PROJ_ROOT, pkgname)
#     tag = ""
#     url="https://github.com/josePereiro/MaxEnt_EColi_paper"
# 
#     if !isdir(pkg_path)
#         cd(PROJ_ROOT)
#         _clone_tag(pkgname, url, tag)
#     end
#     cd(pkg_path)
#     _checkout_tag(tag)
# end

_info("All good")
