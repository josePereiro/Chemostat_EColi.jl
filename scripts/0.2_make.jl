using ProjAssistant
@quickactivate

@time begin
    import Chemostat_EColi
    const ChE = Chemostat_EColi

    import Chemostat_Kayser2005
    const ChK = Chemostat_Kayser2005

    import Chemostat_Nanchen2006
    const ChN = Chemostat_Nanchen2006

    import Chemostat_Folsom2014
    const ChF = Chemostat_Folsom2014

    import Chemostat_InSilico
    const ChIn = Chemostat_InSilico

    using Dates
end

## ------------------------------------------------------
# TODO: add ArgParse capability
JULIA_CMD = "julia"
NTHREADS = max(Sys.CPU_THREADS - 1, 1)
PROJ_ROOT = dirname(projdir(ChE))

## ------------------------------------------------------
# utils
function _info(msg; kwargs...)
    time = now()
    println("\n"^2)
    @info("$(msg) $("-"^45)", time, kwargs...)
    println("\n"^2)
end

function run_script(Proj::Module, arg, args...; flags::Vector = [], nthrs = NTHREADS)
    src = scriptsdir(Proj, arg, args...)
    !isfile(src) && error("script file is missing, file: ", src)
    
    _info("Starting"; src, nthrs)
    cmd = Cmd(string.([JULIA_CMD, "-t$(nthrs)", "--startup-file=no", "--project", src, flags...]))
    run(cmd; wait = true)
    _info("Finished"; src, nthrs)

    return nothing
end

## ------------------------------------------------------
# testing juila
try; run_script(ChE, "0.0_test.jl")
    catch err
    @error("Problem calling julia", JULIA_CMD)
    rethrow(err)
end

# check proj dir
let
    paths = readdir(PROJ_ROOT)
    for proj in [
            "MaxEnt_EColi_paper", 
            "Chemostat_EColi", 
            "Chemostat_Kayser2005", 
            "Chemostat_Nanchen2006", 
            "Chemostat_Folsom2014", 
        ]
        !(proj in paths) && error("proj ", proj, " is missing from the root dir")
    end
    _info("Root dir test passed")
end

## ------------------------------------------------------
# run simulation
let
    # sim2D
    _info("Run simulation")
    run_script(ChIn.Dynamic, ["sim2D"], "1.0_compute_Spaces.jl")
    run_script(ChIn.Dynamic, ["sim2D"], "2.0_run_sim2D.jl")
    run_script(ChIn.Dynamic, ["sim2D"], "3.0_run_MaxEnt.jl")
end

## ------------------------------------------------------
# run EColi packages
let
    _info("Run EColi packages")

    ## ------------------------------------------------------
    # Kayser2005
    run_script(ChK.KayserData, "1_converted_data.jl")
    run_script(ChK.iJR904, "1.0_prepare_gem.jl")
    run_script(ChK.iJR904, "2.0_maxent_ep_variants.jl")
    
    ## ------------------------------------------------------
    # Nanchen2005
    run_script(ChN.iJR904, "1.0_prepare_gem.jl")
    run_script(ChN.iJR904, "2.0_maxent_ep_variants.jl")
    
    ## ------------------------------------------------------
    # Folsom2014
    run_script(ChF.iJR904, "1.0_prepare_gem.jl")
    run_script(ChF.iJR904, "2.0_maxent_ep_variants.jl")

end

## ------------------------------------------------------
# collect
let
    run_script(ChE, ["experiments"], "0.1_LP.jl")
    run_script(ChE, ["experiments"], "0.2_collect_DAT.jl")
end

## ------------------------------------------------------
# make plots
let
    run_script(ChE, ["simulation"], "2.0_sim_plots.jl")
    run_script(ChE, ["experiments"], "1.0_exp_plots.jl")
end

## ------------------------------------------------------
# run latex
let
    tex_proj = joinpath(PROJ_ROOT, "MaxEnt_EColi_paper")
    teximgs_dir = joinpath(tex_proj, "images")
    !isdir(teximgs_dir) && error("image dir not found, path: ", imgs_dir)

    # cp to image dir
    _info("Copying figures"; teximgs_dir)
    for figi in 1:10
        imgname = string("figure_", figi, ".png")
        srcfile = plotsdir(ChE, imgname)
        !isfile(srcfile) && error("image file is missing, file: ", srcfile)
        destfile = joinpath(teximgs_dir, imgname)
        cp(srcfile, destfile; force = true)
        @info("Copied figure file", srcfile, destfile)
    end

    # try to run pdflatex
    _info("Trying pdflatex"; tex_proj)
    tex_mainfile = joinpath(tex_proj, "MaxEnt_EColi.tex")
    !isfile(tex_mainfile) && error("main tex file not found, file: ", tex_mainfile)
    pdflatex_cmd = Cmd(
        ["pdflatex", "-synctex=1", "-interaction=nonstopmode", "-file-line-error", "-recorder", 
        string("-output-directory=", tex_proj), tex_mainfile]
    )
    bibtex_cmd = Cmd(`bibtex  "MaxEnt_EColi"`)
    
    try; 
        cd(tex_proj)
        
        # from vscode recipe
        run(pdflatex_cmd; wait = true)
        run(bibtex_cmd; wait = true)
        run(pdflatex_cmd; wait = true)
        run(bibtex_cmd; wait = true)
        run(pdflatex_cmd; wait = true)

        catch err
            @info("pdflatex/bibtex failed, but the tex project is updated. You must compile it manually. ", tex_proj)
            rethrow(err)
        finally; cd(projdir(ChE))
    end

end

_info("The end")

