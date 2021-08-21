using ProjAssistant
@quickactivate

import Chemostat_EColi
const ChE = Chemostat_EColi
using Dates

## ------------------------------------------------------
# TODO: add ArgParse capability
JULIA_CMD = "julia"
NTHREADS = max(Sys.CPU_THREADS - 1, 1)
PROJ_ROOT = dirname(projdir(ChE))

## ------------------------------------------------------
# utils
function run_script(arg, args...; flags::Vector = [], nthrs = NTHREADS)
    src = scriptsdir(ChE, arg, args...)
    !isfile(src) && error("script file is missing, file: ", src)
    
    println("\n"^3)
    time = now()
    @info("Running $("-"^23)", src, nthrs, time)
    println("\n"^2)

    cmd = Cmd(string.([JULIA_CMD, "-t$(nthrs)", "--startup-file=no", "--project", src, flags...]))
    run(cmd; 
        wait = true
    )

    println("\n"^3)
    time = now()
    @info("Finished $("-"^23)", src, nthrs, time)
    println("\n"^2)

    return nothing
end

## ------------------------------------------------------
# test julia
try; run_script("test.jl")
    catch err
    @error("Problem calling julia", JULIA_CMD)
    rethrow(err)
end

## ------------------------------------------------------
# run simulation

## ------------------------------------------------------
# run EColi packages

## ------------------------------------------------------
# collect
for script in ["0.1_LP.jl", "0.2_collect_DAT.jl"]
    run_script(["experiments"], script)
end

## ------------------------------------------------------
# make plots
run_script(["simulation"], "2.0_sim_plots.jl")
run_script(["experiments"], "1.0_exp_plots.jl")

## ------------------------------------------------------
# update paper