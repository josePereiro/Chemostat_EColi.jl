# Install
let
    dir = pwd()
    @info("Testing", dir)
    
    cmdstr = """\
    git clone --depth 1 --branch main --single-branch \
    https://github.com/josePereiro/Chemostat_EColi.jl Chemostat_EColi && \
    cd Chemostat_EColi && \
    julia --startup-file=no scripts/0.1_install.jl\
    """
    
    run(Cmd(`bash -c $(cmdstr)`); wait = true)
end

# Import
using Chemostat_EColi
using Test

# Tests
@testset "Chemostat_EColi.jl" begin
    @test true
end
