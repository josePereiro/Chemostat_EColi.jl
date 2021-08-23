# Chemostat_EColi

[![Build Status](https://github.com/josePereiro/Chemostat_EColi.jl/workflows/CI/badge.svg)](https://github.com/josePereiro/Chemostat_EColi.jl/actions)

This is the main package for reproducing the results published at ?.

## Installation (macOs or linux)

The installation expect you to have `git` and `julia` (version >= 1.6) installed and on the `PATH` with those names.
For testing that run this:

```bash
git version && echo && julia -e 'println(VERSION)' && echo "Test Passed"
```

For installing all required repositories (this part require internet on the terminal) run the follow command in an empty folder (we'll called `PROJ_ROOT`).

```bash
git clone --depth 1 --branch v0.1.0 --single-branch \
https://github.com/josePereiro/Chemostat_EColi.jl Chemostat_EColi && \
cd Chemostat_EColi && \
julia --startup-file=no scripts/0.1_install.jl
```

To reproduce the results (`WARNING`: this ~ 2 days in a 4 cpu 2019 MacBook Air), run the follow in the `Chemostat_EColi` folder:

```bash
julia --startup-file=no scripts/0.2_make.jl
```

If `bibtex` and `pdflatex` are install and in the `PATH` the make script should reproduce all the figures (at `Chemostat_EColi/plots`) and produce the `.pdf` at `MaxEnt_EColi_paper/MaxEnt_EColi.pdf`.
