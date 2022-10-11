# Chemostat_EColi

<!-- [![Build Status](https://github.com/josePereiro/Chemostat_EColi.jl/workflows/CI/badge.svg)](https://github.com/josePereiro/Chemostat_EColi.jl/actions) -->

[![DOI](https://zenodo.org/badge/355884387.svg)](https://zenodo.org/badge/latestdoi/355884387)

This is the main package for reproducing the results published at (TODO).

## Installation (macOs or linux)

The installation expect you to have `git` and `julia` (version >= 1.6) installed and on the `PATH` with those names.
For testing run:

```bash
git version && echo && julia -e 'println(VERSION)' && echo "Test Passed"
```

For installing all required repositories (this part require internet on the terminal) run the follow command in an empty folder.

```bash
git clone --depth 1 --branch main --single-branch \
https://github.com/josePereiro/Chemostat_EColi.jl Chemostat_EColi && \
cd Chemostat_EColi && \
julia --startup-file=no scripts/0.1_install.jl
```

To reproduce the results (`WARNING`: this takes ~ 2 days in a core i5 cpu 2019 MacBook Air), run the follow command (in the `Chemostat_EColi` folder):

```bash
julia --startup-file=no scripts/0.2_make.jl
```

If `bibtex` and `pdflatex` are installed and in the `PATH` the make script should produce all the figures (contained at `Chemostat_EColi/plots`) and the paper's `.pdf` (at `MaxEnt_EColi_paper/MaxEnt_EColi.pdf`).
