let
    method = ME_MAX_POL
    niders = length(FLX_IDERS)

    stds = map(EXPS) do exp
        DAT[method, :eperr, :flx, FLX_IDERS, exp] ./ abs.(DAT[method, :Fd, :flx, "GLC", exp])
    end
    
    # sGLC
    ids = [:D, :X, :uAC, :uGLC, :uLAC, :uPYR]
    for id in ids
        p = scatter(;xlabel = string(id), ylabel = "std")
        id_vals = [Fd.val(id, exp) for exp in EXPS]
        sidx = sortperm(id_vals)
        for idx in sidx
            x = id_vals[idx]
            ys = stds[idx]
            scatter!(p, fill(x, niders), ys; label = "", color = :black, m = 8)
            scatter!(p, [x], [mean(ys)]; label = "", color = :red, m = (8, :star))
        end
        mysavefig(p, "av_var_vs_"; id)
    end
end    