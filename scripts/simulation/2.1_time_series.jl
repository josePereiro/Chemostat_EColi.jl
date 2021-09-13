## ------------------------------------------------------
function plot_time_series(simid, D, ϵs, cg)

   Dstr = string(round(D; sigdigits = 3))
   tslim = [Inf, -Inf]

   Xp = plot()
   Xlim = [0.0, -Inf]
   
   sgp = plot()
   sglim = [0.0, -Inf]

   max_ϵs_ = maximum(ϵs)
   ls0 = 2.5

   for (ϵi, ϵ) in enumerate(ϵs)
      
      status = get_status(simid, (;D, ϵ, cg))
      (status == UNDONE_SIM_STATUS) && continue

      @info("Doing", D, ϵ, cg, status)

      (ts, Xts) = collect_ts(:Xts, simid, (;D, ϵ, cg))
      (ts, sgts) = collect_ts(:sgts, simid, (;D, ϵ, cg))

      ts = ts ./ maximum(ts) .* 100.0
      tslim = (-5.0, 105.0)

      for (p, panel, ys, ylim, ylabel) in [
         [Xp, "A", Xts, Xlim, _textbf("X")],
         [sgp, "B", sgts, sglim, _textbf("s_g")]
      ]
         ylim .= lims(ys, ylim)
         plot!(p, ts, ys; 
            title = _textbf(panel),
            xlim = tslim, ylim,
            label = "", 
            ylabel,
            color = :black,
            xlabel = _textbf("simulation progress (\\%)"),
            alpha = 0.7,
            lw = ls0 + (2.0 * ls0) * (ϵ / max_ϵs_)
         )
      end
   end
   
   # annotation
   # fontsize = 16
   # text = _textbf("D=", Dstr)
   # annotate!(Xp, [(70.0, last(Xlim)* 0.4, (text, fontsize, :left, :top, :black))])
   # annotate!(sgp, [(70.0, last(sglim)* 0.4, (text, fontsize, :left, :top, :black))])
   
   # ffile = sfig(ChE, [Xp, sgp], 
   #    @fileid, "dyn_time_series", (;D, cg), ".png";
   #    layout = (1, 2)
   # )
   # @info("Done", D, ffile)

   return [Xp, sgp]
end