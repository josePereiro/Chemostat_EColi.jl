## ------------------------------------------------------
function plot_time_series(simid, D, ϵs, cg)

   ϵcolors = palette(:thermal, length(ϵs))

   Dstr = string(round(D; sigdigits = 3 ))
   tslim = [Inf, -Inf]

   Xp = plot()
   Xlim = [0.0, -Inf]
   
   sgp = plot()
   sglim = [0.0, -Inf]

   for (ϵi, ϵ) in enumerate(ϵs)
      thid = threadid()
      
      status = get_status(simid, (;D, ϵ, cg))
      (status == UNDONE_SIM_STATUS) && continue

      @info("Doing", D, ϵ, cg, status)

      (ts, Xts) = collect_ts(:Xts, simid, (;D, ϵ, cg))
      (ts, sgts) = collect_ts(:sgts, simid, (;D, ϵ, cg))

      ts = ts ./ maximum(ts) .* 100.0
      tslim = (-5.0, 105.0)
      
      for (p, ys, ylim, ylabel) in [
         [Xp, Xts, Xlim, "X"],
         [sgp, sgts, sglim, "sg"]
      ]
         ylim .= lims(ys, ylim)
         plot!(p, ts, ys; 
            title = string("D=", Dstr),
            xlim = tslim, ylim,
            label = "", 
            ylabel,
            color = ϵcolors[ϵi], 
            zcolor = ϵ,
            colorbar_title = "ϵ",
            xlabel = "simulation progress (%)",
            lw = 3
         )
      end

   end
   
   ffile = sfig(ChE, [Xp, sgp], 
      @fileid, "dyn_time_series", (;D, cg), ".png";
      layout = (1, 2)
   )
   @info("Done", D, ffile)

   return [Xp, sgp]
end