using ProjAssistant
@quickactivate

# ------------------------------------------------------
@time begin 
   import Chemostat_EColi
   import Chemostat_EColi: lims
   const ChE = Chemostat_EColi
    
   import Chemostat_InSilico
   const Dyn = Chemostat_InSilico.Dynamic
   import Chemostat_InSilico.Dynamic: 
      SimD2, run_simD2!, check_stst, hist, 
      n_ave_conv!, tail_ave,
      Container, vec!, Vi, 
      normalizeP!,
      save_simdat, load_simdat, simdat_file,
      set_status, get_status,
      load_batch, save_batch,
      UNDONE_SIM_STATUS, DEAD_SIM_STATUS,
      EXPLODED_SIM_STATUS, NITERS_SIM_STATUS, 
      STST_SIM_STATUS, collect_ts

   using Plots
   import GR
   !isinteractive() && GR.inline("png")

   using Base.Threads
   using ProgressMeter
   using BenchmarkTools
   using ColorSchemes

end

## ------------------------------------------------------
function pltparams()
   return (;
      titlefont = 16,
      axisfont = 16,
      guidefont = 16,
      colorbar_titlefont = 16,
      xtickfont = 12,
      ytickfont = 12,
      legendfont = 12,
      thickness_scaling = 1.6,
      size = (1220, 940)
   )
end

function _textbf(strs...) 
   strs = replace.(string.(strs), " " => "~")
   string("\$\\textbf{", strs..., "}\$")
end


## ------------------------------------------------------
# Figure 1: chemostat_scheme
let
   
   spaces = lfig(rawdir(ChE, "chemostat_scheme.png"))
   sfig(ChE, spaces, 
      "figure_1", ".png"
   )

end

## ------------------------------------------------------
# Figure 2: 2D_spaces_at_stst
let
   
   spaces = lfig(rawdir(ChE, "2D_spaces_at_stst.png"))
   sfig(ChE, spaces, 
      "figure_2", ".png"
   )

end

## ------------------------------------------------------
# Figure 3: ME_1D_example
let

   example = lfig(rawdir(ChE, "ME_1D_example.png"))
   sfig(ChE, example, 
      "figure_3", ".png"
   )

end

## ------------------------------------------------------
# Figure 4: toy_net_and_time_series
include("2.1_time_series.jl")

let
   params = lglob(Dyn, :dyn, :params, :finite_cg)
   @extract params: Ds ϵs cg simid
   
   # scheme
   toy_net = lfig(rawdir(ChE, "toy_model_scheme.png"))
   
   # time series
   Di = 9
   Xp, sgp = plot_time_series(simid, Ds[Di], ϵs, cg)
   
   # plot
   for p in [Xp, sgp]
      plot!(p; pltparams()...)
   end

   sfig(ChE, [toy_net, Xp, sgp], 
      "figure_4", ".png"; 
      layout = (1, 3)
   )
end

## ------------------------------------------------------
# Figure 5
include("2.2_sim2D_dead_explosion_map.jl")
include("2.3_Pol_Volume_Study.jl")

let
   params = lglob(Dyn, :dyn, :params, :finite_cg)
   @extract params: Ds ϵs cg simid
   
   panel1 = plot_pol_vol_study(simid, Ds, ϵs[1:5:end], cg)
   panel2 = plot_dead_explosion_map()
   
   for p in [panel1, panel2]
      plot!(p; pltparams()...)
   end

   sfig(ChE, [panel1, panel2], 
      "figure_5", ".png", 
      layout = (1,2)
   )

end

## ------------------------------------------------------
# Figure 6
include("2.4_maxent_plots.jl")

let
   params = lglob(Dyn, :dyn, :params, :finite_cg)
   @extract params: Ds ϵs cg simid

   ps = plot_flux_marginals(simid, Ds[8], ϵs[1:3:end], cg)

   for p in ps
      plot!(p; pltparams()...)
   end

   layout = (2, 2)
   ffile = sfig(ChE, ps, 
      "figure_6", ".png";
      layout
   )
   @info("Done", ffile)

end

## ------------------------------------------------------
# Figure 7
include("2.5_dyn_ug_corrs.jl")

let
   
   params = lglob(Dyn, :dyn, :params, :finite_cg)
   @extract params: Ds ϵs cg simid

   ps = _plot_dyn_ug_corrs(simid, Ds, ϵs[1:3:end], cg)

   for p in ps
      plot!(p; pltparams()...)
   end
   
   sfig(ChE, ps, 
      "figure_7", ".png";
      layout = (1, 3)
   )

end 