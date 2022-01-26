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
      hist, n_ave_conv!, Vi, 
      load_simdat, simdat_file, get_status,
      UNDONE_SIM_STATUS, DEAD_SIM_STATUS,
      EXPLODED_SIM_STATUS, NITERS_SIM_STATUS, 
      STST_SIM_STATUS, collect_ts

   using Plots
   using Plots.Measures
   import GR
   !isinteractive() && GR.inline("png")

   using Base.Threads
   using ProgressMeter
   using ColorSchemes
   using LaTeXStrings

end

# ------------------------------------------------------
function _textbf(strs...)
   strs = replace.(string.(strs), " " => "~")
   latexstring(string("\\mathbf{", strs..., "}"))
end

## ------------------------------------------------------
# Collect glucose-limited STST
include("2.1_get_glclim_stst.jl")

simid = :SimD3
params = lglob(Dyn, simid, :params, :finite_cg)
@extract params: Ds ϵs cg
is_glclim_stst = get_glclim_stst()

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
# Figure 4: dyn_time_series
include("2.1_time_series.jl")
let

   # time series
   Di = 8
   Xp, sgp = plot_time_series(simid, Ds[Di], ϵs, cg)

   # plot
   for p in [Xp, sgp]
      plot!(p;
         titlefont = font(27),
         guidefont = font(21),
         xtickfont = font(21),
         ytickfont = font(21),
         thickness_scaling = 1.6,
         size = (1220, 940)
      )
   end

   sfig(ChE, [Xp, sgp],
      "figure_4", ".png"; 
      layout = (1, 2)
   )
end

## ------------------------------------------------------
# Figure 5
include("2.3_Pol_Volume_Study.jl")
include("2.4_glucose_limiting_study.jl")

let
   ϵs_ = ϵs
   panel1 = plot_pol_vol_study(simid, Ds, ϵs_, cg)
   panel2, panel3 = plot_glc_limiting_study(simid, Ds, ϵs_, cg, is_glclim_stst)

   for (title, p) in zip('A':'Z', [panel1, panel2, panel3])
      plot!(p; 
         titlefont = font(22),
         guidefont = font(19),
         xtickfont = font(18),
         ytickfont = font(18),
         legendfont = font(14),
         thickness_scaling = 1.6,
         size = (1220, 940), 
         title
      )
   end

   sfig(ChE, [panel1, panel2, panel3], 
      "figure_5.png"; 
      layout = (1, 3)
   )

end

## ------------------------------------------------------
# Figure 6
include("2.5_flx_corrs.jl")

let

   ps = plot_flx_corrs(simid, Ds, ϵs, cg, is_glclim_stst)

   for p in ps
      plot!(p;
         guidefont = font(26),
         xtickfont = font(21),
         ytickfont = font(21),
         thickness_scaling = 1.6,
         size = (1220, 900),
         rightmargin = 2mm
      )
   end 
   
   @assert length(ps) == 3 * 6
   sfig(ChE, ps, 
      "figure_6", ".png";
      layout = (6, 3)
   )

end 

## ------------------------------------------------------
# Figure 7
include("2.4_plot_marginal.jl")

let
   D_ = Ds[6]
   @show D_
   ϵs_ = ϵs[1:3:end]
   ps = plot_flux_marginals3D(simid, D_, ϵs_, cg, is_glclim_stst)

   plot!(ps[1]; bottom_margin=5mm)
   plot!(ps[4]; bottom_margin=5mm)
   for p in ps
      plot!(p; titlefont = font(22),
         guidefont = font(22),
         xtickfont = font(18),
         ytickfont = font(18),
         legendfont = font(14),
         thickness_scaling = 1.6,
         size = (1220, 940)
      )
   end

   layout = (2, 3)
   ffile = sfig(ChE, ps, 
      "figure_7", ".png";
      layout
   )
   @info("Done", ffile)

end


## ------------------------------------------------------
# Figure 11
include("2.2_sim3D_dead_explosion_map.jl")
let
   p = plot_dead_explosion_map(simid)

   plot!(p; 
      titlefont = font(22),
      guidefont = font(22),
      xtickfont = font(18),
      ytickfont = font(18),
      legendfont = font(14),
      thickness_scaling = 1.6,
      size = (1220, 940)
   )

   ffile = sfig(ChE, p, 
      "figure_11", ".png"
   )
   @info("Done", ffile)

end