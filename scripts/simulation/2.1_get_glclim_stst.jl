function get_glclim_stst()
   
    simid = :SimD3
    params = lglob(Dyn, simid, :params, :finite_cg)
    @extract params: Ds ϵs cg
    
    validator = Dict()
    for (ϵi, ϵ) in enumerate(ϵs)
       for D in Ds
 
          
          ## ------------------------------------------------------
          # Dynamic
          simparams = (;D, ϵ, cg)
          validator[simparams] = false
          status = get_status(simid, simparams)
          (status != STST_SIM_STATUS) && continue
          datfile = Dyn.simdat_file(simid, simparams)
          !isfile(datfile) && continue
          S = ldat(datfile)
          dyn_ug = S.ug_av
          dyn_X = S.X
          
          ## ------------------------------------------------------
          # LP 
          net = S.V.net
          Dyn.fix!(net, :z, D)
          Dyn.ub!(net, :ug, cg * D / S.X)
          Dyn.lb!(net, :ug, 0.0)
          Dyn.reducebox!(net)
          Lug, Uug = Dyn.fva(net, :ug)
          UX = cg * D / Lug
          
          #  glc limited
          abs(Uug - dyn_ug) / abs(Uug + dyn_ug) > 0.1 && continue
 
 
          validator[simparams] = true
             
       end
    end

    return validator
 end