pro fitvth_spec_nov14o4,pid=pid,fid=fid,regid=regid, $
  fiter=fiter,pint=pint,monte=monte,plter=plter,$
  uncert_val=uncert_val,de=de,noplot=noplot

  ; Do a single thermal spectral fit to the November 04 data
  ; By default will fit region L3 for P1 and FPMA over 2.5-7 keV

  ; Optional inputs
  ; pid         -  Pointing of orbit 4 (1 or 2, default 1)
  ; fpid        -  FPM to use ('A' or 'B', default is 'A')
  ; regid       -  Region to fit (0 to 4 for ['D1','D2','L1','L2','L3'], default is 4 'L3')
  ; fiter       -  Energy range to fit over (default is 2.5 to last bin with >10 counts)
  ; pint        -  Ininitial parameters for the fit (default is [42,0.5,1.0])
  ; de          -  Energy binning of the spectrum and response (default 0.2 keV)
  ; monte       -  Do the monte carlo analysis on the errors?
  ; plter       -  Energy range to plot
  ; uncert_val  -  Systematic error to add to ospex via spex_uncert=uncert_val
  ; noplot      -  If keyword set do not to produce plot of histograms (>8.0 to do the plots)

  ; Warning -  Assumes the spectra and response files already exist in out_files/
  ; For more info on OSPEX see http://hesperia.gsfc.nasa.gov/ssw/packages/spex/doc/ospex_explanation.htm

  ; 09-Nov-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Which pointing to use
  if (n_elements(pid) ne 1) then pid=1
  pname='P'+string(pid,format='(i1)')
  ; Which FPM to use
  if (n_elements(fid) ne 1) then fid='A'
  fpmid='FPM'+fid
  ; Which region to do
  if (n_elements(regid) ne 1) then regid='4'
  regnmall=['D1','D2','L1','L2','L3']
  regnm=regnmall[regid]
  ; Starting point of the fit parameters
  if (n_elements(pint) ne 3) then pin=[42, 0.5, 1]
  ; What energy range to plot over ?
  if (n_elements(plter) ne 2) then plter=[1.6,8.5]
  ; Energy bin width of spectrum and response
  if (n_elements(de) ne 1) then de=0.2
  denm='DE'+strmid(string(1000+de*100,format='(i4)'),1,3)
  intnam='INT'+string(1000+100*(de),format='(i4)')

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Load in the spectrum
  restore,file='out_files/specs_'+denm+'_'+pname+fpmid+'.dat'

  spectrum=reform(spc_out.spec_reg[regid,*])
  nspec=n_elements(spectrum)
  sengs=spc_out.eng_edges
  smide=spc_out.engs
  regnm=spc_out.reg_nm[regid]
  dur=spc_out.dur
  lvt=spc_out.lvt_frac
  t1=spc_out.t1
  t2=spc_out.t2

  ; Load in the responses
  fl_rsp='out_files/RSP_'+pname+regnm+fpmid+'_'+intnam+'.dat'
  restore,file=fl_rsp

  ; Assuming saved rsp is just RMF*ARF
  ; so need to /dE to get into (counts/photons)/keV can check via o->getdata(class='spex_drm')
  rsp=rsp/dE

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; start up ospex in noninteractive mode
  set_logenv, 'OSPEX_NOINTERACTIVE', '1
  o = ospex()
  o->set,  spex_autoplot_enable=0, spex_fit_progbar=0, spex_fitcomp_plot_resid=0,spex_fit_manual=0
  o->set, spex_data_source='SPEX_USER_DATA'
  z = o->get(/spex_data_origunits)
  z.data_name = 'NuSTAR'
  o->set, spex_data_origunits=z
  o->set, spex_detectors=fpmid
  o->set, spex_ut_edges=anytim([t1,t2])
  o->set, spex_fit_time_interval=o->get(/spex_ut_edges)
  o->set, spex_ct_edges=edges
  o->set, mcurvefit_itmax=1000
  o->set, mcurvefit_tol=1e-4

  ; Our data so bother with the "use expected errors"
  o->set, spex_error_use_expected=0
  ; Add in 5% uncertainty
  o->set, spex_uncert=0.05
  if (n_elements(uncert_val) ne 0)  then o->set, spex_uncert=uncert_val

  ; Just in counts (input should not be rate or flux)
  spec_in=spectrum
  ; error estimate just Poisson on actual detected counts (again input should be just counts)
  err_in=sqrt(spectrum)
  ; Aactual livetime in seconds
  ltime_in=fltarr(n_elements(spectrum))+(dur*lvt)

  ; Put them and the response into ospex
  o->set, spectrum=spec_in, errors=err_in, livetime=ltime_in, spex_respinfo=rsp

  ; below what number of counts would we consider it to be non-Gaussian (i.e. Poisson) uncertainties?
  gausslim=10.
  if (n_elements(fiter) ne 2) then begin
    ; If no fit range is specified then we fit from 2.5 keV to bins with >10 counts
    ; staring bin in spec_in often 0 so avoid that as well
    maxfite=min(where(spec_in[1:n_elements(spec_in)-1] le gausslim))
    fiter=[2.5,edges[1,maxfite]]
  endif

  ; setup the fit and then do it
  o->set, spex_erange=fiter
  o->set, fit_function='vth'
  o->set, fit_comp_minima= [1e-10, 0.1, 0.5]
  o->set, fit_comp_maxima= [1e10,2.0, 1.5]
  o->set, fit_comp_param=pin
  o->set, fit_comp_free = [1, 1, 0]
  o->dofit, /all
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; get the resulting fitted spectrum model and params
  model = o -> calc_func_components(spex_units='flux',/counts,/all_func,this_interval=0)
  parm=o->get(/spex_summ_params)
  parmerr=o->get(/spex_summ_sigmas)
  chisq=o->get(/spex_summ_chisq)
  tkev=0.08617
  eranfit=o->get(/spex_erange)
  ; The original data in units of count flux
  ; should be the same as spec_in/dE/(lvt*dur)
  obs=o->calc_summ(item='data_count_flux',errors=err,this_interval=0)

  ; Also get it in terms of the photon flux
  pmodel = o -> calc_func_components(spex_units='flux',/photons,/all_func,this_interval=0)
  pobs=o->calc_summ(item='data_photon_flux',errors=perr,this_interval=0)
  id0=where(obs eq 0,nid0)
  if (nid0 gt 0) then begin
    pobs[id0]=0.
    perr[id0]=0.
  endif

  engs=get_edges(model.ct_energy,/mean)
  np=n_elements(engs)

  ; Save all the spectra and fit info out
  fit_out={pname:pname,fpm:fpmid,regid:regnm,$
    eranfit:eranfit,engs:engs,$
    parm:parm,parmerr:parmerr,tmk:parm[1]/tkev,chisq:chisq[0],$
    cnt:spec_in,ecnt:err_in,cnt_mod:model.yvals[*]*dE*lvt*dur,$
    cnt_flx:obs,ecnt_flx:err,cnt_flx_mod:model.yvals[*],$
    ph_flx:pobs,eph_flx:perr,ph_flx_mod:pmodel.yvals[*],$
    dur:dur,lvt_frac:lvt,exp:dur*lvt}

  fname=pname+regnm+fpmid+'_'+intnam
  save,file='out_files/fitvth_'+fname+'.dat',fit_out

  ; The chisq calculated is the reduced chisq, i.e.
  id=where(fit_out.engs ge eranfit[0] and fit_out.engs le eranfit[1])
  aa=(fit_out.cnt_flx-fit_out.cnt_flx_mod)^2/fit_out.ecnt_flx^2
  print,total(aa[id])/(n_elements(id)-2-1.), chisq

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Make a plot of the spectra and fit
  if (keyword_set(noplot) ne 1 and float(!version.release) ge 8.0) then begin
    w=window(dimensions=[400,600],/buffer)
    mrg=0.2
    ylim=[2e-1,2e4]

    flxmin=(fit_out.cnt_flx-fit_out.ecnt_flx) >ylim[0]
    flxmax=(fit_out.cnt_flx+fit_out.ecnt_flx) <ylim[1]

    id=where(fit_out.cnt_flx gt 0.,nid)

    ;    p0=errorplot(fit_out.engs,fit_out.cnt_flx,fit_out.ecnt_flx,/ylog,lines=6,symbol='+',title=pname+' '+fpmid+' '+regnm,$
    ;      yrange=ylim,xrange=plter,ytickunits='Scientific',errorbar_capsize=0,sym_thick=2,errorbar_thick=2,$
    ;      xtitle='',ytitle='count s!U-1!N keV!U-1!N',position=[0.175,0.3,0.975,0.94],/current,xtickformat='(a1)')

    p0=plot(fit_out.engs,fit_out.cnt_flx,/ylog,lines=6,symbol='x',title=pname+' '+fpmid+' '+regnm,$
      yrange=ylim,xrange=plter,ytickunits='Scientific',sym_thick=2,$
      xtitle='',ytitle='count s!U-1!N keV!U-1!N',position=[0.175,0.3,0.975,0.94],/current,xtickformat='(a1)')

    for i=0,nid-1 do !null=plot(fit_out.engs[id[i]]*[1,1],[flxmin[id[i]],flxmax[id[i]]],thick=2,/over,/current)

    plmod=plot(fit_out.engs,fit_out.cnt_flx_mod,color='firebrick',thick=2,/over,/current)
    !null=plot(eranfit[0]*[1,1],ylim,color='gray',lines=2,/over,/current)
    !null=plot(eranfit[1]*[1,1],ylim,color='gray',lines=2,/over,/current)

    resd=(fit_out.cnt_flx-fit_out.cnt_flx_mod)/fit_out.ecnt_flx
    bd=where(finite(resd) ne 1)
    resd[bd]=0
    pres=plot(fit_out.engs,resd,$
      yrange=[-4.5,4.5],xrange=plter,xtit='Energy [keV]',/stair,thick=2,$
      position=[0.175,0.08,0.975,0.28],ytit='(Obs-Mod)/Err',/current)

    !null=text(360,520,string(parm[1]/tkev,format='(f4.1)')+'$\pm$'+$
      string(parmerr[1]/tkev,format='(f4.2)')+' MK ('+string(parm[1],format='(f4.2)')+' keV)',$
      /device,color='firebrick',align=1)
    !null=text(360,500,string(parm[0]*1e3,format='(f4.1)')+'$\pm$'+string(parmerr[0]*1e3,format='(f4.2)')+$
      '$\times$10!U46!N cm!U-3!N',/device,color='firebrick',align=1)
    !null=text(360,480,'$\chi^2=$ '+strcompress(string(chisq,format='(f5.1)'),/rem),/device,color='firebrick',align=1)


    w.save,'figs/fitvth_'+fname+'.pdf',page_size=w.dimensions/100.
    w.close
  endif

end