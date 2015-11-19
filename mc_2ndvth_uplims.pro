pro mc_2ndvth_uplims,regid=regid,de=de,lime=lime

  ; Takes a vth from the ospex fit to the data and then adds in another component to see what is possible
  ; fold through the NuSTAR respose for a given obs setup
  ; most here are region of 120" then different effective area based on where in FoV

  ; Optional inputs
  ; regid       -  Region to fit (0 to 4 for ['D1','D2','L1','L2','L3'], default is 4 'L3')
  ; de          -  Energy binning of the spectrum and response (default 0.2 keV)
  ; lime        -  Energy above which the upper limits are calculated (default: 6keV D1,D2,L1 & 7keV L2, L3)
  ; 13-Nov-2015 IGH
  ; 16-Nov-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Which region to do
  if (n_elements(regid) ne 1) then regid='4'
  regnmall=['D1','D2','L1','L2','L3']
  regnm=regnmall[regid]
  ; Use FPMA for both
  ; Using P2 for D1 and L1
  ; Using P1 for D2, L2 and L3
  if (regid eq 0 or regid eq 2) then pid=2
  if (regid eq 1 or regid eq 3 or regid eq 4) then pid=1
  pname='P'+string(pid,format='(i1)')
  fid='A'
  fpmid='FPM'+fid
  ; Energy bin width of spectrum and response
  if (n_elements(de) ne 1) then de=0.2
  denm='DE'+strmid(string(1000+de*100,format='(i4)'),1,3)
  intnam='INT'+string(1000+100*(de),format='(i4)')
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Load in the fit results for chosen P, FPM and reg
  fname=pname+regnm+fpmid+'_'+intnam
  restore,file='out_files/fitvth_'+fname+'.dat'

  ; Load in the responses
  fl_rsp='out_files/RSP_'+pname+regnm+fpmid+'_'+intnam+'.dat'
  restore,file=fl_rsp

  lvt=fit_out.lvt_frac
  dur=fit_out.dur
  engs=fit_out.engs
  spectrum=fit_out.cnt_flx
  ns=n_elements(spectrum)
  obs=fit_out.cnt
  mod1_counts=fit_out.cnt_mod
  mk2kev=0.08617

  parms=fit_out.parm

  tmks2=5+findgen(8)
  nt2=n_elements(tmks2)
  ems2=dblarr(nt2)
  totsg6=intarr(nt2)

  ; If the energy isn't specified use these defaults
  ; For D1, D2, L1 want the lims of 4 counts >6keV
  ; For L2, L3 want the lims of 4 counts >7keV (as emission brighter and already at that lim)
  if (n_elements(lime) ne 1)  then begin
    if (regid le 2) then lime=6.0
    if (regid gt 2) then lime=7.0
  endif
  
  glime='G'+string(lime, format='(i1)')


  ids=where(engs ge lime)
  idsb=where(engs lt lime and engs ge 2.0)

  clearplot
  !p.charsize=2
  !p.thick=4
  !p.multi=[0,2,1]
  window,0,xsize=800,ysize=500

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;     Find the maximum emission measure for each temperature bin
  for ii=0, nt2-1 do begin
    fmod1 = f_vth( edges, parms )
    cntr=0
    print,tmks2[ii]
    ; What is the limit of counts >6keV that you want to have?
    ; 1\sigma is 1.841
    ; 2\sigma is 3.783
    ; 3\sigma is 6.608
    ; via Gehrels ApJ 1986 Table 1
    c6lim=3.783
    ; Just a large starting value to make sure it loops
    totg6=10
    
    ; Addtional constraint on model 2 EM so that  to be less than model 1
    ph6lim=total(fmod1[idsb])
    ; Just something bigger than above so while loop runs
    phtotl6=2*ph6lim

    ; Stop once the EM of the second f_vth is small enough so that less than c6lim counts
    ; otherwise reduce it and iterate again
    while (totg6 gt c6lim or phtotl6 gt ph6lim) do begin

      if (ii eq 0) then ems0=2.*parms[0] else ems0=ems2[ii-1]

      ems2[ii]=10d^(alog10(ems0)-cntr*0.01)
      parms2=[ems2[ii],tmks2[ii]*mk2kev,1]
      fmod2 = f_vth( edges, parms2 )
      fmod=fmod1+fmod2
      phtotl6=total(fmod2[idsb])

      ; Monte Carlo to get a synthetic NuSTAR spectra for this model and setup
      cntsflx=(rsp*de)#fmod
      syn_spec = cntsflx*de*lvt*dur
      syn_cts = round(total(syn_spec))
      syn_spec_pdf = total(syn_spec, /cum)
      syn_spec_pdf /= max(syn_spec_pdf)
      sim_spec = fltarr(n_elements(syn_spec))
      for kk =0, syn_cts - 1 do sim_spec[min(where(syn_spec_pdf gt randomu(seed)))]++

      ;  How many counts >specific energy in the synthetic spetrum
      totg6=total(sim_spec[ids])

      cntr++

      print,tmks2[ii],ems2[ii]*1d49*1d-46,totg6

      ; Some plot just to see what the model and synthetic count spectra look like
      plot,engs,fmod,/ylog,yrange=[1e-4,2e4],psym=10, title=pname+''+fpmid+''+regnm
      tube_line_colors
      oplot,engs,fmod1,color=1
      oplot,engs,fmod2,color=5
      
      xyouts, 9,1e5,'<6keV: '+ string(ph6lim,format='(i5)'),/data,chars=1.5,align=1,color=1
      xyouts, 9,3e4,'<6keV: '+ string(phtotl6,format='(i5)'),/data,chars=1.5,align=1,color=5
      xyouts, 9,1e4,string(tmks2[ii],format='(i2)')+'MK',/data,chars=1.5,align=1,color=5


      plot,engs,sim_spec,/ylog,yrange=[0.5,1e4],psym=10,ytickf='exp1',ystyle=17
      xyouts, 9,5e3,'>6keV: '+ string(totg6,format='(i3)'),/data,chars=1.5,align=1
 

      wait,0.001
      totsg6[ii]=totg6
    endwhile
  endfor

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !p.multi=0
  window,1,xsize=500,ysize=500
  plot,tmks2,ems2*1d49,/ylog,psym=1,syms=4,$
    xrange=[2,13],yrange=[1d42,1d48],$
    title=pname+''+fpmid+''+regnm,xtitle='T [MK]',ytit='EM [cm!U-3!N]'
  plots,parms[1]/mk2kev,parms[0]*1d49,psym=2


  em1=parms[0]*1d49
  tmk1=parms[1]/mk2kev
  
  uplims={reg:regnm, pname:pname,fpm:fpmid,$
    tmk1:tmk1,em1:em1,tmks2:tmks2,ems2:ems2*1d49}
  
  save,file='out_files/lims_2ndvth_'+fname+'_'+glime+'.dat',uplims

  ;  stop
end