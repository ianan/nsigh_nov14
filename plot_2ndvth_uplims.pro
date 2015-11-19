pro plot_2ndvth_uplims

  ; Plot upper limits to emission from second f_vth produced by mc_2ndvth_uplims
  ;
  ; 17-Nov-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  regnm=['D1','D2','L1','L2','L3']
  nr=n_elements(regnm)
  ;  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;  ; Do D1,D2,L1 >6 keV, L2, L3 >7keV
  ;  ts1b=fltarr(nr)
  ;  es1b=dblarr(nr)
  ;  ts2b=fltarr(nr,8)
  ;  es2b=dblarr(nr,8)
  ;  regb=strarr(nr)
  ;  pnameb=strarr(nr)
  ;
  ;  for i=0,2 do begin
  ;    fname='out_files/lims*'+regnm[i]+'*G6*.dat'
  ;    ff=file_search(fname)
  ;    restore,file=ff
  ;    ts1b[i]=uplims.tmk1
  ;    es1b[i]=uplims.em1
  ;    ts2b[i,*]=uplims.tmks2
  ;    es2b[i,*]=uplims.ems2
  ;    regb[i]=uplims.reg
  ;    pnameb[i]=uplims.pname
  ;  endfor
  ;
  ;  for i=3,4 do begin
  ;    fname='out_files/lims*'+regnm[i]+'*G7*.dat'
  ;    ff=file_search(fname)
  ;    restore,file=ff
  ;    ts1b[i]=uplims.tmk1
  ;    es1b[i]=uplims.em1
  ;    ts2b[i,*]=uplims.tmks2
  ;    es2b[i,*]=uplims.ems2
  ;    regb[i]=uplims.reg
  ;    pnameb[i]=uplims.pname
  ;  endfor
  ;
  ;  @post_outset
  ;  !p.thick=4
  ;  loadct,0,/silent
  ;  tube_line_colors
  ;  !p.multi=0
  ;
  ;  set_plot,'ps'
  ;  device, /encapsulated, /color, /isolatin1, $
  ;    /inches, bits=8, xsize=5, ysize=4,file='figs/all_lims_2ndvth_G6G7.eps'
  ;
  ;  xr=[1.5,13.5]
  ;  plot,ts2b,es2b,/ylog,xtit='T [MK]',ytit='Max(EM!D2!N) [cm!U-3!N]',$
  ;    thick=4,chars=1.1,yrange=[2d41,2d49],xrange=xr,psym=1,/nodata,xstyle=9,position=[0.15,0.15,0.95,0.9]
  ;  plotsym,1,2,thick=4
  ;
  ;  ct=[1,2,5,3,8]
  ;  for i=0, nr-1 do begin
  ;    plots,ts1b[i],es1b[i],color=ct[i],psym=1,syms=1
  ;    dt=ts2b[0,1]-ts2b[0,0]
  ;    for tt=0, n_elements(ts2b[0,*])-1 do oplot,ts2b[i,tt]+[-0.5*dt,0.5*dt],[1,1]*es2b[i,tt],thick=3,color=ct[i]
  ;    oplot,ts2b[i,*],es2b[i,*],psym=8,thick=3,color=ct[i]
  ;    xyouts,10,10d^(48.5-0.5*i),regb[i]+' ('+pnameb[i]+')',color=ct[i],/data,chars=1
  ;  endfor
  ;
  ;  axis, xaxis=1,xtitle='log!D10!N T', xrange=xr,xtickn=string(alog10([2,4,6,8,10,12]*1d6),format='(f4.2)')
  ;
  ;  device,/close
  ;  set_plot, mydevice

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Do all with >6 keV
  ts1=fltarr(nr)
  es1=dblarr(nr)
  ts2=fltarr(nr,8)
  es2=dblarr(nr,8)
  reg=strarr(nr)
  pname=strarr(nr)

  for i=0,nr-1 do begin
    fname='out_files/lims*'+regnm[i]+'*G6*.dat'
    ff=file_search(fname)
    restore,file=ff
    ts1[i]=uplims.tmk1
    es1[i]=uplims.em1
    ts2[i,*]=uplims.tmks2
    es2[i,*]=uplims.ems2
    reg[i]=uplims.reg
    pname[i]=uplims.pname
  endfor

  @post_outset
  !p.thick=4
  loadct,0,/silent
  tube_line_colors
  !p.multi=0

  set_plot,'ps'
  device, /encapsulated, /color, /isolatin1, $
    /inches, bits=8, xsize=5, ysize=4,file='figs/all_lims_2ndvth_G6.eps'

  xr=[1.4,13.]
  ylim=[2d41,9d49]
  plot,ts2,es2,/ylog,xtit='Temperature [MK]',ytit='Emission Measure [cm!U-3!N]',$
    thick=4,chars=1.1,yrange=ylim,xrange=xr,psym=1,/nodata,xstyle=9,position=[0.15,0.15,0.95,0.9]
  plotsym,1,2,thick=4

  ct=[1,2,5,3,8]
  for i=0, nr-1 do begin
    plots,ts1[i],es1[i],color=ct[i],psym=1,syms=1
    dt=ts2[0,1]-ts2[0,0]
    for tt=0, n_elements(ts2[0,*])-1 do oplot,ts2[i,tt]+[-0.5*dt,0.5*dt],[1,1]*es2[i,tt],thick=3,color=ct[i]
    oplot,ts2[i,*],es2[i,*],psym=8,thick=3,color=ct[i]
    xyouts,10,10d^(49.0-0.5*i),reg[i]+' ('+pname[i]+')',color=ct[i],/data,chars=1
  endfor

  xyouts, 2.5,10d^42.4,'IsoTh',color=0,/data,chars=0.8,align=0
  xyouts, 2.5,10d^42.0,'Up Lim',color=0,/data,chars=0.8,align=0

  plots,2,10d^42.5,psym=1,syms=0.8
  plots,2,10d^42.1,psym=8,syms=0.6
  oplot,[1.85,2.15],[10d^42.1,10d^42.1]

  axis, xaxis=1,xtitle='log!D10!N T', xrange=xr,xtickn=string(alog10([2,4,6,8,10,12]*1d6),format='(f4.2)')

  device,/close
  set_plot, mydevice

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Plot in comparison to other small FoV observations

  set_plot,'ps'
  device, /encapsulated, /color, /isolatin1, $
    /inches, bits=8, xsize=5, ysize=4,file='figs/all_lims_2ndvth_G6_compare.eps'


  plot,ts2,es2,/ylog,xtit='Temperature [MK]',ytit='Emission Measure [cm!U-3!N]',$
    thick=4,chars=1.1,yrange=ylim,xrange=xr,psym=1,/nodata,xstyle=9,position=[0.15,0.15,0.95,0.9]
  plotsym,1,2,thick=4

  ct=replicate(120,5);[1,2,5,3,8]
  for i=0, nr-1 do begin
    plots,ts1[i],es1[i],color=ct[i],psym=1,syms=1
    dt=ts2[0,1]-ts2[0,0]
    for tt=0, n_elements(ts2[0,*])-1 do oplot,ts2[i,tt]+[-0.5*dt,0.5*dt],[1,1]*es2[i,tt],thick=3,color=ct[i]
    oplot,ts2[i,*],es2[i,*],psym=8,thick=3,color=ct[i]
    ;      xyouts,10,10d^(48.5-0.5*i),reg[i]+' ('+pname[i]+')',color=ct[i],/data,chars=1
  endfor

  ; Del Zanna & Mason 2014
  ; min of the EM loci curve (min EM< max T of formation)
  ; but if steep DEM spectra then over estimate
  tz1=10d^6.9/1d6;10d^6.7/1d6
  emz1=10d^27.3*(15*7.25d7)^2.0;10d21*(10d^6.725-10d^6.675)*(15*7.25d7)^2.0
  plots,tz1,emz1,color=1,psym=6,thick=3

  ; actual DEM values
  tz=10d^6.7/1d6
  emz=10d21*(10d^6.725-10d^6.675)*(15*7.25d7)^2.0
  plots,tz,emz,color=1,psym=8,thick=3
  oplot,tz+[-.25,.25],[emz,emz],color=1,thick=3

  tz=10d^[6.45,6.5]/1d6
  emz=[10d^23.7*(10d^6.475-10d^6.425),10d^21.2*(10d^6.525-10d^6.475)]*(15*7.25d7)^2.0
  oplot,tz,emz,color=1,psym=-7,thick=3,lines=2

  ; Skylark/KAP for Parkinson 1975 and Del Zanna & Mason 2014
  tz=10d^[6.4,6.5,6.6]/1d6
  ; Assuming DlogT=0.05 and area is 180"x180" for conversion from cm-5 K-1 to cm-3
  emz=[10d^21*(10d^6.425-10d^6.375),10d^21*(10d^6.525-10d^6.475),$
    10d^20.5*(10d^6.625-10d^6.575)]*(180*7.25d7)^2.0
  oplot,tz,emz,color=5,psym=-7,thick=3,lines=2

  ; brosius et al. 2014 (EUNIS)
  ; again max of the T
  tb=[1.6,8.9]
  emb=[2.79d46,2.12d45]
  plots,tb,emb,color=8,psym=6,thick=3

  xyouts, 12,10d^49.0,'NuSTAR',color=ct[0],/data,chars=0.9,align=1
  xyouts, 12,10d^48.5,'SMM/FCS (DZM)',color=1,/data,chars=0.9,align=1
  xyouts, 12,10d^48.0,'Skylark/KAP (P+DZM)',color=5,/data,chars=0.9,align=1
  xyouts, 12,10d^47.5,'EUNIS (BDR)',color=8,/data,chars=0.9,align=1

  xyouts, 2.5,10d^43.2,'G(T!Dmax!N)',color=0,/data,chars=0.8,align=0
  xyouts, 2.5,10d^42.8,'DEM',color=0,/data,chars=0.8,align=0
  xyouts, 2.5,10d^42.4,'IsoTh',color=0,/data,chars=0.8,align=0
  xyouts, 2.5,10d^42.0,'Up Lim',color=0,/data,chars=0.8,align=0

  plots,2,10d^43.3,psym=6,syms=0.8
  oplot,[1.74,2.21],10d^42.9*[1,1],psym=-7,syms=0.8,lines=2
  plots,2,10d^42.5,psym=1,syms=0.8
  plots,2,10d^42.1,psym=8,syms=0.6
  oplot,[1.85,2.15],[10d^42.1,10d^42.1]


  axis, xaxis=1,xtitle='log!D10!N T', xrange=xr,xtickn=string(alog10([2,4,6,8,10,12]*1d6),format='(f4.2)')

  device,/close
  set_plot, mydevice
  
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Plot in comparison to Full Disc observations

  ; Load in the fraction each region is of the fulldisk flux in SDO/AIA/94A FeXVIII
  restore,file='out_files/roi_no14_ffper.dat'
  fracs=ff_out.prcntfe18/100.

  es2fd=es2
  for i=0,nr-1 do es2fd[i,*]=es2[i,*]/fracs[i]

  set_plot,'ps'
  device, /encapsulated, /color, /isolatin1, $
    /inches, bits=8, xsize=5, ysize=4,file='figs/all_lims_2ndvth_G6_comparefd.eps'

  plot,ts2,es2,/ylog,xtit='Temperature [MK]',ytit='Emission Measure [cm!U-3!N]',$
    thick=4,chars=1.1,yrange=ylim,xrange=xr,psym=1,/nodata,xstyle=9,position=[0.15,0.15,0.95,0.9]
  plotsym,1,2,thick=4

  ct=replicate(7,5);[1,2,5,3,8]
  for i=0, nr-1 do begin
    plots,ts1[i],es1[i]/fracs[i],color=ct[i],psym=1,syms=1
    dt=ts2[0,1]-ts2[0,0]
    for tt=0, n_elements(ts2[0,*])-1 do oplot,ts2[i,tt]+[-0.5*dt,0.5*dt],[1,1]*es2fd[i,tt],thick=3,color=ct[i]
    oplot,ts2[i,*],es2fd[i,*],psym=8,thick=3,color=ct[i]
  endfor

  ; caspi et al. 2015
  tc=[2.9,11]
  emc=[4.9d48,1.4d46]
  plots,tc[0],emc[0],color=2,psym=1;,syms=2
  plots,tc[1],emc[1],color=2,psym=1;,syms=2

  ; caspi DEM fit over 1.82 to 10.7MK, with index=6.31, norm 5e43 cm-3 keV-1 at 23.2MK (2kev)
  mk2kev=0.08617
  tcd=[1.82,10.7]
  nts=10
  tcds=tcd[0]+findgen(nts)*(tcd[1]-tcd[0])/(nts-1.0)
  tcks=tcds*mk2kev
  dem_c=5d43*(2/tcks)^(6.31)
  print,total(dem_c*(tcks[1]-tcks[0]))
  oplot,tcds,dem_c*(tcks[1]-tcks[0]),lines=2,color=2
  oplot,[max(tcds),min(tcds)],minmax(dem_c*(tcks[1]-tcks[0])),color=2,psym=7

  ; miceli et al. 2015
  ti=[2.73,6.6]
  emi=[63d46,0.027d46]
  plots,ti[0],emi[0],color=4,psym=1;,syms=2
  plots,ti[1],emi[1],color=4,psym=1;,syms=2

  ; McTiernan 2009
  tm=[6,8,11]
  emm=[1d46,1d47,1d48]

  oplot,tm[0:1],emm[1]*[1,1],color=11,lines=1,psym=-1
  oplot,tm[0:1],emm[0]*[1,1],color=11,lines=1,psym=-1
  oplot,tm[0]*[1,1],emm[0:1],color=11,lines=1,psym=-1
  oplot,tm[1]*[1,1],emm[0:1],color=11,lines=1,psym=-1

  xyouts, 12,10d^49.0,'NuSTAR FD Approx',color=ct[0],/data,chars=0.9,align=1
  xyouts, 12,10d^48.5,'CORONAS-Photon/SphinX (M)',color=4,/data,chars=0.9,align=1
  xyouts, 12,10d^48.0,'SDO/EVE/X123 (C)',color=2,/data,chars=0.9,align=1
  xyouts, 12,10d^47.5,'RHESSI (McT)',color=11,/data,chars=0.9,align=1

  xyouts, 2.5,10d^42.8,'DEM',color=0,/data,chars=0.8,align=0
  xyouts, 2.5,10d^42.4,'IsoTh',color=0,/data,chars=0.8,align=0
  xyouts, 2.5,10d^42.0,'Up Lim',color=0,/data,chars=0.8,align=0

  oplot,[1.74,2.21],10d^42.9*[1,1],psym=-7,syms=0.8,lines=2
  plots,2,10d^42.5,psym=1,syms=0.8
  plots,2,10d^42.1,psym=8,syms=0.6
  oplot,[1.85,2.15],[10d^42.1,10d^42.1]

  axis, xaxis=1,xtitle='log!D10!N T', xrange=xr,xtickn=string(alog10([2,4,6,8,10,12]*1d6),format='(f4.2)')

  device,/close
  set_plot, mydevice

  stop
end