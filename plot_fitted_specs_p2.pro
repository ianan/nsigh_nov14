pro plot_fitted_specs_p2,plter=plter,de=de

  ; PLot showing the regions selected in P2 and the resulting spectra and fits for FPMA and FPMB
  ; Assumes make_map_nov14o4.pro and fitvth_spec_nov14o4.pro has already been run
  ;
  ;  Will also need tvim.pro from the esrg package
  ;
  ; Optional inputs:
  ;   de        - Energy binning of spectrum (default 0.2keV)
  ;   plter       -  Energy range to plot
  ;
  ; For colortable 74 need IDL >8.0
  ;
  ; 09-Nov-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  @post_outset
  pid=2
  pname='P'+string(pid,format='(i1)')
  grdid=0
  grdname='G'+string(grdid,format='(i1)')
  if (n_elements(plter) ne 2) then plter=[2.3,7.9]

  if (n_elements(de) eq 0) then de=0.2
  denm='DE'+strmid(string(1000+de*100,format='(i4)'),1,3)
  intnam='INT'+string(1000+100*(de),format='(i4)')

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; The new regions of interest
  ; in S/C arcsec
  ; loads in the roi, xcsp, ycsp,rnm,wid
  restore,file='out_files/roi_no14.dat'
  nr=n_elements(rnm)

  !p.multi=[0,3,1]
  loadct,74,/silent
  reverse_ct
  tvlct,r,g,b,/get
  r[0]=0
  g[0]=0
  b[0]=0
  r[1]=255
  g[1]=255
  b[1]=255
  tvlct,r,g,b

  ; For P2 only want to do it for D1 and L1
  idds=[0,2]
  nr=n_elements(idds)

  set_plot,'ps'
  device, /encapsulated, /color, /isolatin1,/inches, $
    bits=8, xsize=12, ysize=3.5,file='figs/fittedSpecsall_'+pname+'_'+intnam+'.eps'
  !p.charsize=2.0
  dfac=1.
  dfac0=1.

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Load in the maps, plot and then overplot the selected regions
  fits2map,'out_files/maps_'+pname[0]+grdname+'_EG2_FPMA.fits',mma_in
  fits2map,'out_files/maps_'+pname[0]+grdname+'_EG2_FPMB.fits',mmb_in
  sub_map,mma_in,mma,xrange=[500,1500],yrange=[-150,850]
  sub_map,mmb_in,mmb,xrange=[500,1500],yrange=[-150,850]

  mall=[mma.data,mmb.data]
  nx=n_elements(mma.data[*,0])
  xa0=mma.xc-0.5*mma.dx*nx
  ya0=mma.yc-0.5*mma.dy*nx

  xcspa=(xcsp-xa0)/mma.dx
  xcspb=nx+(xcsp-xa0)/mma.dx
  ycspa=(ycsp-ya0)/mma.dx
  wida=wid/mma.dx

  tvim,alog10(gauss_smooth(mall,3)),range=[-2,1],$
    pcharsize=0.01,position=[0.01,0.2,0.35,0.8],/noframe
  for xx=0,nr-1 do box_igh,xcspa[idds[xx]],ycspa[idds[xx]],wida,color=0,thick=2
  for xx=0,nr-1 do box_igh,xcspb[idds[xx]],ycspa[idds[xx]],wida,color=0,thick=2
  for xx=0,nr-1 do xyouts,xcspa[idds[xx]]+0.5*wida+5,ycspa[idds[xx]]-0.5*wida-5,rnm[idds[xx]],color=0,chars=0.8
  for xx=0,nr-1 do xyouts,xcspb[idds[xx]]+0.5*wida+5,ycspa[idds[xx]]-0.5*wida-5,rnm[idds[xx]],color=0,chars=0.8


  tube_line_colors
  xyouts, 0.5*nx,nx,'P1 FPMA',align=0.5,/data,color=1,chars=1
  xyouts, 1.5*nx,nx,'P1 FPMB',align=0.5,/data,color=2,chars=1

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Label with the livetime information - so need to load in fits to get it
  restore,file='out_files/fitvth_'+pname+rnm[0]+'FPMA_'+intnam+'.dat'
  xyouts, 0.5*nx,-0.01*nx,strcompress(string(fit_out.exp,format='(f5.2)'),/rem)+'s ('+$
    strcompress(string(fit_out.lvt_frac*100.,format='(f5.2)')+'%)',/rem),align=0.5,/data,color=1,chars=1

  restore,file='out_files/fitvth_'+pname+rnm[0]+'FPMB_'+intnam+'.dat'
  xyouts, 1.5*nx,-0.01*nx,strcompress(string(fit_out.exp,format='(f5.2)'),/rem)+'s ('+$
    strcompress(string(fit_out.lvt_frac*100.,format='(f5.2)')+'%)',/rem),align=0.5,/data,color=2,chars=1

  tmk=0.08617
  ylim=[3e-1,2e4]*dfac
  xlim=plter

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; For each region plot the FPMA and FPMB spectra and fits
  for xx=0,nr-1 do begin
    i=idds[xx]
    restore,file='out_files/fitvth_'+pname+rnm[i]+'FPMA_'+intnam+'.dat'

    plot,[1,1],[1,1],xtit='Energy [keV]',xrange=xlim,/ylog,$
      yrange=ylim,ytit='counts keV!U-1!N s!U-1!N',ytickf='exp1',$
      ;        yrange=ylim,ytit='counts s!U-1!N cm!U-2!N keV!U-1!N',ytickf='exp1',$
      title=pname+' '+rnm[i]
    oplot,[1,10],[10.,10]/de/fit_out.exp ,lines=1,color=150
    oploterr,fit_out.engs,fit_out.cnt_flx,fit_out.ecnt_flx,color=1,psym=10,thick=5,/nohat,errthick=3,errcolor=1
    ;oploterr,fit_out.engs,fit_out.cnt_flx,fit_out.ecnt_flx,color=1,psym=1,thick=4,/nohat,errthick=4,errcolor=1
    oplot,fit_out.engs,fit_out.cnt_flx_mod,color=1,thick=5


    xyouts,plter[0]+0.95*(plter[1]-plter[0]),10d^3.9,'A: '+string(fit_out.tmk,format='(f4.1)')+' MK, '$
      +string(fit_out.parm[0]*1e3,format='(f4.1)')+' x10!U46!N cm!U-3!N',color=1,chars=1.,/data,align=1

    oplot,fit_out.eranfit[0]*[1,1],ylim,line=2,thick=2,color=150
    oplot,fit_out.eranfit[1]*[1,1],[ylim[0],2e3],line=2,thick=2,color=150

    restore,file='out_files/fitvth_'+pname+rnm[i]+'FPMB_'+intnam+'.dat'

    oploterr,fit_out.engs,fit_out.cnt_flx,fit_out.ecnt_flx,color=2,psym=10,thick=5,/nohat,errthick=3,errcolor=2
    ;oploterr,fit_out.engs,fit_out.cnt_flx,fit_out.ecnt_flx,color=2,psym=1,thick=4,/nohat,errthick=4,errcolor=2
    oplot,fit_out.engs,fit_out.cnt_flx_mod,color=2,thick=5

    xyouts,plter[0]+0.95*(plter[1]-plter[0]),10d^3.5,'B: '+string(fit_out.tmk,format='(f4.1)')+' MK, '$
      +string(fit_out.parm[0]*1e3,format='(f4.1)')+' x10!U46!N cm!U-3!N',color=2,chars=1.,/data,align=1


  endfor

  device,/close
  set_plot, mydevice


end