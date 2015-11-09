pro plot_aiasxi_maps_nov14o4

  ; Produce comparison plots of the NuSTAR maps and SDO/AIA and GEOS/SXI
  ; FYI< the SXI and AIA maps aren't imcluded here and script just for
  ; reference as to how the figure is produced
  ;
  ; For some of the colour tables need IDL >8.0
  ;
  ; 09-Nov-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;
  ; Location of SDO/AIA mean map
  aiadir='~/data/aia/ev_01112014/'
  ; Location of GOES/SXI mean map
  sxidir='~/data/sxi/'

  grdid=0
  grdname='G'+string(grdid,format='(i1)')

  pids=[1,2,3,4]
  pname='P'+string(pids,format='(i1)')
  fits2map,'out_files/maps_'+pname[0]+grdname+'_E2_4_FPMA.fits',ma124
  fits2map,'out_files/maps_'+pname[1]+grdname+'_E2_4_FPMA.fits',ma224
  ma24=ma124
  ma24.data=ma124.data+ma224.data
  ma24.id='LC FPMA G0 2-4 keV'
  fits2map,'out_files/maps_'+pname[0]+grdname+'_E2_4_FPMB.fits',mb124
  fits2map,'out_files/maps_'+pname[1]+grdname+'_E2_4_FPMB.fits',mb224
  mb24=mb124
  mb24.data=mb124.data+mb224.data
  mb24.id='LC FPMB G0 2-4 keV'

  fits2map,'out_files/maps_'+pname[0]+grdname+'_E4_6_FPMA.fits',ma146
  fits2map,'out_files/maps_'+pname[1]+grdname+'_E4_6_FPMA.fits',ma246
  ma46=ma146
  ma46.data=ma146.data+ma246.data
  ma46.id='LC FPMA G0 4-6 keV'
  fits2map,'out_files/maps_'+pname[0]+grdname+'_E4_6_FPMB.fits',mb146
  fits2map,'out_files/maps_'+pname[1]+grdname+'_E4_6_FPMB.fits',mb246
  mb46=mb146
  mb46.data=mb146.data+mb246.data
  mb46.id='LC FPMB G0 4-6 keV'

  sr=4
  mb24s=mb24
  mb24s.data=gauss_smooth(mb24.data,sr)
  mb46s=mb46
  mb46s.data=gauss_smooth(mb46.data,sr)

  ma24s=ma24
  ma24s.data=gauss_smooth(ma24.data,sr)
  ma46s=ma46
  ma46s.data=gauss_smooth(ma46.data,sr)

  mab24=ma24
  mab24.data=gauss_smooth(ma24.data+mb24.data,sr)
  mab24.id='LC FPMA+FPMB G0 2-4keV'

  rmab24=ma24
  rmab24.data=ma24.data+mb24.data
  rmab24.id='LC FPMA+FPMB G0 2-4keV'
  rmab24=grid_map(rmab24,dx=12.5,dy=12.5)


  mab46=ma46
  mab46.data=gauss_smooth(ma46.data+mb46.data,sr)
  mab46.id='LC FPMA+FPMB G0 4-6keV'

  rmab46=ma46
  rmab46.data=ma46.data+mb46.data
  rmab46.id='LC FPMA+FPMB G0 4-6keV'
  rmab46=grid_map(rmab46,dx=12.5,dy=12.5)

  xr=[350,1350]
  yr=[-650,350]

  mar24=rebin_map(ma24,100,100)
  mbr24=rebin_map(mb24,100,100)
  mar46=rebin_map(ma46,100,100)
  mbr46=rebin_map(mb46,100,100)

  mabr24=mar24
  mabr24.data=mar24.data+mar24.data
  mabr46=mar46
  mabr46.data=mar46.data+mar46.data

  mrat=mabr24

  bd24=where(mabr24.data lt 5e-1)
  bd46=where(mabr46.data lt 4e-2)

  dmab24=mabr24.data
  dmab24[bd24]=0.
  dmab46=mabr46.data
  dmab46[bd46]=0.

  mrat.data=dmab46/dmab24

  @post_outset

  d=fltarr(2,6)
  d[*,0]=[2,42]
  d[*,1]=[4,1e3]
  d[*,2]=[1e2,3e3]
  d[*,3]=[1e2,3e3]
  d[*,4]=[1e2,3e3]
  d[*,5]=[4,100]
  wavn=['94','131','171','193','211','335']

  fits2map,aiadir+'P1_'+wavn[0]+'_mn.fits',m94
  fits2map,aiadir+'P1_'+wavn[5]+'_mn.fits',m335
  fits2map,aiadir+'P1_'+wavn[1]+'_mn.fits',m131
  fits2map,aiadir+'P1_FE18_mn.fits',mfe18


  fits2map,sxidir+'P1_SXI_Be12_mn.fits',ms
  ms.xc=ms.xc-45
  ms.yc=ms.yc+15


  set_plot,'ps'
  device, /encapsulated, /color, / isolatin1,/inches, $
    bits=8, xsize=12, ysize=6,file='figs/Orb4_8maps_FPMAB.eps'

  !p.multi=[0,4,2]
  !p.charsize=1.4

  xr=[350,1200]
  yr=[-550,300]

  ; just make a black background for the whole image
  loadct,0,/silent
  restore,file='out_files/roi_no14.dat'
  nr=n_elements(rnm)

  loadct,62;56
  reverse_ct
  tvlct,r,g,b,/get
  r[0]=0
  b[0]=0
  g[0]=0
  r[255]=255
  b[255]=255
  g[255]=255
  tvlct,r,g,b
  plot_map,mab24,/log,dmin=1.,dmax=18,$
    xrange=xr,yrange=yr,bottom=0,/limb,grid_spacing=20,gthick=2,lcolor=255,gcolor=255,$
    tit='NuSTAR 2-4 keV ('+ string(12.4/4.,format='(f3.1)')+'-'+string(12.4/2.,format='(f3.1)')+STRING(197B)+')'
  for i=0,nr-1 do box_igh,xcsp[i],ycsp[i],wid,color=255,thick=2

  loadct,64;49
  reverse_ct
  tvlct,r,g,b,/get
  r[0]=0
  b[0]=0
  g[0]=0
  r[255]=255
  b[255]=255
  g[255]=255
  tvlct,r,g,b
  plot_map,mab46,/log,dmin=0.03,dmax=0.4,$
    xrange=xr,yrange=yr,bottom=0,/limb,grid_spacing=20,gthick=2,lcolor=255,gcolor=255,$
    tit='NuSTAR 4-6 keV ('+string(12.4/6.,format='(f3.1)')+'-'+string(12.4/4.,format='(f3.1)')+STRING(197B)+')'
  for i=0,nr-1 do box_igh,xcsp[i],ycsp[i],wid,color=255,thick=2

  ;
  loadct,39,/silent
  tvlct,r,g,b,/get
  r[254]=255
  b[254]=255
  g[254]=255
  r[255]=0
  b[255]=0
  g[255]=0
  tvlct,r,g,b
  mrat=mab24
  bd24=where(mab24.data lt 1.0)
  bd46=where(mab46.data lt 0.03)
  ;
  dmab24=mab24.data
  dmab24[bd24]=0.
  dmab46=mab46.data
  dmab46[bd46]=0.

  mmnn=-1.0
  mmxx=-0.92

  mrat.data=(dmab46-dmab24)/(dmab46+dmab24);dmab46/dmab24
  sub_map,mrat,smrat,xrange=xr,yrange=yr
  plot_map,smrat,dmin=mmnn,dmax=mmxx,title='NuSTAR (I!D46!N-I!D24!N)/(I!D46!N+I!D24!N)',$
    xrange=xr,yrange=yr,bottom=35,/limb,grid_spacing=20,gthick=2,lcolor=254,gcolor=254
  for i=0,nr-1 do box_igh,xcsp[i],ycsp[i],wid,color=254,thick=2
  plot_map_cb_igh,[mmnn,mmxx],35,255-35,position=[0.57,0.925,0.66,0.94],color=254,plchars=1.1,$
    format='(f5.2)',major=4

  for i=0,nr-1 do xyouts,xcsp[i]+65,ycsp[i]-70,rnm[i],color=254,chars=0.8

  loadct,63,/silent
  reverse_ct
  tvlct,r,g,b,/get
  r[0]=0
  b[0]=0
  g[0]=0
  tvlct,r,g,b
  plot_map,mfe18,/log,color=0,dmin=2,dmax=25,$
    xrange=xr,yrange=yr,grid_spacing=20,lcolor=255,gcolor=255,/limb,tit='SDO/AIA '+wavn[0]+ STRING(197B)+'/FeXVIII'
  loadct,0,/silent
  for i=0,nr-1 do box_igh,xcsp[i],ycsp[i],wid,color=255,thick=2


  loadct,60,/silent
  reverse_ct
  tvlct,r,g,b,/get
  r[0]=0
  b[0]=0
  g[0]=0
  tvlct,r,g,b
  plot_map,ms,/log,/limb,dmin=72,dmax=2400,bot=0,$
    xrange=xr,yrange=yr,grid_spacing=20,lcolor=255,gcolor=255,tit='GOES/SXI 10-20'+STRING(197B)
  loadct,0,/silent
  for i=0,nr-1 do box_igh,xcsp[i],ycsp[i],wid,color=255,thick=2

  i=1
  aia_lct,wave=wavn[i],/load
  plot_map,m131,/log,color=0,dmin=d[0,i],dmax=d[1,i],$
    xrange=xr,yrange=yr,grid_spacing=20,lcolor=255,gcolor=255,/limb,/no_fill,tit='SDO/AIA '+wavn[i]+ STRING(197B)
  loadct,0,/silent
  for i=0,nr-1 do box_igh,xcsp[i],ycsp[i],wid,color=255,thick=2

  i=5
  aia_lct,wave=wavn[i],/load
  plot_map,m335,/log,color=0,dmin=d[0,i],dmax=d[1,i],$
    xrange=xr,yrange=yr,grid_spacing=20,lcolor=255,gcolor=255,/limb,tit='SDO/AIA '+wavn[i]+ STRING(197B)
  loadct,0,/silent
  for i=0,nr-1 do box_igh,xcsp[i],ycsp[i],wid,color=255,thick=2
  i=0
  aia_lct,wave=wavn[i],/load
  plot_map,m94,/log,color=0,dmin=d[0,i],dmax=d[1,i],$
    xrange=xr,yrange=yr,grid_spacing=20,lcolor=255,gcolor=255,/limb,tit='SDO/AIA '+wavn[i]+ STRING(197B)
  loadct,0,/silent
  for i=0,nr-1 do box_igh,xcsp[i],ycsp[i],wid,color=255,thick=2

  device,/close
  set_plot, mydevice

;  ffn='figs/Orb4_8maps_FPMAB'
;  spawn,'ps2pdf -dPDFSETTINGS=/prepress -sEPSCrop '+ffn+'.eps ' +ffn+'.pdf'

end