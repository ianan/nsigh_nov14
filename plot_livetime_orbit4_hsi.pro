pro plot_livetime_orbit4_hsi,maindir=maindir

  ; Produce a nice time profile lightcurve of the NuSTAR livetime, RHESSI and GOES fluxes
  ; 09-Nov-2015
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Where is the NuSTAR data located
  if (n_elements(maindir) ne 1) then maindir='~/data/ns_data/obs2_bg2/

  t1='01-Nov-2014 21:30'
  t2='01-Nov-2014 22:35'

  ; Names of the subdirectories per pointing of orbit 4
  sdirs=['20012001_Sol_14305_AR2192_1/','20012002_Sol_14305_AR2192_2/',$
    '20012003_Sol_14305_AR2192_3/','20012004_Sol_14305_AR2192_4/']
  dnms=['20012001001','20012002001', '20012003001','20012004002']

  ddname=dnms[0]
  subdir=sdirs[0]
  hka_file=maindir+subdir+ddname+'/hk/nu'+ddname+'A_fpm.hk'
  hka1 = mrdfits(hka_file, 1, hkahdr)
  gd=where(hka1.livetime lt 4.5e-3 and hka1.livetime gt 0.)
  hka1=hka1[gd]
  hkatims1=anytim(hka1.time+anytim('01-Jan-2010'),/yoh,/trunc)
  ddname=dnms[1]
  subdir=sdirs[1]
  hka_file=maindir+subdir+ddname+'/hk/nu'+ddname+'A_fpm.hk'
  hka2 = mrdfits(hka_file, 1, hkahdr)
  gd=where(hka2.livetime lt 5e-2 and hka2.livetime gt 0.)
  hka2=hka2[gd]
  hkatims2=anytim(hka2.time+anytim('01-Jan-2010'),/yoh,/trunc)
  ddname=dnms[2]
  subdir=sdirs[2]
  hka_file=maindir+subdir+ddname+'/hk/nu'+ddname+'A_fpm.hk'
  hka3 = mrdfits(hka_file, 1, hkahdr)
  gd=where(hka3.livetime lt 5e-2 and hka3.livetime gt 0.)
  hka3=hka3[gd]
  hkatims3=anytim(hka3.time+anytim('01-Jan-2010'),/yoh,/trunc)
  ddname=dnms[3]
  subdir=sdirs[3]
  hka_file=maindir+subdir+ddname+'/hk/nu'+ddname+'A_fpm.hk'
  hka4 = mrdfits(hka_file, 1, hkahdr)
  gd=where(hka4.livetime lt 5e-2 and hka4.livetime gt 0.)
  hka4=hka4[gd]
  hkatims4=anytim(hka4.time+anytim('01-Jan-2010'),/yoh,/trunc)

  mint=strarr(4)
  mint[0]=anytim(mean(anytim(hkatims1)),/yoh)
  mint[1]=anytim(mean(anytim(hkatims2)),/yoh)
  mint[2]=anytim(mean(anytim(hkatims3)),/yoh)
  mint[3]=anytim(mean(anytim(hkatims4)),/yoh)


  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Setup the postscript plotting
  @post_outset
  !p.multi=[0,1,3]
  !p.charsize=1.7

  set_plot,'ps'
  device, /encapsulated, /color, /isolatin1,/inches, $
    bits=8, xsize=5, ysize=5,file='figs/Orb4_Livetime_hsi.eps'

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Plot the NuSTAR livetim
  utplot,hkatims1,hka1.livetime,/ylog,yrange=[1e-1,9],/nodata,$
    ytitle='NuSTAR Livetime %',position=[0.15,0.67,0.95,0.97],xtit='',xtickf='(a1)',timer=[t1,t2]
  tube_line_colors
  outplot,hkatims1,hka1.livetime*100,color=5,thick=5
  outplot,hkatims2,hka2.livetime*100,color=4,thick=5
  outplot,hkatims3,hka3.livetime*100,color=8,thick=5
  outplot,hkatims4,hka4.livetime*100,color=12,thick=5

  evt_grid,mint,labpos=0.13,labels='P'+['1','2','3','4'],labcolor=[5,4,8,12],/noarrow,/noline,/data,labsize=1.1

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Get the RHESSI data if not already got it
  if (file_test('out_files/hsi_data_01112014.dat') eq 0) then begin
    obj = hsi_spectrum()
    obj-> set, obs_time_interval= [' 1-Nov-2014 21:30:00.000', ' 1-Nov-2014 23:00:00.000']
    obj-> set, seg_index_mask= [1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    obj-> set, sp_data_unit= 'Rate'
    obj-> set, sp_energy_binning= [4, 8]
    obj-> set, sp_time_interval= 4
    obj-> set, sum_flag= 1
    rdata = reform(obj->getdata())
    rtims=anytim(anytim('1-Nov-2014 21:30:00')+4*findgen(n_elements(rdata)),/yoh,/trunc)
    save,file='out_files/hsi_data_01112014.dat',rtims,rdata
  endif
  restore,file='out_files/hsi_data_01112014.dat'

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Get the RHESSI obs summary if not already got it
  if (file_test('out_files/hsi_flags_01112014.dat') eq 0) then begin
    obj =hsi_obs_summary(obs_time= [' 1-Nov-2014 21:30:00.000', ' 1-Nov-2014 23:00:00.000'])
    rflags = obj -> getdata(class='flag')
    rftims = anytim(obj -> getaxis(/ut, class='flag'),/yoh,/trunc)
    rinfo = obj -> get(/info, class='flag')
    rflag_changes = obj -> changes()
    save,file='out_files/hsi_flags_01112014.dat',rflags,rftims,rinfo,rflag_changes
  endif
  restore,file='out_files/hsi_flags_01112014.dat'

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Plot the RHESSI data and flags
  saas=rflag_changes.saa_flag.start_times[0]
  saae=rflag_changes.saa_flag.end_times[0]

  eclps=rflag_changes.eclipse_flag.start_times[1]
  eclpe=rflag_changes.eclipse_flag.end_times[1]

  utplot,rtims,rdata,yrange=[0.01,49],$
    ytitle='RHESSI Rate D168',position=[0.15,0.4,0.95,0.65],xtit='',xtickf='(a1)',$
    timer=[t1,t2],thick=4,psym=10,/nodata
  outplot,rtims,rdata,thick=4,psym=10,color=2
  evt_grid,'01-Nov-14 21:34',labpos=[5],labels=['4-8 keV'],/data,$
    labsize=1.,/labonly,labcolor=2,align=0,/noarrow

  outplot,anytim([saas,saae],/yoh), [42,42],color=10,thick=3,lines=2
  outplot,anytim([eclps,eclpe],/yoh),[40,40],color=10,thick=3
  evt_grid,'01-Nov-14 21:50',labpos=[43],labels=['SAA'],/data,$
    labsize=0.8,/labonly,labcolor=10,align=1,/noarrow
  evt_grid,'01-Nov-14 22:12',labpos=[41],labels=['NIGHT'],/data,$
    labsize=0.8,/labonly,labcolor=10,align=1,/noarrow

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Get the GOES data if not already got it
  if (file_test('out_files/goes_data_01112014.dat') eq 0) then begin
      a = ogoes()
      a->set,tstart=anytim(t1,/yoh),tend=anytim(t2,/yoh),/three
      glow=a->getdata(/low)
      ghigh=a->getdata(/high)
      gtim = a->getdata(/times)
      gutbase = a->get(/utbase)
      gtime=anytim(anytim(gutbase)+gtim,/yoh,/trunc)
      save,file='out_files/goes_data_01112014.dat',glow,gtime,ghigh
  endif
  restore,file='out_files/goes_data_01112014.dat'

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Plot the GOES data
  utplot,gtime,glow*1d7,ytitle='!3GOES [x10!U-7!N Wm!U-2!N]',$
    /nodata,yrange=[3.5,4.4],timer=[t1,t2],position=[0.15,0.1,0.95,0.38]
  ;    outplot,gtime,ghigh*1d9,color=2,thick=4
  outplot,gtime,glow*1d7,color=1,thick=4
  evt_grid,'01-Nov-14 21:34',labpos=[3.6],labels=['1-8'+string(197b)],/data,$
    labsize=1.,/labonly,labcolor=1,align=0,/noarrow

  device,/close
  set_plot, mydevice

end
