pro make_gradehist_nov14o4,pid=pid,maindir=maindir,noplot=noplot,fid=fid

  ;  Make, save out and plot (default) histogram of event grades for
  ;  whole FoV and each individual regions within the P1 and P2 pointings from Nov 2014
  ;
  ;  Start by filtering out the "bad" pixels in FPMA and only CHU23 (P1,P2) (P3 CHU23, P4 CHU13)
  ;  Shift P1 so matches P2 pointing
  ;
  ;         Optional inputs:
  ;         pid       - 1=P1 (default) or 2=P2 ?
  ;         fid       - A (default) or B?
  ;         maindir   - Where Nov data is kept - maindir of the ftp structed dirs
  ;         noplot    - Keyword not to produce plot of histograms
  ;                   - Need >8.0 to do the plots (was done on v8.5)
  ;
  ;  Assumes that make_roi_nov14.pro has already been run, so roi_nov14.dat exists in same directory
  ;
  ; 05-11-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (n_elements(pid) ne 1) then pid=1
  pname='P'+string(pid,format='(i1)')

  if (n_elements(fid) ne 1) then fid='A'

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Where I've put the sunpos corrected nupipeline outoput file
  if (n_elements(maindir) ne 1) then maindir='~/data/ns_data/obs2_bg2/

  ; Names of the subdirectories per pointing of orbit 4
  ; ONly really need P1 and P2 here but might as well keep in the P
  sdirs=['20012001_Sol_14305_AR2192_1/','20012002_Sol_14305_AR2192_2/',$
    '20012003_Sol_14305_AR2192_3/','20012004_Sol_14305_AR2192_4/']
  subdir=sdirs[pid-1]
  dnms=['20012001001','20012002001', '20012003001','20012004002']
  ddname=dnms[pid-1]
  chm=[13,13,13,10]
  chmn=['CHU23','CHU23','CHU23','CHU13']
  chumask=chm[pid-1]
  chunam=chmn[pid-1]

  xcs=[690,710,710,605]
  ycs=[-19,-19,-19,46]
  xc=xcs[1]; using P2 as refernce map centre
  yc=ycs[1]; using P2 as refernce map centre

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Get the *_cl_sunpost.evt file

  cl_file=maindir+subdir+ddname+'/event_cl/nu'+ddname+fid+'06_cl_sunpos.evt'
  evt = mrdfits(cl_file, 1,evth)

  ; Before doing anything else need to filter out the "bad" pixels in FPMA
  ; these were the ones BG had previously identified - caused the "hard knots" in the data
  ; https://github.com/NuSTAR/nustar_solar/blob/master/solar_mosaic_20141211/combine_events.pro

  if (fid eq 'A') then begin

    use = bytarr(n_elements(evt)) + 1
    thisdet = where(evt.det_id eq 2)
    badones = where(evt[thisdet].rawx eq 16 and evt[thisdet].rawy eq 5, nbad)
    if nbad gt 0 then use[thisdet[badones]]=0
    badones = where(evt[thisdet].rawx eq 24 and evt[thisdet].rawy eq 22, nbad)
    if nbad gt 0 then use[thisdet[badones]]=0

    thisdet = where(evt.det_id eq 3)
    badones = where(evt[thisdet].rawx eq 22 and evt[thisdet].rawy eq 1, nbad)
    if nbad gt 0 then use[thisdet[badones]]=0
    badones = where(evt[thisdet].rawx eq 15 and evt[thisdet].rawy eq 3, nbad)
    if nbad gt 0 then use[thisdet[badones]]=0
    badones = where(evt[thisdet].rawx eq 0 and evt[thisdet].rawy eq 15, nbad)
    if nbad gt 0 then use[thisdet[badones]]=0

    evt=evt[where(use)]
  endif


  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Doing it in just one CHU so can filter out other CHU as well
  chufile = file_search(maindir+subdir+ddname+'/hk/', '*chu123.fits')
  for chunum= 1, 3 do begin
    chu = mrdfits(chufile, chunum)
    maxres = 20 ;; [arcsec] maximum solution residual
    qind=1 ; From KKM code...
    if chunum eq 1 then begin
      mask = (chu.valid EQ 1 AND $          ;; Valid solution from CHU
        chu.residual LT maxres AND $  ;; CHU solution has low residuals
        chu.starsfail LT chu.objects AND $ ;; Tracking enough objects
        chu.(qind)(3) NE 1)*chunum^2       ;; Not the "default" solution
    endif else begin
      mask += (chu.valid EQ 1 AND $            ;; Valid solution from CHU
        chu.residual LT maxres AND $    ;; CHU solution has low residuals
        chu.starsfail LT chu.objects AND $ ;; Tracking enough objects
        chu.(qind)(3) NE 1)*chunum^2       ;; Not the "default" solution
    endelse
  endfor

  ; make time binning of chus to evt data
  chu_comb = round(interpol(mask, chu.time, evt.time))

  ; filter out bad CHUs
  ida2=where(chu_comb eq chumask)
  evt=evt[ida2]
  ;  engs=1.6+0.04*evt.pi

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Setup the pixel and binning sizes
  ; Get the same values if using evtah or evtbh
  ttype = where(stregex(evth, "TTYPE", /boolean))
  xt = where(stregex(evth[ttype], 'X', /boolean))
  xpos = (strsplit( (strsplit(evth[ttype[max(xt)]], ' ', /extract))[0], 'E', /extract))[1]
  npix = sxpar(evth, 'TLMAX'+xpos)
  pix_size = abs(sxpar(evth,'TCDLT'+xpos))

  centerx = round(xc / pix_size) + npix * 0.5
  centery = round(yc / pix_size) + npix * 0.5
  im_size = 1037. / pix_size
  im_width = round(im_size * 2.)
  
  ; grade binning
  ; 0 to 32
  ngrd=33

  ; Shifts relative to P2 pointing position - as that matches AIA
  ; Rememeber in this one the X is still opposite direction to S/C X (so need to *-1)
  xshf=-1*round((xcs[1]-xcs[pid-1])/pix_size)
  yshf=round((ycs[1]-ycs[pid-1])/pix_size)
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Do the histograms of the grades per pixel
  ims=intarr(ngrd,1688,1688)

  evtx=evt.x-xshf
  evty=evt.y-yshf

  for i=0, ngrd-1 do begin
    ide=where(evt.grade eq i,nid)
  ;  print,i
    if (nid gt 1) then begin
      ; this data still  has the x opposite direction to standard solar coords
      pixinds = (npix - evtx[ide]) + evty[ide] * npix
      im_hist = histogram(pixinds, min = 0, max = npix*npix-1, binsize = 1)
      im = reform(im_hist, npix, npix)
      im= im[(centerx-im_width):(centerx+im_width-1), (centery-im_width):(centery+im_width-1)]
      ims[i,*,*] = im
    endif
    if (nid eq 1) then begin
      pixinds = (npix - evtx[ide]) + evty[ide] * npix
      im_hist = lonarr(npix*npix)
      im_hist[pixinds]=1
      im = reform(im_hist, npix, npix)
      im= im[(centerx-im_width):(centerx+im_width-1), (centery-im_width):(centery+im_width-1)]
      ims[i,*,*] = im
    endif
  endfor

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Can take an even smaller region which covers all 4 pointings
  npp=n_elements(ims[0,0,*])
  subx=[300,1200]
  suby=[550,1450]

  ims=ims[*,subx[0]:subx[1],suby[0]:suby[1]]

  pxs=pix_size
  x0=xc-npp*0.5*pxs+pxs*subx[0]
  y0=yc-npp*0.5*pxs+pxs*suby[0]

  newxc=x0+0.5*n_elements(ims[0,*,0])*pxs
  newyc=y0+0.5*n_elements(ims[0,0,*])*pxs

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Histogram the grades
  grades=indgen(ngrd)
  grd_hist_all=reform(total(total(ims,2),2))
  g2124_g0_all=total(ims[21:24,*,*])/(total(ims[0,*,*])*1.0)

  ; loads in the roi, xcsp, ycsp,rnm,wid
  restore,file='out_files/roi_no14.dat'
  nr=n_elements(rnm)
  grd_hist_reg=fltarr(nr,ngrd)
  g2124_g0_reg=fltarr(nr)

  ; Convert region coords from S/C arcsec to pixels in ims
  xcs=round((xcsp-x0)/pxs)
  ycs=round((ycsp-y0)/pxs)
  iwid=round(wid/pxs)

  ; do all for P1 but only D1 and L1 for P2
  if (pid eq 1) then ids=indgen(5)
  if (pid eq 2) then ids=[0,2]

  nids=n_elements(ids)

  ; Histograms per each region
  print,pname+' FPM'+fid
  print, 'FoV',total(ims),total(ims[0,*,*]),total(ims[21:24,*,*]),' ',$
    string(total(ims[21:24,*,*])/total(ims[0,*,*]),format='(e9.3)')
  for i=0, nids-1 do begin
    r=ids[i]
    xr=xcs[r]+0.5*[-iwid,iwid]
    yr=ycs[r]+0.5*[-iwid,iwid]
    gtemp=ims[*,xr[0]:xr[1],yr[0]:yr[1]]
    grd_hist_reg[r,*]=reform(total(total(gtemp,2),2))  
    g2124_g0_reg[r]=total(gtemp[21:24,*,*])/total(gtemp[0,*,*])
    
    print,rnm[r],total(gtemp),total(gtemp[0,*,*]),total(gtemp[21:24,*,*]),' ',$
      string(total(gtemp[21:24,*,*])/total(gtemp[0,*,*]),format='(e9.3)')
    
  endfor
  ; Save out the information
  
  grd_out={pname:pname,fpm:'FPM'+fid,$
    grd_hist_all:grd_hist_all,g2124_g0_all:g2124_g0_all,$
    reg_nm:rnm,grd_hist_reg:grd_hist_reg,g2124_g0_reg:g2124_g0_reg}
    
  save,file='out_files/grds_'+pname+'FPM'+fid+'.dat',grd_out

  ; plot it out if have new enough idl and no /noplot
  if (keyword_set(noplot) ne 1 and float(!version.release) ge 8.0) then begin
    w=window(dimensions=[900,450],/buffer)
    ;    p0=plot(grades,grd_hist_all/total(grd_hist_all),/histogram,/ylog,title='FoV',$
    ;      yrange=[1e-5,2],xrange=[-0.5,32.5],thick=2,xtitle='Grade',ytitle='Frequency',layout=[3,2,1],/current)
    ;for i=0, nr-1 do !null=plot(grades,grd_hist_reg[i,*]/total(grd_hist_reg[i,*]),/histogram,/ylog,title=rnm[i],$
    ;  yrange=[1e-5,2],xrange=[-0.5,32.5],thick=2,xtitle='Grade',ytitle='Frequency',layout=[3,2,2+i],/current)
    mrg=0.2
    p0=barplot(grades,grd_hist_all,/ylog,axis_style=1,$
      title='FoV, $G_{21-24}/G_0$='+string(g2124_g0_all*100.,format='(f6.3)')+'%',$
      yrange=[0.5,1e5],xrange=[-3,35],xtitle='Grade',ytitle='Number',layout=[3,2,1],/current,margin=mrg)
    yaxis=axis('Y',location='right',title='Freq',coord_transform=[0,total(grd_hist_all)],$
      target=p0,tickunits='Scientific',color='grey')
    for i=0, nr-1 do begin
      zz=barplot(grades,grd_hist_reg[i,*],/ylog,axis_style=1,margin=mrg,$
        title=rnm[i]+', $G_{21-24}/G_0$='+string(g2124_g0_reg[i]*100.,format='(f6.3)')+'%',$
        yrange=[0.5,1e5],xrange=[-3,35],xtitle='Grade',ytitle='Number',layout=[3,2,2+i],/current)
      if (total(grd_hist_reg[i,*]) gt 0.0) then $
        yaxis=axis('Y',location='right',title='Freq',coord_transform=[0,1.0/total(grd_hist_reg[i,*])],$
        target=zz,tickunits='Scientific',color='grey')

    endfor
    tt=text(0,0,pname+' FPM'+fid,/device,target=w)

    w.save,'figs/grade_hist_'+pname+'FPM'+fid+'.pdf',page_size=w.dimensions/100.
    w.close

    ;    stop
  endif

end