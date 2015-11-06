pro make_specs_nov14o4_de,pid=pid,maindir=maindir,fid=fid,de=de,noplot=noplot

  ;  Make the spectra (energy,x,y) for orbit 4 data in P1 or P2 over the whole FoV and
  ;  each of the chosen regions using an energy binning of de between 1.6 and 10 keV
  ;
  ;  Starts by filtering out out "bad" pixels in FPMA and only the one CHU time CHU23 (P1,P2) and then
  ;  only those in Grade 0 (to make sure have only single pixel hit events). We also get the livetime
  ;  correction and shift pointing so all match P2 before saving it all out
  ;
  ;         Optional inputs:
  ;         pid       - 1=P1 (default) or 2=P2 ?
  ;         fid       - A (default) or B?
  ;         de        - Energy binning of spectrum (default 0.2keV)
  ;         maindir   - Where Nov data is kept - maindir of the ftp structed dirs
  ;         noplot    - Keyword not to produce plot of histograms
  ;                   - Need >8.0 to do the plots (was done on v8.5)
  ;
  ; Assumes that make_roi_nov14.pro has already been run, so roi_nov14.dat exists in out_files/
  ;
  ; 06-11-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (n_elements(pid) ne 1) then pid=1
  pname='P'+string(pid,format='(i1)')

  if (n_elements(de) ne 1) then de=0.2
  denm='DE'+strmid(string(1000+de*100,format='(i4)'),1,3)

  if (n_elements(fid) ne 1) then fid='A'

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Where I've put the sunpos corrected nupipeline outoput file
  if (n_elements(maindir) ne 1) then maindir='~/data/ns_data/obs2_bg2/

  ; Names of the subdirectories per pointing of orbit 4
  sdirs=['20012001_Sol_14305_AR2192_1/','20012002_Sol_14305_AR2192_2/',$
    '20012003_Sol_14305_AR2192_3/','20012004_Sol_14305_AR2192_4/']
  subdir=sdirs[pid-1]
  dnms=['20012001001','20012002001', '20012003001','20012004002']
  ddname=dnms[pid-1]
  chm=[13,13,13,10]
  chmn=['CHU23','CHU23','CHU23','CHU13']
  chumask=chm[pid-1]
  chunam=chmn[pid-1]

  ; These were manally found
  ; The original P2 one from BG
  ; Assume P3 the same as same CHU as P2 and not obvious features for matching to AIA
  ; P1 from IGH found in comparison to AIA 94/335
  ; P4 from IGH found in comparison to AIA 94/335 but from previous orbits (as no features in this orbit)
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

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Filter out just the grade 0 events
  evt=evt[where(evt.grade eq 0)]

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Doing it in just one CHU so can filter out other CHU as well
  ; Based on https://github.com/NuSTAR/nustar_solar/blob/master/solar_mosaic_20141211/solar_mosaic_hk.pro
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
  chu_filt=where(chu_comb eq chumask,nchu)
  evt=evt[chu_filt]

  ; What fraction of the "good" grade 0 events were in our chosen CHU combination ?
  inchu_frac=nchu/(n_elements(chu_comb)*1.0)

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Get the livetime info
  hk_file=maindir+subdir+ddname+'/hk/nu'+ddname+fid+'_fpm.hk'
  hk = mrdfits(hk_file, 1, hkhdr)
  ; NuSTAR times are seconds from 01-Jan-2010 so convert to which anytim format
  hktims=anytim(hk.time+anytim('01-Jan-2010'))

  ; Assuming duration of observation is continous (etc no major jumps in CHU)
  t1=anytim(min(evt.time)+anytim('01-Jan-2010'),/yoh)
  t2=anytim(max(evt.time)+anytim('01-Jan-2010'),/yoh)
  dur=max(evt.time)-min(evt.time)

  lvid=where(hktims ge anytim(t1) and hktims lt anytim(t2))
  lvtcor=mean(hk[lvid].livetime)

  ; Can also get this information from the header file BUT again assumes continuous and gaps
  ; Exposure=livetime=duration*livetime fraction
  exploc=where(stregex(evth, "EXPOSURE",/boolean))
  expos=abs((strsplit(evth[exploc],/extract))[1])
  expos_dl=dur*lvtcor

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

  ; energy binning for output spectrum/datacube
  maxe=10.
  mine=1.6
  eid=mine+de*indgen(2+(maxe-mine)/de)
  mide=get_edges(eid,/mean)
  nengs=n_elements(mide)

  ims=intarr(nengs,1688,1688)

  ; Shifts relative to P2 pointing position - as that matches AIA
  ; Rememeber in this one the X is still opposite direction to S/C X (so need to *-1)
  xshf=-1*round((xcs[1]-xcs[pid-1])/pix_size)
  yshf=round((ycs[1]-ycs[pid-1])/pix_size)
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Create the data cube
  evtengs=1.6+0.04*evt.pi
  evtx=evt.x-xshf
  evty=evt.y-yshf

  for i=0, nengs-1 do begin
    ide=where(evtengs ge eid[i] and evtengs lt eid[i+1],nid)
    print,eid[i],' --- ',eid[i+1]
    if (nid gt 1) then begin
      ; this data still  has the x opposite direction to standard solar coords
      pixinds = (npix - evtx[ide]) + evty[ide] * npix
      im_hist = histogram(pixinds, min = 0, max = npix*npix-1, binsize = 1)
      im = reform(im_hist, npix, npix)
      im= im[(centerx-im_width):(centerx+im_width-1), (centery-im_width):(centery+im_width-1)]
      ims[i,*,*] = im
    endif
    if (nid eq 1) then begin
      ; this data still  has the x opposite direction to standard solar coords
      pixinds = (npix - evtx[ide]) + evty[ide] * npix
      im_hist = lonarr(npix*npix)
      im_hist[pixinds]=1
      im = reform(im_hist, npix, npix)
      im= im[(centerx-im_width):(centerx+im_width-1), (centery-im_width):(centery+im_width-1)]
      ims[i,*,*] = im
    endif
  endfor


  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Save all this out
  npp=n_elements(ims[0,0,*])
  ; Can take an even smaller region which covers all 4 pointings
  ; file output should then "only" be about 0.6GB not 2.6GB!
  subx=[300,1200]
  suby=[550,1450]
  ims=ims[*,subx[0]:subx[1],suby[0]:suby[1]]

  pxs=pix_size
  x0=xc-npp*0.5*pxs+pxs*subx[0]
  y0=yc-npp*0.5*pxs+pxs*suby[0]

  newxc=x0+0.5*n_elements(ims[0,*,0])*pxs
  newyc=y0+0.5*n_elements(ims[0,0,*])*pxs

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Make spectra out of the data cube

  spec_all=reform(total(total(ims,2),2))

  ; loads in the roi, xcsp, ycsp,rnm,wid
  restore,file='out_files/roi_no14.dat'
  nr=n_elements(rnm)
  spec_reg=fltarr(nr,nengs)

  ; Convert region coords from S/C arcsec to pixels in ims
  xcs=round((xcsp-x0)/pxs)
  ycs=round((ycsp-y0)/pxs)
  iwid=round(wid/pxs)

  ; do all for P1 but only D1 and L1 for P2
  if (pid eq 1) then ids=indgen(5)
  if (pid eq 2) then ids=[0,2]

  nids=n_elements(ids)

  ; Spectrum per each region
  print,pname+' FPM'+fid
  for i=0, nids-1 do begin
    r=ids[i]
    xr=xcs[r]+0.5*[-iwid,iwid]
    yr=ycs[r]+0.5*[-iwid,iwid]
    spec_reg[r,*]=reform(total(total(ims[*,xr[0]:xr[1],yr[0]:yr[1]],2),2))

  endfor
  ; Save out the information

  spc_out={pname:pname,fpm:'FPM'+fid,chu:chunam,inchu_frac:inchu_frac,$
    t1:t1,t2:t2,dur:dur,lvt_frac:lvtcor,expos:expos,expos_dl:expos_dl,$
    engs:mide,eng_edges:eid,spec_all:spec_all,$
    reg_nm:rnm,spec_reg:spec_reg}

  save,file='out_files/specs_'+denm+'_'+pname+'FPM'+fid+'.dat',spc_out

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Plot it out if have new enough idl and no /noplot
  if (keyword_set(noplot) ne 1 and float(!version.release) ge 8.0) then begin
    w=window(dimensions=[900,600],/buffer)
    mrg=0.2
    p0=plot(mide,spec_all,/ylog,title='FoV',/stair, margin=mrg,$
      yrange=[0.5,1e5],xrange=[1.0,10.0],xtitle='Energy [keV]',ytitle='NuSTAR Counts',layout=[3,2,1],/current)
    for i=0, nr-1 do zz=plot(mide,spec_reg[i,*],/ylog,/stair,margin=mrg, title=rnm[i],$
      yrange=[0.5,1e5],xrange=[1.0,10.0],xtitle='Energy [keV]',ytitle='NuSTAR Counts',layout=[3,2,2+i],/current)
    tt=text(0,2,pname+' FPM'+fid+' '+$
      string(expos_dl,format='(f6.2)')+'s ('+string(dur,format='(f6.1)')+'s, '+$
      string(lvtcor*100,format='(f5.2)')+'%)',/device,target=w)

    w.save,'figs/specs_'+pname+'FPM'+fid+'.pdf',page_size=w.dimensions/100.
    w.close
  endif
  ;stop
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


end