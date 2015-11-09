pro make_map_nov14o4,pid=pid,maindir=maindir,$
  grdid=grdid,erang=erang,plot=plot

  ;  Make the map in a given energy range for orbit 4 data in FPMA and FPMB
  ;         Filter out "bad" pixels in FPMA
  ;         Filter out so only CHU23 (P1,P2,P3) or CHU13 (P4)
  ;         Filter out events other than GRADE 0 - should minimise pileup
  ;         Get the livetime correction per FPM
  ;         Shift pointing so all match P2
  ;         Make the map
  ;         Plot the map
  ;
  ;         Optional inputs:
  ;         pid     1,2,3,4 for P1,P2,P3,P4
  ;         grdid   0,1,2 for Event grade 0, all or 21-24
  ;         erang   Energy range covered by the map, if only 1 elements then >erang
  ;                  (default is >2 keV)
  ;         plot    Want to plot the maps? (default no)
  ;
  ;         For non-IGH use need to change
  ;         maindir - where Nov data is kept - maindir of the ftp structed dirs
  ;
  ;
  ; 09-Nov-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (n_elements(pid) ne 1) then pid=1
  pname='P'+string(pid,format='(i1)')
  if (n_elements(grdid) ne 1) then grdid=0
  grdname='G'+string(grdid,format='(i1)')
  ; Only want single pixel hit: grade 0
  if (grdid eq 0) then grdrn=[0,0]
  ; Want all the grades: grade =>0
  if (grdid eq 1) then grdrn=[0,32]
  ; Want just the corner second pixel hits: grade 21 to 24
  if (grdid eq 2) then grdrn=[21,24]

  if (n_elements(erang) lt 1 or n_elements(erang) gt 2) then erang =2
  if (n_elements(erang) eq 1) then begin
    emin=erang
    emax=100.
    eid='>'+strcompress(string(erang,format='(i2)'),/rem)+' keV'
    enme='EG'+strcompress(string(erang,format='(i2)'),/rem)
  endif
  if (n_elements(erang) eq 2) then begin
    emin=erang[0]
    emax=erang[1]
    eid=strcompress(string(erang[0],format='(i2)'),/rem)+'-'+$
      strcompress(string(erang[1],format='(i2)'),/rem)+' keV'
    enme='E'+strcompress(string(erang[0],format='(i2)'),/rem)+'_'+$
      strcompress(string(erang[1],format='(i2)'),/rem)
  endif

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
  ; Get the *_cl_sunpost.evt files

  cla_file=maindir+subdir+ddname+'/event_cl/nu'+ddname+'A06_cl_sunpos.evt'
  evta = mrdfits(cla_file, 1,evtah)

  clb_file=maindir+subdir+ddname+'/event_cl/nu'+ddname+'B06_cl_sunpos.evt'
  evtb = mrdfits(clb_file, 1,evtbh)

  ; Before doing anything else need to filter out the "bad" pixels in FPMA
  ; these were the ones BG had previously identified - caused the "hard knots" in the data
  ; https://github.com/NuSTAR/nustar_solar/blob/master/solar_mosaic_20141211/combine_events.pro

  use = bytarr(n_elements(evta)) + 1
  thisdet = where(evta.det_id eq 2)
  badones = where(evta[thisdet].rawx eq 16 and evta[thisdet].rawy eq 5, nbad)
  if nbad gt 0 then use[thisdet[badones]]=0
  badones = where(evta[thisdet].rawx eq 24 and evta[thisdet].rawy eq 22, nbad)
  if nbad gt 0 then use[thisdet[badones]]=0

  thisdet = where(evta.det_id eq 3)
  badones = where(evta[thisdet].rawx eq 22 and evta[thisdet].rawy eq 1, nbad)
  if nbad gt 0 then use[thisdet[badones]]=0
  badones = where(evta[thisdet].rawx eq 15 and evta[thisdet].rawy eq 3, nbad)
  if nbad gt 0 then use[thisdet[badones]]=0
  badones = where(evta[thisdet].rawx eq 0 and evta[thisdet].rawy eq 15, nbad)
  if nbad gt 0 then use[thisdet[badones]]=0

  evta=evta[where(use)]

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
  achu_comb = round(interpol(mask, chu.time, evta.time))
  bchu_comb = round(interpol(mask, chu.time, evtb.time))

  ; filter out bad CHUs and the requested grades
  ida2=where(achu_comb eq chumask and evta.grade ge grdrn[0] and evta.grade le grdrn[1])
  evta=evta[ida2]
  a_engs=1.6+0.04*evta.pi

  idb2=where(bchu_comb eq chumask and evtb.grade ge grdrn[0] and evtb.grade le grdrn[1])
  evtb=evtb[idb2]
  b_engs=1.6+0.04*evtb.pi

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Only want the counts in energy range
  ida2=where(a_engs ge emin and a_engs lt emax)
  evta=evta[ida2]
  a_engs=1.6+0.04*evta.pi

  idb2=where(b_engs ge emin and b_engs lt emax)
  evtb=evtb[idb2]
  b_engs=1.6+0.04*evtb.pi

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Get the livetime info
  hka_file=maindir+subdir+ddname+'/hk/nu'+ddname+'A_fpm.hk'
  hka = mrdfits(hka_file, 1, hkahdr)
  hkatims=anytim(hka.time+anytim('01-Jan-2010'))

  t1a=anytim(min(evta.time)+anytim('01-Jan-2010'),/yoh,/trunc)
  t2a=anytim(max(evta.time)+anytim('01-Jan-2010'),/yoh,/trunc)

  lvida=where(hkatims ge anytim(t1a) and hkatims lt anytim(t2a))
  lvtcora=mean(hka[lvida].livetime)

  hkb_file=maindir+subdir+ddname+'/hk/nu'+ddname+'B_fpm.hk'
  hkb = mrdfits(hkb_file, 1, hkbhdr)
  hkbtims=anytim(hkb.time+anytim('01-Jan-2010'))

  t1b=anytim(min(evtb.time)+anytim('01-Jan-2010'),/yoh,/trunc)
  t2b=anytim(max(evtb.time)+anytim('01-Jan-2010'),/yoh,/trunc)

  lvidb=where(hkbtims ge anytim(t1b) and hkbtims lt anytim(t2b))
  lvtcorb=mean(hkb[lvidb].livetime)

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Setup the pixel and binning sizes
  ; Get the same values if using evtah or evtbh
  ttype = where(stregex(evtah, "TTYPE", /boolean))
  xt = where(stregex(evtah[ttype], 'X', /boolean))
  xpos = (strsplit( (strsplit(evtah[ttype[max(xt)]], ' ', /extract))[0], 'E', /extract))[1]
  npix = sxpar(evtah, 'TLMAX'+xpos)
  pix_size = abs(sxpar(evtah,'TCDLT'+xpos))

  centerx = round(xc / pix_size) + npix * 0.5
  centery = round(yc / pix_size) + npix * 0.5
  im_size = 1037. / pix_size
  im_width = round(im_size * 2.)
  ims=intarr(1688,1688)

  ; Shifts relative to P2 pointing position - as that matches AIA
  ; Rememeber in this one the X is still opposite direction to S/C X (so need to *-1)
  xshf=-1*round((xcs[1]-xcs[pid-1])/pix_size)
  yshf=round((ycs[1]-ycs[pid-1])/pix_size)
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Do FPMA
  engs=a_engs
  evtx=evta.x-xshf
  evty=evta.y-yshf

  ; this data still  has the x opposite direction to standard solar coords
  pixinds = (npix - evtx) + evty * npix
  im_hist = histogram(pixinds, min = 0, max = npix*npix-1, binsize = 1)
  im = reform(im_hist, npix, npix)
  ims= im[(centerx-im_width):(centerx+im_width-1), (centery-im_width):(centery+im_width-1)]

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  npp=n_elements(ims[0,*])
  ; Can take an even smaller region which covers all 4 pointings
  ; file output should then "only" be about 0.6GB not 2.6GB!
  subx=[300,1200]
  suby=[550,1450]
  ims=ims[subx[0]:subx[1],suby[0]:suby[1]]

  pxs=pix_size
  x0=xc-npp*0.5*pxs+pxs*subx[0]
  y0=yc-npp*0.5*pxs+pxs*suby[0]

  newxc=x0+0.5*n_elements(ims[*,0])*pxs
  newyc=y0+0.5*n_elements(ims[0,*])*pxs

  dur=anytim(t2a)-anytim(t1a)
  time=t1a
  ang = pb0r(t1a,/arcsec,l0=l0)

  ima2_lvt=ims/(float(lvtcora)*dur)

  mapa=make_map(ima2_lvt,dx=pxs,dy=pxs,xc=newxc,yc=newyc,$
    time=time,dur=dur,id='LC FPMA '+grdname+' '+eid+string(lvtcora*100,format='(f6.2)')+'%',$
    l0=l0,b0=ang[1],rsun=ang[2])

  map2fits,mapa,'out_files/maps_'+pname+grdname+'_'+enme+'_FPMA.fits'

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Do FPMB
  engs=b_engs
  evtx=evtb.x-xshf
  evty=evtb.y-yshf

  ; this data still  has the x opposite direction to standard solar coords
  pixinds = (npix - evtx) + evty * npix
  im_hist = histogram(pixinds, min = 0, max = npix*npix-1, binsize = 1)
  im = reform(im_hist, npix, npix)
  ims= im[(centerx-im_width):(centerx+im_width-1), (centery-im_width):(centery+im_width-1)]

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  npp=n_elements(ims[0,*])
  ; Can take an even smaller region which covers all 4 pointings
  ; file output should then "only" be about 0.6GB not 2.6GB!
  subx=[300,1200]
  suby=[550,1450]
  ims=ims[subx[0]:subx[1],suby[0]:suby[1]]

  pxs=pix_size
  x0=xc-npp*0.5*pxs+pxs*subx[0]
  y0=yc-npp*0.5*pxs+pxs*suby[0]

  newxc=x0+0.5*n_elements(ims[*,0])*pxs
  newyc=y0+0.5*n_elements(ims[0,*])*pxs

  dur=anytim(t2b)-anytim(t1b)
  time=t1b
  ang = pb0r(t1b,/arcsec,l0=l0)

  imb2_lvt=ims/(float(lvtcorb)*dur)

  mapb=make_map(imb2_lvt,dx=pxs,dy=pxs,xc=newxc,yc=newyc,$
    time=time,dur=dur,id='LC FPMB '+grdname+' '+eid+string(lvtcorb*100,format='(f6.2)')+'%',$
    l0=l0,b0=ang[1],rsun=ang[2])

  map2fits,mapb,'out_files/maps_'+pname+grdname+'_'+enme+'_FPMB.fits'
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  if keyword_set(plot) then begin
    loadct,39,/silent
    !p.multi=[0,2,1]
    plot_map,mapa,/log,chars=1.5,tit=mapa.id,/limb,grid_spacing=15
    plot_map,mapb,/log,chars=1.5,tit=mapb.id,/limb,grid_spacing=15
  endif

  ;  stop

end