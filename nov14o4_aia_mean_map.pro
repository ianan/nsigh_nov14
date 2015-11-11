pro nov14o4_aia_mean_map,aid=aid, sub_all=sub_all

  ; Produce a mean map for the specified SDO/AIA channel for comparison to the NuSTAR P1 data
  ; Preps the *.fits files before averaging them over the time of P1 NuSTAR obs

  ;  Optional inputs:
  ;     aid   - Which AIA channel? (0-5, for 94A,131A,171A,193A,211A,335A, default 0 or 94A)
  ;     sub_all   - Save out all the submaps to a single fits file (default no as BIG file)
  ;
  ; 11-Nov-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; start/end of P1
  t1=' 1-Nov-2014 21:34:49'
  t2=anytim(anytim(t1)+790,/yoh,/trunc)

  ; which SDO/AIA channel do we want?
  wavn=['94','131','171','193','211','335']
  if (n_elements(aid) eq 0) then aid=0
  wvn=wavn[aid]

  ; Where the AIA data is stored and the mean map will be outputted
  aiadir='~/data/aia/ev_01112014/'
  ff=file_search(aiadir+'*.'+wvn+'A_*.fits')
  idt=strpos(ff[0],'01T')

  ftims=anytim('01-Nov-2014 '+strmid(ff,idt+3,2)+':'+strmid(ff,idt+6,2)+':'+strmid(ff,idt+9,2))

  timer=[t1,t2]

  gd=where(ftims ge anytim(t1) and ftims le anytim(t2),ngd)
  ; Sub-map region that covers AR of interest in P1 (and P2)
  xr=[350,1350]
  yr=[-650,350]

  for i=0, ngd-1 do begin
    print,i+1,' -- ',ngd
    ; Load in and prep the file
    aia_prep,ff[gd[i]],-1,ind,data
    ; Convert data and header to map
    index2map,ind,data*1.0,map_temp
    ; Extract the submap region
    sub_map,map_temp,smap_temp,xrange=xr,yrange=yr

    ; get rid of negatives (introduced by subtracting spikes ?)
    map_temp.data=map_temp.data >0.
    smap_temp.data=smap_temp.data >0.
    if (i eq 0) then fmp=map_temp else fmp.data=fmp.data+map_temp.data
    if (i eq 0) then smp=smap_temp else smp.data=smp.data+smap_temp.data
    if (i eq 0) then smps=smap_temp else smps=[smps,smap_temp]

  endfor

  ; Saves out all sub_maps to a single fits file - BIG
  if keyword_set(sub_all) then map2fits,smps,aiadir+'P1_'+wvn+'.fits'

  ; Average the sub_map and save it to a fits file
  ; Assumes all are of the same exposure - you should manually check smps.dur
  smp.data=smp.data/(ngd*1.0)
  smp.id='SDO/AIA '+wvn+'A, Mean over '+smps[0].time+' to '+anytim(smps[ngd-1].time,/time,/yoh)
  map2fits,smp,aiadir+'P1_'+wvn+'_mn.fits'

  ; Average the whole map and save it to a fits file
  ; Assumes all are of the same exposure - you should manually check smps.dur
  fmp.data=fmp.data/(ngd*1.0)
  fmp.id='SDO/AIA '+wvn+'A, Mean over '+smps[0].time+' to '+ anytim(smps[ngd-1].time,/time,/yoh)
  map2fits,fmp,aiadir+'P1_'+wvn+'_mn_full.fits'

;  clearplot
;  aia_lct,wave=wvn,/load
;  plot_map,smp,/log,title=smp.id
end