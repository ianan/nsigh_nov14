pro make_reg_file,x,y,rad,time,id
  ; Take the sun centred x,y in arsec and create a ds9 like *.reg file
  ; for a circular region (of radius rad) at that RA and dec
  ;
  ; Needs nustar_sunpoint.pro to run
  ;
  ; 09-Nov-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  pp=nustar_sunpoint([x,y],time=time)
  ras=pp[0]
  dcs=pp[1]
  radec,ras,dcs,rahr,ramin,rasec,dchr,dcmin,dcsec

  if (rasec ge 10.) then begin
    rash=string(rahr,format='(i3)')+':'+string(ramin,format='(i3)')+':'+string(rasec,format='(f6.3)')
  endif else begin
    rash=string(rahr,format='(i3)')+':'+string(ramin,format='(i3)')+':0'+string(rasec,format='(f6.3)')
  endelse

  if (dcsec ge 10.) then begin
    dcsh=string(dchr,format='(i3)')+':'+string(dcmin,format='(i3)')+':'+string(dcsec,format='(f6.3)')
  endif else begin
    dcsh=string(dchr,format='(i3)')+':'+string(dcmin,format='(i3)')+':0'+string(dcsec,format='(f6.3)')
  endelse

  rash=strcompress(rash,/rem)
  dcsh=strcompress(dcsh,/rem)
  rads=string(rad,format='(f7.2)')+'"'

  fileout='out_files/circ_'+id+'.reg'
  openw, lun, /get_lun, fileout
  printf, lun, '# Region file format: DS9 version 4.1'
  printf, lun, 'global color=white dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
  printf, lun, 'fk5'
  printf, lun, 'circle('+rash+','+dcsh+','+rads+')'
  close, lun
  free_lun, lun

end

pro make_reg_file_multi,x,y,rad,time,id
  ; Take the sun centred x,y in arsec and create a ds9 like *.reg file
  ; for a circular region (of radius rad) at that RA and dec
  ; puts all the regions in one file - useful for loading/comparison in ds9
  ;
  ; Needs nustar_sunpoint.pro to run
  ;
  ; 09-Nov-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  nr=n_elements(x)

  fileout='out_files/circ_'+id+'.reg'
  openw, lun, /get_lun, fileout
  printf, lun, '# Region file format: DS9 version 4.1'
  printf, lun, 'global color=white dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
  printf, lun, 'fk5'

  for ii=0,nr-1 do begin


    pp=nustar_sunpoint([x[ii],y[ii]],time=time)
    ras=pp[0]
    dcs=pp[1]
    radec,ras,dcs,rahr,ramin,rasec,dchr,dcmin,dcsec

    if (rasec ge 10.) then begin
      rash=string(rahr,format='(i3)')+':'+string(ramin,format='(i3)')+':'+string(rasec,format='(f6.3)')
    endif else begin
      rash=string(rahr,format='(i3)')+':'+string(ramin,format='(i3)')+':0'+string(rasec,format='(f6.3)')
    endelse

    if (dcsec ge 10.) then begin
      dcsh=string(dchr,format='(i3)')+':'+string(dcmin,format='(i3)')+':'+string(dcsec,format='(f6.3)')
    endif else begin
      dcsh=string(dchr,format='(i3)')+':'+string(dcmin,format='(i3)')+':0'+string(dcsec,format='(f6.3)')
    endelse

    rash=strcompress(rash,/rem)
    dcsh=strcompress(dcsh,/rem)
    rads=string(rad,format='(f7.2)')+'"'

    printf, lun, 'circle('+rash+','+dcsh+','+rads+')'


  endfor

  close, lun
  free_lun, lun

end

pro make_reg_file_multi_sq,x,y,wid,ang,time,id
  ; Take the sun centred x,y in arsec and create a ds9 like *.reg file
  ; for a square (of width wid) at that RA and dec
  ; puts all the regions in one file - useful for loading/comparison in ds9
  ;
  ; Note - nuproducts/numkarf/psf doesn't supports this shape, just for comparison in ds9
  ;
  ; Needs nustar_sunpoint.pro to run
  ;
  ; 09-Nov-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  nr=n_elements(x)

  fileout='out_files/sq_'+id+'.reg'
  openw, lun, /get_lun, fileout
  printf, lun, '# Region file format: DS9 version 4.1'
  printf, lun, 'global color=white dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
  printf, lun, 'fk5'

  for ii=0,nr-1 do begin


    pp=nustar_sunpoint([x[ii],y[ii]],time=time)
    ras=pp[0]
    dcs=pp[1]
    radec,ras,dcs,rahr,ramin,rasec,dchr,dcmin,dcsec

    if (rasec ge 10.) then begin
      rash=string(rahr,format='(i3)')+':'+string(ramin,format='(i3)')+':'+string(rasec,format='(f6.3)')
    endif else begin
      rash=string(rahr,format='(i3)')+':'+string(ramin,format='(i3)')+':0'+string(rasec,format='(f6.3)')
    endelse

    if (dcsec ge 10.) then begin
      dcsh=string(dchr,format='(i3)')+':'+string(dcmin,format='(i3)')+':'+string(dcsec,format='(f6.3)')
    endif else begin
      dcsh=string(dchr,format='(i3)')+':'+string(dcmin,format='(i3)')+':0'+string(dcsec,format='(f6.3)')
    endelse

    rash=strcompress(rash,/rem)
    dcsh=strcompress(dcsh,/rem)
    wids=string(wid,format='(f7.2)')+'"'
    angs=strcompress(string(ang,format='(f5.2)'),/rem)

    printf, lun, 'box('+rash+','+dcsh+','+wids+','+wids+','+angs+')'


  endfor

  close, lun
  free_lun, lun

end


pro get_spectra_reg

  ; For the regions selected in Sun centred arcsec what are the RA/Dec postion of the region
  ; for ds9 and nuproducts. Note that the output should be checked/maybe manualled shifted
  ; in ds9 vs the original *.evt file (the RA/Dec one, not the sunpos corrected one)
  ;
  ; 09-Nov-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Load in the regions
  restore,file='out_files/roi_no14.dat'
  nr=n_elements(rnm)

  ; match the area of box to circle
  rad=wid/sqrt(!pi)

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; For P1
  id='P1'
  time=anytim(mean(anytim(['01-Nov-14 21:34:49','01-Nov-14 21:47:59'])),/yoh,/trunc)
  shftx=-70.
  shfty=-10.

  ; Make a *.reg text file for each region
  for rr=0,nr-1 do make_reg_file,xcsp[rr]+shftx,ycsp[rr]+shfty,rad,time,id+rnm[rr]

  ; Make sinlge *.reg text file containing multiple regions
  ; useful for testing in ds9 - not sure if nupipeline can handle it though
  make_reg_file_multi,xcsp+shftx,ycsp+shfty,rad,time,id

  ; Make sinlge *.reg text file containing multiple regions
  ; useful for testing in ds9 - not sure if nupipeline can handle it though
  ; this is the box one and mostly just for testing
  make_reg_file_multi_sq,xcsp+shftx,ycsp+shfty,wid,23.45,time,id

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; For P2
  id='P2'
  time=anytim(mean(anytim(['01-Nov-14 21:48:52','01-Nov-14 22:00:59'])),/yoh,/trunc)
  shftx=-90.
  shfty=-10.

  ; Make a *.reg text file for each region
  for rr=0,nr-1 do make_reg_file,xcsp[rr]+shftx,ycsp[rr]+shfty,rad,time,id+rnm[rr]

  ; Make sinlge *.reg text file containing multiple regions
  ; useful for testing in ds9 - not sure if nupipeline can handle it though
  make_reg_file_multi,xcsp+shftx,ycsp+shfty,rad,time,id

  ; Make sinlge *.reg text file containing multiple regions
  ; useful for testing in ds9 - not sure if nupipeline can handle it though
  ; this is the box one and mostly just for testing
  make_reg_file_multi_sq,xcsp+shftx,ycsp+shfty,wid,23.45,time,id

end