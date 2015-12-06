;+
; NAME:
;     NUSTAR_SUNPOINT
; PURPOSE:
;     NuSTAR pointing for solar observations, epoch J2000
; CATEGORY:
; CALLING SEQUENCE:
;     point = nustar_sunpoint()            (returns solar [RA, dec])
;     point = nustar_sunpoint(time=time)   (returns solar [RA, dec] at time)
;     point = nustar_sunpoint([925.,0.],time=time) points near W limb
;     point = nustar_sunpoint([0.,-925.],time=time,note='S Pole') points and prints
; INPUTS:
;     [x,y] heliographic angular distances, arc s, default [0,0]
;     time, ANYTIM time (default now)
; OPTIONAL (KEYWORD) INPUT PARAMETERS:
;     epoch_in, year of input coordinates (default current time in years)
;     epoch_out, year of output coordinates (default 2000)
; OUTPUTS:
;     basic pointing coordinates (RA, dec) in degrees, for epoch 2000
; CALLS: JULDAY, SUNPOS, PB0R
; RESTRICTIONS:
;     small-angle approximation; needs more sines and cosines
; MODIFICATION HISTORY:
;     written 16-Jul-2014 (HSH)
;     18-Jul-14 fixed bugs found by LG and JCMO (HSH)
;     27-Jul-14, fixed bug found by AM and set output to J2000 by default (HSH)
;     31-Oct-14, added print statement (HSH)
;     03-Dec-15, corrected HSH's wrong epoch as fraction of year (IGH)
;-

function nustar_sunpoint, offset, time=time, epoch_from=epoch_from, $
  epoch_to=epoch_to, qdebug=qdebug, note=note

  default, time, !stime
  default, offset, [0.,0.]
  default, epoch_to, 2000.
  text = anytim(time,/ext)
  sunpos,julday(text[5],text[4],text[6],text[0],text[1],text[2]),ra,dec
;  epoch = 1.*text[6]+text[5]/12.+text[4]/365.
  year=string(text[6],format='(i4)')
  epoch=year+(anytim(time)-anytim('01-Jan-'+year))/(anytim('31-Dec-'+year+' 24:00')-anytim('01-Jan-'+year))
  
  if keyword_set(epoch_from) then epoch = epoch_from
  pang = (pb0r(time))[0]/!radeg
  matrix = [[cos(pang),sin(pang)],[-sin(pang),cos(pang)]]
  xy = offset/3600.
  xy1 = matrix#xy
  dec1 = dec + xy1[1]
  ra1 = ra - xy1[0]/cos(dec/!radeg)
  ra = ra1 & dec = dec1
  ; FK5 defined on J2000

  precess,ra,dec,epoch,epoch_to
  if n_elements(note) eq 1 then begin
    format = '(a20,2f10.3,a20)'
    print,time,ra,dec,note,format=format
  endif

  if keyword_set(qdebug) then stop

  return, [ra,dec]
end
