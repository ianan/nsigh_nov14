function pile_convol, spectrum

  ;The piled-up spectrum (pilespec) is initialized as zero.
  ;Then you step through the spectrum one bin at a time.
  ;To the piled-up spectrum you add the counts in the current bin
  ;times the whole spectrum above it shifted up in energy
  ;by the energy of the current bin.  You use only the spectrum
  ;above the current bin to avoid double counting.  This is a sort
  ;of auto-convolution, I think.

  ;Note this only works when the energy scale starts at E=0.  I can
  ;try to write version where the scale starts at nonzero energy if
  ;it's needed.
  
  ; from D. Smith 

  pilespec = spectrum*0.
  piler = spectrum
  pilee = spectrum

  tpilee = pilee/total(pilee)
  tshift = total( /cum, piler)
  nume = n_elements(spectrum)

  for n=1l,n_elements(spectrum)-2 do begin

    shiftspec = piler[0:n-1]
    if ( tshift[n-1] EQ 0.) then normshift = shiftspec $
    else         normshift =  shiftspec*tpilee[n]

    ;This is the version by Richard Schwartz:
    pilespec[n+1] = pilespec[n+1:*]+ normshift

    ;This is the version by me, slower but more explicit in the indices.
    ;Does the same thing.
    ;    top = (2*n+1) < nume
    ;   pilespec[n+1:top-1] = pilespec[n+1:top-1]+ normshift[0:top-1-(n+1)]

  end

  return, pilespec
end