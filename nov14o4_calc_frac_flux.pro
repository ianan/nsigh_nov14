pro nov14o4_calc_frac_flux

  ; What fraction of the AIA 94 and FeXVIII flus is in each of the selected regions
  ; and then saves it out (Note AIA mean maps not included here)
  ;
  ; 11-Nov-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Load in the SDO/AIA mean maps
  aiadir='~/data/aia/ev_01112014/'

  fits2map,aiadir+'P1_FE18_mn_full.fits',mfe18
  fits2map,aiadir+'P1_94_mn_full.fits',m94

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Loads in the roi, xcsp, ycsp,rnm,wid
  restore,file='out_files/roi_no14.dat'
  nr=n_elements(rnm)

  flux94_all=total(m94.data > 0.)
  flux94_reg=fltarr(nr)

  fluxfe18_all=total(mfe18.data > 0.)
  fluxfe18_reg=fltarr(nr)

  for i=0, nr-1 do begin
    xrang=xcsp[i]+0.5*wid*[-1,1]
    yrang=ycsp[i]+0.5*wid*[-1,1]

    sub_map,m94,sm94,xrange=xrang,yrange=yrang
    flux94_reg[i]=total(sm94.data > 0.)
    sub_map,mfe18,smfe18,xrange=xrang,yrange=yrang
    fluxfe18_reg[i]=total(smfe18.data > 0.)

    print,rnm[i],100*flux94_reg[i]/flux94_all,100*fluxfe18_reg[i]/fluxfe18_all
  endfor

  ff_out={reg_id:rnm,$
    flux94_all:flux94_all,flux94_reg:flux94_reg,$
    fluxfe18_all:fluxfe18_all,fluxfe18_reg:fluxfe18_reg,$
    prcnt94:100*flux94_reg/flux94_all, prcntfe18:100*fluxfe18_reg/fluxfe18_all}

  save,file='out_files/roi_no14_ffper.dat',ff_out
  
end