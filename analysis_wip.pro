pro analysis_wip

  ; firstly make sure the *.dat with the regions of interest co-ords are saved out
  make_roi_nov14

  ; Then make the grade histograms per fov and region
  ; need IDL >8.0 for the plots to work
  make_gradehist_nov14o4,pid=1,fid='A'
  make_gradehist_nov14o4,pid=1,fid='B'
  make_gradehist_nov14o4,pid=2,fid='A'
  make_gradehist_nov14o4,pid=2,fid='B'

  ; Then make the spectra per fov and regions
  ; again need IDL >8.0 for the plots to work
  make_specs_nov14o4_de,pid=1,fid='A'
  make_specs_nov14o4_de,pid=1,fid='B'
  make_specs_nov14o4_de,pid=2,fid='A'
  make_specs_nov14o4_de,pid=2,fid='B'

  ; Then test to see how much of the spectrum is piled-up
  specs_pileup_test,pid=1,fid='A'
  specs_pileup_test,pid=1,fid='B'
  specs_pileup_test,pid=2,fid='A'
  specs_pileup_test,pid=2,fid='B'

  ; Make the srm for ospex
  regnm=['D1','D2','L1','L2','L3']
  for i=0, n_elements(regnm)-1 do make_ns_srm,pid=1,fid='A',regnm=regnm[i]
  for i=0, n_elements(regnm)-1 do make_ns_srm,pid=1,fid='B',regnm=regnm[i]
  make_ns_srm,pid=2,fid='A',regnm='D1'
  make_ns_srm,pid=2,fid='A',regnm='L1'
  make_ns_srm,pid=2,fid='B',regnm='D1'
  make_ns_srm,pid=2,fid='B',regnm='L1'

  ; Fit an isothermal model to the spectrum
  for i=0, 4 do fitvth_spec_nov14o4,pid=1,fid='A',regid=i
  for i=0, 4 do fitvth_spec_nov14o4,pid=1,fid='B',regid=i
  fitvth_spec_nov14o4,pid=2,fid='A',regid=0
  fitvth_spec_nov14o4,pid=2,fid='A',regid=2
  fitvth_spec_nov14o4,pid=2,fid='B',regid=0
  fitvth_spec_nov14o4,pid=2,fid='B',regid=2
  
  stop
end