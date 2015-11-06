pro grade_analysis_wip

  ; firstly make sure the *.dat with the regions of interest co-ords are saved out
  make_roi_nov14

  ; Then make the grade histograms per fov and region
  ; need IDL >8.0 for the plots to work
  make_gradehist_nov14o4,pid=1,fid='A'
  make_gradehist_nov14o4,pid=1,fid='B'
  make_gradehist_nov14o4,pid=2,fid='A'
  make_gradehist_nov14o4,pid=2,fid='B'


  stop
end