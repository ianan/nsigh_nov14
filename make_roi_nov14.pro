pro make_roi_nov14

  ; Save out the ROI info for the Nov 2014 data
  ; Centres in S/C arcsec of boxes with widthes of wid
  ; 02-Jul-2015 IGH
  
  ; Add in some background/ghost-ray regions
  ; 19-Oct-2015 IGH 

  xcsp=[820,660,1005,1010,970];,790,760,1300,1250]
  ycsp=[105,-260,145,-90,-240];,-100,-300,-100,-300]

  rnm=['D1','D2','L1','L2','L3'];,'B1','B2','B3','B4']
  wid=120.
  
  save,file='out_files/roi_no14.dat',xcsp,ycsp,rnm,wid
  
  xcsp=[820,660,1005,1010,970,790,790,1300,1250]
  ycsp=[105,-260,145,-90,-240,-100,-290,-100,-300]

  rnm=['D1','D2','L1','L2','L3','B1','B2','B3','B4']
  wid=120.

  save,file='out_files/roi_no14_bg.dat',xcsp,ycsp,rnm,wid
end
