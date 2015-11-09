pro box_igh,xc,yc,wid,color=color,thick=thick,linestyle=linestyle
  if (n_elements(color) ne 1) then color=0
  if (n_elements(thick) ne 1) then thick=1
  if (n_elements(linestyle) ne 1) then linestyle=0
  xr=xc+0.5*[-wid,wid]
  yr=yc+0.5*[-wid,wid]
  oplot,xr,yr[0]*[1,1],color=color,thick=thick,linestyle=linestyle
  oplot,xr,yr[1]*[1,1],color=color,thick=thick,linestyle=linestyle
  oplot,xr[0]*[1,1],yr,color=color,thick=thick,linestyle=linestyle
  oplot,xr[1]*[1,1],yr,color=color,thick=thick,linestyle=linestyle

end