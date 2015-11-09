pro make_ns_srm,pid=pid,fid=fid,regnm=regnm,de=de,nupdir=nupdir

  ; Make the SRM for OSPEX from the NuSTAR RMF and ARF files
  ;
  ; Optional inputs
  ; pid     -  Pointing of orbit 4 (1 or 2, default 1)
  ; fpmid   -  FPM to use ('A' or 'B', default is 'A')
  ; regnm   -  Region to fit ['D1','D2','L1','L2','L3'], default is'L3')
  ; de      -  Energy bin width (default 0.2 keV from 1.6 to 10 keV)
  ; nupdir  -  Where the output of the nupipeline ARF and RMF are stored
  ;
  ; 09-Nov-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Which pointing to use
  if (n_elements(pid) ne 1) then pid=1
  pname='P'+string(pid,format='(i1)')
  ; Which FPM to use
  if (n_elements(fid) ne 1) then fid='A'
  fpmid='FPM'+fid
  ; Which region to do
  if (n_elements(regnm) ne 1) then regnm='L3'
  ; Where are the nupipeline outputs stored
  ; These would have been produced by something like
  ; nuproducts indir=./ instrument=FPMA steminputs=nu20012001001 outdir=./P1L3A extended=no runmkarf=yes runmkrmf=yes infile=nu20012001001A06_cl.evt pilow=0 pihigh=85 bkgextract=no srcregionfile=circ_P1L3.reg attfile=./nu20012001001_att.fits hkfile=./nu20012001001A_fpm.hk 
  ; The *.reg files were made via get_spectra_reg.pro (after checking in ds9 with the original *.evt files)
   if (n_elements(nupdir) ne 1) then nupdir='../rsp/'+pname+regnm+fpmid+'/'

  if (n_elements(de) ne 1) then de=0.2
  maxe=10.
  mine=1.6
  eid=mine+de*indgen(2+(maxe-mine)/de)
  newe=get_edges(eid,/mean)

  intnam='INT'+string(1000+100*de,format='(i4)')
  nspec=n_elements(newe)
  sengs=[newe-0.5*de,newe[nspec-1]+0.5*de]

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ; These are the outputs from nupipeline
  flrmf=file_search(nupdir,'*.rmf')
  flarf=file_search(nupdir,'*.arf')

  ; The idl procedures to make these are mostly at
  ; https://lost-contact.mit.edu/afs/physics.wisc.edu/home/craigm/lib/idl/spectral/
  ; vcol2arr.pro and pointer_value.pro are at
  ; https://lost-contact.mit.edu/afs/physics.wisc.edu/home/craigm/lib/idl/util/
  rmfread, flrmf, rmf_str, compressed='UNCOMPRESSED'
  rmf = transpose(*(rmf_str.data))
  edges = *(rmf_str.ebins)
  fxbopen, unit, flarf, 'SPECRESP', hh, errmsg=err
  fxbreadm, unit, ['ENERG_LO','ENERG_HI','SPECRESP'], elo, ehi, arf, errmsg=err
  fxbclose, unit
  rsp=rmf
  for rr = 0, n_elements(rmf[0,*])-1 do rsp[*,rr]=rmf[*,rr]*arf[rr]

  nmx=min(where(edges[1,*] gt sengs[nspec]))

  edges=edges[*,0:nmx]
  rsp=rsp[0:nmx,0:nmx]

  medges=get_edges(edges,/mean)

  new_edges=fltarr(2,nspec)
  new_edges[0,*]=sengs[0:nspec-1]
  new_edges[1,*]=sengs[1:nspec]
  new_rsp_t=fltarr(nspec,nmx+1)
  for yy=0,nmx do new_rsp_t[*,yy]=interpol(rsp[*,yy],medges,newe)
  new_rsp=fltarr(nspec,nspec)
  for xx=0, nspec-1 do new_rsp[xx,*]=interpol(new_rsp_t[xx,*],medges,newe)
  edges=new_edges

  rsp=(newe[1]-newe[0])*new_rsp/(medges[1]-medges[0])
  rsp=rsp

  fl_rsp='out_files/RSP_'+pname+regnm+fpmid+'_'+intnam+'.dat'

  save,file=fl_rsp,edges,rsp

end
