pro specs_pileup_test,pid=pid,fid=fid,de=de,noplot=noplot

  ; Test how much of a NuSTAR spectrum is affected by pile-up Nov 14 P1 or P2
  ;
  ;  Starts by filtering out out "bad" pixels in FPMA and only the one CHU time CHU23 (P1,P2) and then
  ;  only those in Grade 0 (to make sure have only single pixel hit events). We also get the livetime
  ;  correction and shift pointing so all match P2 before saving it all out
  ;
  ;         Optional inputs:
  ;         pid       - 1=P1 (default) or 2=P2 ?
  ;         fid       - A (default) or B?
  ;         de        - Energy binning of spectrum (default 0.2keV)
  ;         noplot    - Keyword not to produce plot of histograms
  ;                   - Need >8.0 to do the plots (was done on v8.5)
  ;
  ; Assumes that make_roi_nov14.pro has already been run, so roi_nov14.dat exists in out_files/
  ; Assumes that make_gradehist_nov14o4 has been run so that the required  grds*.dat is in out_files/
  ; Assumes that make_specs_nov14o4_de has been run so that the required spec*.dat is in out_files/
  ;
  ; Still confused by the factor of 2(ish)
  ;
  ; 06-11-2015 IGH
  ; 13-11-2015 IGH - Added in the fitted f_vth model for comparison
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (n_elements(pid) ne 1) then pid=1
  pname='P'+string(pid,format='(i1)')

  if (n_elements(de) ne 1) then de=0.2
  denm='DE'+strmid(string(1000+de*100,format='(i4)'),1,3)
  intnam='INT'+string(1000+100*(de),format='(i4)')

  if (n_elements(fid) ne 1) then fid='A'
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Load in the spectra structure, spc_out
  restore,file='out_files/specs_'+denm+'_'+pname+'FPM'+fid+'.dat'
  ; Load in the grade structure, grd_out
  restore,file='out_files/grds_'+pname+'FPM'+fid+'.dat'

  nr=n_elements(spc_out.reg_nm)
  nengs=n_elements(spc_out.engs)
  specs=spc_out.spec_reg
  mide=spc_out.engs
  engs=spc_out.eng_edges
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Load in the fitted spectra

  spec_mod=fltarr(nr,nengs)

  if (pid eq 1) then ids=indgen(5)
  if (pid eq 2) then ids=[0,2]
  nids=n_elements(ids)

  for xx=0, nids-1 do begin
    i=ids[xx]
    restore,file='out_files/fitvth_'+pname+spc_out.reg_nm[i]+'FPM'+fid+'_'+intnam+'.dat'
    spec_mod[i,*]=fit_out.cnt_mod

  endfor

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Work out the pile-up factor as=0.25*total(grade 21-24)/total(grade 0)
  puf=0.25*grd_out.g2124_g0_reg

  ; Calculate the pile-up using different methods
  pl_specs=fltarr(nr,nengs)
  pl2_specs=fltarr(nr,nengs)
  plcvl_specs=fltarr(nr,nengs)

  for rr=0,nr-1 do begin
    ; David's method but needs the spectrum to start at 0 so expand spectrum down
    ; num bins from zero
    num_zero=round(engs[0]/(engs[1]-engs[0]))
    spec0=[replicate(0,num_zero),reform(specs[rr,*])]
    de=(engs[1]-engs[0])
    mide0=[findgen(num_zero)*de+0.5*de,mide]

    ; Calculate the double shifted spectrum
    pdu_spec=interpol(spec0,mide0,mide0/2.0)
    ; The NuSTAR pile-up fraction (PUF) is the fraction of events when another one actually occured
    ; so the pile-up spectrum is the double shifted one * PUF (?? since double shifting the single pixel spectrum ??)
    ; This version is about *2 of pl2_specs and plcvl_specs
    pl_specs[rr,*]=pdu_spec[num_zero:n_elements(pdu_spec)-1]*puf[rr]

    ; David's version of the double shifted spectrum also normalises by total(spec)/total(shifted)=0.5
    ; So assuming the pile-up fraction is the fraction of counts when two occured ???
    pdu_temp=pdu_spec
    pdu_temp*=total(specs[rr,*])/total(pdu_spec)*puf[rr]
    pl2_specs[rr,*]=pdu_temp[num_zero:n_elements(pdu_spec)-1]

    ; David's convolution version also normalises by total(spec)/total(convoled)=2.3
    ; So assuming the pile-up fraction is the fraction of counts when two occured ???
    plcvl_temp =pile_convol(spec0)
    plcvl_temp *= total(spec0)/total(plcvl_temp)*puf[rr]
    plcvl_specs[rr,*]=plcvl_temp[num_zero:n_elements(plcvl_temp)-1]

  endfor

  ; Don't want to save anything for D2, L2, or L3 with P2
  if (pid eq 2) then specs[[1,3,4],*]=0.

  ; How many counts do we have piled-up above some chosen energies?
  id5=where(engs ge 5.)
  id6=where(engs ge 6.)
  id7=where(engs ge 7.)

  cnt_obs=intarr(3,5)
  cnt_mod=fltarr(3,5)
  cnt_pile=fltarr(3,5)
  cnt_pile2=fltarr(3,5)
  cnt_pilecvl=fltarr(3,5)

  for i=0, nr-1 do begin

    cnt_obs[0,i]=total(specs[i,id5])
    cnt_obs[1,i]=total(specs[i,id6])
    cnt_obs[2,i]=total(specs[i,id7])

    cnt_mod[0,i]=total(spec_mod[i,id5])
    cnt_mod[1,i]=total(spec_mod[i,id6])
    cnt_mod[2,i]=total(spec_mod[i,id7])

    cnt_pile[0,i]=total(pl_specs[i,id5])
    cnt_pile[1,i]=total(pl_specs[i,id6])
    cnt_pile[2,i]=total(pl_specs[i,id7])

    cnt_pile2[0,i]=total(pl2_specs[i,id5])
    cnt_pile2[1,i]=total(pl2_specs[i,id6])
    cnt_pile2[2,i]=total(pl2_specs[i,id7])

    cnt_pilecvl[0,i]=total(plcvl_specs[i,id5])
    cnt_pilecvl[1,i]=total(plcvl_specs[i,id6])
    cnt_pilecvl[2,i]=total(plcvl_specs[i,id7])
  endfor

  print,pname,' FPM',fid
  print, 'Observed Counts'

  for xx=0, nids-1 do begin
    i=ids[xx]
    print,spc_out.reg_nm[i],' ',string(cnt_obs[*,i],format='(i5)')
  endfor

  print, 'Model Counts'
  for xx=0, nids-1 do begin
    i=ids[xx]
    print,spc_out.reg_nm[i],' ',string(cnt_mod[*,i],format='(f6.3)')
  endfor

  print, 'Observed - Convolved Counts'
  for xx=0, nids-1 do begin
    i=ids[xx]
    print,spc_out.reg_nm[i],' ',string(cnt_obs[*,i]-cnt_pilecvl[*,i],format='(f6.3)')
  endfor

  print,'Confidence to reject null hypothesis'
  ; null hypothesis = f-vth/isothermal is all there is in the observation
  ; 1\sigma=0.8413, 2\sigma=0.9772, 3\sigma=0.9987
  ; So need value >0.9987 to be 3\sigma confident we have something more than just f_vth
  for xx=0, nids-1 do begin
    i=ids[xx]
    print,spc_out.reg_nm[i],' ',$
      string(1-igamma(cnt_obs[*,i]-cnt_pilecvl[*,i],cnt_mod[*,i]),format='(f6.3)')
  endfor

  ;  print, 'Double shifted - IGH Scaling'
  ;  for i=0, nr-1 do print,spc_out.reg_nm[i],' ',string(cnt_pile[*,i],format='(f6.3)')
  ;
  ;  print, 'Double shifted - DS Scaling'
  ;  for i=0, nr-1 do print,spc_out.reg_nm[i],' ',string(cnt_pile2[*,i],format='(f6.3)')

  ;  print, 'Convolved - DS Scaling'
  ;  for i=0, nr-1 do print,spc_out.reg_nm[i],' ',string(cnt_pilecvl[*,i],format='(f6.3)')

  ; Plot it out if have new enough idl and no /noplot
  if (keyword_set(noplot) ne 1 and float(!version.release) ge 8.0) then begin
    !except=0
    w=window(dimensions=[900,600],/buffer)
    mrg=0.2

    for i=0, nr-1 do begin
      zz=plot(mide,specs[i,*],/ylog,/stair,margin=mrg,ytickunits='scientific', $
        title=spc_out.reg_nm[i]+', $G_{21-24}/G_0$='+string(grd_out.g2124_g0_reg[i]*100.,format='(f6.3)')+'%',$
        yrange=[1e-2,1e4],xrange=[1.0,10.0],xtitle='Energy [keV]',ytitle='NuSTAR Counts',layout=[3,2,2+i],/current)
      ; for some reason /stair (hist centred on x) goes weird so need to do /hist (hist starts on x)
      !null=plot(mide-de*0.5,pl_specs[i,*],/hist,color='firebrick',/over)
      !null=plot(mide-de*0.5,pl2_specs[i,*],/hist,color='forest green',/over)
      !null=plot(mide-de*0.5,plcvl_specs[i,*],/hist,color='medium blue',/over)
    endfor

    tt=text(0,2,pname+' FPM'+fid+' '+$
      string(spc_out.expos_dl,format='(f6.2)')+'s ('+string(spc_out.dur,format='(f6.1)')+'s, '+$
      string(spc_out.lvt_frac*100,format='(f5.2)')+'%)',/device,target=w)

    !null=text(50,500,'$C(E)$',/device,target=w)
    !null=text(50,475,'$P\times C(E/2)$',color='firebrick',/device,target=w)
    !null=text(50,450,'$P\times C(E/2)\times \Sigma C(E)/\Sigma C(E/2)$',color='forest green',/device,target=w)
    !null=text(50,425,'$P\times CNVL\times \Sigma C(E)/\Sigma CNVL$',color='medium blue', /device,target=w)

    w.save,'figs/pileup_'+pname+'FPM'+fid+'.pdf',page_size=w.dimensions/100.
    w.close
    !except=1
  endif

end