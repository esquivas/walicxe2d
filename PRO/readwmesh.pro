FUNCTION READWMESH,np,nx,ny,nlevs,path,nout,xmax,ymax,nbrootx,nbrooty
;READS MRSH FROM THE BINARY OUTPUT FROM THE 2D VERSION OF WALIXCE

;Important parameters that should coincide with those in
;parameters.f90
;np           number of processors
;nx(y)        block size in x(y)
;nlevs        # of levels
;x(y)max      # box size in x(y)
;nbrootx(y)   # of parent blocks in x(y) directions
;nout         # of output to read
nbmin   = lonarr(np)
nbmax   = lonarr(np)
nbleafs = lonarr(np)

dx=fltarr(nlevs)
dy=fltarr(nlevs)

for i=0,nlevs-1 do begin
    dx[i]=(xmax/nbrootx)/nx/(2.^I)
    dy[i]=(ymax/nbrooty)/ny/(2.^I)
endfor

print,dx
print,'***'
print,dx

;------------------------------------------------------
for rank=0,np-1 do begin
    OPENR, Unit, path+'/pointers' $
      +string(rank,'(i3.3)')+'.'  $
      +string(nout,'(i3.3)')+'.bin' $
      ,/get_lun,/swap_if_big_endian
    stream  = assoc(Unit,lonarr(1))
    nbmin  [rank] = (stream[0])[0]
    nbmax  [rank] = (stream[1])[0]
    nbleafs[rank] = (stream[2])[0]
    close,/all
endfor 
;------------------------------------------------------
mesh=fltarr(fix(total(nbleafs)),6) ; x0, and y0, level,NP
ll=0l
;------------------------------------------------------
;   READ THE MESH
for rank=0,np-1 do begin
    OPENR, Unit, path+'/pointers' $
      +string(rank,'(i3.3)') $
      +'.'+string(nout,'(i3.3)')+'.bin' $
      ,/get_lun,/swap_if_big_endian
    stream  = assoc(Unit,lonarr(1))
    nbmin  [rank] = (stream[0])[0]
    nbmax  [rank] = (stream[1])[0]
    nbleafs[rank] = (stream[2])[0]
    lp     = lonarr(nbmax[rank]-nbmin[rank]+1,8)
    leaf   = lonarr(nbmax[rank]-nbmin[rank]+1  )
    icoord = lonarr(nbmax[rank]-nbmin[rank]+1,2)
    if ((stream[3])[0] gt nlevs) then  k=4l else k=3l

    for i=0,7 do begin
        for nb=0,nbmax[rank]-nbmin[rank] do begin
            lp[nb,i]=stream[k]
            k=k+1
        endfor
    endfor
    for nb=0,nbmax[rank]-nbmin[rank] do begin
         ;    keep in the leaf variable the level of each block
        if (stream[k] ne 0) then leaf[nb]=lp[nb,0]
        k=k+1
    endfor
    for i=0,1 do begin
        for nb=0,nbmax[rank]-nbmin[rank] do begin
            icoord[nb,i]=stream[k]
            k=k+1
        endfor
    endfor
    
    close,/all
    
    for nb=0,nbmax[rank]-nbmin[rank] do begin
        if (leaf[nb] ne 0) then begin 
            i=leaf[nb]-1
            mesh[ll,0]= icoord[nb,0]    *dx[i]
            mesh[ll,1]=(icoord[nb,0]+nx)*dx[i]
            mesh[ll,2]= icoord[nb,1]    *dy[i]
            mesh[ll,3]=(icoord[nb,1]+ny)*dy[i]
            mesh[ll,4]= leaf[nb]
            mesh[ll,5]= rank
            ll=ll+1
        endif
    endfor
endfor

RETURN, MESH


end
