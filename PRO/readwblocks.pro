FUNCTION READWBLOCKS,np,nx,ny,nlevs,path,nout,xmax,ymax,nbrootx,nbrooty,mesh,neq,ieq, double=double
;READS THE BINARY OUTPUT FROM THE 2D VERSION OF WALIXCE

;Important parameters that should coincide with those in
;parameters.f90
;np           number of processors
;nx(y)        block size in x(y)
;neq          # of equations
;nlevs        # of levels
;x(y)max      # box size in x(y)
;nbrootx(y)   # of parent blocks in x(y) directions

nbmin   = lonarr(np)
nbmax   = lonarr(np)
nbleafs = lonarr(np)

dx=fltarr(nlevs)
dy=fltarr(nlevs)

for i=0,nlevs-1 do begin
    dx[i]=(xmax/nbrootx)/nx/(2.^I)
    dy[i]=(ymax/nbrooty)/ny/(2.^I)
endfor

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

;------------------------------------------------------
;   READ THE BLOCKS
if keyword_set( double ) then begin
    blocks=dblarr(fix(total(nbleafs)),neq,nx+4,ny+4)
endif else begin
    blocks=fltarr(fix(total(nbleafs)),neq,nx+4,ny+4)
endelse

ll=0l
for rank=0,np-1 do begin
    OPENR, Unit, path+'/blocks' $
      +string(rank,'(i3.3)') $
      +'.'+string(nout,'(i3.3)')+'.bin' $
      ,/get_lun,/swap_if_big_endian

    if keyword_set( double ) then begin
            stream  = assoc(Unit,dblarr(1,neq,nx+4,ny+4))
    endif else begin
        stream  = assoc(Unit,fltarr(1,neq,nx+4,ny+4))
    endelse

    for NB=0,nbleafs[rank]-1 do begin
        blocks[ll,*,*,*]=stream[nb]
        ll=ll+1
    endfor
    
    close,/all

endfor
;------------------------------------------------------
if keyword_set( double ) then begin
    map=dblarr(2^(nlevs-1)*nx*nbrootx,2^(nlevs-1)*ny)
endif else begin
    map=fltarr(2^(nlevs-1)*nx*nbrootx,2^(nlevs-1)*ny)    
endelse

for nb=0,(size(mesh))[1]-1 do begin
    x0=fix(mesh[nb,0]/dx[nlevs-1])
    y0=fix(mesh[nb,2]/dy[nlevs-1])
    if (mesh[nb,4] eq nlevs) then begin
        x1=x0+nx-1
        y1=y0+ny-1
        map[x0:x1,y0:y1]=blocks[nb,ieq,2:nx+1,2:ny+1]
    endif else begin
        fac=2^(nlevs-mesh[nb,4])
        x1=x0+fac*nx-1
        y1=y0+fac*nx-1
        map[x0:x1,y0:y1]= $
          congrid(reform(blocks[nb,ieq,2:nx+1,2:ny+1]),fac*nx,fac*ny,/center)
    endelse
endfor

return, map

end
