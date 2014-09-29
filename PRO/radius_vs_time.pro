!P.charthick=2.
!P.thick=1.
!X.thick=2.
!Y.thick=2.
;Important parameters that should coincide with those in
;parameters.f90
np=8
nx=16
ny=16
neq=6
nlevs=7
xmax=1.
ymax=.5
nbrootx=2
nbrooty=1
path='./R1.0pc-DMC-ion/BIN'
outpref='./EPS/R0.4pc-dmc-rad-vs-t-ionized'
double=1
iffile=1
rad=fltarr(101)
imax=nbrootx*nx*2^(nlevs-1)
xphys=100.
x=findgen(101)+1

;------------------------------------------------------
FOR nout=0,100 DO BEGIN

mesh=readwmesh   (np,nx,ny,nlevs,path,nout,xmax,ymax,nbrootx,1)
rho =readwblocks (np,nx,ny,nlevs,path,nout,xmax,ymax,nbrootx,1,mesh,neq,0,double=double)

i=(size(rho))[1]-1
j=(size(rho))[2]-1
while (rho[i,j] lt 1.5 ) do begin
    i=i-1
    j=j-1
    RAD[NOUT]=Sqrt(float(i-1024)^2+float(j)^2)
endwhile

ENDFOR

rad=rad*xphys/float(imax)



;-----------------------------------------------------------------------
;   dimensiones (en centimetros) de la grafica
fig_size_x=10.
fig_size_y=9.
;   definimos un dx, dy para usar sistema de coordenadas normalizado
dx_fig= 1./fig_size_x
dy_fig= 1./fig_size_y

file_out=outpref+'.eps'
;-----------------------------------------------------------------------
if iffile eq 1 then begin
    set_plot,'PS'
    device,filename=file_out,/CMYK, $
      /encapsulated,xsize=fig_size_x,ysize=fig_size_y $
      ,/color,bits=24
endif
;-----------------------------------------------------------------------

plot,x,rad,/ylog,/xlog,xs=1,ys=1,xrange=[1,100],yrange=[5,50],psym=-1 $
  ,xtitle='t [kyr]', ytitle='R [pc]',symsize=.5

y=6.*x^(2./5.)
y2= 10*x^(1./4.)
oplot,x,y,linestyle=1
oplot,x,y2,linestyle=2

;-----------------------------------------------------------------------
if iffile eq 1 then begin
    device,/close
    set_plot,'X'
;    spawn, 'epstopdf '+file_out+'&'  ; cooment
endif
;-----------------------------------------------------------------------

;FOR nout=41,50 DO BEGIN

;mesh=readwmesh   (np,nx,ny,nlevs,path,nout,xmax,ymax,nbrootx,1)
;rho =readwblocks (np,nx,ny,nlevs,path,nout,xmax,ymax,nbrootx,1,mesh,neq,0,double=double)
;
;i=0
;while (rho[i,1000] lt 1.5 ) do begin
;    i=i+1
;    RAD[NOUT]=sqrt( (imax/2-i-1)^2+(1000.)^2 )
;endwhile
;ENDFOR


end
