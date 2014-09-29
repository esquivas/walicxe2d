!P.charthick=2.
!P.thick=1.
!X.thick=2.
!Y.thick=2.
;Important parameters that should coincide with those in
;parameters.f90
np=16
nx=16
ny=16
neq=6
nlevs=7
nbmin   = intarr(np)
nbmax   = intarr(np)
nbleafs = intarr(np)
xmax=1.0
ymax=0.25
xphys=3.
yphys=3.*.25
nbrootx=4
nbrooty=1
path='/datos_io/esquivel/W2D-test/BIN/'
double=1
outpref='./EPS/R1.0pc-DMC-ion-'
iffile=0
zoom=0
small=0.
minrho=.1    ;0.01
maxrho=1E4  ;1. 
minT=1e3
maxT=1e6
Cv=1.5
Tempsc= 1.E4*(Cv+1.)/Cv
;------------------------------------------------------
FOR nout=15,15 DO BEGIN

mesh=readwmesh   (np,nx,ny,nlevs,path,nout,xmax,ymax,nbrootx,1)

rho  =readwblocks (np,nx,ny,nlevs,path,nout,xmax,ymax,nbrootx,1,mesh,neq,0,double=double)
rhoN =readwblocks (np,nx,ny,nlevs,path,nout,xmax,ymax,nbrootx,1,mesh,neq,4,double=double)

; for the temperature
etot=readwblocks(np,nx,ny,nlevs,path,nout,xmax,ymax,nbrootx,1,mesh,neq,3,double=double)

px  =readwblocks(np,nx,ny,nlevs,path,nout,xmax,ymax,nbrootx,1,mesh,neq,1,double=double)

py  =readwblocks(np,nx,ny,nlevs,path,nout,xmax,ymax,nbrootx,1,mesh,neq,2,double=double)

pgas=(etot-0.5*(px*px+py*py)/rho)/ Cv  > 1e-15
rho = rho                              > 1e-15

dentot=(2.*rho-rhon)  ;Cooling
;dentot=rho    ; adiabatic y dmc
temp=(pgas/dentot)*Tempsc

print,'nout=',nout
print,minmax(rho)
print,minmax(temp)
print,minmax(pgas)

;px=0
;py=0
;pgas=0


;scale=1.3*1.66e-24
;rho=rho*scale

 rho=rho>minrho
 rho=rho<maxrho
 rho[0,0]=minrho
 rho[0,1]=maxrho

 temp=temp > minT
 temp=temp < maxT
 temp[0,0]= minT
 temp[0,1]= maxT

temp=rotate(temp,7)
;f=rotate(f,7)
;------------------------------------------------------
;------------------------------------------------------
;   make the plots
;------------------------------------------------------
if iffile eq 0 then begin

!P.Multi=[0,1,3]

xax=(findgen((size(rho))[1])+1)/((size(rho))[1])*xphys
yax=(findgen((size(rho))[2])+1)/((size(rho))[2])*yphys

plot,[0,xmax],[0,ymax],/nodata,xs=4,ys=4,/isotropic,ticklen=0 $
  ,title='t='+string(nout*1000,format='(i6)')+' yr'
for nb=0,(size(mesh))[1]-1 do begin
    plots,[mesh[nb,0],mesh[nb,0]],[mesh[nb,2],mesh[nb,3]]
    plots,[mesh[nb,1],mesh[nb,1]],[mesh[nb,2],mesh[nb,3]]
    plots,[mesh[nb,0],mesh[nb,1]],[mesh[nb,2],mesh[nb,2]]
    plots,[mesh[nb,0],mesh[nb,1]],[mesh[nb,3],mesh[nb,3]]
endfor


loadct,1

imagescbar,rho,xax,yax,/multi,/bar,/aspect,/noerase,background=0 $
  ,color=255,/log,min=minrho,max=maxrho,charsize=2.5,/force

;loadct,9
;imagescbar,f,/multi,/bar,/aspect,/noerase,background=0 $
;  ,color=255,charsize=2.5

loadct,3
imagescbar,temp,xax,yax,/multi,/bar,/aspect,/noerase,background=0 $
  ,color=255,/log,min=mint,max=maxt,charsize=2.5,/force

if zoom eq 1 then begin
wait,3
!P.Multi=[0,1,2]
loadct,1
imagescbar,rho [800:1200,0:350],/multi,/bar,/aspect,/noerase,background=0 $
  ,color=255,/log,min=minrho,max=maxrho,charsize=2.5
loadct,3
imagescbar,temp[800:1200,674:1023],/multi,/bar,/aspect,/noerase,background=0 $
  ,color=255,/log,min=mint,max=maxt,charsize=2.5
endif

endif

;-----------------------------------------------------------------------
;   dimensiones (en centimetros) de la grafica
fig_size_x=20.
fig_size_y=26.3
;   definimos un dx, dy para usar sistema de coordenadas normalizado
dx_fig= 1./fig_size_x
dy_fig= 1./fig_size_y

file_out=outpref+string(nout,format='(i3.3)')+'.eps'
;-----------------------------------------------------------------------
if iffile eq 1 then begin
    set_plot,'PS'
    device,filename=file_out,/CMYK, $
      /encapsulated,xsize=fig_size_x,ysize=fig_size_y $
      ,/color,bits=24
endif
;-----------------------------------------------------------------------
if small eq 1 then begin
    rho=congrid(rho,1024,512)
    temp=congrid(temp,1024,512)
    rho[0,0]=minrho
    rho[0,1]=maxrho
    temp[0,0]= minT
    temp[0,1]= maxT
endif

;primero la malla
;aux=fltarr(4,1)
if iffile eq 1 then begin

loadct,1
;tvscl,alog10(aux),xsize=16,ysize=4.,1.,9.5,/centimeters
plot,[0,xmax],[0,ymax],/nodata,xs=4,ys=4,/isotropic,ticklen=0 $
  ,position=[1.*dx_fig,17.5*dy_fig,17.*dx_fig,25.5*dy_fig],/normal $
 ,title='t='+string(nout,format='(i3)')+' 000 yr',charsize=1.3

for nb=0,(size(mesh))[1]-1 do begin
    plots,[mesh[nb,0],mesh[nb,0]],[mesh[nb,2],mesh[nb,3]]
    plots,[mesh[nb,1],mesh[nb,1]],[mesh[nb,2],mesh[nb,3]]
    plots,[mesh[nb,0],mesh[nb,1]],[mesh[nb,2],mesh[nb,2]]
    plots,[mesh[nb,0],mesh[nb,1]],[mesh[nb,3],mesh[nb,3]]
endfor

col=0

loadct,1
;la densidad
tvscl,alog10(rho),xsize=16,ysize=8.,1.,9.,/centimeters

range_bar=[minrho,maxrho]
colorbar,ncolors=256 ,charthick=2.0,range=range_bar $
  ,/right,/vertical,charsize=1. $
  ,/YLOG,format='(e10.1)',color=col $
  ,pos=[17.3*dx_fig,9.*dy_fig,17.7*dx_fig,16.*dy_fig] $
  ,yticks=0

xyouts,17.2*dx_fig,16.3*dy_fig,textoidl('n_H [cm^{-3}]'),/normal,charthick=2.0

xyouts,17.2*dx_fig,8.3*dy_fig,textoidl(' T  [ K ]'),/normal,charthick=2.0 

loadct,3
;la temperatura
tvscl,alog10(temp),xsize=16,ysize=8.,1.,1.,/centimeters

range_bar=[minT,maxT]
colorbar,ncolors=256 ,charthick=2.0,range=range_bar $
  ,/right,/vertical,charsize=1. $
  ,/YLOG,format='(e10.1)',color=col $
  ,pos=[17.3*dx_fig,1.*dy_fig,17.7*dx_fig,8.*dy_fig] $
  ,yticks=0

endif

;-----------------------------------------------------------------------
if iffile eq 1 then begin
    device,/close
    set_plot,'X'
    spawn, 'epstopdf '+file_out+'&'  ; cooment
endif
;-----------------------------------------------------------------------
!P.multi=0
;rhosc=1.3
;psc=5./3.*1.66E-24*8.3145E7*1E4

;
;imagescbar,pgas[950:1050,0:100]*psc,/aspect,/bar,charsize=2.




ENDFOR

end
