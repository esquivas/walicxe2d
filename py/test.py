import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from walicxe2d_utils import *

path = '/Users/esquivel/Desktop/diable/data/esquivel/Walixce-2D/LE-JET/'
runs = ['Tau40AD','Tau40ISO','Tau40DMC','Tau40AD-2']

nproc = 16
nx    = 32
ny    = 32
nlevs = 7
neqs  = 5
xmax  = 4.#6e17
ymax  = 1.#1.5e17
nbx   = 4
nby   = 1
nxtot = 4.*nx*2**(nlevs-1)
dx    = xmax/nxtot

firstrun = False

rhomin= 10. 
rhomax= 400.

tau = 40.*3.156e7
v0  = 200e5

plt.ion()

for ii in range(np.size(runs)):

	if firstrun:
		pos  = np.zeros(200,dtype = 'd')
		time = np.zeros(200,dtype = 'd')
		flag = True

		for nout in range(0,200):

			#  read the mesh and get a density (neq = 0), map
			mesh = read_mesh(nproc,nx,ny,nlevs,nout,xmax,ymax,nbx,nby,path=path+runs[ii]+'/BIN/')
			rho = read_map(nproc, nx, ny, nlevs, neqs, 0, nout, xmax, ymax, nbx, nby, path=path+runs[ii]+'/BIN/')

			plt.figure(ii+2, figsize=(10,2.5))
			plt.clf()
			plt.axes().set_aspect('equal')

			#   density plot
			plt.imshow(rho,origin='lower',interpolation='none',cmap='plasma',
			  norm=LogNorm(),extent=[0.,xmax,0.,ymax], vmin = rhomin, vmax=rhomax)
			plt.title(runs[ii])
			
			#   add the mesh as an overlay
			for nb in range(mesh.shape[0]):
			  plt.plot([ mesh[nb,0],mesh[nb,0] ],[ mesh[nb,1],mesh[nb,3] ],color='dimgray',linestyle='solid',alpha=0.25, linewidth=1.)
			  plt.plot([ mesh[nb,0],mesh[nb,2] ],[ mesh[nb,1],mesh[nb,1] ],color='dimgray',linestyle='solid',alpha=0.25, linewidth=1.)
			  plt.plot([ mesh[nb,2],mesh[nb,0] ],[ mesh[nb,3],mesh[nb,3] ],color='dimgray',linestyle='solid',alpha=0.25, linewidth=1.)
			  plt.plot([ mesh[nb,2],mesh[nb,2] ],[ mesh[nb,1],mesh[nb,3] ],color='dimgray',linestyle='solid',alpha=0.25, linewidth=1.)
			
			#.  add colorbar
			plt.colorbar(shrink=0.8)

			#.  save to pdf
			#plt.savefig('fig'+str(nout).zfill(3)+'.pdf', transparent=True, bbox_inches='tight')

			Lx = rho.shape[1]
			i = Lx-1
			while flag:
				if np.any(rho[:,i] > 110.):
					print 'position of the head', nout, i
					print  'time:', 5.*nout ,'pos:', i*dx*6e17/4
					pos[nout] = i*dx*6e17/4
					time[nout] = 5.*3.156e7*nout
					flag = False
				i = i - 1
				if (i==0): flag = False
			flag= True

		np.savez('Rvst-'+runs[ii]+'.npz',pos, time)

	else:
		plt.figure(0)
		if (ii == 0) : plt.clf()
		data = np.load('Rvst-'+runs[ii]+'.npz')
		time = data['arr_1']/tau
		pos  = data['arr_0']/v0/tau
		plt.plot(time,pos,label=runs[ii], linewidth = 2)
		plt.xlabel(r'$t/\tau$')
		plt.ylabel(r'$x_h/v_0\,\tau$')
		vel  = np.zeros(200,dtype = 'd')
		for jj in range(1,200):
				vel[jj]=(pos[jj]-pos[jj-1])/.125
		plt.figure(1)
		if (ii == 0) : plt.clf()
		plt.plot(time,vel,label=runs[ii], linewidth = 2)
		plt.xlabel(r'$t/\tau$')
		plt.ylabel(r'$v_h/v_0$')

		plt.figure(0)  ;  plt.legend(loc='lower right')
		plt.figure(1)  ;  plt.legend(loc='lower right')
