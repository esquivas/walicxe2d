#!/usr/bin/python
import numpy as np
import struct
import scipy.ndimage


'''
  Reads mesh from Walicxe
'''
def read_mesh(nproc, nx, ny, nlevs, nout, xmax, ymax, 
              nbrootx, nbrooty, path='./',verbose = False):
  nbmin   = np.empty(nproc,dtype='i4')
  nbmax   = np.empty(nproc,dtype='i4')
  nbleafs = np.empty(nproc,dtype='i4')
  dx      = np.empty(nlevs,dtype='d')
  dy      = np.empty(nlevs,dtype='d')

  for nl in range(nlevs):
    dx[nl] = xmax/nbrootx/nx/(2.**nl)
    dy[nl] = ymax/nbrooty/ny/(2.**nl)

  for rank in range(nproc):
    file_in = path+'pointers'+str(rank).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
    if verbose: 
      print 'to read:', file_in
    f = open(file_in,'rb')
    nbmin[rank], nbmax[rank], nbleafs[rank], nblocks = struct.unpack('4i',f.read(16))
    f.close()

  ll = 0
  mesh = np.empty([np.sum(nbleafs),6],dtype='d')
  for rank in range(nproc):
    file_in = path+'pointers'+str(rank).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
    f = open(file_in,'rb')
    nbmin[rank], nbmax[rank], nbleafs[rank], nblocks = struct.unpack('4i',f.read(16))
    lp     = np.empty( [ (nbmax[rank] - nbmin[rank]+1)  , 8 ], dtype='i' )
    leaf   = np.empty( [ (nbmax[rank] - nbmin[rank]+1)      ], dtype='i' )
    coords = np.empty( [ (nbmax[rank] - nbmin[rank]+1)  ,2  ], dtype='i' )

    #  read lp's
    for i in range(8):
      #print i
      for nb in range(nbmax[rank]-nbmin[rank]+1):
        lp[nb,i] = struct.unpack('1i',f.read(4))[0]
    #  read leafs
    for nb in range(nbmax[rank]-nbmin[rank]+1):
      aux = struct.unpack('1i',f.read(4))[0]
      if (aux == -1 ):
        leaf[nb] = 1
      else:
        leaf[nb] = 0
    #  read icoords
    for i in range(2):
      for nb in range(nbmax[rank]-nbmin[rank]+1):
        coords[nb,i] = struct.unpack('1i',f.read(4))[0]
    
    f.close()

    for  nb in range(nbmax[rank]-nbmin[rank]+1):
      if leaf[nb] == 1 :
        lev = lp[nb,0]-1
        mesh[ll,0] = ( coords[nb,0]    ) * dx[lev]
        mesh[ll,1] = ( coords[nb,1]    ) * dy[lev]
        mesh[ll,2] = ( coords[nb,0]+nx ) * dx[lev]
        mesh[ll,3] = ( coords[nb,1]+ny ) * dy[lev]
        mesh[ll,4] = lp[nb,0]
        mesh[ll,5] = rank
        ll = ll + 1
    
  return(mesh)


'''
  Reads the NEQ(+1) conserved variable from walicxe
'''
def read_blocks(nproc, nx, ny, nlevs, neqs, nout, xmax, ymax, 
              nbrootx, nbrooty, path='./',verbose = False, double=True):
  nbmin   = np.empty(nproc,dtype='i4')
  nbmax   = np.empty(nproc,dtype='i4')
  nbleafs = np.empty(nproc,dtype='i4')
  dx      = np.empty(nlevs,dtype='d')
  dy      = np.empty(nlevs,dtype='d')
  for nl in range(nlevs):
    dx[nl] = xmax/nbrootx/nx/(2.**nl)
    dy[nl] = ymax/nbrooty/ny/(2.**nl)
  
  for rank in range(nproc):
    file_in = path+'pointers'+str(rank).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
    if verbose: 
      print 'to read:', file_in
    f = open(file_in,'rb')
    nbmin[rank], nbmax[rank], nbleafs[rank], nblocks = struct.unpack('4i',f.read(16))
    f.close()
  
  #  declare arrays to store the block
  if double:
    blocks = np.zeros([ny+4,nx+4,neqs,nbleafs.sum()],dtype='d')
  else:
    blocks = np.zeros([ny+4,nx+4,neqs,nbleafs.sum()],dtype='f')
  #  read all blocks
  ll = 0
  for rank in range(nproc):
    file_in = path+'blocks'+str(rank).zfill(3)+'.'+str(nout).zfill(3)+'.bin'
    f = open(file_in,'rb')
    for nb in range(nbleafs[rank]):
      for i in range(ny+4):
        for j in range(nx+4):
          for neq in range(neqs):
            if double:
              blocks[j,i,neq,ll] = struct.unpack('d',f.read(8))[0]
            else:
              blocks[j,i,neq,ll] = struct.unpack('f',f.read(4))[0]
      ll = ll + 1
    f.close()
  return(blocks)

'''
    constructs a 2D map of eq NEQ(+1), interpolated to max resolution
'''
def read_map(nproc, nx, ny, nlevs, neqs, neq, nout, xmax, ymax, 
            nbrootx, nbrooty, path='./',verbose = False, double=True):
  nxtot = nbrootx*nx*2**(nlevs-1)
  nytot = nbrooty*ny*2**(nlevs-1)
  dx    = xmax/nxtot
  dy    = ymax/nytot
  if double: 
    map = np.zeros([nxtot,nytot],dtype='d')
  else:
    map = np.zeros([nxtot,nytot],dtype='f')

  mesh = read_mesh(nproc,nx,ny,nlevs,nout,xmax,ymax,nbrootx,nbrooty,path=path)
  blocks = read_blocks(nproc,nx,ny,nlevs,neqs,nout,xmax,ymax,nbrootx,nbrooty,path=path)

  for nb in range(mesh.shape[0]):
    x0 = np.int(mesh[nb,0]/dx)
    y0 = np.int(mesh[nb,1]/dy)
    if (np.int(mesh[nb,4]) == nlevs):
      x1 = x0+nx
      y1 = y0+ny
      map[x0:x1,y0:y1] = blocks[2:ny+2,2:nx+2,neq,nb]
    else:
      factor = 2**(nlevs-np.int(mesh[nb,4]))
      x1 = x0 + factor*nx
      y1 = y0 + factor*ny
      A =  blocks[2:ny+2,2:nx+2,neq,nb]
      B = scipy.ndimage.zoom(A, factor, order=1)
      map[x0:x1,y0:y1] = B
  return(map.T)




