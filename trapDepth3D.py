# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 18:26:57 2020

@author: lizaz
"""
from scipy.signal import argrelextrema, argrelmax, argrelmin, find_peaks
import matplotlib.pyplot as plt, numpy as np, scipy.constants as ct
from electropy.electrode_layouts.ring_electrodes import tworing_electrode
np.set_printoptions(precision=2) 

def CreatePotentialMesh(s, nx=5, ny=7, nz=9, x = [-2,2], y = [-2,2], z = [0.1,4.1]):
    
    xmin, xmax = x
    ymin, ymax = y
    zmin, zmax = z
    
    x = np.linspace(xmin,xmax,nx) 
    y = np.linspace(ymin,ymax,ny) 
    z = np.linspace(zmin,zmax,nz)    
    
    yy,xx,zz = np.meshgrid(y,x,z)
        
    p = np.zeros([nx,ny,nz])
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                p[i,j,k] = s.potential([x[i],y[j],z[k]])[0]
    
    return xx, yy, zz, p

def trapDepth3D(xx, yy, zz, p, plotting=False, mode = 'J', q=1, m=1, o=1, l=1, u=1):

    x = xx[:,0,0]
    y = yy[0,:,0]
    z = zz[0,0,:]

    nx = len(x)
    ny = len(y)
    nz = len(z)
    
    i,j,k= argrelmin(p)
    
    try:    
        idx_min = np.argmin(p[i,j,k])
   
    except ValueError:
        return 0,0,0
        
        pass

    x_min = x[i[idx_min]]
    y_min = y[j[idx_min]]
    z_min = z[k[idx_min]]

    p_resh =  p.reshape(nx*ny*nz)
    x_resh = xx.reshape(nx*ny*nz)
    y_resh = yy.reshape(nx*ny*nz)
    z_resh = zz.reshape(nx*ny*nz)
    idx = find_peaks(p_resh)[0]

    idx_max = np.argmin(p_resh[idx])
    x_max = x_resh[idx[idx_max]]
    y_max = y_resh[idx[idx_max]]
    z_max = z_resh[idx[idx_max]]
     
    if plotting == True:
         
     plt.figure()
     plt.contourf(yy[int((nx-1)/2),:,:],zz[int((nx-1)/2),:,:],p[int((nx-1)/2),:,:],50)
     plt.plot(y_max,z_max,'o')

     plt.figure()
     plt.contourf(xx[:,int((ny-1)/2),:],zz[:,int((ny-1)/2),:],p[:,int((ny-1)/2),:],50)
     plt.plot(x_max,z_max,'o')
     fig, axs = plt.subplots(1,3)
     
     axs[0].plot(z,p[int((nx-1)/2),int((ny-1)/2),:]/u**2)
     axs[0].plot(z_min,p[int((nx-1)/2),int((ny-1)/2),k[idx_min]]/u**2,'o')
     axs[0].plot(z_max,p[int((nx-1)/2),int((ny-1)/2), np.where(z==z_max)]/u**2,'o')
     axs[0].set_ylim(0,0.05)          
       
     axs[1].plot(y,p[int((nx-1)/2),:,int((nz-1)/2)]/u**2)
     axs[1].plot(y_min,p[int((nx-1)/2),j[idx_min],int((nz-1)/2)]/u**2,'o')

     axs[2].plot(x,p[:,int((ny-1)/2),int((nz-1)/2)]/u**2)
     axs[2].plot(x_min,p[i[idx_min],int((ny-1)/2),int((nz-1)/2)]/u**2,'o')

    D = (p[np.where(x==x_max),np.where(y==y_max),np.where(z==z_max)] - p[np.where(x==x_min),np.where(y==y_min),np.where(z==z_min)])[0][0]*(q ** 2) / (4 * m * o ** 2 * l ** 2)
    delta_xyz = np.sqrt((x_max-x_min)**2+(y_max-y_min)**2+(z_max-z_min)**2)
    xyz_min = [x_min, y_min, z_min]
    
    if D<0:
        
        return 0,0,0
    
    else:
    
        if D/ct.k > (3/2*ct.k*300/ct.elementary_charge):
       
            if mode == 'J':
                return D, delta_xyz, xyz_min #output in Joules
            else:
                return D/ct.k, delta_xyz, xyz_min #output in Kelvins

        else:
            print('The particle cannot be trapped.')
            return D/ct.k, delta_xyz, xyz_min 
 
