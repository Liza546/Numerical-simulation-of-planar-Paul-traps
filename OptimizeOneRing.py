# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 14:44:06 2020

@author: lizaz
"""
import matplotlib.pyplot as plt, numpy as np, scipy.constants as ct
from electropy.electrode_layouts.ring_electrodes import onering_electrode
from electropy.postprocessing.trapDepth3D import CreatePotentialMesh, trapDepth3D
np.set_printoptions(precision=2) 



def OptimizeOneRing(r, R , q=1, m=1, o=1, l=1, u=1):
    Optimum = -1
    trapDepth = []
    U = np.zeros([len(r),len(R)])
    for i in range(len(r)):
        # print('i=',i)
        for j in range(len(R)):
            # print('j=',j)
            if r[i] < R[j]:
                s = onering_electrode(r=r[i], R=R[j],segments=100)
                xx,yy,zz,p = CreatePotentialMesh(s, nx=5, ny=7, nz=9, x = [-0.4,0.4], y = [-0.4,0.4], z = [0.1,3])
                tD, delta_xyz, xyz_min = trapDepth3D(xx,yy,zz,p, plotting = False, mode='J', q=q, m=m, o=o, l=l, u=u)
                U[i][j] = tD
              
                if Optimum < U[i][j]:
                    Optimum = U[i][j]
                    r_opt = r[i]
                    R_opt = R[j]
                    
            else: 
                U[i][j] = 0
               
    return U, Optimum, r_opt, R_opt      

# <codecell>
l=1e-3
R_in = np.linspace(0.01,1,10)
R_out = np.linspace(0.01,2,19)
U, Optimum, r_opt, R_opt = OptimizeOneRing(r = R_in, R = R_out, q=10*ct.elementary_charge, m=3.77e-16, o=10*1e3, l=1e-3, u=1e3 )
rr_outer, rr_inner = np.meshgrid(R_out*l,R_in*l)
# <codecell>

plt.contourf(rr_outer,rr_inner,U,levels=20)
plt.plot(R_opt*l,r_opt*l, 'bx')
cbar=plt.colorbar()
cbar.set_label('Trap depth [J]')
plt.xlabel('r_outer [mm]')
plt.ylabel('r_inner [mm]')
plt.title('Optimized one-ring trap (res=1um, l=1e-3)')  