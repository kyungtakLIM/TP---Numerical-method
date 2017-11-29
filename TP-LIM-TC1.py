
# coding: utf-8

# In[94]:

import numpy as np
import matplotlib.pyplot as plt
import csv
import time
import numpy.linalg as la

from pylab import *
from scipy.constants import e, c, epsilon_0, m_e
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3


# In[95]:

# mesh

Nx = 100
Ny = 100
dx = 1.0/Nx
dy = 1.0/Ny
print("mesh")
print("Nx="),Nx,("Ny="),Ny


# In[96]:

# ADD HERE THE DEFINITION OF X and XM 	
# X = cell center / Xm = grid points
Xm = [ dx * i for i in range(Nx)]
Ym = [ dx * j for j in range(Ny)]
#print (Xm)

X = []
Y = []
for i in range(Nx-1):
    X.append( (Xm[i] + Xm[i+1]) / 2 )
    Y.append( (Ym[i] + Ym[i+1]) / 2 )


# In[97]:

# Potentials on the anode and cathode
Ua = 100       #        in V (right side)
Uc = 10        # (left side)

eps0 = 8.85E-14   # F/cm


# In[98]:

# dielectric layers
epsdiel=1    # default =1 (no dielectric)

nxdg = Nx/10                    # dielectric anode side
nxdd = Nx/10                    #dielectric cathode side  


# In[99]:

# definition of epsilon
# Size of system is 50, in two end points [ 0 , nxdg ] and [Nx - nxdd ,Nx ] we use relative permittivity
# in the plasma regime, we utilise permittivity vacuum
def E(i,j):
    if i < nxdg or i > Nx - nxdd:
        return eps0*epsdiel
    else:
        return eps0


# In[100]:

# Coefficients

def Ve(i,j):
    Ei=E(i,j)
    EI=E(i+1,j)
    if i == Nx-1:    # Boundary condition right side
        return 2*Ei/(dx*dx)
    else:
        return 2*Ei*EI/(dx*dx*(Ei+EI))

def Vo(i,j):
    Ei=E(i,j)
    EI=E(i-1,j)
    if i == 0:
        return 2*Ei/(dx*dx)   # Boundary condition left side
    else:
        return 2*Ei*EI/(dx*dx*(Ei+EI))

def Vn(i,j):
    Ej=E(i,j)
    EJ=E(i,j+1)
    if j == Ny-1:
        return 0
    else:
        return 2*Ej*EJ/(dy*dy*(Ej+EJ))

def Vs(i,j):
    Ej=E(i,j)
    EJ=E(i,j-1)
    if j == 0:
        return 0
    else:
        return 2*Ej*EJ/(dy*dy*(Ej+EJ))

def Vc(i,j):
    return Ve(i,j)+Vo(i,j)+Vn(i,j)+Vs(i,j)

def rho(i,j):
    if i == 0 :
        return -2*E(i,j)/(dx*dx)*Uc
    else:
        if i == Nx-1:
            return -2*E(i,j)/(dx*dx)*Ua
        else:
            return 0

def V(i,j):
    if i == -1 or i == Nx or j == -1 or j == Ny:
        return 0
    else:
        return Vact[i,j]


# In[101]:

t2=time.clock()

Vact = np.zeros((Nx,Ny))

# CODE HERE THE GAUSSE-SEIDEL AND/OR SOR


N_iteration = 2000


for k in range(N_iteration):
    for i in range(Nx-1):
        for j in range(Ny-1):
            Vact[i,j] = ( Ve(i,j)*V(i+1,j) + Vo(i,j) * V(i-1,j) + Vn(i,j) * V(i,j+1) + Vs(i,j) * V(i,j-1) - rho(i,j) ) / Vc(i,j)
# END OF THE GAUSSE-SEIDEL AND/OR SOR

t2=time.clock()-t2


# In[102]:

# FIGURES
fig, ax = plt.subplots(dpi = 150)
contour_plot = plt.contourf(np.linspace(dx/2. , 1.-dx/2., Nx), np.linspace( dx/2. , 1.-dx/2. , Ny), np.transpose(Vact), 100, extend='both');
plt.colorbar(contour_plot);
plt.xlabel("x [cm]");
plt.ylabel("y [cm]");
plt.grid()
plt.savefig("V_xy_%g.png" %N_iteration);
plt.show();

fig, ax = plt.subplots(dpi = 150)
plt.plot(np.linspace(dx/2.,1.-dx/2., Nx),Vact[:,Ny/2],'ro');
plt.xlabel("x [cm]");
plt.ylabel("V [V]");
plt.grid()
plt.savefig("V_x_%g.png" %N_iteration);
plt.show();
plt.clf();


# In[103]:

# WRITE HERE THE DERIVATION OF THE ELECTRIC FIELD AND PLOT IT
Ex_field = np.zeros((Nx,Ny))
Ey_field = np.zeros((Nx,Ny))
for i in range(Nx):
    for j in range(Ny):
        Ex_field[i,j] = -2 * E(i+1,j) * (V(i+1,j) - V(i,j)) / (dx * (E(i,j) + E(i+1,j) ))

for i in range(Nx):
    for j in range(Ny):
        Ey_field[i,j] = -2 * E(i,j+1) * (V(i,j+1) - V(i,j)) / (dy * (E(i,j) + E(i,j+1) ))


# In[104]:

fig, ax = plt.subplots(dpi = 150)
ax.imshow(Ex_field)
plt.xlabel("x [cm]");
plt.ylabel("y [cm]");
plt.title("Electricfield in $x$-direction")
plt.savefig("Ex_field_%g.png" %N_iteration);
plt.show()

fig, ax = plt.subplots(dpi = 150)
plt.xlabel("x [cm]");
plt.ylabel("y [cm]");
plt.title("Electricfield in $y$-direction")
plt.savefig("Ey_field_%g.png" %N_iteration);
plt.imshow(Ey_field)
plt.show()


# In[ ]:



