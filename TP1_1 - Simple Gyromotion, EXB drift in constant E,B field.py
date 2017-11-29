
# coding: utf-8

# In[14]:

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy.constants import e, c, epsilon_0, m_e
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
# This command to import some physicals constants


# In[15]:

"""
This program is to calculate the method particle in cell.
The goal is to determine the motion of a charged particle of mass m and
charge q in a given electric ( E ) andmagnetic ( B )
eld.

INPUT :  electric field, magnetic field
OUTPUT :  Velocity, position, Acceleration

Date of coding : 17 / 11 / 2017
Modification of coding : 17 / 11 / 2017
Programmer : lim alice
"""


# In[16]:

N = 10 **3
B0 = 1             # Magnetic field
E0 = 1             # Electric field
m_e = 1
m_p = 3
ee = -1
ep = 1


# ### Simple Gyro-motion

# In[17]:

x = np.zeros((N,3))
v = np.zeros((N,3))
v[0,:] = np.array([1,1,1])
B = np.array([0, 0, B0])
E = np.array([0, 0, 0])
dt = ( 0.1 * m_e ) / (np.abs(ee) * B0)
time = np.linspace(0, dt*N, N)


# In[18]:

for i in range(N-1):
    x[i+1,:] = x[i,:] + v[i,:] * dt
    v_minus = v[i,:] + (ee/m_e) * E * (dt / 2)
    t = (np.abs(ee) * dt) / (2 * m_e) * B
    v_prime = v_minus + np.cross(v_minus, t)
    s = 2 / ( 1 + np.abs(t) **2) * t
    v_plus = v_minus + np.cross(v_prime,s)
    v[i+1,:] = v_plus + (ee/m_e) * E * ( dt / 2)


# In[19]:

fig, ax = plt.subplots(dpi = 100)
plt.plot(x[:,0],x[:,1])
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.grid(True)
plt.title('Derivation')
plt.show()


# In[20]:

fig, ax = plt.subplots(dpi = 100)
fig.patch.set_facecolor('white')
ax = fig.add_subplot(111, projection='3d')
ax.plot(x[:,0], x[:,1], x[:,2], label='parametric curve')
ax.legend()
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()


# ### E X B drift

# In[21]:

N = 10 **3
B0 = 3             # Magnetic field
E0 = 1             # Electric field
m_e = 1
m_p = 3
ee = -1
ep = 1


# In[22]:

x = np.zeros((N,3))
v = np.zeros((N,3))
v[0,:] = np.array([1,1,1])
B = np.array([0, 0, B0])
E = np.array([0, 0, 0])
dt = ( 0.1 * m_e ) / (np.abs(ee) * B0)
time = np.linspace(0, dt*N, N)


# In[23]:

mass = [m_e, m_p]
charge = [ee, ep]
fig, ax = plt.subplots(dpi = 100)
for mass, charge in zip(mass,charge):
    x = np.zeros((N,3))
    v = np.zeros((N,3))
    v[0,:] = np.array([1,1,1])
    B = np.array([B0, 0, 0])
    E = np.array([0, E0, 0])
    dt = ( 0.1 * mass ) / (np.abs(charge) * B0)
    time = np.linspace(0, dt*N, N)
    for i in range(N-1):
        x[i+1,:] = x[i,:] + v[i,:] * dt
        v_minus = v[i,:] + (charge/mass) * E * (dt / 2)
        t = (charge * dt) / (2 * mass) * B
        v_prime = v_minus + np.cross(v_minus, t)
        s = 2 / ( 1 + np.abs(t) **2) * t
        v_plus = v_minus + np.cross(v_prime,s)
        v[i+1,:] = v_plus + (charge/mass) * E * ( dt / 2)
        
    plt.plot(x[:,0],x[:,1], label='charge = %g' %charge)
    plt.xlabel('x-axis')
    plt.ylabel('y-axis')
    plt.grid(True)
    plt.title('E X B drift')
plt.legend(loc=2, prop={'size': 6})
plt.show()


# In[26]:

fig, ax = plt.subplots(dpi = 100)
fig.patch.set_facecolor('white')
ax = fig.add_subplot(111, projection='3d')
ax.plot(x[:,0], x[:,1], x[:,2], label='E X B drift in 3D')
ax.legend()
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()


# In[ ]:



