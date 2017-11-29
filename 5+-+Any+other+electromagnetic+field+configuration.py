
# coding: utf-8

# In[2]:

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy.constants import e, c, epsilon_0, m_e
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
# This command to import some physicals constants


# In[3]:

N = 10 ** 4
B0 = 1           # Magnetic field
E0 = 1          # Electric field
m_e = 1
m_p = 3
ee = -1
ep = 1


# In[4]:

x = np.zeros((N,3))
v = np.zeros((N,3))
v[0,:] = np.array([1,0,0])
#E = np.zeros((N,3))
B = np.array([0, 0, B0])
E = np.array([E0, 0, 0])
dt = ( 0.1 * m_e ) / (np.abs(ee) * B0)
time = np.linspace(0, dt*N)


# In[6]:

mass = [m_e, m_p]
charge = [ee, ep]
fig, ax = plt.subplots(dpi = 100)
alpha = 3
for mass, charge in zip(mass,charge):
    x = np.zeros((N,3))
    x[0,:] = np.array([1,1,0])
    v = np.zeros((N,3))
    B = np.zeros((N,3))
    v[0,:] = np.array([1,0,0])
    E = np.array([0, 0, 0])
    B[0,:] = np.array([0, 0, 0])
    dt = ( 0.1 * mass ) / (np.abs(charge) * B0)
    time = np.linspace(0, dt*N, N)
    
    for i in range(N-1):
        x[i+1,:] = x[i,:] + v[i,:] * dt
        #E[i+1,0] = E0 * time[i+1]
        B[i+1,:] = np.array([0,0,1+alpha * np.sqrt(x[i,0]**2 + x[i,1]**2)])
        v_minus = v[i,:] + (charge/mass) * E * (dt / 2)
        t = (charge * dt) / (2 * mass) * B[i]
        v_prime = v_minus + np.cross(v_minus, t)
        s = 2 / ( 1 + np.abs(t) **2) * t
        v_plus = v_minus + np.cross(v_prime,s)
        v[i+1,:] = v_plus + (charge/mass) * E * ( dt / 2)
        
    plt.plot(x[:,0],x[:,1], label='charge = %g' %charge)
    plt.plot(B[:,0],B[:,1], label='magnetic')
    plt.xlabel('x-axis')
    plt.ylabel('y-axis')
    plt.grid(True)
    plt.title('Magnetic mirror confinement')
plt.legend(loc=2, prop={'size': 6})
plt.show()


# In[48]:

print(x[:,1])


# In[31]:

mass = [m_e, m_p]
charge = [ee, ep]
fig, ax = plt.subplots(dpi = 100)
alpha = 3
for mass, charge in zip(mass,charge):
    x = np.zeros((N,3))
    x[0,:] = np.array([0,0,0])
    v = np.zeros((N,3))
    B = np.zeros((N,3))
    v[0,:] = np.array([1,0,0])
    E = np.array([0, 0, 0])
    B[0,:] = np.array([0, 0, 0])
    dt = ( 0.1 * mass ) / (np.abs(charge) * B0)
    time = np.linspace(0, dt*N, N)
    
    for i in range(N-1):
        x[i+1,:] = x[i,:] + v[i,:] * dt
        #E[i+1,0] = E0 * time[i+1]
        B[i+1,:] = np.array([0,0,x[i,0]])
        v_minus = v[i,:] + (charge/mass) * E * (dt / 2)
        t = (charge * dt) / (2 * mass) * B[i]
        v_prime = v_minus + np.cross(v_minus, t)
        s = 2 / ( 1 + np.abs(t) **2) * t
        v_plus = v_minus + np.cross(v_prime,s)
        v[i+1,:] = v_plus + (charge/mass) * E * ( dt / 2)
        
    plt.plot(x[:,0],x[:,1], label='charge = %g' %charge)
    plt.xlabel('x-axis')
    plt.ylabel('y-axis')
    plt.grid(True)
    plt.title('Magnetic mirror confinement')
plt.legend(loc=2, prop={'size': 6})
plt.show()


# In[32]:

mass = [m_e, m_p]
charge = [ee, ep]
fig, ax = plt.subplots(dpi = 100)
fig.patch.set_facecolor('white')
ax = fig.add_subplot(111, projection='3d')
alpha = 3
for mass, charge in zip(mass,charge):
    x = np.zeros((N,3))
    x[0,:] = np.array([1,1,0])
    v = np.zeros((N,3))
    B = np.zeros((N,3))
    v[0,:] = np.array([1,0,1])
    E = np.array([0, 0, 0])
    B[0,:] = np.array([0, 0, 0])
    dt = ( 0.1 * mass ) / (np.abs(charge) * B0)
    time = np.linspace(0, dt*N, N)
    
    for i in range(N-1):
        x[i+1,:] = x[i,:] + v[i,:] * dt
        #E[i+1,0] = E0 * time[i+1]
        B[i+1,:] = np.array([0,0,1- alpha * np.sqrt(x[i,0]**2 + x[i,1]**2)])
        v_minus = v[i,:] + (charge/mass) * E * (dt / 2)
        t = (charge * dt) / (2 * mass) * B[i]
        v_prime = v_minus + np.cross(v_minus, t)
        s = 2 / ( 1 + np.abs(t) **2) * t
        v_plus = v_minus + np.cross(v_prime,s)
        v[i+1,:] = v_plus + (charge/mass) * E * ( dt / 2)
        
    ax.plot(x[:,0], x[:,1], x[:,2], label='mass = %g, charge = %g' %(mass, charge))
ax.legend(loc=2, prop={'size': 6})
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()


# In[ ]:



