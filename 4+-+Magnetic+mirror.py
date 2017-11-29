
# coding: utf-8

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy.constants import e, c, epsilon_0, m_e
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
# This command to import some physicals constants


# In[7]:

N = 2*10 ** 3
B0 = 1           # Magnetic field
E0 = 3          # Electric field
m_e = 1
m_p = 3
ep = 1
ee = -1


# In[8]:

mass = [m_e]
charge = [ee]
fig, ax = plt.subplots(dpi = 150)
Length = 30 # mirror point
for mass, charge in zip(mass,charge):
    x = np.zeros((N,3))
    x[0,:] = np.array([-3*Length,-3*Length,-3*Length])
    #x[0,:] = np.array([0,0,0])
    v = np.zeros((N,3)) 
    B = np.zeros((N,3))
    v[0,:] = np.array([1,1,0])
    E = np.array([0, 0, 0])
    B[0,:] = np.array([B0, 0, 0])
    dt = ( 0.1 * mass ) / (np.abs(charge) * B0)
    time = np.linspace(0, dt*N, N)
    
    for i in range(N-1):
        x[i+1,:] = x[i,:] + v[i,:] * dt
        B[i+1,0] = B0 * ( 1 + ( x[i+1,0] / Length ) ** 2 )
        v_minus = v[i,:] + (charge/mass) * E * (dt / 2)
        t = (charge * dt) / (2 * mass) * B[i+1]
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


# In[9]:

mass = [m_e, m_p]
charge = [ee, ep]
fig, ax = plt.subplots(dpi = 100)
fig.patch.set_facecolor('white')
ax = fig.add_subplot(111, projection='3d')
ax.plot(x[:,0], x[:,1], x[:,2])
ax.legend(loc=2, prop={'size': 6})
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()


# In[10]:

mass = [m_e]
charge = [ee]
fig, ax = plt.subplots(dpi = 150)
Length = 30 # mirror point
for mass, charge in zip(mass,charge):
    x = np.zeros((N,3))
    x[0,:] = np.array([-3*Length,-3*Length,-3*Length])
    #x[0,:] = np.array([0,0,0])
    v = np.zeros((N,3)) 
    B = np.zeros((N,3))
    v[0,:] = np.array([1,1,0])
    E = np.array([0, 0, 0])
    B[0,:] = np.array([B0, 0, 0])
    dt = ( 0.1 * mass ) / (np.abs(charge) * B0)
    time = np.linspace(0, dt*N, N)
    
    for i in range(N-1):
        x[i+1,:] = x[i,:] + v[i,:] * dt
        B[i+1,0] = 2*B0 * ( 1 + ( x[i+1,0] / Length ) ** 2 )
        v_minus = v[i,:] + (charge/mass) * E * (dt / 2)
        t = (charge * dt) / (2 * mass) * B[i+1]
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


# In[11]:

mass = [m_e, m_p]
charge = [ee, ep]
fig, ax = plt.subplots(dpi = 100)
fig.patch.set_facecolor('white')
ax = fig.add_subplot(111, projection='3d')
ax.plot(x[:,0], x[:,1], x[:,2])
ax.legend(loc=2, prop={'size': 6})
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()


# In[14]:

mass = [m_e]
charge = [ee]
fig, ax = plt.subplots(dpi = 150)
Length = 30 # mirror point
for mass, charge in zip(mass,charge):
    x = np.zeros((N,3))
    x[0,:] = np.array([-3*Length,-3*Length,-3*Length])
    #x[0,:] = np.array([0,0,0])
    v = np.zeros((N,3)) 
    B = np.zeros((N,3))
    v[0,:] = np.array([1,1,0])
    E = np.array([0, 0, 0])
    B[0,:] = np.array([B0, 0, 0])
    dt = ( 0.1 * mass ) / (np.abs(charge) * B0)
    time = np.linspace(0, dt*N, N)
    
    for i in range(N-1):
        x[i+1,:] = x[i,:] + v[i,:] * dt
        B[i+1,0] = 0.5*B0 * ( 1 + ( x[i+1,0] / Length ) ** 2 )
        v_minus = v[i,:] + (charge/mass) * E * (dt / 2)
        t = (charge * dt) / (2 * mass) * B[i+1]
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


# In[15]:

mass = [m_e, m_p]
charge = [ee, ep]
fig, ax = plt.subplots(dpi = 100)
fig.patch.set_facecolor('white')
ax = fig.add_subplot(111, projection='3d')
ax.plot(x[:,0], x[:,1], x[:,2])
ax.legend(loc=2, prop={'size': 6})
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()


# In[5]:

mass = [m_e, m_p]
charge = [ee, ep]
fig, ax = plt.subplots(dpi = 100)
fig.patch.set_facecolor('white')
ax = fig.add_subplot(111, projection='3d')
ax.plot(x[:,0], x[:,1], x[:,2])
ax.legend(loc=2, prop={'size': 6})
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()


# In[ ]:



