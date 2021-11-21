# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 20:07:49 2021

@author: Beni
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 23:17:37 2021

@author: Beni
"""

#background equations 

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


omr=0
omc=0
omb=0

def background(s,N):
    x=s[0]
    y=s[1]
  

    lam = 1
    Q = 0.07
    qs = 2


    dxdN=-((1/2)*x*(3*qs*x**2-3-3*y**2+omr)
           +2.45/(2*qs)*(lam*y**2-Q*(1-qs*x**2-y**2-omr-omb)))
   
    dydN=-(1/2)*y*(3*qs*x**2-x*lam*2.45-3*y**2+omr+3)
   
   

    return [dxdN,dydN]


N=np.linspace(0,-5,1000000)

s0=[-0.595,0.0475]
s1=odeint(background,s0,N) 
x1=s1[:,0]
y1=s1[:,1]

s0=[-0.535,0.0575]
s2=odeint(background,s0,N) 
x2=s2[:,0]
y2=s2[:,1]

s0=[-0.31,0.0575]
s3=odeint(background,s0,N) 
x3=s3[:,0]
y3=s3[:,1]

s0=[-0.115,0.0575]
s4=odeint(background,s0,N) 
x4=s4[:,0]
y4=s4[:,1]

s0=[0.075,0.0625]
s5=odeint(background,s0,N) 
x5=s5[:,0]
y5=s5[:,1]

s0=[0.255,0.0575]
s6=odeint(background,s0,N) 
x6=s6[:,0]
y6=s6[:,1]

s0=[0.565,0.0375]
s7=odeint(background,s0,N) 
x7=s7[:,0]
y7=s7[:,1]

   


plt.xlim([-1, 1])
plt.ylim([0,1.2])
plt.xlabel("x")
plt.ylabel("$y$")


plt.plot(x1,y1)
plt.plot(x2,y2)
plt.plot(x3,y3)
plt.plot(x4,y4)
plt.plot(x5,y5)
plt.plot(x6,y6)
plt.plot(x7,y7)

