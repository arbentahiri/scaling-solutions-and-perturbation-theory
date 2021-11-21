# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 21:55:39 2021

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
    ps = 15
    r=0.01

    dxdN=-x *((ps ** 3* 
       x ** 7* (-6 *qs *x + 2.45* lam* y ** 2 + 
         2.45* Q *(-1 + omb + omr + qs* x**2 + y**2)) + 
      ps ** 2 *r *x ** 4* 
       y ** 2* (omr *(-2 + 2.45* (lam - 7* Q) *x) + 
         x *(30 *qs* x + 
            2.45* lam* (-3 + omb + 3 *qs *x ** 2 - 2 *y ** 2) - 
            2.45* Q *(-5 + 7 *omb + 5 *qs *x ** 2 + 5* y ** 2))) + 
      ps *r ** 2*x ** 2* 
       y ** 4* (omr *(8 + 3 *2.45* (-2 *lam + 5 *Q) *x) + 
         x *(-42 *qs* x + 
            2.45*lam*(6 - 6 *omb - 6 *qs *x ** 2 + y ** 2) + 
            2.45* Q *(-7 + 15 *omb + 7* qs *x ** 2 + 7 *y ** 2))) - 
      3 *r ** 3* 
       y ** 6*(omr *(2 + 3 *2.45* (-lam + Q) *x) + 
         x *(-6 *qs* x + 2.45* (-lam + Q) *qs *x ** 2 + 
            2.45* (lam - 3 *lam* omb + 
               Q *(-1 + 3 *omb + y ** 2)))))/(2 *ps ** 3*qs *x ** 8 - 
      6 *r ** 3* 
       y ** 6* (-1 + 3 *omb + 3 *omr + 2 *qs* x ** 2 + y ** 2) - 
      2 *ps ** 2*r *x ** 4*
       y ** 2* (-3 + omb + omr + 8 *qs *x ** 2 + 3 *y ** 2) + 
      2 *ps *r ** 2* x ** 2* 
       y ** 4* (-6 + 6 *omb + 6 *omr + 13 *qs* x ** 2 + 6 *y ** 2)) + 
   1/(2 *ps *x ** 2 - 
     6 *r *y ** 2) *(6* r *y ** 4 - 
      3* r* y ** 
        2*(2 + 3 *omb + 4* omr + 4 *qs *x ** 2 + 
         3 *(1 - omb - 
            omr - (-2* r *
              y ** 2 + (qs *x ** 2 + y ** 2) *(ps* x ** 2 - 
                r *y ** 2))/(ps *x ** 2 - 3 *r *y ** 2))) + 
      ps *x ** 
        2* (3 *omb + 4 *omr + 6 *qs* x ** 2 + 
         3 *(1 - omb - 
            omr - (-2* r *
              y ** 2 + (qs *x ** 2 + y ** 2) *(ps *x ** 2 - 
                r *y ** 2))/(ps *x ** 2 - 3 *r *y ** 2)))))
                    
    dydN=y *(1/(2 *ps* x ** 2 - 
     6 *r *y ** 2)*(-6 *r* y ** 4 + 
      3 *r* y ** 
        2* (2 + 3* omb + 4 *omr + 4 *qs* x ** 2 + 
         3* (1 - omb - 
            omr - (-2 *r* 
              y ** 2 + (qs *x ** 2 + y ** 2) *(ps* x ** 2 - 
                r *y ** 2))/(ps* x**2 - 3* r* y**2))) - 
      ps *x ** 
        2* (3* omb + 4 *omr + 6 *qs* x ** 2 + 
         3 *(1 - omb - 
            omr - (-2 *r *
              y ** 2 + (qs* x ** 2 + y ** 2) *(ps *x ** 2 - 
                r *y ** 2))/(ps *x ** 2 - 3 *r* y ** 2)))) + 
   1/2 *lam* x *2.45)
    
    dombdN=-omb*(3*qs*x**2-3*y**2+omr)
    
    domrdN=-omr*(3*qs*x**2-1-3*y**2+omr)
    return [dxdN,dydN]



N=np.linspace(0,-5,1000000)

s0=[3.75455599e-02,0.02]
s1=odeint(background,s0,N) 
x1=s1[:,0]
y1=s1[:,1]

s0=[-0.535,0.05575]
s2=odeint(background,s0,N) 
x2=s2[:,0]
y2=s2[:,1]

s0=[0.35,0.0575]
s3=odeint(background,s0,N) 
x3=s3[:,0]
y3=s3[:,1]

s0=[0.115,0.0575]
s4=odeint(background,s0,N) 
x4=s4[:,0]
y4=s4[:,1]

s0=[0.175,0.0625]
s5=odeint(background,s0,N) 
x5=s5[:,0]
y5=s5[:,1]

s0=[-0.255,0.0575]
s6=odeint(background,s0,N) 
x6=s6[:,0]
y6=s6[:,1]

s0=[0.565,0.0375]
s7=odeint(background,s0,N) 
x7=s7[:,0]
y7=s7[:,1]
   
s0=[-0.655,0.075]
s8=odeint(background,s0,N) 
x8=s8[:,0]
y8=s8[:,1]


plt.xlim([-1, 1])
plt.ylim([0,1.2])
plt.xlabel("$x$")
plt.ylabel("$y$")


plt.plot(x1,y1)
plt.plot(x2,y2)
plt.plot(x3,y3)
plt.plot(x4,y4)
plt.plot(x5,y5)
plt.plot(x6,y6)
plt.plot(x7,y7)
plt.plot(x8,y8)