# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 05:33:02 2021

@author: Beni
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math

#x0=3.99082578e-05
#y0=6.68179879e-04
#omb0=9.04827903e-02
#omr0=0.04490344

x0=1e-13
y0=1/9*x0
omb0=3.6e-6
omr0=0.999962



def background(s,N):
    x=s[0]
    y=s[1]
    omb=s[2]
    omr=s[3]
  

    lam = 1.0
    Q = 0.07
    qs = 2
    ps=15
    r=-0.01

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

    return [dxdN,dydN,dombdN,domrdN]


s0=[x0, y0, omb0,omr0]
#N=np.linspace(5,-5,1000000)
N=np.linspace(18.23435,-5,1000000)
s=odeint(background,s0,N)    


x=s[:,0]
y=s[:,1]
omb=s[:,2]
omr=s[:,3]

qs=2
r=-0.01
ps=15

omphi=((ps*x**2-r*y**2)*(qs*x**2+y**2)-2*r*y**2)/(ps*x**2-3*r*y**2)
omc=1-omphi-omr-omb
wphi=(qs*x**2-y**2)/omphi
eps=-3/2*(omc+omb+4/3*omr+omphi*(1+wphi))
weffe=-1-2/3*eps



plt.xlim([-2,18])
plt.ylim([-1.2,1.2])


plt.plot(N,x)
plt.plot(N,y)
plt.plot(N,wphi)
plt.plot(N,weffe)
plt.xlabel("ln(1+z)")


#plt.xlim([-2,18])
#plt.ylim([1e-8,6])
#plt.xlabel("ln(1+z)")

#plt.semilogy(N,omphi)
#plt.semilogy(N,omc)
#plt.semilogy(N,omb)
#plt.semilogy(N,omr)


print(omphi)
print(omc)
print(omb)
print(omr)
print(wphi)