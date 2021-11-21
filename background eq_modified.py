# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 04:13:02 2021

@author: Beni
"""

#background equations 

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt



x0=1e-13
y0=1/10*x0
omb0=5.8e-6
omr0=0.999962

def background(s,N):
    x=s[0]
    y=s[1]
    omb=s[2]
    omr=s[3]
    
    
    lam = 1
    Q = 0.07
    qs=2
    s=-0.1
   

   
    dxdN=-x *(1/(2 *(x + 2* s *y))* (-3* s *y ** 3 + 
      x *(3 *omb + 4 *omr + 6 *qs *x ** 2 + 
         3 *(1 - omb - omr - (
            qs* x ** 2* (x + s* y) + y *(s + x *y + s *y ** 2))/(
            x + 2 *s *y))) + 
      s *y* (3 + 6 *omb + 8 *omr + 9 *qs *x ** 2 + 
         6* (1 - omb - omr - (
            qs* x ** 2 *(x + s* y) + y *(s + x *y + s *y ** 2))/(
            x + 2 *s *y)))) + (x* (s *
         y *(x + 2 *s *y) ** 
          2* (1 - omb - omr - (
           qs *x ** 2*(x + s* y) + y *(s + x *y + s *y ** 2))/(
           x + 2 *s *y)) *2.45 - 
        1/x*(x + s* y)* (6 *qs* x ** 2*(x + s* y) *(2* x + 3 *s* y) - 
           2* s* y *(x* (-3 + 3 *omb + 4 *omr + 3 *y ** 2 + 
                 3* (1 - omb - omr - (
                    qs* x ** 2*(x + s* y) + 
                    y *(s + x *y + s *y ** 2))/(x + 2 *s *y))) + 
              s *y *(-3 + 6 *omb + 8* omr + 3 *y ** 2 + 
                 6 *(1 - omb - omr - (
                    qs *x ** 2* (x + s *y) + 
                    y *(s + x *y + s *y ** 2))/(x + 2 *s *y)))) + 
           x *(y* (-2* lam* x ** 2* y - 
                 4 *lam *(s ** 2)* (y ** 3) + 
                 s *x* (-1 + qs *x ** 2 + y ** 2 - 
                    6 *lam *y ** 2)) + 
              2 *Q *(x + 2 *s *y) ** 
                2*(1 - omb - omr - (
                 qs *x ** 2* (x + s *y) + y *(s + x *y + s* y ** 2))/(
                 x + 2 *s* y))) *2.45)))/(2 *
       x *(x + s *y) *(s *y *(-1 + y ** 2) + 
         qs *x *(2 *x ** 2 + 7 *s* x *y + 4 *(s ** 2) *(y ** 2))) - 
      2 *s* y *(x + 2 *s *y) ** 
        2* (1 - omb - omr - (
         qs *x ** 2*(x + s* y) + y *(s + x *y + s *y ** 2))/(
         x + 2 *s *y))))
               
               
    dydN=-(1/2)*y*(3*qs*x**2-x*lam*2.45-3*y**2+omr+3)
    
    dombdN=-omb*(3*qs*x**2-3*y**2+omr)
    
    domrdN=-omr*(3*qs*x**2-1-3*y**2+omr)

    return [dxdN,dydN,dombdN,domrdN]


s0=[x0, y0, omb0,omr0]
N=np.linspace(18.23435,-5,1000000)
s=odeint(background,s0,N)    


x=s[:,0]
y=s[:,1]
omb=s[:,2]
omr=s[:,3]

qs=2
s=-0.1
omphi=((1+s*y*x**-1)*(qs*x**2+y**2)+s*y*x**-1)/(1+2*s*y*x**-1)
omc=1-omphi-omr-omb
weffe=-1+2*qs*x**2+omb+omc+4/3*omr
wphi=(qs*x**2-y**2)/omphi



plt.xlim([-4,18.23435])
plt.ylim([1e-8,6])


plt.semilogy(N,omb)
plt.semilogy(N,omr)
plt.semilogy(N,omc)
plt.semilogy(N,omphi)

print(omphi)
