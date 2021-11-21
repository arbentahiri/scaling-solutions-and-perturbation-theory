#background equations 

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math




def background(s,N):
    x=s[0]
    y=s[1]
    omb=s[2]
    omr=s[3]
  

    lam = 1
    Q = 0.07
    qs = 2


    dxdN=-((1/2)*x*(3*qs*x**2-3-3*y**2+omr)
           +2.45/(2*qs)*(lam*y**2-Q*(1-qs*x**2-y**2-omr-omb)))
   
    dydN=-(1/2)*y*(3*qs*x**2-x*lam*2.45-3*y**2+omr+3)
   
    dombdN=-omb*(3*qs*x**2-3*y**2+omr)
    domrdN=-omr*(3*qs*x**2-1-3*y**2+omr)

    return [dxdN,dydN,dombdN,domrdN]


s0=[1e-13,1e-14,5.8e-6,0.999962]
#s0=[-2.23828677e-02,7.27951208e-04,1.46981786e-01,0.04527429]
N=np.linspace(18.23435,-5,1000000)
#N=np.linspace(5,-5,1000000)
s=odeint(background,s0,N)    

x=s[:,0]
y=s[:,1]
omb=s[:,2]
omr=s[:,3]

lam = 1
Q = 0.07
qs = 2



omc=1-qs*x**2-y**2-omr-omb
omphi=qs*x**2+y**2
wphi=(qs*x**2-y**2)/(qs*x**2+y**2)
weffe=-1+2*qs*x**2+omb+omc+4/3*omr
ephi=-3+2.45/(2*qs*x)*(lam*y**2-Q*omc)


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

print(omb)
print(omr)
print(x)
print(y)
