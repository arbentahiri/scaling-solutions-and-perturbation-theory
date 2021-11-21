# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 04:04:18 2021

@author: Beni
"""

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
    X1=s[4]
    X2=s[5]
    Y1=s[6]
    Y2=s[7]
  

    lam = 1
    Q = 0.04
    qs = 1
    m=2
    beta=0
    k=21900
    H0=73
    cs2=1-m*beta/qs
    
    omphi=qs*x**2+y**2
    
    eps=-3*qs*x**2-3/2*(1-qs*x**2-y**2-omb-omr)-3/2*omb-2*omr
    
    ephi=-3+2.45/(2*qs*x)*(lam*y**2-Q*(1-omb-omr-omphi))
    
    v=(2*m*beta*(1+2*beta)*(5+eps+2*ephi)*x**2+(2+eps+2.45*Q*x)*(1+(2-m)*beta)*(1-omb-omr-omphi))/(2*m*beta*(1+2*beta)*x**2+(1+(2-m)*beta)*(1-omb-omr-omphi)) 
    
    r1=2*Q*(3*Q*(1-omb-omr-omphi)+2.45*m*beta*x*(2+ephi+2.45*Q*x))/(3*(1-omb-omr-omphi)*(1+(2-m)*beta))
    r2=(2*m*beta*(1+2*beta)*x**2)/((1-omb-omr-omphi)*(1+(2-m)*beta))
    Gcc= (1+r1)/(1+r2)
    Gcb=1/(1+r2)
    


    dxdN=((1/2)*x*(3*qs*x**2-3-3*y**2+omr)
           +2.45/(2*qs)*(lam*y**2-Q*(1-qs*x**2-y**2-omr-omb)))
   
    dydN=(1/2)*y*(3*qs*x**2-x*lam*2.45-3*y**2+omr+3)
   
    dombdN=omb*(3*qs*x**2-3*y**2+omr)
   
    domrdN=omr*(3*qs*x**2-1-3*y**2+omr)
       
    dX1dN=X2
    
    dX2dN=-v*X2+3/2*((1-omb-omr-omphi)*X1*Gcc+omb*Y1*Gcb)
    
    dY1dN=Y1
    
    dY2dN=-(2+eps)*Y2+3/2*((1-omb-omr-omphi)*X1+omb*Y1)

    return [dxdN,dydN,dombdN,domrdN,dX1dN,dX2dN,dY1dN,dY2dN]


s0=[-2.23828677e-02,7.27951208e-04,1.46981786e-01,0.04527429,np.exp(-5),np.exp(-5),np.exp(-5),np.exp(-5)]
N=np.linspace(-5,5,1000000)
s=odeint(background,s0,N)    

x=s[:,0]
y=s[:,1]
omb=s[:,2]
omr=s[:,3]
X1=s[:,4]
X2=s[:,5]
Y1=s[:,6]
Y2=s[:,7]

lam = 1
Q = 0.04
qs = 1
m=2
beta=0

cs2=1-m*beta/qs
k=21900
H0=73



omphi=qs*x**2+y**2
omc=1-omphi-omr-omb
wphi=(qs*x**2-y**2)/(qs*x**2+y**2)
weffe=-1+2*qs*x**2+omb+omc+4/3*omr
ephi=-3+2.45/(2*qs*x)*(lam*y**2-Q*omc)

eps=-3*qs*x**2-3/2*(1-qs*x**2-y**2-omb-omr)-3/2*omb-2*omr
    
ephi=-3+2.45/(2*qs*x)*(lam*y**2-Q*(1-omb-omr-omphi))
    
v=(2*m*beta*(1+2*beta)*(5+eps+2*ephi)*x**2+(2+eps+2.45*Q*x)*(1+(2-m)*beta)*(1-omb-omr-omphi))/(2*m*beta*(1+2*beta)*x**2+(1+(2-m)*beta)*(1-omb-omr-omphi)) 

r1=2*Q*(3*Q*(1-omb-omr-omphi)+2.45*m*beta*x*(2+ephi+2.45*Q*x))/(3*(1-omb-omr-omphi)*(1+(2-m)*beta))
r2=(2*m*beta*(1+2*beta)*x**2)/((1-omb-omr-omphi)*(1+(2-m)*beta))
Gcc= (1+r1)/(1+r2)
Gcc1= (1+r1)/(1+r2)
Gcb=1/(1+r2)

xprime=1/2*x*(3*qs*x**2-3-3*y**2+omr)+2.45/(2*qs)*(lam*y**2-Q*(1-qs*x**2-y**2-omr-omb))
yprime=1/2*y*(3*qs*x**2-x*lam*2.45-3*y**2+omr+3)
ombprime=omb*(3*qs*x**2-3*y**2+omr)
omrprime=omr*(3*qs*x**2-1-3*y**2+omr)
omphiprime=2*(qs*x*xprime+y*yprime)
omcprime=-omphiprime-ombprime-omrprime
X2prime=-v*X2+3/2*((1-omb-omr-omphi)*X1*Gcc+omb*Y1*Gcb)
Y2prime=-(2+eps)*Y2+3/2*((1-omb-omr-omphi)*X1+omb*Y1)

K=np.exp(-N)*k/H0*(np.exp(-4*N)*omr+np.exp(-3*N)*omb+np.exp(-3*N)*(1-omphi-omb-omr)+omphi)**-1/2
dKdN=(-np.exp(-N)*k/H0*(np.exp(-4*N)*omr+np.exp(-3*N)*(1-omphi-omr)+omphi)**(-1/2))-np.exp(-N)*k/(2*H0)*((np.exp(-4*N)*omr+np.exp(-3*N)*(1-omphi-omr)+omphi)**-3/2)*(-4*np.exp(-4*N)*omr+np.exp(-4*N)*omrprime-3*np.exp(-3*N)*(1-omphi-omr)+np.exp(-3*N)*(-omrprime-2*qs*x*xprime-2*y*yprime)+2*qs*x*xprime+2*y*yprime)
r1=2*Q*(3*Q*omc+2.45*m*beta*x*(2+ephi+2.45*Q*x))/(3*omc*(1+(2-m)*beta))
r2=(2*m*beta*(1+2*beta)*x**2)/(omc*(1+(2-m)*beta))
phiN=1/(qs*cs2*K**2)*(m*beta*X2-2.45*Q*omc*X1/(2*x))
dphindN=-2/(qs*cs2*K**3)*dKdN*(m*beta*X2-2.45*Q*omc*X1/(2*x))+1/(qs*cs2*K**2)*(m*beta*X2prime-2.45*Q*omcprime*X1/(2*x)+2.45*Q*(1-omb-omphi-omr)*xprime*X1/(2*x**2)-2.45*Q*(1-omb-omphi-omr)*X2/(2*x))

dm=(X1+2.45*Q*x*phiN)*omc/(omb+omc)+Y1*omb/(omb+omc)

dmprime=(X2+2.45*Q*(xprime*phiN+x*dphindN))*omc/(omb+omc)+(X1+2.45*Q*x*phiN)*(omcprime/(omc+omb)-omc*(ombprime+omcprime)/(omb+omc)**2)+Y2*omb/(omc+omb)-Y1*(ombprime+omcprime)*omb/(omb+omc)**2+Y1*ombprime/(omc+omb)


sigma80=0.811
dm0=0.81720115
fsigma8=sigma80*dmprime/dm0
#print(dm)

##############################################################
##############################################################
##############################################################

def background(s,N):
    x=s[0]
    y=s[1]
    omb=s[2]
    omr=s[3]
    X1=s[4]
    X2=s[5]
    Y1=s[6]
    Y2=s[7]
  

    lam = 1
    Q = 0.04
    qs = 3
    m=2
    beta=1
    k=21900
    H0=73
    cs2=1-m*beta/qs
    
    omphi=qs*x**2+y**2
    
    eps=-3*qs*x**2-3/2*(1-qs*x**2-y**2-omb-omr)-3/2*omb-2*omr
    
    ephi=-3+2.45/(2*qs*x)*(lam*y**2-Q*(1-omb-omr-omphi))
    
    v=(2*m*beta*(1+2*beta)*(5+eps+2*ephi)*x**2+(2+eps+2.45*Q*x)*(1+(2-m)*beta)*(1-omb-omr-omphi))/(2*m*beta*(1+2*beta)*x**2+(1+(2-m)*beta)*(1-omb-omr-omphi)) 
    
    r1=2*Q*(3*Q*(1-omb-omr-omphi)+2.45*m*beta*x*(2+ephi+2.45*Q*x))/(3*(1-omb-omr-omphi)*(1+(2-m)*beta))
    r2=(2*m*beta*(1+2*beta)*x**2)/((1-omb-omr-omphi)*(1+(2-m)*beta))
    Gcc= (1+r1)/(1+r2)
    Gcb=1/(1+r2)
    


    dxdN=((1/2)*x*(3*qs*x**2-3-3*y**2+omr)
           +2.45/(2*qs)*(lam*y**2-Q*(1-qs*x**2-y**2-omr-omb)))
   
    dydN=(1/2)*y*(3*qs*x**2-x*lam*2.45-3*y**2+omr+3)
   
    dombdN=omb*(3*qs*x**2-3*y**2+omr)
   
    domrdN=omr*(3*qs*x**2-1-3*y**2+omr)
       
    dX1dN=X2
    
    dX2dN=-v*X2+3/2*((1-omb-omr-omphi)*X1*Gcc+omb*Y1*Gcb)
    
    dY1dN=Y1
    
    dY2dN=-(2+eps)*Y2+3/2*((1-omb-omr-omphi)*X1+omb*Y1)

    return [dxdN,dydN,dombdN,domrdN,dX1dN,dX2dN,dY1dN,dY2dN]


s0=[-2.23828677e-02,7.27951208e-04,1.46981786e-01,0.04527429,np.exp(-5),np.exp(-5),np.exp(-5),np.exp(-5)]
N=np.linspace(-5,5,1000000)
s=odeint(background,s0,N)    

x=s[:,0]
y=s[:,1]
omb=s[:,2]
omr=s[:,3]
X1=s[:,4]
X2=s[:,5]
Y1=s[:,6]
Y2=s[:,7]

lam = 1
Q = 0.04
qs = 3
m=2
beta=1

cs2=1-m*beta/qs
k=21900
H0=73



omphi=qs*x**2+y**2
omc=1-omphi-omr-omb
wphi=(qs*x**2-y**2)/(qs*x**2+y**2)
weffe=-1+2*qs*x**2+omb+omc+4/3*omr
ephi=-3+2.45/(2*qs*x)*(lam*y**2-Q*omc)

eps=-3*qs*x**2-3/2*(1-qs*x**2-y**2-omb-omr)-3/2*omb-2*omr
    
ephi=-3+2.45/(2*qs*x)*(lam*y**2-Q*(1-omb-omr-omphi))
    
v=(2*m*beta*(1+2*beta)*(5+eps+2*ephi)*x**2+(2+eps+2.45*Q*x)*(1+(2-m)*beta)*(1-omb-omr-omphi))/(2*m*beta*(1+2*beta)*x**2+(1+(2-m)*beta)*(1-omb-omr-omphi)) 

r1=2*Q*(3*Q*(1-omb-omr-omphi)+2.45*m*beta*x*(2+ephi+2.45*Q*x))/(3*(1-omb-omr-omphi)*(1+(2-m)*beta))
r2=(2*m*beta*(1+2*beta)*x**2)/((1-omb-omr-omphi)*(1+(2-m)*beta))
Gcc= (1+r1)/(1+r2)
Gcc2= (1+r1)/(1+r2)
Gcb=1/(1+r2)

xprime=1/2*x*(3*qs*x**2-3-3*y**2+omr)+2.45/(2*qs)*(lam*y**2-Q*(1-qs*x**2-y**2-omr-omb))
yprime=1/2*y*(3*qs*x**2-x*lam*2.45-3*y**2+omr+3)
ombprime=omb*(3*qs*x**2-3*y**2+omr)
omrprime=omr*(3*qs*x**2-1-3*y**2+omr)
omphiprime=2*(qs*x*xprime+y*yprime)
omcprime=-omphiprime-ombprime-omrprime
X2prime=-v*X2+3/2*((1-omb-omr-omphi)*X1*Gcc+omb*Y1*Gcb)
Y2prime=-(2+eps)*Y2+3/2*((1-omb-omr-omphi)*X1+omb*Y1)

K=np.exp(-N)*k/H0*(np.exp(-4*N)*omr+np.exp(-3*N)*omb+np.exp(-3*N)*(1-omphi-omb-omr)+omphi)**-1/2
dKdN=(-np.exp(-N)*k/H0*(np.exp(-4*N)*omr+np.exp(-3*N)*(1-omphi-omr)+omphi)**(-1/2))-np.exp(-N)*k/(2*H0)*((np.exp(-4*N)*omr+np.exp(-3*N)*(1-omphi-omr)+omphi)**-3/2)*(-4*np.exp(-4*N)*omr+np.exp(-4*N)*omrprime-3*np.exp(-3*N)*(1-omphi-omr)+np.exp(-3*N)*(-omrprime-2*qs*x*xprime-2*y*yprime)+2*qs*x*xprime+2*y*yprime)
r1=2*Q*(3*Q*omc+2.45*m*beta*x*(2+ephi+2.45*Q*x))/(3*omc*(1+(2-m)*beta))
r2=(2*m*beta*(1+2*beta)*x**2)/(omc*(1+(2-m)*beta))
phiN=1/(qs*cs2*K**2)*(m*beta*X2-2.45*Q*omc*X1/(2*x))
dphindN=-2/(qs*cs2*K**3)*dKdN*(m*beta*X2-2.45*Q*omc*X1/(2*x))+1/(qs*cs2*K**2)*(m*beta*X2prime-2.45*Q*omcprime*X1/(2*x)+2.45*Q*(1-omb-omphi-omr)*xprime*X1/(2*x**2)-2.45*Q*(1-omb-omphi-omr)*X2/(2*x))

dm1=(X1+2.45*Q*x*phiN)*omc/(omb+omc)+Y1*omb/(omb+omc)

dmprime1=(X2+2.45*Q*(xprime*phiN+x*dphindN))*omc/(omb+omc)+(X1+2.45*Q*x*phiN)*(omcprime/(omc+omb)-omc*(ombprime+omcprime)/(omb+omc)**2)+Y2*omb/(omc+omb)-Y1*(ombprime+omcprime)*omb/(omb+omc)**2+Y1*ombprime/(omc+omb)


sigma80=0.811
dm01=0.81601742
fsigma81=sigma80*dmprime1/dm01
#print(dm1)

############################################################
############################################################
############################################################

def background(s,N):
    x=s[0]
    y=s[1]
    omb=s[2]
    omr=s[3]
    X1=s[4]
    X2=s[5]
    Y1=s[6]
    Y2=s[7]
  

    lam = 1
    Q = 0.02
    qs = 2
    m=3
    beta=0.5
    k=21900
    H0=73
    cs2=1-m*beta/qs
    
    omphi=qs*x**2+y**2
    
    eps=-3*qs*x**2-3/2*(1-qs*x**2-y**2-omb-omr)-3/2*omb-2*omr
    
    ephi=-3+2.45/(2*qs*x)*(lam*y**2-Q*(1-omb-omr-omphi))
    
    v=(2*m*beta*(1+2*beta)*(5+eps+2*ephi)*x**2+(2+eps+2.45*Q*x)*(1+(2-m)*beta)*(1-omb-omr-omphi))/(2*m*beta*(1+2*beta)*x**2+(1+(2-m)*beta)*(1-omb-omr-omphi)) 
    
    r1=2*Q*(3*Q*(1-omb-omr-omphi)+2.45*m*beta*x*(2+ephi+2.45*Q*x))/(3*(1-omb-omr-omphi)*(1+(2-m)*beta))
    r2=(2*m*beta*(1+2*beta)*x**2)/((1-omb-omr-omphi)*(1+(2-m)*beta))
    Gcc= (1+r1)/(1+r2)
    Gcb=1/(1+r2)
    


    dxdN=((1/2)*x*(3*qs*x**2-3-3*y**2+omr)
           +2.45/(2*qs)*(lam*y**2-Q*(1-qs*x**2-y**2-omr-omb)))
   
    dydN=(1/2)*y*(3*qs*x**2-x*lam*2.45-3*y**2+omr+3)
   
    dombdN=omb*(3*qs*x**2-3*y**2+omr)
   
    domrdN=omr*(3*qs*x**2-1-3*y**2+omr)
       
    dX1dN=X2
    
    dX2dN=-v*X2+3/2*((1-omb-omr-omphi)*X1*Gcc+omb*Y1*Gcb)
    
    dY1dN=Y1
    
    dY2dN=-(2+eps)*Y2+3/2*((1-omb-omr-omphi)*X1+omb*Y1)

    return [dxdN,dydN,dombdN,domrdN,dX1dN,dX2dN,dY1dN,dY2dN]


s0=[-2.23828677e-02,7.27951208e-04,1.46981786e-01,0.04527429,np.exp(-5),np.exp(-5),np.exp(-5),np.exp(-5)]
N=np.linspace(-5,5,1000000)
s=odeint(background,s0,N)    

x=s[:,0]
y=s[:,1]
omb=s[:,2]
omr=s[:,3]
X1=s[:,4]
X2=s[:,5]
Y1=s[:,6]
Y2=s[:,7]

lam = 1
Q = 0.02
qs = 2
m=3
beta=0.5

cs2=1-m*beta/qs
k=21900
H0=73



omphi=qs*x**2+y**2
omc=1-omphi-omr-omb
wphi=(qs*x**2-y**2)/(qs*x**2+y**2)
weffe=-1+2*qs*x**2+omb+omc+4/3*omr
ephi=-3+2.45/(2*qs*x)*(lam*y**2-Q*omc)

eps=-3*qs*x**2-3/2*(1-qs*x**2-y**2-omb-omr)-3/2*omb-2*omr
    
ephi=-3+2.45/(2*qs*x)*(lam*y**2-Q*(1-omb-omr-omphi))
    
v=(2*m*beta*(1+2*beta)*(5+eps+2*ephi)*x**2+(2+eps+2.45*Q*x)*(1+(2-m)*beta)*(1-omb-omr-omphi))/(2*m*beta*(1+2*beta)*x**2+(1+(2-m)*beta)*(1-omb-omr-omphi)) 

r1=2*Q*(3*Q*(1-omb-omr-omphi)+2.45*m*beta*x*(2+ephi+2.45*Q*x))/(3*(1-omb-omr-omphi)*(1+(2-m)*beta))
r2=(2*m*beta*(1+2*beta)*x**2)/((1-omb-omr-omphi)*(1+(2-m)*beta))
Gcc= (1+r1)/(1+r2)
Gcc3= (1+r1)/(1+r2)
Gcb=1/(1+r2)

xprime=1/2*x*(3*qs*x**2-3-3*y**2+omr)+2.45/(2*qs)*(lam*y**2-Q*(1-qs*x**2-y**2-omr-omb))
yprime=1/2*y*(3*qs*x**2-x*lam*2.45-3*y**2+omr+3)
ombprime=omb*(3*qs*x**2-3*y**2+omr)
omrprime=omr*(3*qs*x**2-1-3*y**2+omr)
omphiprime=2*(qs*x*xprime+y*yprime)
omcprime=-omphiprime-ombprime-omrprime
X2prime=-v*X2+3/2*((1-omb-omr-omphi)*X1*Gcc+omb*Y1*Gcb)
Y2prime=-(2+eps)*Y2+3/2*((1-omb-omr-omphi)*X1+omb*Y1)

K=np.exp(-N)*k/H0*(np.exp(-4*N)*omr+np.exp(-3*N)*omb+np.exp(-3*N)*(1-omphi-omb-omr)+omphi)**-1/2
dKdN=(-np.exp(-N)*k/H0*(np.exp(-4*N)*omr+np.exp(-3*N)*(1-omphi-omr)+omphi)**(-1/2))-np.exp(-N)*k/(2*H0)*((np.exp(-4*N)*omr+np.exp(-3*N)*(1-omphi-omr)+omphi)**-3/2)*(-4*np.exp(-4*N)*omr+np.exp(-4*N)*omrprime-3*np.exp(-3*N)*(1-omphi-omr)+np.exp(-3*N)*(-omrprime-2*qs*x*xprime-2*y*yprime)+2*qs*x*xprime+2*y*yprime)
r1=2*Q*(3*Q*omc+2.45*m*beta*x*(2+ephi+2.45*Q*x))/(3*omc*(1+(2-m)*beta))
r2=(2*m*beta*(1+2*beta)*x**2)/(omc*(1+(2-m)*beta))
phiN=1/(qs*cs2*K**2)*(m*beta*X2-2.45*Q*omc*X1/(2*x))
dphindN=-2/(qs*cs2*K**3)*dKdN*(m*beta*X2-2.45*Q*omc*X1/(2*x))+1/(qs*cs2*K**2)*(m*beta*X2prime-2.45*Q*omcprime*X1/(2*x)+2.45*Q*(1-omb-omphi-omr)*xprime*X1/(2*x**2)-2.45*Q*(1-omb-omphi-omr)*X2/(2*x))

dm2=(X1+2.45*Q*x*phiN)*omc/(omb+omc)+Y1*omb/(omb+omc)

dmprime2=(X2+2.45*Q*(xprime*phiN+x*dphindN))*omc/(omb+omc)+(X1+2.45*Q*x*phiN)*(omcprime/(omc+omb)-omc*(ombprime+omcprime)/(omb+omc)**2)+Y2*omb/(omc+omb)-Y1*(ombprime+omcprime)*omb/(omb+omc)**2+Y1*ombprime/(omc+omb)


sigma80=0.811
dm02=0.80310072
fsigma82=sigma80*dmprime2/dm02

#print(dm2)

########################################################################
########################################################################
########################################################################
plt.xlim([-1,5])
plt.ylim([0,1.4])
plt.xlabel("z")
plt.ylabel("$G_{cc}/G$")
plt.plot(np.exp(-N)-1,Gcc1)
plt.plot(np.exp(-N)-1,Gcc2)
plt.plot(np.exp(-N)-1,Gcc3)

#plt.xlim([0,2])
#plt.ylim([0.2,0.5])
#plt.xlabel("z")
#plt.ylabel("$f\sigma_8$")
#plt.plot(np.exp(-N)-1,fsigma8)
#plt.plot(np.exp(-N)-1,fsigma81)
#plt.plot(np.exp(-N)-1,fsigma82)

#plt.xlim([2,-18])
#plt.ylim([1e-8,6])
#plt.semilogy(N,omphi)
#plt.semilogy(N,omc)
#plt.semilogy(N,omb)
#plt.semilogy(N,omr)

