# 1D Burgers eqn for flow with a step function as IC
#Lax-Wendroff scheme

import numpy
from matplotlib import pyplot

#domain def
nx=101
tmax=2
xmax=4
sigma=0.5
sigma_an=2


dx=xmax/(nx-1)
dt=sigma*dx
nt=int((tmax/dt)+1)


dt_an = dx/0.5
nt_an = int((tmax/dt_an)+1)


x=numpy.linspace(0,xmax,nx)

#IC
u=numpy.zeros(nx)
u[0:int(2/dx)+1]=1

u_an=numpy.zeros(nx)
u_an[0:int(2/dx)+1]=1

#plotting IC
pyplot.figure(figsize=(11,7),dpi=300)
pyplot.plot(x,u,label='IC')


#Numerical solution
for n in range(nt):
    un=u.copy()
    Fn=(un**2)/2
    An=un
    
    u[1:-1]=un[1:-1]-(sigma/2)*(Fn[2:]-Fn[0:-2])+((sigma**2)/4)*(((An[2:]+An[1:-1])*(Fn[2:]-Fn[1:-1]))-(An[0:-2]+An[1:-1])*(Fn[1:-1]-Fn[0:-2]))
    
    
    #BC
    u[0]=1
    u[-1]=0
    
    
#Analytical Sol
for n in range(nt_an):
    un_an=u_an.copy()
    u_an[1:]=un_an[0:-1]
    
    #BC
    u_an[0]=1
    u_an[-1]=0
   
#plotting results

pyplot.plot(x,u,label='LW solution')
pyplot.plot(x,u_an,'--',label='Analytical solution')
pyplot.legend()
pyplot.xlabel('x [m]')
pyplot.ylabel('u [m/s]')
pyplot.title('Lax-Wendroff scheme')

