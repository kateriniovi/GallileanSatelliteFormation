# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
Code Created by Aikaterini-Niovi Triantafyllaki with the help of Akbar Selemani for the purposes of the course Planetary 
Systems, Aix-Marseille University M2, 2021/22
"""


import numpy as np
import matplotlib.pyplot as plt



time = 3600*24*365.25*1000 # convertion time to s

tau = 10*time #  accretion time s #SI
delta_z = 10e3 # Thickness of surface layer m #SI
Rf = 2000e3 # final radius if no mass lost m #SI
taulist=[0,tau,5000]
rhoB = 1800 #background density kg m-3 #SI
rhorock = 3500 # rock density kg m-3 #SI
rhoice = 920 # ice density kg m-3 #SI
Mf = rhoB*(4/3)*np.pi*Rf*Rf*Rf # final mass  kg #SI

G = 6.67408e-11 # universal constante of the gravitation: m3 kg-1 s-2 #SI
sigma = 5.670374419e-8 # Stefan-Boltzman constant: kg s-3 K-4 #SI
kB = 1.38064852e-23 # Botzman constant:  m2 kg s-2 K-1 #SI

T0 = 273.16 # temperature that water melts, (0 C), #SI
T_back = 200 #backround temperature K #SI
P0 = 611.73 # pressure at T0 kg m-1 s-2 #SI

Mg = 0.01802 # Water molecular mass kg mol-1 #SI? is mol si?
Lv = 2.5e6 # Latent heat J kg-1 #SI
cp = 1300 # Specific heat J kg-1 K-1 #SI

R = 8.314 # perfect gas cst J mol-1 K-1 #SI
R_spe = R/Mg # perfect gas cst devided by the water molecular gas gives us specific gas cst J kg-1 K-1 #SI



def RK4(dxdt,x,t,h,args=()):


    k1 = dxdt(x,t,*args)
    k2 = dxdt(x+k1*h/2,t+h/2,*args)
    k3 = dxdt(x+k2*h/2,t+h/2,*args)
    k4 = dxdt(x+k3*h,t+h,*args)

    x +=  h*(k1+2*k2+2*k3+k4)/6
    
    return x

#def RK4 (f,y,t, args=()):
#    n = len(t)
#    y = np.zeros(n)
#    y[0] = y0
#    for i in range(n - 1):
#        h = t[i+1] - t[i]
#        k1 = f(y[i], t[i], *args)
#        k2 = f(y[i] + k1 * h / 2., t[i] + h / 2., *args)
#        k3 = f(y[i] + k2 * h / 2., t[i] + h / 2., *args)
#        k4 = f(y[i] + k3 * h, t[i] + h, *args)
#        y[i+1] = y[i] + (h / 6.) * (k1 + 2*k2 + 2*k3 + k4)
#    return y


def dMdt(Ms,t,Ml):

    if Ms + Ml >= Mf:
        dMdt = 0
    else:
        dMdt = 7.15*(Ms/Mf)**(2/3)*((Mf-Ms)/tau)

    return dMdt


def Phi(Ml,t,u,r,rho_s,T,Ms):


    Phi = 0

    if T>=T0:
        Phi = 4*np.pi*r*r*rho_s*u
        
    E_escape = Phi*(Lv+G*Ms/r)
    E_excess = G*Ms*dMdt(Ms,t,Ml)/r - 4*np.pi*r*r*sigma*(T**4-T_back**4) - cp*dMdt(Ms,t,Ml)*(T-T_back)

    if E_escape > abs(E_excess):
        Phi = abs(E_excess) / (Lv + G*Ms/r)

    return Phi

def dTdt(T,t,r,E_excess,E_escape,rho):
   
    
    dTdt = E_excess - E_escape
    dTdt /= cp*rho*(4/3)*np.pi*(r*r*r - (r - min(r,delta_z))**3)    
    return dTdt

def computation(h,t,Ms,Ml,T,r,u,rho):


    ### Satellite mass computation ###
    Ms = RK4(dMdt,Ms,t,h,args=(Ml,)) # Mass satellite: kg

    ### surfaoce velocity calculation as seen by the paper of
    Ps = P0*np.exp((Lv/R_spe/T0)*(1 - T0/T)) # Suface pressure: Pa
    u0 = R_spe*T # sound velocity in square: (m s-1)**2
    rc = G*Ms/(2*u0) # radius of the satellite + atmosphere: m
    rho_s = Ps / u0 # surface density kg m-3
    r = ((3*Ms)/(4*np.pi*rho))**(1/3) # Surface radius: m

    # Surface velocity: m s-1
    
    if T<T0:  
        u=0 # no vertical velocity
    else:
        u = np.sqrt(u0)*(rc*rc/(r*r))*np.exp(-1/2 + 2*(1-rc/r))
        if u>np.sqrt(u0):
            u = np.sqrt(u0)
       
    ### Mass lost function calculation ###
    Ml = RK4(Phi,Ml,t,h,args=(u,r,rho_s,T,Ms)) # Mass lost: kg

    ### density computation ###
    f_i = Ml/(Ms + Ml)
    rho = rhoB +  f_i*rhoice   
    
    ### Temperature calculation: K ###

    # Energy to form satellite
    E_excess = G*Ms*dMdt(Ms,t,Ml)/r - 4*np.pi*r*r*sigma*(T**4-T_back**4) - cp*dMdt(Ms,t,Ml)*(T-T_back)
    # Energy of gas escape
    E_escape = Phi(Ml,t,u,r,rho_s,T,Ms)*(Lv + G*Ms/r)

      
    T = RK4(dTdt,T,t,h,args=(r,E_excess,E_escape,rho)) # Surface temperature: K 


        
    
                    
    return Ms,Ml,T,r,rc,u,rho

if __name__ == '__main__':
    
   

    ## Initial condition ##
    N = 5000 # number of point
    h = tau/N # time step
    t = 0.0 #0.0 # initial time

    Ms = 0.001
    Ml = 0.001
    T = T_back
    r = 0.001
    u = 0.0
    rc = 0.0
    rho = rhoB


    t += h # next step
    
    #### begin while ####
    n = 0

    t += h 
    n += 1
    tl=[]
    MsMf=[]
    MlMf=[]
    Tl=[]
    rl=[]
    rcl=[]
    rrcl=[]
    ul=[]
    rhol=[]
    Msl=[]
    Mfl=[]
    Mll=[]
    #l stands for list

    for i in  range(len(taulist)):
        while t<tau:
        
            Ms,Ml,T,r,rc,u,rho = computation(h,t,Ms,Ml,T,r,u,rho)

        # save on out file

            t += h 
            n += 1
            tl.append(t/time)
            MsMf.append(Ms/Mf)
            MlMf.append(Ml/Mf)
            Tl.append(T)
            rl.append(r/1000)
            rcl.append(rc/1000)
            rrcl.append(rc/r)
            ul.append(u)
            rhol.append(rho)
            Msl.append(Ms)
            Mfl.append(Mf)
            Mll.append(Ml)

    ax1 = plt.subplot(2,2,1)

    ax1.plot(tl,Tl, color='black')
    ax1.set_xlabel(r'time ($kyears$)')
    ax1.set_ylabel(r'Temperature ($K$)')
    ax1.set_ylim(195,300)
    
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.plot(tl,rl,color='red')  
    ax2.set_xlabel(r'time ($kyears$)')
    ax2.set_ylabel(r'radius ($km$)', labelpad=-200, fontsize=8) #this from the right side of the plot
    ax2.set_ylim(0,1900)
    ax2.tick_params(axis='y', labelcolor=color)

    
    ax3 = plt.subplot(2,2,2)
    ax3.plot(tl,rrcl, color='black')
    ax3.set_xlabel(r'time ($times$)')
    ax3.set_ylabel(r'$r_c/r$')
    
    ax4 = ax3.twinx()
    ax4.semilogy(tl,ul, color='red') #log scale? how--> yscale("log") #changed it in append
    ax4.tick_params(axis='y', labelcolor=color)
    
    plt.subplot(2,2,3)
        # plot of Ms and Ml
    plt.plot(tl,MsMf,color='blue',label='$M_S$')
    plt.plot(tl,MlMf,color='orange',label='$M_L$')
    plt.legend()
    plt.xlim(0,10)
    plt.ylim(0.0,0.85)

    plt.xlabel(r'time ($times$)')
    plt.ylabel(r'$M/M_f$')


    # plot of density
    plt.subplot(2,2,4)
    plt.plot(tl,rhol,color='black')

    plt.xlim(0,10)

    plt.xlabel(r'time ($times$)')
    plt.ylabel(r'density ($kg$ $m^{-3}$)')
    
    plt.savefig('Bierson.png',format='png')
    plt.show()
