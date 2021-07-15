import numpy as np
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import *

def rex(Z,V):
    #print(Z)
    k1a=(7.056*10**10)*np.exp(-(125000/(8.314*Z[6])))
    k2a=(1.224*10**12)*np.exp(-(145000/(8.314*Z[6])))
    Cto=float(CTo.get())
    Ft=Z[0]+Z[1]+Z[2]+Z[3]+Z[4]+Z[5]
    To=660
    Ca=Cto*(Z[0]/Ft)*(To/Z[6])
    Cb=Cto*(Z[1]/Ft)*(To/Z[6])
    Cc=Cto*(Z[2]/Ft)*(To/Z[6])
    Cd=Cto*(Z[3]/Ft)*(To/Z[6])
    Ce=Cto*(Z[4]/Ft)*(To/Z[6])
    Cf=Cto*(Z[5]/Ft)*(To/Z[6])
    r1a=(-k1a*Ca)
    r2a=(-k2a*Ca)
    #print(r2a)
    
    dAdV=r1a+r2a
    dBdV=3.5*r1a +5.5*r2a
    dCdV=-r1a
    dDdV=-4*r1a-5*r2a
    dEdV=-2*r2a
    dFdV=-2*r2a
    H1=-271900-(59.68*(Z[6]-298))
    H2=126500+(101.202*(Z[6]-298))
    k=Adiabatic.get()
    if (k=="Yes"):
        dGdV=(r1a*H1+r2a*H2)/(20.056*Z[0]+29.526*Z[1]-72.015*Z[2]+33.933*Z[3]+29.556*Z[5]+27.437*Z[4])
    else:
        dGdV=((r1a*H1+r2a*H2)-4000*(Z[6]-373))/(20.056*Z[0]+29.526*Z[1]-72.015*Z[2]+33.933*Z[3]+29.556*Z[5]+27.437*Z[4])
        
    return [dAdV,dBdV,dCdV,dDdV,dEdV,dFdV,dGdV]





def solve(event):
    V=np.linspace(0,100000,num=10000)
    #W=V*900
    Zo=[float(Fao.get()),float(Fbo.get()),0,0,0,0,660]
    Conc=odeint(rex,Zo,V)
    #Conc=solve_ivp(rex,(0,2),Zo,t_eval=V)
    #print(len(Conc))
    Fa=Conc[:,0]
    Fb=Conc[:,1]
    Fc=Conc[:,2]
    Fd=Conc[:,3]
    Fe=Conc[:,4]
    Ff=Conc[:,5]
    T=Conc[:,6]




    X=(float(Fao.get())-Fa)/float(Fao.get())
    selec=(7.056*10**10/1.224*10**12)*np.exp((145000-12500)/(8.314*T))
    plt.rcParams["figure.figsize"] = [11, 11]
    fig ,axs=plt.subplots(3,2)
    fig.suptitle("Imporant Graphs")

    axs[0,0].plot(V,Fd,'-',label = "Fd")
    axs[0,0].plot(V,Fa,'--',label = "Fa")
    axs[0,0].plot(V,Fc,label = "Fc")
    axs[0,0].plot(V,Fe,label = "Fe")
    axs[0,0].plot(V,Ff,label = "Ff")
    axs[0,0].set(xlabel="W")
    axs[0,0].legend()
    axs[0,1].plot(V,T)
    axs[0,1].set(xlabel="W",ylabel="T")
    axs[1,0].plot(V,X)
    axs[1,0].set(xlabel="W",ylabel="X")
    axs[1,1].plot(T,X)
    axs[1,1].set(xlabel="T",ylabel="X")
    axs[2,0].plot(T,selec)
    axs[2,0].set(xlabel="T",ylabel="Inst. Selectivity")
    axs[2,1].plot(V,Fb)
    axs[2,1].set(xlabel="W",ylabel="Fb")
    plt.show()


global screen
screen=Tk()
screen.geometry("600x600")
screen.title("MALEIC ANHYDRIDE")
Label(screen, text="MALEIC ANHYDRIDE PRODUCTION PLANT",bg="red",width="500",height="4").pack()
Label(text="").pack()
Label(text="").pack()
    
Fao=StringVar()
Fbo=StringVar()
CTo=StringVar()
Adiabatic=StringVar()
label1 = tk.Label(screen, text="Fao(n-butane)")
label1.place(x="5", y="100")
label2 = tk.Label(screen, text="Fbo(oxygen)")
label2.place(x="10", y="150")
label3 = tk.Label(screen, text="Cto(Tot. conc.)")
label3.place(x="5", y="200")
label3 = tk.Label(screen, text="Adiabatic")
label3.place(x="5", y="250")
Entry(screen,textvariable=Fao,).place(x="100",y="100")
Entry(screen,textvariable=Fbo).place(x="100",y="150")
Entry(screen,textvariable=CTo).place(x="100",y="200")
Entry(screen,textvariable=Adiabatic).place(x="100",y="250")
Label(screen,text="").pack()
button_1=Button(text="Solve",width="30")
button_1.bind("<Button-1>",solve)
button_1.place(x="100",y="300")

screen.mainloop()
    

    
    


"""   
def rex(Z,V):
    k1a=(7.056*10**10)*math.exp(4000*(1/300-1/Z[3]))
    k2a=(1.224*10**12)*math.exp(9000*(1/300-1/Z[3]))
    Cto=48
    Ft=Z[0]+Z[1]+Z[2]
    To=423
    Ca=Cto*(Z[0]/Ft)*(To/Z[3])
    Cb=Cto*(Z[1]/Ft)*(To/Z[3])
    Cc=Cto*(Z[2]/Ft)*(To/Z[3])
    r1a=-k1a*Ca
    r2a=-k2a*Ca*Ca
    
    dAdV=r1a+r2a
    dBdV=-r1a
    dCdV=-r2a/2
    dDdV=(4000*(373-Z[3])+(-r1a)*20000+(-r2a)*60000)/(90*Z[0]+90*Z[1]+180*Z[2])
    return [dAdV,dBdV,dCdV,dDdV]
"""









    
    
    


