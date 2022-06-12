import numpy as np
import math
from Emanuel.RK4 import *
fm1 = 10 #Flujo masico
fm2 = 10 # Flujo masico
epsilon = 10
gB = 10 #Grados Brix
NBf=10
Dintcz = 10
SepBf = 10
kTJ = 10
KTtb = 10 #Conductividad termica
Dexttb = 10
TeA = 10
TeJ = 10


wtb = 10 #Espesor del tubos

Cpj = 4.187*(1 - 0.006*gB) #Capacidad calorifica especifica del jugo.

#T = 273.15 #Grados Kelvin.

def Cpa(T):
    return 7.2951406*(10**(-8))*(T**3) - 7.431358269*(10**(-5))*(T**2) + 0.02620653*T + 1.0156045
#Capacidad calorifica del agua.
def pj(T):
    return 1179.7 - 0.354*T  #Densidad del jugo.
def pA(T):
    return -0.004323394923*((T-273.15)**2) - 0.04038343824*(T-273.15) + 1000.807908 #Densidad del agua.

Ntb = 10 #Número de tubos.
Ltb = 10 #Longitud de los tubos.
Dinttb = 10 #Diametro interno de los tubos.
Pitb = 10


Vj = Ntb*Ltb*((np.pi*(Dinttb**2))/4) #Volumen del jugo.

Lcz = 10 #Longitud de la coraza.
Dintcz = 10 #Diametro interno de la coraza.
Dexttb = 10 #Diametro externo de los tubos.

VA = Ntb*Lcz*((np.pi*(Dintcz**2))/4) - Ntb*Ltb*((np.pi*(Dexttb**2))/4)

def Mj(T):
    return pj(T)*Vj #Masa del jugo.

def MA(T):
    return pA(T)*VA #Masa del agua.

Mtb = 10 #Masa del tubo.
Mcz = 10 #Masa de la coraza.
kTA = 10

def Hj(T):
    return 5.46829932*(10**(-4))*((T-273.15)**2) + 3.694618471*(T-273.15)-0.633984 #Entalpia del jugo.
Deqcz = 4((((Pitb**2)*math.sqrt(3))/4) - np.pi*(Dexttb**2))/((np.pi*Dexttb)/2)

def muA(T):
    return (0.0002471250997*(T-273.15)**2)-0.03226252622(T-273.15)+1.545370351
def muJ(T):
    return (0.047)/(T-273.15)

def NReA(T):
    return (pA(T)*Deqcz*VA)/(muA(T))
AFtb = (np.pi*(Dinttb**2)/(4))

def velJ(T):
    return (fm1/pj(T))/(AFtb)

def NReJ(T):
    return (pj(T)*Vj*Dinttb)/muJ(T)

def NPrA(T):
    return (Cpa(T)*muA(T))/(kTA)

def hA(T):
    return (kTJ*Deqcz)*(0.36*(NReA(T))**(0.55)*(NPrA(T))**(0.333))

def HA(T):
    return 4.180358722*(T-273.15) + 0.3506658477 #Entalpia del agua.
AHtbJ = np.pi*Dinttb*Ltb
AHtbA = np.pi*Dexttb*Ltb
AH = (1/2)*(AHtbJ+ AHtbA)

def deltaT1(T4):
    return TeJ - T4

def deltaT2(T2):
    return T2 - TeA

def deltaTML(T2,T4):
    return (deltaT1(T4)-deltaT2(T2))/(math.log(deltaT1(T4)/deltaT2(T2)))

AHtb = (AHtbJ + AHtbA)/2

def NPrj(T):
    return (Cpj*muJ(T))/(kTJ)

def hJ(T):
    return (kTJ/Dinttb)*(0.0243*(NReJ(T))**(0.8)*(NPrj(T))*0.333)

def U(T):
    return 1/(AH*((1/(hJ(T)*AHtbJ))+(wtb/(KTtb*AHtb))+(1/(hA(T)*AHtbA))))

def dQIdeC(T,T2,T4):
    return U(T)*AH*deltaTML(T2,T4)

##dQtbA_A = hA*AHtbA*(TtbA-TA); 

def dQtbA_A(T,T2,T4):
    return dQIdeC(T,T2,T4)






##dQJ_tbJ = hJ*AHtbJ*(TJ-TtbJ);
def dQJ_tbJ(T,T2,T4):
    return dQIdeC(T,T2,T4)

def dQtbJ_tbA(T,T2,T4):
    return dQJ_tbJ(T,T2,T4) - dQtbA_A(T,T2,T4)


dQA_czA = 0

dQczA_czat = 0 #Presuntamente todo esto es igual a 0.

dQtbJ_tbA = 0

def NNu(T):
    return (hJ(T)*Dinttb)/(kTJ)


Septb = Pitb - Dexttb
AFcz = (Dintcz*Septb*SepBf)/Pitb

def fD(T):
    return (-2*math.log10((epsilon/Dinttb)-(5.02/NReJ(T))*(math.log10((epsilon/(Dinttb*3.71)+(14.5/NReJ(T)))))))**(-2) #que es epsilon?

def hf1_2(T):
    return fD(T)*((Ltb*(velJ(T))**2)/(Dinttb*2))

def hfcz(T):
    return fD(T)*(Dintcz*(NBf+1))/(Deqcz)

def vAcz(T):
    return ((fm1/pA(T))/AFcz)


############################################################
T1 = 10
def fT2(t, T2):
    return (1/(Cpj*Mj(T2)))*(fm1*Hj(T1)-fm2*Hj(T2)-dQJ_tbJ(T2,T2,T4))
