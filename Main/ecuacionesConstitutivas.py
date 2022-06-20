
import numpy as np
from RK4 import *
fm1 = 1250 #Flujo masico
fm3 = 3500 # Flujo masico
epsilon = 0              
gB = 20 #Grados Brix
NBf = 500 #Número de Bafles
Dintcz = 0.023 #Diametro interno de la coraza
SepBf = 0.241 #Separación entre Bafles
kTJ = .73 #Conductividad termica del jugo
KTtb = 16.3 #Conductividad termica de los tubos #####
Dexttb = 0.019 #Diametro externo de los tubos
TeA = 293.15 #Temperatura de entrada del agua
TeJ = 323.15 #Temperatura de entrada del jugo
Cpj = 4.187*(1 - 0.006*gB) #Capacidad calorifica especifica del jugo
Cpa =  lambda T4: 7.2951406e-8*(T4**3) - 7.431358269e-5*(T4**2) + 0.02620653*T4 + 1.0156045
#Capacidad calorifica del agua.
pj =  lambda T2: 1179.7 - 0.354*T2  #Densidad del jugo.
pA = lambda T4: -0.004323394923*((T4-273.15)**2) - 0.04038343824*(T4-273.15) + 1000.807908
 #Densidad del agua.
Ntb = 532 #Número de tubos.
Ltb = 4 #Longitud de los tubos.
Dinttb = 0.016 #Diametro interno de los tubos.
Pitb = 0.02381 #Separacion entre centros de tubos contiguos
wtb = Dexttb - Dinttb #Espesor del tubos
Vj = Ntb*Ltb*((np.pi*(Dinttb**2))/4) #Volumen del jugo.
Lcz = 3 #Longitud de la coraza.


VA = Ntb*Lcz*((np.pi*(Dintcz**2))/4) - Ntb*Ltb*((np.pi*(Dexttb**2))/4) #volumen del agua

Mj = lambda T2: pj(T2)*Vj*(1/0.1) #Masa del jugo.

MA = lambda T4 : pA(T4)*VA*(1/0.01)  #Masa del agua.

Mtb = 1.120*Ltb #Masa del tubo.
Mcz = 2 #Masa de la coraza.
kTA = 0.58 #Confuctividad termica del agua.

Hj = lambda T2:5.46829932e-4*((T2-273.15)**2) + 3.694618471*(T2-273.15)-0.633984 #Entalpia del jugo.

Deqcz = 4*((((Pitb**2)*np.sqrt(3))/4) - (np.pi*(Dexttb**2)/8))/((np.pi*Dexttb)/2)
#Diametro equivalente de la coraza

muA  = lambda T4: (0.0002471250997*(T4-273.15)**2) - 0.03226252622*(T4-273.15) + 1.545370351 
#viscosidad del agua

muJ  = lambda T2: (0.047)/(T2-273.15) #Viscosidad del jugo

velA = lambda T4: (fm3/pA(T4))/(AFcz) #Velocidad del agua

NReA = lambda T4: (pA(T4)*Deqcz*velA(T4))/(muA(T4)) #Número de Reynolds a las condiciones del agua

AFtb = (np.pi*(Dinttb**2)/(4)) #Area de flujo del tubo

velJ = lambda T2: (fm1/pj(T2))/(AFtb) #Velocidad del jugo

NReJ = lambda T2: (pj(T2)*velJ(T2)*Dinttb)/muJ(T2) #Número de Reynolds a las condiciones del jugo

NPrA = lambda T4: (Cpa(T4)*muA(T4))/(kTA) #Número de Prandtl a las condiciones del agua

hA = lambda T4:(kTA/Deqcz)*(0.36*(NReA(T4)**(0.55))*(NPrA(T4))**(1./3)) 
#Correlación de McAdams a las condiciones del agua

HA =  lambda T4: 4.180358722*(T4-273.15) + 0.3506658477  if ((T4-273.15) <= 50) else (8.64*8.314/18)*647.15*((T4/647.15-0.4221)/(1-647.15))**0.053 #Entalpia del agua.
#Entalpía del agua

AHtbJ = np.pi*Dinttb*Ltb #Area interna del tubo en contacto con el jugo

AHtbA = np.pi*Dexttb*Ltb #Area externa del tubo en contacto con el agua

AH = (1./2)*(AHtbJ+ AHtbA) #Area de intercambio del tubo

deltaT1 = lambda T4: TeJ - T4 #diferencia temperatura entrada

deltaT2 = lambda T2: T2 - TeA #diferencia temperatura salida

deltaTML =  lambda T2,T4:(deltaT1(T4)-deltaT2(T2))/(np.log(deltaT1(T4)/deltaT2(T2)))
#Flujo en contracorriente

AHtb = (AHtbJ + AHtbA)/2 #Area de intercambio del tubo

NPrj = lambda T2: (Cpj*muJ(T2))/(kTJ) #Número de Prandtl a las condiciones del jugo

hJ = lambda T2:(kTJ/Dinttb)*(0.0243*(NReJ(T2)**(4./5))*(NPrj(T2)**(1./3)))
#Correlación de McAdams a las condiciones del jugo.

U =lambda T2,T4: 1/(AH*((1/(hJ(T2)*AHtbJ))+(wtb/(KTtb*AHtb))+(1/(hA(T4)*AHtbA))))
#Coeficiente global de transferencia de calor

dQIdeC = lambda T2,T4:U(T2,T4)*AH*deltaTML(T2,T4)
#Flujo de calor del intercambiador de calor

dQtbA_A = lambda T2,T4: dQIdeC(T2,T4)
#Flujo de calor del tubo que está en contacto con el agua hacia el agua

dQJ_tbJ =lambda T2,T4: dQIdeC(T2,T4)
#Flujo de calor del jugo hacia el tubo que está en contacto con el jugo

dQA_czA = 0
#Flujo de calor del agua hacia la coraza que est

dQczA_czat = 0 #Presuntamente todo esto es igual a 0.

dQtbJ_tbA = 0

NNu = lambda T2:(hJ(T2)*Dinttb)/(kTJ)

Septb = Pitb - Dexttb

AFcz = (Dintcz*Septb*SepBf)/Pitb

fD = lambda T2:(-2*np.log10((epsilon/Dinttb)/3.71-(5.02/NReJ(T2))*(np.log10((epsilon/(Dinttb*3.71)+(14.5/NReJ(T2)))))))**(-2) #que es epsilon?

hf1_2= lambda T2,T4: fD(T2)*((Ltb*(velJ(T2)**2))/(Dinttb*2));

hfcz= lambda T2: fD(T2)*(Dintcz*(NBf+1))/(Deqcz);

vAcz =lambda T4: ((fm1/pA(T4))/AFcz)


T1 = 323.15
T3 = 293.15
CPtb = 155
CPcz = 26

