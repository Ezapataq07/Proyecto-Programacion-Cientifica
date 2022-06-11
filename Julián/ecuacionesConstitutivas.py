import numpy as np
import math
fm1 = 10 #Flujo masico
fm3 = 10 # Flujo masico

gB = 10 #Grados Brix

Cpj = 4.187*(1 - 0.006*gB) #Capacidad calorifica especifica del jugo.

T = 273.15 #Grados Kelvin.

Cpa = 7.2951406*(10**(-8))*(T**3) - 7.431358269*(10**(-5))*(T**2) + 0.02620653*T + 1.0156045
#Capacidad calorifica del agua.

pj = 1179.7 - 0.354*T  #Densidad del jugo.

pA = -0.004323394923*((T-273.15)**2) - 0.04038343824*(T-273.15) + 1000.807908 #Densidad del agua.

Ntb = 10 #Número de tubos.
Ltb = 10 #Longitud de los tubos.
Dinttb = 10 #Diametro interno de los tubos.

Vj = Ntb*Ltb*((np.pi*(Dinttb**2))/4) #Volumen del jugo.

Lcz = 10 #Longitud de la coraza.
Dintcz = 10 #Diametro interno de la coraza.
Dexttb = 10 #Diametro externo de los tubos.

VA = Ntb*Lcz*((np.pi*(Dintcz**2))/4) - Ntb*Ltb*((np.pi*(Dexttb**2))/4)

Mj = pj*Vj #Masa del jugo.

MA = pA*VA #Masa del agua.

Mtb = 10 #Masa del tubo.
Mcz = 10 #Masa de la coraza.

Hj = 5.46829932*(10**(-4))*((T-273.15)**2) + 3.694618471*(T-273.15)-0.633984 #Entalpia del jugo.

HA = 4.180358722*(T-273.15) + 0.3506658477 #Entalpia del agua.

dQtbJ_tbA = dQJ_tbJ - dQtbA_A

dQczA_czat = dQA_czA - dQczat_at #Presuntamente todo esto es igual a 0.

dQJ_tbJ = hj*AHtbJ*(TJ-TtbJ);  dQJ_tbJ = dQIdeC;

dQtbJ_tbA = (KTtb/wtb)*AHtb*(TtbJ-TtbA)

dQtbA_A = hA*AHttbA*(TtbA-TA); dQtbA_A = dQIdeC;

dQA_czA = hA*AHczA*(TA - TczA)

dQA_czA = (kTcz/wcz)*AHcz*(TA-Tczat)

dQIdeC = u*AH*deltaTML

deltaTML = (deltaT1-deltaT2)/(math.log(deltaT1/deltaT2))

deltaT1 = TeJ - T4

Tej = 10

deltaT2 = T2 - TeA

TeA = 10

U = 1/(AH*((1/(hJ*AHtbJ))+(wtb/(KTtb*AHtb))+(1/(hA*AHtbA))))

AH = (1/2)*(AHtb+ AHtBA)

AHtbJ = np.pi*Dinttb*Ltb
AHtbA = np.pi*Dexttb*Ltb

Dexttb = 10

AHtb = (AHtbJ + AHtbA)/2

KTtb = 10 #Conductividad termica

wtb = 10 #Espesor del tubos

NNu = 0.0243*NRe08*NPr0333*(mu/mu**w)**0.14

NNu = (hJ*Dinttb)/(kTJ)

hJ = (kTJ/Dinttb)*(0.0243*NRe08*NPr0333)

kTJ = 10

NReJ = (pJ*vJ*Dinttb)/muJ

vJ = (dmJ/pJ)/(AFtb)

AFtb = (np.pi*(Dinttb**2)/(4))

muJ = (0.047)/(T-273.15)

NPr = (CPJ*muJ)/(kTJ)

hA = (kTJ*Deqcz)*(0.36*NReA055*NPrA0333)

Deqcz = 4((((Pitb**2)*math.sqrt(3))/4) - np.pi*(Dexttb**2))/((np.pi*Dexttb)/2)

Pitb = 10

AFcz = (Dintcz*Septb*SepBf)/Pitb

Dintcz = 10

SepBf = 10

Septb = Pitb - Dexttb
muA = (0.0002471250997*(T-273.15)**2)-0.03226252622(T-273.15)+1.545370351;
hf1_2=fD*((Ltb*(vjtb)**2)/(Dinttb*2));
fD = (-2*math.log10((epsilon/Dinttb)-(5.02/NReJ)*(math.log10((epsilon/(Dinttb*3.71)+(14.5/NReJ))))))**(-2); #que es epsilon?
hfcz=fD*((Dintcz*(NBf+1))/(Deqcz*(muA/muA**w)**0.14));
NBf=10;
vAcz = ((fm1/pA)/AFcz);