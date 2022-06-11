import numpy as np
import math
fm1 = 10 #Flujo masico
fm3 = 10 # Flujo masico
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
T2=10
T4=10

wtb = 10 #Espesor del tubos

Cpj = 4.187*(1 - 0.006*gB) #Capacidad calorifica especifica del jugo.

T = 273.15 #Grados Kelvin.

Cpa = 7.2951406*(10**(-8))*(T**3) - 7.431358269*(10**(-5))*(T**2) + 0.02620653*T + 1.0156045
#Capacidad calorifica del agua.

pj = 1179.7 - 0.354*T  #Densidad del jugo.

pA = -0.004323394923*((T-273.15)**2) - 0.04038343824*(T-273.15) + 1000.807908 #Densidad del agua.

Ntb = 10 #Número de tubos.
Ltb = 10 #Longitud de los tubos.
Dinttb = 10 #Diametro interno de los tubos.
Pitb = 10


Vj = Ntb*Ltb*((np.pi*(Dinttb**2))/4) #Volumen del jugo.

Lcz = 10 #Longitud de la coraza.
Dintcz = 10 #Diametro interno de la coraza.
Dexttb = 10 #Diametro externo de los tubos.

VA = Ntb*Lcz*((np.pi*(Dintcz**2))/4) - Ntb*Ltb*((np.pi*(Dexttb**2))/4)

Mj = pj*Vj #Masa del jugo.

MA = pA*VA #Masa del agua.

Mtb = 10 #Masa del tubo.
Mcz = 10 #Masa de la coraza.
kTA = 10

Hj = 5.46829932*(10**(-4))*((T-273.15)**2) + 3.694618471*(T-273.15)-0.633984 #Entalpia del jugo.
Deqcz = 4((((Pitb**2)*math.sqrt(3))/4) - np.pi*(Dexttb**2))/((np.pi*Dexttb)/2)
muA = (0.0002471250997*(T-273.15)**2)-0.03226252622(T-273.15)+1.545370351;
muJ = (0.047)/(T-273.15)
NReA = (pA*Deqcz*VA)/(muA)
AFtb = (np.pi*(Dinttb**2)/(4))
velJ = (fm1/pj)/(AFtb)

NReJ = (pj*Vj*Dinttb)/muJ
NPrA = (Cpa*muA)/(kTA)
hA = (kTJ*Deqcz)*(0.36*(NReA)**(0.55)*NPrA**(0.333))
HA = 4.180358722*(T-273.15) + 0.3506658477 #Entalpia del agua.
AHtbJ = np.pi*Dinttb*Ltb
AHtbA = np.pi*Dexttb*Ltb
AH = (1/2)*(AHtbJ+ AHtbA)
deltaT1 = TeJ - T4
deltaT2 = T2 - TeA

deltaTML = (deltaT1-deltaT2)/(math.log(deltaT1/deltaT2))

dQIdeC = U*AH*deltaTML
##dQtbA_A = hA*AHtbA*(TtbA-TA); 
dQtbA_A = dQIdeC;
NPrj = (Cpj*muJ)/(kTJ)
hJ = (kTJ/Dinttb)*(0.0243*(NReJ)**(0.8)*NPrj*0.333)
AHtb = (AHtbJ + AHtbA)/2

U = 1/(AH*((1/(hJ*AHtbJ))+(wtb/(KTtb*AHtb))+(1/(hA*AHtbA))))
##dQJ_tbJ = hJ*AHtbJ*(TJ-TtbJ);
dQJ_tbJ = dQIdeC
dQtbJ_tbA = dQJ_tbJ - dQtbA_A

dQA_czA = 0

dQczA_czat = 0 #Presuntamente todo esto es igual a 0.

dQtbJ_tbA = 0

NNu = (hJ*Dinttb)/(kTJ)

Septb = Pitb - Dexttb
AFcz = (Dintcz*Septb*SepBf)/Pitb
fD = (-2*math.log10((epsilon/Dinttb)-(5.02/NReJ)*(math.log10((epsilon/(Dinttb*3.71)+(14.5/NReJ))))))**(-2) #que es epsilon?
hf1_2=fD*((Ltb*(velJ)**2)/(Dinttb*2));
hfcz=fD*(Dintcz*(NBf+1))/(Deqcz);
vAcz = ((fm1/pA)/AFcz)