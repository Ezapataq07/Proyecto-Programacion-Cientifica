
import numpy as np
from RK4 import *
from ecuacionesConstitutivas import *

# La siguiente función es una función vectorial que contiene las ecuaciones diferenciales a resolver
def F(t, W):
    dT2 = (0.1/(Cpj*Mj(W[0])))*(fm1*Hj(T1)-fm1*Hj(W[0])-dQJ_tbJ(W[0],W[2]))
    dTtb = (1/(CPtb*Mtb))*(dQtbJ_tbA)
    dT4 = (0.01/(Cpa(W[2])*MA(W[2])))*(fm3*HA(T3)-fm3*HA(W[2])+dQtbA_A(W[0],W[2])-dQA_czA)
    dTcz = (1/(CPcz*Mcz))*dQczA_czat
    return np.array([dT2,dTtb,dT4,dTcz])

# Se definen las condiciones iniciales
W0 = np.array([T1-0.05,T1-0.05,T3+0.1,T3+0.1])

# Se resuelve el sistema por medio del método de RK4
t, Sol = RK4_Vec(F, W0, 0, 400, 20000)

# Se convierte las temperaturas de K a °C
Sol = (Sol - 273.15) 

# Se grafican las temperaturas de salida del jugo y del agua respectivamente
Plot_Component(0, t, Sol)
Plot_Component(2, t, Sol)