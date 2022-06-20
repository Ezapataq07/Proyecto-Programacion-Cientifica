
import numpy as np
from RK4 import *
from ecuacionesConstitutivas import *

def F(t, W):
    f0 = (1/(Cpj*Mj(W[0])))*(fm1*Hj(T1)-fm1*Hj(W[0])-dQJ_tbJ(W[0],W[2]))
    f1 = (1/(CPtb*Mtb))*(dQtbJ_tbA)
    f2 = (1/(Cpa(W[2])*MA(W[2])))*(fm3*HA(T3)-fm3*HA(W[2])+dQtbA_A(W[0],W[2])-dQA_czA)
    f3 = (1/(CPcz*Mcz))*dQczA_czat
    return np.array([f0,f1,f2,f3])

W0 = np.array([323.10,100,293.25,100])

t, Sol = RK4_Vec(F, W0, 0, 400, 20000)

Sol = (Sol - 273.15) 

Plot_Component(0, t, Sol)
Plot_Component(2, t, Sol)