from ecuacionesConstitutivas import *
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def F(A,t):
    T2, Ttb, T4, Tcz = A
    return [(1/(Cpj*Mj(T2)))*(fm1*Hj(T1)-fm1*Hj(T2)-dQJ_tbJ(T2,T4)),
            (1/(CPtb*Mtb))*(dQtbJ_tbA),
            (1/(Cpa(T4)*MA(T4)))*(fm3*HA(T3)-fm3*HA(T4)+dQtbA_A(T2,T4)-dQA_czA),
            (1/(CPcz*Mcz))*dQczA_czat]


W0 = [323.10,100,293.25,100]
t = np.linspace(0,400,1000)

sol = odeint(F, W0, t)

fig, axes = plt.subplots()
axes.plot(t, sol[:,0], 'r')
#axes.plot(t, sol[:,2], 'r')
plt.grid()
plt.show()