""" En este programa se resuelven las ecuaciones diferenciales del problema utilizando una librería de SciPy
"""
from ecuacionesConstitutivas import *
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# La siguiente función es una función vectorial que contiene las ecuaciones diferenciales a resolver
def F(A,t):
    T2, Ttb, T4, Tcz = A
    return [(1/(Cpj*Mj(T2)))*(fm1*Hj(T1)-fm1*Hj(T2)-dQJ_tbJ(T2,T4)),
            (1/(CPtb*Mtb))*(dQtbJ_tbA),
            (1/(Cpa(T4)*MA(T4)))*(fm3*HA(T3)-fm3*HA(T4)+dQtbA_A(T2,T4)-dQA_czA),
            (1/(CPcz*Mcz))*dQczA_czat]


# Se definen las condiciones iniciales
W0 = np.array([T1-0.05,T1-0.05,T3+0.1,T3+0.1])
# Se define el intervalo de tiempo a solucionar
t = np.linspace(0,400,1000)

# se hace la solución del sistema utilizando la función odeint de SciPy
sol = odeint(F, W0, t)

# Se grafican las temperaturas de salida del jugo y del agua respectivamente
fig, axes = plt.subplots()
axes.plot(t, sol[:,0], 'r')
#axes.plot(t, sol[:,2], 'r')
plt.grid()
plt.show()