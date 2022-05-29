"""(
El siguiente programa implementa el método numérico Runge Kutta 4 para solucionar ecuaciones difereciales del tipo
PVI:
        y'(t) = f(t, y(t))     para        a<=t<=b
        y(a) = y0

Dicho método es un algoritmo iterativo en el cuál se divide el dominio en N parte, y haciendo h=(b-a)/N se tienen
las siguientes definiciones:

    y_i+1 = y_i + (h/6)(k_1 + 2k_2 + 2k_3 + k_4), con:

    k_1 = f(t_i, y_i)
    k_2 = f(t_i + h/2, y_i + (h/2)k_1)
    k_3 = f(t_i + h/2, y_i + (h/2)k_2)
    k_4 = f(t_i+1, y_i + hk_3)
"""
#Para ejemplo inicial y comprobación de la implementación se trabaja el siguiente problema:
#           y'(t) = 4(t^3)(arctan(2y(t))),   1<=t<=2
#           y(1) = 2


import numpy as np

def f(t, y):
    return 4*t**3*np.arctan(2*y)

def RK4(f, y0, a, b, N):
    h=(b-a)/(N-1)
    y=np.zeros(N)
    k1=np.zeros(N)
    k2=np.zeros(N)
    k3=np.zeros(N)
    k4=np.zeros(N)

    t=np.linspace(a,b,N)
    y[0]=y0

    for i in range(N-1):
        k1=f(t[i],y[i])
        k2=f(t[i]+h/2,y[i]+k1*h/2)
        k3=f(t[i]+h/2,y[i]+k2*h/2)
        k4=f(t[i+1],y[i]+h*k3)

        y[i+1]=y[i]+(h/6)*(k1+2*k2+2*k3+k4)
    
    print(t)
    print(y)
    return y

RK4(f, 2, 1, 2, 4)




