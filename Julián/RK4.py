"""
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
# Para ejemplo inicial y comprobación de la implementación se trabaja el siguiente problema:
#           y'(t) = 4(t^3)(arctan(2y(t))),   1<=t<=2
#           y(1) = 2



import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

layout = go.Layout(
    title="",    
    plot_bgcolor="#FFFFFF",
    hovermode="x",
    hoverdistance=100, # Distance to show hover label of data point
    spikedistance=1000, # Distance to show spike
    xaxis=dict(
        title="time",
        linecolor="#BCCCDC",
        showspikes=True, # Show spike line for X-axis
        # Format spike
        spikethickness=2,
        spikedash="dot",
        spikecolor="#999999",
        spikemode="across",
    ),
    yaxis=dict(
        title="",
        linecolor="#BCCCDC"
    )
)

def f(t, y):
    return 4*t**3*np.arctan(2*y)

def RK4(f, y0, a, b, N):
    N=N+1
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
    
    fig = go.Figure(data=[go.Scatter(x=t,y=y)],layout=layout)
    fig.show()
    return t, y

#RK4(f, 2, 1, 2, 1000)

""" Implementación de RK4 Vectorial 
"""



def RK4_Vec(F, W0, a, b, N):
    N=N+1
    h=(b-a)/(N-1)
    W = np.zeros(N, dtype=np.ndarray)
    W[0]=W0
    t = np.linspace(a,b,N)
    
    for i in range(N-1):
        k1=F(t[i],W[i])
        k2=F(t[i]+h/2,W[i]+k1*h/2)
        k3=F(t[i]+h/2,W[i]+k2*h/2)
        k4=F(t[i+1],W[i]+h*k3)

        W[i+1]=W[i]+(h/6)*(k1+2*k2+2*k3+k4)
    return t, W

def F(t, W):
    f0 = W[1]
    f1 = W[2]
    f2 = np.log(t**2+4)-np.sin(W[1])-np.exp(t)*W[0]
    return np.array([f0,f1,f2])

W0 = np.array([1,0,8])



def Plot_Component(component, t, W):
    y = np.zeros(len(W))
    for i in range(len(W)):
        y[i] = W[i][component]
    fig = go.Figure(data=[go.Scatter(x=t,y=y)],layout=layout)
    fig.show()

##t, Sol = RK4_Vec(F, W0, 3, 5, 100)
##Plot_Component(0,t,Sol)

