import pandas as pd
import numpy as np

##Para el sistema de ecuaciones diferenciales
from scipy.integrate import odeint 
np.random.seed(1)
#Para generar las muestras 
from pyDOE import lhs
from scipy.stats.distributions import norm

#Para la visualización
import seaborn as sns
import matplotlib.pyplot as plt 

from RK4 import *
from ecuacionesConstitutivas import *

# Se define el modelo general del Intercambiador de Calor
def IdeC_Model(t, W0):
  # La siguiente función es una función vectorial que contiene las ecuaciones diferenciales a resolver
  def F(A,t):
      T2, Ttb, T4, Tcz = A
      return [(1/(Cpj*Mj(T2)))*(fm1*Hj(T1)-fm1*Hj(T2)-dQJ_tbJ(T2,T4)),
              (1/(CPtb*Mtb))*(dQtbJ_tbA),
              (1/(Cpa(T4)*MA(T4)))*(fm3*HA(T3)-fm3*HA(T4)+dQtbA_A(T2,T4)-dQA_czA),
              (1/(CPcz*Mcz))*dQczA_czat]
  # Se realiza la solución del sistema usando una librería de SciPy
  sol = odeint(F, W0, t)
  return sol.T

#GENERAR LAS MUESTRAS:

  #PASO 1: Generación de las muestras aleatorias. 
# Se especifica el número de muestras
sample_num = 500

# Se toman muestras aleatorias para la temperatura del jugo a la entrada (T1) y para la temperatura del agua a la entrada (T3), suponiendo que 
# tienen una distribución normal
# La media de la distribución se toman los valores de temperatura especificados en el ejemplo
# Se toma una desviación estándar de 2 grados
T1_samples = np.random.normal(T1,2,sample_num)
T3_samples = np.random.normal(T3,2,sample_num)


  # Se pueden visualizar las muestras transformadas 
df_samples = pd.DataFrame({r"$T_1$":T1_samples,
                          r"$T_3$":T3_samples})

plt.rcParams["axes.labelsize"] = 15
sns.jointplot(data=df_samples, x=r"$T_1$", y=r"$T_3$", height=5);

plt.show()
#SIMULACIÓN MONTECARLO: Para este paso se crea un ciclo para predecir la evolución de S, I y R, para cada muestra de 𝛽 y 𝛾.

# Monte Carlo
# Se define el intervalo temporal
t = np.linspace(0,400,2000)
# Se crean las matrices que contienen las diferentes soluciones
T2_s, Ttb_s, T4_s, Tcz_s = (np.zeros((t.shape[0], sample_num)) for i in range(4))

# En el siguiente ciclo se soluciona el modelo para las diferentes muestras de flujos másicos
for i in range(sample_num):
    W0 = np.array([T1_samples[i]-0.05,100,T3_samples[i]+0.1,100])
    T2, Ttb, T4, Tcz = IdeC_Model(t,W0)
    
    T2_s[:,i] = T2
    Ttb_s[:,i] = Ttb
    T4_s[:,i] = T4
    Tcz_s[:,i] = Tcz


# Se extraen las variables de interés, en este caso la Temperatura mínima del Jugo a la salida
min_T2 = T2_s.min(axis=0)  
min_time = t[T2_s.argmin(axis=0)]

#Visualización de la incertidumbre: Para presentar la incertidumbre en la salida se presentan histogramas y gráficas de dispersión.
#El paquete Seaborn contiene la función jointplot, la cual permite mostrar simultáneamente las dos cosas

# Visualizar las incertidumbre en las salidas
outputs = pd.DataFrame({"Min $T_2$":min_T2,
                          "Occurrence Time (s)":min_time})

plt.rcParams["axes.labelsize"] = 15
h = sns.jointplot(data=outputs, x="Occurrence Time (s)", y="Min $T_2$", 
                  # xlim=(5,165), ylim=(-20,1120),
                  marginal_kws=dict(bins=20),
                  alpha=0.6,
                  height=6)
"""De la figura de arriba, se puede ver cómo las salidas están variando según los parámetros de entrada inciertos. Sin embargo en la animación a
continuación se puede observar con mayor claridad cómo varían las variables de salida en términos de las diferentes combinaciones de gamma y beta"""



# Se extraen las variables de interés, en este caso la temperatura máxima del Agua a la salida
max_T4 = T4_s.max(axis=0)  
max_time = t[T4_s.argmax(axis=0)]

#Visualización de la incertidumbre: Para presentar la incertidumbre en la salida se presentan histogramas y gráficas de dispersión.
#El paquete Seaborn contiene la función jointplot, la cual permite mostrar simultáneamente las dos cosas

# Visualizar las incertidumbre en las salidas
outputs = pd.DataFrame({"Max $T_4$":max_T4,
                          "Occurrence Time (s)":max_time})

plt.rcParams["axes.labelsize"] = 15
h = sns.jointplot(data=outputs, x="Occurrence Time (s)", y="Max $T_4$", 
                  # xlim=(5,165), ylim=(-20,1120),
                  marginal_kws=dict(bins=20),
                  alpha=0.6,
                  height=6)
"""De la figura de arriba, se puede ver cómo las salidas están variando según los parámetros de entrada inciertos. Sin embargo en la animación a
continuación se puede observar con mayor claridad cómo varían las variables de salida en términos de las diferentes combinaciones de gamma y beta"""
plt.show()