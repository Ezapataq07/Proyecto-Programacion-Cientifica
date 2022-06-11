"""Método MonteCarlo para el cálculo de incertidumbres: Ejemplo de aplicación
Caso de estudio: En este caso se modela la propagación de de una enfermedad en la población. En donde 
S(t): número de individuos quienes son susceptibles a adquirir la enfermedad pero todavía no están infectados, 
I(t): número de individuos infectados 
R(t): número de individuos quienes se han recuperado y son inmunes a la enfermedad.
Así, se describe la evolución temporal de S(t), I(t) y R(t) con un sistema de ecuaciones diferenciales ordinarias:

dS/dt= -(beta)SI/N
dI/dt= (beta)SI/N - aI
dR/dt= (gamma)I

Donde gamma representa la tasa de recuperación y beta denota la tasa de infección y N es el número total de individuos
en la población.

Por lo general la estimación de beta y gamma contiene incertidumbre. Para este ejemplo, se asume que gamma y beta ya han sido
determinados y que siguen una distribución de probabilidad normal bivariada:

Mean values (beta,gamma)=(0.22,0.1), Variance values (beta,gamma)=(2e-4, 1e-4) y covariance value= 4e-5

Y se asume que I(0)=8, R(0)=0 y N=1000

El algoritmo a seguir es el siguiente: 
1. Generar un número grande de muestras de gamma y beta de su distribución de probabilidad
2. Para cada muestra insertar gamma y beta en el sistema de ecuaciones diferenciales y resolverlo para obtener 
las salidas de interés
3. Basados en las predicciones se estiman las cantidades estadísticas"""

#Se importan los paquetes necesarios:

##Para análisis de datos
import pandas as pd
import numpy as np

##Para el sistema de ecuaciones diferenciales
from scipy.integrate import odeint 

#Para generar las muestras 
from pyDOE import lhs
from scipy.stats.distributions import norm

#Para la visualización
import seaborn as sns
import matplotlib.pyplot as plt 

#Creación de animaciones 
from matplotlib import animation
from IPython.display import HTML
from celluloid import Camera

def SIR_model(beta, gamma, t, N, I0, R0):
#La función resuelve el sistema de ecuaciones para simular la propagación de la infección
    
    S0 = N - I0 - R0
    
    def deriv(y, t, N, beta, gamma):
        S, I, R = y
        dSdt = -beta * S * I / N
        dIdt = beta * S * I / N - gamma * I
        dRdt = gamma * I
        return dSdt, dIdt, dRdt
    
    ret = odeint(deriv, [S0, I0, R0], t, args=(N, beta, gamma))
    
    return ret.T 

""" La función tiene como argumentos la tasa de infección 𝛽, y la tasa de recuperación 𝛾, una malla de puntos de tiempo t (en días) para calcular 
la evolución de la pandemia, el tamaño de la población N, y los casos iniciales infectados y recuperados I0 y R0, respectivamente."""

#GENERAR LAS MUESTRAS:

  #PASO 1: Generación de las muestras aleatorias. Para esto se usa el método Latin Hypercube Sampling, este método se usa desde el paquete pyDOE.
  #Después de crear las muestras, se deben transformar ya que las muestras generadas con pyDOE.lhs están distribuidas uniformemente entre (0,1). 
  #Por esto, se deben transformar para que queden acordes a la distribución normal bivariada.

  # Se especifica el número de muestras
sample_num = 1000

  # Se usa LHS para generar las muestras 
uni_samples = lhs(n=2, samples=sample_num, criterion='maximin') 
"""lhs tiene 3 criterios: número de parámetros (integer), número de muestras a generar (integer), 
criterion es el criterio con el que se separan los puntos (string)"""

"""El array de 2D tiene 1000 filas y 2 columnas, cada columna contiene 1000 valores aleatorios sacados de una distribución uniforme U(0,1)"""

  # Transformando las muestras a la distribución normal estandar 
std_norm_samples = np.zeros_like(uni_samples)
for i in range(2):
    std_norm_samples[:,i] = norm(loc=0, scale=1).ppf(uni_samples[:,i])

  #PASO 2: Convertir la distribución normal estandar en la distribución normal bivariada asumida inicialmente

  # Se especifica la distribución de los parámetros
mean_val = np.array([0.22,0.1]) 
cov_val = np.array([[2e-4, 4e-5],[4e-5,1e-4]]) #matriz de covarianza (columna 1 beta, columna dos gamma)

  # Se transforma en la distribución especificada
L = np.linalg.cholesky(cov_val) #teorema 
target_samples = mean_val + std_norm_samples.dot(L.T) #L es una matriz triangular obtenida con cholesky
beta_samples, gamma_samples = target_samples[:,0], target_samples[:,1] 

  # Se pueden visualizar las muestras transformadas 
df_samples = pd.DataFrame({r"$\beta$":beta_samples,
                          r"$\gamma$":gamma_samples})

plt.rcParams["axes.labelsize"] = 15
sns.jointplot(data=df_samples, x=r"$\beta$", y=r"$\gamma$", height=5);
#SIMULACIÓN MONTECARLO: Para este paso se crea un ciclo para predecir la evolución de S, I y R, para cada muestra de 𝛽 y 𝛾.

# Monte Carlo
t = np.arange(0, 101, 2)
susceptible, infection, recovery = (np.zeros((t.shape[0], sample_num)) for i in range(3))

for i in range(sample_num):
    S, I, R= SIR_model(beta_samples[i], gamma_samples[i],
                       t, I0=8, N=1000, R0=0)
    
    susceptible[:,i] = S
    infection[:,i] = I
    recovery[:,i] = R

# Se extraen las variables de interés 
max_infection = infection.max(axis=0)  
max_time = t[infection.argmax(axis=0)]

#Visualización de la incertidumbre: Para presentar la incertidumbre en la salida se presentan histogramas y gráficas de dispersión.
#El paquete Seaborn contiene la función jointplot, la cual permite mostrar simultáneamente las dos cosas

# Visualizar las incertidumbre en las salidas
outputs = pd.DataFrame({"Max Infection":max_infection,
                          "Occurrence Day":max_time})

plt.rcParams["axes.labelsize"] = 15
h = sns.jointplot(data=outputs, x="Occurrence Day", y="Max Infection", 
                  # xlim=(5,165), ylim=(-20,1120),
                  marginal_kws=dict(bins=12),
                  alpha=0.6,
                  height=6)
"""De la figura de arriba, se puede ver cómo las salidas están variando según los parámetros de entrada inciertos. Sin embargo en la animación a
continuación se puede observar con mayor claridad cómo varían las variables de salida en términos de las diferentes combinaciones de gamma y beta"""