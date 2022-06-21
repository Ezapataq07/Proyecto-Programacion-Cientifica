Este conjunto de módulos contiene la implementación del modelo de Intercambiador de Calor presentado en "", se solucionan las ecuaciones
diferenciales importantes y se implementa el método Montecarlo para caracteriar el comportamiento de las variables de salida a cambios en las
entradas.

        *   RK4.py contiene la implementación del método Runge Kutta 4 vectorial 
        *   ecuacionesConstitutivas.py contiene la definición de todas las constantes y funciones necesarias para la solución del modelo
        *   solucion.py contiene la solución del modelo mediante el método RK4
        *   solucionSciPy.py contiene la solucion del modelo mediante una librería de solucionSciPy
        *   IncertidumbreTemperatura.py contiene la implementación del método Montecarlo para evaluar las temperaturas de salida según 
                                        diferentes Temperaturas de entradas
        *   IncertidumbreFlujo.py contiene la implementación del método Montecarlo para evaluar las temperaturas de salida según
                                  diferentes flujos de entradas

Para comprobar las implementaciones basta con ejecutar los últimos 4 archivos mencionados previamente, en ese orden