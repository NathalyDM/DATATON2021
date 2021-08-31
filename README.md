# DATATON2021
Emplear  las bases del Perú sobre los casos COVID, centros de vacunación, sospechosos, personas vacunadas, fallecidos, donaciones y presupuestos COVID ejecutados en las distintas regiones para determinar variables de interés que puedan ser tomadas en cuenta al momento de proponer programas para prevenir la tercera ola.

## Analisis de Data

### Análisis de Variables en R 
A continuación se mostrará la logica para el análisis de dependencia de variables. En primer lugar se crea un nuevo slot en la data que contenga 0 si para ese en ese mes una región haya tenido una mayor cantidad de vacunados y 1 de lo contrario. 

```bash
vacunados_dep_mes$statu<-1
for (i in 1:nrow(vacunados_dep_mes)){
    if (vacunados_dep_mes$Vacunados[i]>29733){
      vacunados_dep_mes$statu[i]<-0
    }
}
```

Posteriormente se realiza un análisis con survival : 

```
covariates <- c("car_por_sup", "temperatura", "avance")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(Mes, statu)~', x)))
                        
univ_models <- lapply( univ_formulas, function(x){coxph(x, data =  vacunados_dep_mes)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                          x <- summary(x)
                          p.value<-signif(x$wald["pvalue"], digits=2)
                          wald.test<-signif(x$wald["test"], digits=2)
                          beta<-signif(x$coef[1], digits=2);#coeficient beta
                          HR <-signif(x$coef[2], digits=2);#exp(beta)
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          HR <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                          res<-c(beta, HR, wald.test, p.value)
                          names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                        "p.value")
                          return(res)
                          #return(exp(cbind(coef(x),confint(x))))
                         })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)
```
Y finalmente se realiza un analisis de Cox los resultados de esta función pueden ser empleados posteriormente. 

```
res.cox <- coxph(Surv(Mes, statu) ~ avance+temperatura, data =  vacunados_dep_mes)
summary(res.cox)
```


## Implementación de Modelo SEIRD EN Python

## Implementación de Modelo SEIRD EN Python

### Librerías empleadas

```bash
Python 
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = "notebook"
%matplotlib inline
plt.style.use('ggplot')
import kaleido
```

### Construcción del Modelo en Python
El modelo fue construido guíandose de el Modelo SEIR del estudio "SEIR and Regression Model based COVID-19 outbreak predictions in India" (Pandey & Gaurav, 2020). 


![1_-TCB--QLTcyjQu9Mi30Ulg](https://user-images.githubusercontent.com/40121093/131426121-6dc7134f-1d8a-4c09-8eb9-ef8d24b80556.png)

- Beta: Ratio de infección
- Sigma: Ratio de incubación
- Gamma: Ratio de recuperación
- Mu : Tasa de motalidad

Las ecuaciones diferenciales que modelan dicho modelo son: 

![graphic-4](https://user-images.githubusercontent.com/40121093/131426124-690af8d3-f193-4537-b0cd-61327a9d5451.gif)

`Funciones Matemáticas`:

```markdown
def ode_model(z, t, beta, sigma, gamma, mu):
    """
    Reference https://www.idmod.org/docs/hiv/model-seir.html
    """
    S, E, I, R, D = z
    N = S + E + I + R + D
    dSdt = -beta*S*I/N
    dEdt = beta*S*I/N - sigma*E
    dIdt = sigma*E - gamma*I - mu*I
    dRdt = gamma*I
    dDdt = mu*I
    return [dSdt, dEdt, dIdt, dRdt, dDdt]


def ode_solver(t, initial_conditions, params):
    initE, initI, initR, initN, initD = initial_conditions
    beta, sigma, gamma, mu = params
    initS = initN - (initE + initI + initR + initD)
    res = odeint(ode_model, [initS, initE, initI, initR, initD], t, args=(beta, sigma, gamma, mu))
    return res

```

Entry format is described as follows:



## Adding Steps

Use the `@step` decorator. The main loop passes in the master file list and [jinja2][j2] environment.

## Referencias 
Pandey, Gaurav. (2020). SEIR and Regression Model-based COVID-19 outbreak predictions in India (Preprint). 10.2196/preprints.19406. 

