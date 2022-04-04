import pandas as pd
from matplotlib import pyplot as plt

def getLimiarBeta(hasVaccine, r):
    mi = 1/70
    sigmaP = 0.73 # sigmaP = 0.9
    omegaP = 1/4
    omegaM = 1/10
    rho = 0.05
    gama = 0.0 # gama = 0.75
    deltaM = 1/12
    deltaP = 1/12
    if hasVaccine:
        sigmaP = 0.9
        gama = 0.75

    nominator = (mi+(sigmaP*omegaP)+(1-sigmaP)*omegaM)
    denominator = (rho*(1-gama)*(omegaP*sigmaP*(mi+deltaM)+r*omegaM*(1-sigmaP)*(mi+deltaP)))

    return (nominator/denominator)

deltas = []

for i in range(30):
    delta = (getLimiarBeta(True, i) + getLimiarBeta(False, i))/getLimiarBeta(False, i)
    deltas.append(delta)

plt.plot(range(30), deltas)
plt.legend(['Delta', 'R'])
plt.show()