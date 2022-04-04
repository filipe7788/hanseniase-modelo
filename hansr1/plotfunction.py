import pandas as pd
from matplotlib import pyplot as plt

df = pd.read_csv('sr1.txt', header=None, delimiter=r"\s+")
Pauci = df.iloc[:, 4].tolist()
beta = df.iloc[:, 0].tolist()

df2 = pd.read_csv('sr2.txt', header=None, delimiter=r"\s+")
PauciVac = df2.iloc[:, 4].tolist()
betaVac = df2.iloc[:, 0].tolist()

plt.plot(betaVac, Pauci, label='Paucibacilar')
plt.plot(betaVac, PauciVac, label='Paucibacilar com vacina')
plt.legend()
plt.show()