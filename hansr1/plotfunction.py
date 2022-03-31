import pandas as pd
from matplotlib import pyplot as plt

df = pd.read_csv('sr1.txt', header=None, delimiter=r"\s+")
betaP = df.iloc[:, 4].tolist()
betaM = df.iloc[:, 5].tolist()
beta = df.iloc[:, 0].tolist()

plt.plot(beta, betaP, label='betaP')
plt.plot(beta, betaM, label='betaM')
plt.legend()
plt.show()