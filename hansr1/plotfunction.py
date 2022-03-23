import pandas as pd

df = pd.read_csv("teste.dat", delim_whitespace=True, index_col=0, on_bad_lines='skip')
df.plot(subplots=True)