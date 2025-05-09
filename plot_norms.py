# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv("compare_norms.csv", sep=",", header=0, index_col=None)

# %%
data.groupby(['language', 'method', 'ic', 'boundary', 'dt_out']).count()

# %%
for i in data.loc[(data['language'] == "Julia") & (data['method'] == 'NMG') & (
        data["ic"] == '25p') & (data['boundary'] == 'neumann')]['pathname']:
    print(i)
# %%
