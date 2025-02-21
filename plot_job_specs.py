# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# %%
data = pd.read_csv("Job_specs.csv", sep=",", header=0, index_col=None)

data = data.loc[data["timesteps"] == 2000]
sns.catplot(
    data,
    kind="bar",
    y="time (secs)",
    x="solver",
    hue="dt",
    col="language",
    height=4,
    aspect=0.5,
)
plt.show()
# %%
