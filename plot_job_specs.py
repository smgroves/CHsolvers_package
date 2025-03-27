# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['pdf.fonttype'] = 'TrueType'
# %%
data = pd.read_csv("Job_specs.csv", sep=",", header=0, index_col=None)

data = data.loc[data["timesteps"] == 2000]
# %%
data = data.loc[data['solver'] != "FD"]
# %%
g = sns.catplot(
    data.loc[data["nx"] == 128],
    kind="bar",
    y="time (secs)",
    x="solver",
    hue="bc",
    col="language",
    height=4,
    aspect=0.5,
    # log=True
    palette="Set2"
)
g.figure.get_axes()[0].set_yscale('log')
# plt.show()
plt.savefig("./output/compare_runtime_neumann_periodic.pdf")
# %%
sns.lineplot(
    data,
    y="time (secs)",
    x="nx",
    hue="solver",
    style="language",
)
plt.show()

# %%
