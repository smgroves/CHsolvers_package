# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['pdf.fonttype'] = 'TrueType'
# %%
data = pd.read_csv("Job_specs.csv", sep=",", header=0, index_col=None)
data_rivanna = pd.read_csv("Job_specs_rivanna.csv",
                           sep=",", header=0, index_col=None)

data['comp'] = "local"
data_rivanna['comp'] = "rivanna"

# %%
for i, r in data_rivanna.loc[data_rivanna['t_iter'] == 20].iterrows():
    new_r = r.copy()
    new_r['t_iter'] = new_r['t_iter'] * 100
    new_r['elapsed_time(s)'] = new_r['elapsed_time(s)'] * 100
    new_r['comp'] = "rivanna_extrap"
    data_rivanna.loc[data_rivanna.shape[0]+1] = new_r

# %%
all_data = pd.concat([data, data_rivanna], join='outer', ignore_index=True)
# %%
all_data = all_data.loc[all_data["t_iter"] == 2000]
all_data = all_data.loc[all_data['solver'] != "FD"]

# %%
all_data = all_data.loc[~((all_data["language"] == "Julia")
                          & (all_data["comp"] == "local"))]
# %%
bc = "periodic"
GridSize = 512
g = sns.catplot(
    all_data.loc[(all_data["GridSize"] == GridSize) &
                 (all_data["boundary"] == bc)],
    kind="bar",
    y="elapsed_time(s)",
    x="solver",
    # hue="bc",
    col="language",
    height=4,
    aspect=0.5,
    col_order=["Python", "MATLAB", "Julia"],
    # log=True
    # palette="Set2"
)
plt.ylim(1, 1e6)
g.figure.get_axes()[0].set_yscale('log')
plt.suptitle(
    f"Elapsed Time for {GridSize}x{GridSize} Grid Size, 2000 Iterations, {bc.capitalize()} BC", y=1.)
g.set_axis_labels("Solver", "Elapsed Time (log[sec])")
plt.tight_layout()
# plt.show()
plt.savefig(f"./output/compare_runtime_{bc}_{GridSize}.pdf")
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
