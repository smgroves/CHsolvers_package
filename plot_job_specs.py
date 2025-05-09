# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['pdf.fonttype'] = 'TrueType'
# %%
# data = pd.read_csv("Job_specs.csv", sep=",", header=0, index_col=None)
data_rivanna = pd.read_csv("Job_specs_rivanna.csv",
                           sep=",", header=0, index_col=None)

# data['comp'] = "local"
data_rivanna['comp'] = "rivanna"

# %%
for length in [20, 200]:
    for i, r in data_rivanna.loc[data_rivanna['t_iter'] == length].iterrows():
        new_r = r.copy()
        new_r['t_iter'] = new_r['t_iter'] * 2000//length
        new_r['elapsed_time(s)'] = new_r['elapsed_time(s)'] * 2000//length
        new_r['comp'] = "rivanna_extrap"
        data_rivanna.loc[data_rivanna.shape[0]+1] = new_r

# %%
# all_data = pd.concat([data, data_rivanna], join='outer', ignore_index=True)
all_data = data_rivanna.copy()
# %%
all_data = all_data.loc[all_data["t_iter"] == 2000]
all_data = all_data.loc[all_data['method'] != "FD"]
all_data = all_data.loc[all_data['comp'] != "local"]

# %% Check if all expected cases ran on Rivanna
all_data = all_data.sort_values(
    by=['method', 'language', 'GridSize', 'print', 'boundary', 'IC_percent_+1'])
for language in all_data['language'].unique():
    for method in all_data['method'].unique():
        for GridSize in all_data['GridSize'].unique():
            for bc in all_data['boundary'].unique():
                for IC in all_data['IC_percent_+1'].unique():
                    for print_data in all_data['print'].unique():
                        # print(bc, method, GridSize, language, IC)
                        if all_data.loc[(all_data["GridSize"] == GridSize) &
                                        (all_data["boundary"] == bc) &
                                        (all_data["method"] == method) &
                                        (all_data["language"] == language) &
                                        (all_data["IC_percent_+1"] == IC) &
                                        (all_data["print"] == print_data)].shape[0] == 0:
                            print(
                                f"Missing {language} {method} {GridSize} {bc} {IC} {print_data}")
# (all_data.groupby(['language', 'method', 'GridSize', 'print', 'boundary',
#  'IC_percent_+1']).count()['date']).to_csv("Job_specs_Rivanna_counts.csv")
# %%
    for method in all_data['method'].unique():
        for GridSize in all_data['GridSize'].unique():
            for bc in all_data['boundary'].unique():
                for IC in all_data['IC_percent_+1'].unique():
                    for print_data in all_data['print'].unique():
                        # print(bc, method, GridSize, language, IC)
                        if all_data.loc[(all_data["GridSize"] == GridSize) &
                                        (all_data["boundary"] == bc) &
                                        (all_data["method"] == method) &
                                        (all_data["language"] == "Python") &
                                        (all_data["IC_percent_+1"] == IC) &
                                        (all_data["print"] == print_data)].shape[0] == 0:
                            print(
                                f"{GridSize} {bc} {method} {print_data} {IC}")
# %%
# all_data = all_data.loc[~((all_data["language"] == "Julia")
#   & (all_data["comp"] == "local"))]
# %%
bc = "neumann"
GridSize = 128
g = sns.catplot(
    all_data.loc[(all_data["GridSize"] == GridSize) &
                 (all_data["boundary"] == bc) &
                 (all_data["print"] == False)],
    kind="bar",
    y="elapsed_time(s)",
    x="method",
    # hue="bc",
    col="language",
    height=4,
    aspect=0.4,
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
plt.savefig(f"./output/compare_runtime_{bc}_{GridSize}_no_print.pdf")
# %%
# %%
bc = "periodic"
GridSize = 128
g = sns.catplot(
    all_data.loc[(all_data["GridSize"] == GridSize) &
                 (all_data["boundary"] == bc) &
                 (all_data["print"] == False)],
    kind="bar",
    y="elapsed_time(s)",
    x="method",
    # hue="bc",
    col="language",
    height=4,
    aspect=0.4,
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
plt.savefig(f"./output/compare_runtime_{bc}_{GridSize}_no_print.pdf")

# %% compare printing to no print
bc = "neumann"
GridSize = 512
g = sns.catplot(
    all_data.loc[(all_data["GridSize"] == GridSize) &
                 (all_data["boundary"] == bc)],
    kind="bar",
    y="elapsed_time(s)",
    x="method",
    hue="print",
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
plt.savefig(f"./output/compare_runtime_{bc}_{GridSize}_printing_effect.pdf")

# %% compare grid sizes
bc = "neumann"
GridSizes = [512, 256, 128, 64]
print_data = True
fig, ax = plt.subplots(1, 3, figsize=(9, 5))
for i, language in enumerate(["Python", "MATLAB", "Julia"]):
    for gr in GridSizes:
        g = sns.barplot(
            data=all_data.loc[(all_data["GridSize"] == gr) &
                              (all_data["boundary"] == bc) &
                              (all_data["print"] == print_data) &
                              (all_data["language"] == language)],
            y="elapsed_time(s)",
            x="method",
            alpha=1,
            log=True,
            label=f"{gr}x{gr}",
            ax=ax[i],
        )
    ax[i].set_ylim(1, 1e7)
    ax[i].set_title(language)
    # plt.set_yscale('log')
plt.suptitle(
    f"Elapsed Time for 2000 Iterations, {bc.capitalize()} BC, print to file: {print_data}", y=1.)
plt.tight_layout()
plt.legend()
# plt.show()
plt.savefig(f"./output/compare_runtime_{bc}_compare_gridsize_with_print.pdf")
# %%
bc = "neumann"
GridSize = 512
g = sns.swarmplot(
    all_data.loc[(all_data["GridSize"] == GridSize) &
                 (all_data["boundary"] == bc)],
    y="elapsed_time(s)",
    hue="language",
    x="print",

)
g.figure.get_axes()[0].set_yscale('log')

plt.show()

# %%
g = sns.lineplot(
    all_data,
    y="elapsed_time(s)",
    x="GridSize",
    hue="method",
    style="language",
)
g.figure.get_axes()[0].set_yscale('log')
g.figure.get_axes()[0].set_xscale('log')

plt.show()

# %%

bc = "periodic"
g = sns.catplot(
    all_data.loc[
        (all_data["boundary"] == bc) &
        (all_data['language'] == "Julia")
    ],
    kind="bar",
    y="mem_allocated(MB)",
    x="solver",
    # # hue="bc",
    col="GridSize",
    # height=4,
    aspect=0.35,
)
# plt.ylim(1, 1e7)
g.figure.get_axes()[0].set_yscale('log')

plt.suptitle(
    f"Memory Allocated for 2000 Iterations, {bc.capitalize()} BC, Julia", y=1.)
g.set_axis_labels("Solver", "Memory Allocated (MB)")
# plt.tight_layout()
# plt.show()
plt.savefig(f"./output/compare_memalloc_{bc}_Julia.pdf")


# %%
