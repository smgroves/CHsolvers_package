# %%
from statsmodels.stats.diagnostic import lilliefors
from scipy import stats
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np
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
# %% FIGURE 2B
# bc = "neumann"
GridSize = 128
g = sns.catplot(
    all_data.loc[(all_data["GridSize"] == GridSize) &
                 #  (all_data["boundary"] == bc) &
                 (all_data["print"] == False)],
    kind="bar",
    y="elapsed_time(s)",
    x="method",
    hue="boundary",
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
    f"Elapsed Time for {GridSize}x{GridSize} Grid Size, 2000 Iterations", y=1.)
g.set_axis_labels("Solver", "Elapsed Time (log[sec])")
plt.tight_layout()
plt.show()
# plt.savefig(f"./output/both_bc_compare_runtime_{GridSize}_no_print.pdf")

# %% FIGURE 2B and 2C REDONE: gridsize 128, no print, both BC, all languages, NMG  or SAV only
GridSize = 128
method = "NMG"
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
g = sns.barplot(
    all_data.loc[(all_data["GridSize"] == GridSize) &
                 (all_data["method"] == method) &
                 (all_data["print"] == False)],
    y="elapsed_time(s)",
    x="language",
    hue="boundary",
    hue_order=["periodic", "neumann"],
    order=["Python", "MATLAB", "Julia"],
    errorbar=None,
    # log=True
    # palette="Set2"
)
sns.swarmplot(
    all_data.loc[(all_data["GridSize"] == GridSize) &
                 (all_data["method"] == method) &
                 (all_data["print"] == False)],
    y="elapsed_time(s)",
    x="language",
    hue="boundary",
    hue_order=["periodic", "neumann"],
    order=["Python", "MATLAB", "Julia"],
    dodge=True,
    # log=True
    # palette="Set2"
)
plt.ylim(1, 1e6)
g.figure.get_axes()[0].set_yscale('log')
plt.suptitle(
    f"Elapsed Time for {GridSize}x{GridSize} Grid Size, 2000 Iterations, {method}", y=1.)
plt.xlabel("Solver")
plt.ylabel("Elapsed Time (log[sec])")
plt.tight_layout()
# plt.show()
plt.savefig(
    f"./output/both_bc_compare_runtime_{GridSize}_{method}_no_print_Figure_2BC.pdf")
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


# %% FIGURE 2A: IC and endpoints for periodic and neumann

indir = "/Users/smgroves/Documents/GitHub/CHsolvers_package/output/output_Rivanna_Julia"
phi_name_MG_periodic = "NMG_Julia_2000_dt_5.50e-06_Nx_128_periodic_dtout_10_phi.csv"
phi_MG_p = np.genfromtxt(
    f"{indir}/{phi_name_MG_periodic}",
    delimiter=",",
)
phi_MG_p = phi_MG_p.reshape(-1, 128, 128).transpose(1, 2, 0)

phi_name_MG_neumann = "NMG_Julia_2000_dt_5.50e-06_Nx_128_neumann_dtout_10_phi.csv"
phi_MG_n = np.genfromtxt(
    f"{indir}/{phi_name_MG_neumann}",
    delimiter=",",
)
phi_MG_n = phi_MG_n.reshape(-1, 128, 128).transpose(1, 2, 0)


# %% save individual plots

# timepoints = [0, 10, 20, 100, 1000, 2000]
timepoints = [0, 200]
dt_out = 10
dt = 5.5e-6
for bc in ["periodic", 'neumann']:
    for timepoint in timepoints:
        normalize_phis = mcolors.TwoSlopeNorm(vcenter=0, vmin=-1, vmax=1)
        if bc == "periodic":
            phi_MG = phi_MG_p
            title = "Periodic"
        else:
            phi_MG = phi_MG_n
            title = "Neumann"

        s = sns.heatmap(
            phi_MG[:, :, timepoint],
            square=True,
            cmap=cm.RdBu_r,
            norm=normalize_phis,
            cbar=False,
            linewidths=0.0,
        )
        plt.xticks(ticks=[], labels=[])
        plt.yticks(ticks=[], labels=[])
        plt.title(f"Time= {timepoint*dt*dt_out}, {title} BC")
        plt.tight_layout()
        plt.savefig(
            f"{indir}/MG_Julia_2000_dt_5.5e-6_t_{title}_{timepoint*dt*dt_out:.2e}.pdf",
            bbox_inches="tight",
            pad_inches=0,
            dpi=300,
        )
        plt.close()

# %% Statistical test comparisons
GridSize = 128

df = all_data.loc[(all_data["GridSize"] == GridSize) &
                  (all_data["print"] == False)].reset_index()


# %%
# Log transform all the simulations in 2B, subtract the mean of each group, then qq plot the mean-centered data.


def log_transform_and_center(df):
    df['log_elapsed_time'] = np.log10(df['elapsed_time(s)'])
    grouped = df.groupby(['language', 'method', 'boundary'])
    df['mean_log_elapsed_time'] = grouped['log_elapsed_time'].transform('mean')
    df['centered_log_elapsed_time'] = df['log_elapsed_time'] - \
        df['mean_log_elapsed_time']
    return df


df = log_transform_and_center(df)
# %%
# Create a QQ plot for each pair across method
a = df.loc[df['method'] == 'NMG', 'centered_log_elapsed_time']
b = df.loc[df['method'] == 'SAV', 'centered_log_elapsed_time']
percs = np.linspace(0, 100, 21)
qn_a = np.percentile(a, percs)
qn_b = np.percentile(b, percs)

plt.plot(qn_a, qn_b, ls="", marker="o")
plt.plot([-.2, .2], [-.2, .2], color='red', ls='--')
plt.xlabel("Mean-Centered Log-transformed NMG Elapsed Time Quantiles")
plt.ylabel("Mean Centered Log-transformed SAV Elapsed Time Quantiles")
plt.title("QQ Plot of NMG vs SAV Elapsed Time")

# %%

# 2. Lillifors Test
lilliefors_test = lilliefors(a.values - b.values)
print("\nLilliefors Test:")
print("Statistic:", lilliefors_test[0])
print("P-value:", lilliefors_test[1])

# 3. Jarque-Bera Test
jarque_bera_test = stats.jarque_bera(a.values - b.values)
print("\nJarque-Bera Test:")
print("Statistic:", jarque_bera_test.statistic)
print("P-value:", jarque_bera_test.pvalue)

# 4. Anderson-Darling Test
anderson_darling_test = stats.anderson(a.values - b.values)
print("\nAnderson-Darling Test:")
print("Statistic:", anderson_darling_test.statistic)
print("Critical values:", anderson_darling_test.critical_values)
print("Significance levels:", anderson_darling_test.significance_level)
# %% t test for NMG vs SAV for Julia, neumann boundary, grid size 128
aa = df.loc[(df['method'] == 'NMG') & (df['language'] == "Julia")
            & (df['boundary'] == "neumann"), 'log_elapsed_time']
bb = df.loc[(df['method'] == 'SAV') & (df['language'] == "Julia")
            & (df['boundary'] == "neumann"), 'log_elapsed_time']
# %%
# Figure 2B t test
stats.ttest_ind(aa, bb)
# TtestResult(statistic=np.float64(1.7237699829831481e-15), pvalue=np.float64(0.9999999999999987), df=np.float64(4.0))

# %%
