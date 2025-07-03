# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
data = pd.read_csv("compare_norms_v3.csv", sep=",", header=0, index_col=None)


# %%
# want to group plots so that each plot shows a different ic and boundary and GridSize
# language, method, dt_out should be overlaid
# x axis = timepoint
# y axis = norm
# want to plot all norms on the same plot
# subset = data.loc[data['language'] != "Python"]
for ic in sorted(data["ic"].unique()):
    fig, ax = plt.subplots(nrows=len(data["boundary"].unique()),
                           ncols=len(data["GridSize"].unique()),
                           figsize=(16, 8))
    for i, boundary in enumerate(sorted(data["boundary"].unique())):
        for j, grid_size in enumerate(sorted(data["GridSize"].unique())):
            subset = data[
                (data["ic"] == ic)
                # & (data["language"] != "Python")
                & (data["boundary"] == boundary)
                & (data["GridSize"] == grid_size)
                & (data["norm"] != 0)
            ]
            sns.lineplot(
                data=subset,
                x="timepoint",
                y="norm",
                hue="language",
                style="method",
                palette="tab10",
                markers=True,
                dashes=False,
                ax=ax[i, j],
            )
            ax[i, j].set_title(
                f"ic={ic}, bc={boundary}, nx={grid_size}")
            ax[i, j].set_xlabel("Timepoint")
            ax[i, j].set_ylim(0, np.max(subset["norm"]+5))
            ax[i, j].set_ylabel("Norm")
    plt.legend(title="Language and Method")
    plt.suptitle(f"Norms for ic={ic}", fontsize=16)
    plt.tight_layout()

    plt.show()

# %%
l2 = data = pd.read_csv("L2_error_python.csv", sep=",",
                        header=0, index_col=None)
print(l2.head())

# for ic in sorted(l2["ic"].unique()):
# ic accidentally got smushed together with dt_out so we want 200050, 200075, and 200025
# for ic in sorted(l2["dt_out"].unique()):
#     if ic in [200050, 200075, 200025]:
#         continue
#     else:
#         print(ic)
fig, ax = plt.subplots(nrows=len(l2["boundary"].unique()),
                       ncols=len(l2["GridSize"].unique()),
                       figsize=(16, 8))
for i, boundary in enumerate(sorted(l2["boundary"].unique())):
    for j, grid_size in enumerate(sorted(l2["GridSize"].unique())):
        for ic in sorted(l2["ic"].unique()):

            subset = l2[
                (l2["ic"] == ic)  # get rid of the first 4 characters
                # & (data["language"] != "Python")
                & (l2["boundary"] == boundary)
                & (l2["GridSize"] == grid_size)
                & (l2["dt_out"] == 1)
            ]
            print(subset.head())
            sns.scatterplot(
                data=subset,
                x="timepoint",
                y="l2_error",
                markers="o",
                edgecolor=None,
                # dashes=False,
                label=f"{ic}",
                ax=ax[i, j],
            )
            ax[i, j].set_title(
                f"bc={boundary}, nx={grid_size}")
            ax[i, j].set_xlabel("Timepoint")
            # ax[i, j].set_ylim(0, np.max(subset["l2_error"]))
            ax[i, j].set_ylabel("L2(NMG- SAV)")
            plt.legend()
plt.suptitle(
    f"L2 Error between NMG and SAV", fontsize=16)
plt.tight_layout()

plt.show()

# %%
julia = data = pd.read_csv("L2_error_Julia.csv", sep=",",
                           header=0, index_col=None)
print(julia.head())

# for ic in sorted(l2["ic"].unique()):
# ic accidentally got smushed together with dt_out so we want 200050, 200075, and 200025
# for ic in sorted(l2["dt_out"].unique()):
#     if ic in [200050, 200075, 200025]:
#         continue
#     else:
#         print(ic)
fig, ax = plt.subplots(nrows=len(julia["boundary"].unique()),
                       ncols=len(julia["GridSize"].unique()),
                       figsize=(16, 8))
for i, boundary in enumerate(sorted(julia["boundary"].unique())):
    for j, grid_size in enumerate(sorted(julia["GridSize"].unique())):
        for ic in sorted(julia["ic"].unique()):

            subset = julia[
                (julia["ic"] == ic)  # get rid of the first 4 characters
                # & (data["language"] != "Python")
                & (julia["boundary"] == boundary)
                & (julia["GridSize"] == grid_size)
                # & (julia["dt_out"] == 10)
            ]
            sns.scatterplot(
                data=subset,
                x="timepoint",
                y="l2_error",
                markers="o",
                edgecolor=None,
                # dashes=False,
                label=f"{ic}",
                ax=ax[i, j],
            )
            ax[i, j].set_title(
                f"bc={boundary}, nx={grid_size}")
            ax[i, j].set_xlabel("Timepoint")
            # ax[i, j].set_ylim(0, np.max(subset["l2_error"]))
            ax[i, j].set_ylabel("L2(NMG- SAV)")
            plt.legend()
plt.suptitle(
    f"L2 Error between NMG and SAV", fontsize=16)
plt.tight_layout()

plt.show()

# %%
