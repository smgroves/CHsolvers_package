# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
data = pd.read_csv("compare_norms_v2.csv", sep=",", header=0, index_col=None)


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
