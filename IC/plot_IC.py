# Figure S1G
# %%
import matplotlib.font_manager as fm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib import font_manager
from matplotlib import rcParams

# for f in fm.findSystemFonts(fontpaths=None, fontext='ttf'):
#     if 'arial' in f.lower():
#         print(f)


# arial_path = "/System/Library/Fonts/Supplemental/Arial.ttf"  # e.g., from Step 1
# arial_font = font_manager.FontProperties(fname=arial_path)
# rcParams["font.family"] = arial_font.get_name()
plt.rcParams['pdf.use14corefonts'] = True

# %%
indir = "/Users/smgroves/Documents/GitHub/CHsolvers_package/IC"
phi_name = "initial_phi_128_smooth_n_relax_4_50p_neumann_v2.csv"

phi = np.genfromtxt(
    f"{indir}/{phi_name}",
    delimiter=",",
)

normalize_phis = mcolors.TwoSlopeNorm(vcenter=0, vmin=-1, vmax=1)
s = sns.heatmap(
    phi,
    square=True,
    cmap=cm.RdBu_r,
    norm=normalize_phis,
    cbar=False,
    linewidths=0.0,
)
plt.xticks(ticks=[], labels=[])
plt.yticks(ticks=[], labels=[])
# plt.title(f"Time= {timepoint*dt}")
plt.tight_layout()
plt.savefig(
    f"{indir}/{phi_name}.png",
    bbox_inches="tight",
    pad_inches=0,
    dpi=300,
)
plt.close()
# plt.show()

# %%

# %%

# %%
