# %%
from matplotlib.lines import Line2D
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re
from matplotlib import font_manager
from matplotlib import rcParams

# for f in fm.findSystemFonts(fontpaths=None, fontext='ttf'):
#     if 'arial' in f.lower():
#         print(f)


# arial_path = "/System/Library/Fonts/Supplemental/Arial.ttf"  # e.g., from Step 1
# arial_font = font_manager.FontProperties(fname=arial_path)
# rcParams["font.family"] = arial_font.get_name()
plt.rcParams['pdf.use14corefonts'] = True


def extract_info(s):
    match = re.search(r'_Nx_(\d+)_([a-z]+)_dtout_(\d+p)', s)
    if match:
        return int(match.group(1)), match.group(2), match.group(3)
    return None, None, None


def convert_to_long(df, dt, dt_out, timepoints):
    df[['Nx', 'boundary', 'IC']] = df['compared_files'].apply(
        # use first file only
        lambda x: pd.Series(extract_info(x.split(',')[0]))
    )

    value_vars = [col for col in df.columns if col in [
        str(i) for i in range(timepoints)]]
    df_long = df.melt(id_vars=['Nx', 'boundary', 'IC'],
                      value_vars=value_vars,
                      var_name='timepoint',
                      value_name='value')
    df_long['timepoint'] = (df_long['timepoint'].astype(float)-1) * dt_out*dt

    return df_long


# %% COMPARING 5.5e-6 to 5.5e-7

dt = 5.5e-6
dt_out = 10
timepoints = 201

df = pd.read_csv(f'./output/rmse/rmse_NMG_SAV_Julia_dt_{dt}.csv')
df_large_dt = convert_to_long(df, dt, dt_out, timepoints)

dt = 5.5e-7
dt_out = 10
timepoints = 2001

df = pd.read_csv(f'./output/rmse/rmse_NMG_SAV_Julia_dt_{dt}_v2.csv')
df2 = pd.read_csv(f'./output/rmse/rmse_NMG_SAV_Julia_dt_{dt}_v3.csv')

df = pd.concat([df, df2], ignore_index=True)

df_small_dt = convert_to_long(df, dt, dt_out, timepoints)
df_small_dt = df_small_dt[df_small_dt['timepoint'].isin(
    df_large_dt['timepoint'].unique())]  # filter Nx = 1


# %%

for bc in ['periodic', 'neumann']:
    fig, ax = plt.subplots(figsize=(4, 3))
    # for nx in df_small_dt['Nx'].unique():

    # for ic in df_small_dt['IC'].unique():
    lines = ['--', '-']
    dt_labels = ['5.5e-7', '5.5e-6']

    for i, (d, line) in enumerate(zip([df_small_dt, df_large_dt], lines)):
        # df_tmp = d[(d['boundary'] == bc) & (d['IC'] == ic)]
        df_tmp = d[(d['boundary'] == bc) & (d['Nx'] == 512)]

        sns.lineplot(
            data=df_tmp,
            x='timepoint',
            y='value',
            hue='IC',
            # palette=sns.color_palette("Greys")[1:5],
            hue_order=['1075p', '1050p', '1025p'],
            linestyle=line,
            ax=ax,
        )

    ax.set_title(f'RMSE(NMG-SAV) Over Time for BC: {bc}')
    ax.set_xlabel('Time')
    ax.set_ylabel('RMSE')
    ax.set_ylim(0, .5)
    ax.set_xlim(0, 0.011)
    # Grab Seaborn-generated legend for hue='Nx'
    handles1, labels1 = ax.get_legend_handles_labels()

    # Remove automatic legend
    ax.legend_.remove()

    # First legend (for Nx)
    legend1 = ax.legend(handles1, labels1, title="IC", loc="upper right")
    ax.add_artist(legend1)

    # # Second legend (for dt line styles)
    # dt_handles = [
    #     Line2D([0], [0], linestyle=style,
    #            color='black', label=f'dt = {lbl}')
    #     for style, lbl in zip(lines, dt_labels)
    # ]
    # legend2 = ax.legend(dt_handles, [f'dt = {lbl}' for lbl in dt_labels],
    #                     title="dt", loc="upper left", bbox_to_anchor=(1.02, 1))
    # ax.add_artist(legend2)

    plt.tight_layout()
    # plt.show()
    plt.savefig(
        f"./output/rmse/rmse_NMG_SAV_Julia_dt_comparison_color_IC_{bc}_512.pdf",
    )


# %% code consistency
dt = 5.5e-6
dt_out = 1
timepoints = 21
# Comparing MATLAB and Python implementations
df = pd.read_csv(f'./output/rmse/rmse_NMG_SAV_Python.csv')
df_long = convert_to_long(df, dt, dt_out, timepoints)

# %%

for bc in df_small_dt['boundary'].unique():
    plt.figure(figsize=(10, 6))
    df_tmp = df_long[(df_long['boundary'] == bc)]
    sns.lineplot(data=df_tmp, x='timepoint', y='value', hue='IC',
                 style='Nx', markers=False)
    plt.title(f'Python RMSE(NMG-SAV) Over Time for BC: {bc}')
    plt.xlabel('Time')
    plt.ylabel('RMSE')
    plt.legend(title='Nx')
    plt.tight_layout()
    plt.show()

# %%

# %%
