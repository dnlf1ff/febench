import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica','Arial']
# plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.size'] = 13
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.titlesize'] = 15
plt.rcParams['axes.titlelocation'] = 'right'
plt.rcParams['axes.titlepad'] = 5

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['xtick.labelbottom'] = True    
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.loc'] = 'lower right'
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'

plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.size'] = 2.5
plt.rcParams['xtick.minor.width'] = 0.3

plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.minor.size'] = 4
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.minor.width'] = 0.3
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['legend.facecolor'] = 'white'

plt.rcParams['figure.figsize'] = 3.4, 2.5
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['figure.subplot.bottom'] = 0.18
plt.rcParams['figure.subplot.left'] = 0.1
plt.rcParams['figure.subplot.right'] = 0.97
plt.rcParams['figure.subplot.top'] = 0.98
plt.rcParams['figure.subplot.wspace'] = 0.03
plt.rcParams['figure.subplot.hspace'] = 0.03


solutes= ['Co', 'Cr', 'Cu', 'Mn', 'Mo', 'Nb', 'Ni', 'Ti', 'V']
dft = df['DFT'].tolist()


solute_color_dict = {
    'Co': '#95bcff',
    'Cr': '#efa366',
    'Cu': '#03bc25',
    'Mn': '#ffa9ec',
    'Mo': '#b1ffba',
    'Nb': '#adadad',
    'Ni': '#e2e2e2',
    'Ti': '#00d4df',
    'V': '#fee900',
}

color_dict= {
    'ompa-omat': '#02baae',
    'ompa-mpa': '#72b8b3',
    'omni-omat': '#0073ff',
    'omni-mpa': '#8cbbf5',
    'omni-matpes_pbe': '#0047b3',
    'ORB-omat': '#d48f3f',
    'ORB-mpa': '#d9ba98',
    'eSEN-omat':'#e8d803',
    'eSEN-oam': '#f5ea74',
    'DPA-mp': '#d19adf',
    'DPA-omat': '#e829f6',
    'UMA-omat': '#ff0000',
    'UMA-omc': '#ff7f7f',
    'MACE-omat': '#00b300',
    'MACE-mpa': '#66ff66',
    'grace-omat': '#8b4513',
    'grace-oam': '#deb887',
    'nequip-oam': '#a52a2a',
}


def tm_3panel(mlips, df, filename):
    fig, ax = plt.subplots(1,3 , figsize=(9,3), sharey=True)
    ax = ax.flatten()
    lims = [-0.5, 0.35]
    plt.rcParams['xtick.major.size'] = 7
    ticks = [-0.25, 0, 0.25]
    plt.rcParams['ytick.major.size'] = 7

    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    for i, mlip in enumerate(mlips):
        mlip_values = df[mlip].tolist()
        ax[i].axhline(0, color='grey', linestyle='-', linewidth=1, zorder=0)
        ax[i].plot([-5,5], [-5,5], color='grey', linestyle='--', linewidth=1, zorder=0)
        ax[i].axvline(0, color='grey', linestyle='-', linewidth=1, zorder=0)
        # ax.set_title(mlip)
        ax[i].set_box_aspect(1)
        ax[i].set_xlim(lims)
        ax[i].set_ylim(lims)
        ax[i].set_xticks(ticks)
        ax[i].set_yticks(ticks)
        ax[i].text(0.05, 0.95, mlip, transform=ax[i].transAxes, fontsize=17, fontweight='bold',verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', edgecolor='none',alpha=0.8))
        for j, solute in enumerate(solutes):
            ax[i].scatter(dft[j],mlip_values[j], color='k', marker='x', alpha=1, zorder=2)        
            ax[i].scatter(dft[j], mlip_values[j], color=solute_color_dict[solute], s=220, edgecolors='k', linewidth=1, alpha=0.7)
        
        if i == 0:
            ax[i].set_ylabel("MLIP", fontsize=15, labelpad=3)
        else:
            ax[i].get_yaxis().set_visible(False)

        ax[i].set_xlabel("DFT", fontsize=15, labelpad=5)
    fig.tight_layout()
    plt.savefig(f'{filename}.png', dpi=300, bbox_inches='tight', transparent=True)
    plt.show()


def tm_legend():
    fig, ax = plt.subplots(figsize=(2,1.2))
    for j, solute in enumerate(solutes):
        ax.scatter(1,1, color = solute_color_dict[solute], s=120, edgecolors='k', linewidth=1, alpha=0.7, label=solute)
        ax.scatter(1,1, color='k', marker='x', alpha=1, zorder=2)
    ax.axis('off')
    ax.set_xlim(0,0.1)
    ax.set_ylim(0,0.1)
    ax.legend(loc='center', ncol=3, frameon=True, fontsize='medium', markerscale=1.2, borderpad = 0.4, labelspacing=0.5, handlelength=1.2, handletextpad=0.5, columnspacing=1.3,)
    fig.tight_layout()
    fig.savefig('tm_legend.png', dpi=300, bbox_inches='tight', transparent=True)

def plot_tm_grid(label):
    df = pd.read_csv(f'{label}/tm_aug_{label}.csv')
    dft = df['dft'].tolist()
    fig = plt.figure(figsize=(10,14))
    # 3 main rows (4 models each), plus 1 row for grace
    gs = gridspec.GridSpec(5, 4, figure=fig)
    lims = [-0.55, 0.55]

    # First 12 models in top 3 rows
    for i, m in enumerate(models[:-6]):
        print(m)
        row, col = divmod(i, 4)
        ax = fig.add_subplot(gs[row, col])
        mlip_values = df[m].tolist()
        ax.text(0.05, 0.95, m, transform=ax.transAxes, fontsize=15, fontweight='bold',verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', edgecolor='none',alpha=0.8))
        ax.set_box_aspect(1)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.axhline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        ax.axvline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        for j, solute in enumerate(solutes):
            ax.scatter(dft[j],mlip_values[j], color='k', marker='x', alpha=1, zorder=3)        
            ax.scatter(dft[j], mlip_values[j], color=solute_color_dict[solute], s=200, edgecolors='k', linewidth=1, alpha=0.7, zorder=2)
            
        ax.plot([-0.6,0.6],[-0.6,0.6],"k--",lw=1, zorder=1, alpha=0.8)
        ax.set_xticks([]);
        if i % 4 == 0:
            ax.set_ylabel("MLIP", fontsize=17, labelpad=3)
        else:
            ax.set_yticks([])

    # Last two models on their own row
    for k, m in enumerate(models[-6:-3]):
        ax = fig.add_subplot(gs[3, k])
        mlip_values = df[m].tolist()
        ax.text(0.05, 0.95, m, transform=ax.transAxes, fontsize=15, fontweight='bold',verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', edgecolor='none',alpha=0.8))
        ax.set_box_aspect(1)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.axhline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        ax.axvline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        for j, solute in enumerate(solutes):
            ax.scatter(dft[j],mlip_values[j], color='k', marker='x', alpha=1, zorder=3)        
            ax.scatter(dft[j], mlip_values[j], color=solute_color_dict[solute], s=200, edgecolors='k', linewidth=1, alpha=0.7, zorder=2)
        ax.plot([-0.6,0.6],[-0.6,0.6],"k--",lw=1, zorder=1, alpha=0.8)
        # ax.set_xlabel("DFT", fontsize=17, labelpad=3)
        ax.set_xticks([])
        if k == 0:
            ax.set_ylabel("MLIP", fontsize=17, labelpad=3)
        else:
            ax.set_yticks([])
        
    for l, m in enumerate(models[-3:]):
        ax = fig.add_subplot(gs[4, l])
        mlip_values = df[m].tolist()
        ax.text(0.05, 0.95, m, transform=ax.transAxes, fontsize=15, fontweight='bold',verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', edgecolor='none',alpha=0.8))
        ax.set_box_aspect(1)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.axhline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        ax.axvline(0, color='grey', linestyle='--', linewidth=1, zorder=0)
        for j, solute in enumerate(solutes):
            ax.scatter(dft[j],mlip_values[j], color='k', marker='x', alpha=1, zorder=3)        
            ax.scatter(dft[j], mlip_values[j], color=solute_color_dict[solute], s=200, edgecolors='k', linewidth=1, alpha=0.7, zorder=2)
        ax.plot([-0.6,0.6],[-0.6,0.6],"k--",lw=1, zorder=1, alpha=0.8)
        ax.set_xlabel("DFT", fontsize=17, labelpad=3)

        if l == 0:
            ax.set_ylabel("MLIP", fontsize=17, labelpad=3)
        else:
            ax.set_yticks([])

    plt.savefig(f'{label}.png', dpi=300, bbox_inches='tight', transparent=True)


    # plt.tight_layout()
    plt.show()

# for label in label_list:
#     plot_tm_grid(label)

def tm_elwise(mlips, df, filename=None):
    lims = [-0.5, 0.35]
    fig, axes = plt.subplots(4,5, figsize=(10,14), sharex=True, sharey=True)
    dft = df['DFT'].tolist()
    for i, (ax, mlip) in enumerate(zip(axes.flatten(), mlip_list)):
        mlip_values = df[mlip].tolist()
        ax.axhline(0, color='grey', linestyle='-', linewidth=1, zorder=0)
        ax.plot([-5,5], [-5,5], color='grey', linestyle='--', linewidth=1, zorder=0)
        ax.axvline(0, color='grey', linestyle='-', linewidth=1, zorder=0)
        # ax.set_title(mlip)
        ax.set_box_aspect(1)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.text(0.05, 0.95, mlip, transform=ax.transAxes, fontsize=15, fontweight='bold',verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', edgecolor='none',alpha=0.8))
        for j, solute in enumerate(solutes):
            ax.scatter(dft[j],mlip_values[j], color='k', marker='x', alpha=1, zorder=2)        
            ax.scatter(dft[j], mlip_values[j], color=solute_color_dict[solute], s=200, edgecolors='k', linewidth=1, alpha=0.7)
        
        if i % 4 == 0:
            ax.set_ylabel("MLIP", fontsize=17, labelpad=3)
        else:
            ax.get_yaxis().set_visible(False)
        if i < 12:
            ax.get_xaxis().set_visible(False)
        else:
            ax.set_xlabel("DFT", fontsize=17, labelpad=5)
                
    for ax in axes.flatten()[len(mlip_list):]:
        ax.axis("off")
    # fig.tight_layout()

    plt.savefig('tm_el_transparent.png', dpi=300, bbox_inches='tight', transparent=True)
    plt.show()

def tm_1panel(mlip, df, filename=None):
    if filename == None:
        filename = mlip
    fig, ax = plt.subplots(figsize=(3,3))
    mlips = ['UMA-omat']
    lims = [-0.5, 0.35]
    plt.rcParams['xtick.major.size'] = 7
    ticks = [-0.25, 0, 0.25]
    plt.rcParams['ytick.major.size'] = 7
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    # ax = ax.flatten()
    mlip_values = df[mlip].tolist()
    ax.axhline(0, color='grey', linestyle='-', linewidth=1, zorder=0)
    ax.plot([-5,5], [-5,5], color='grey', linestyle='--', linewidth=1, zorder=0)
    ax.axvline(0, color='grey', linestyle='-', linewidth=1, zorder=0)
    # ax.set_title(mlip)
    ax.set_box_aspect(1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.text(0.05, 0.95, mlip, transform=ax.transAxes, fontsize=17, fontweight='bold',verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', edgecolor='none',alpha=0.8))
    for j, solute in enumerate(solutes):
        ax.scatter(dft[j],mlip_values[j], color='k', marker='x', alpha=1, zorder=2)        
        ax.scatter(dft[j], mlip_values[j], color=solute_color_dict[solute], s=220, edgecolors='k', linewidth=1, alpha=0.7)
    ax.set_ylabel("MLIP", fontsize=15, labelpad=3)
    ax.set_xlabel("DFT", fontsize=15, labelpad=5)
    fig.tight_layout()
    plt.savefig(f'{filename}.png', dpi=300, bbox_inches='tight', transparent=True)
    plt.show()


def plot_mae():
    plt.figure(figsize=(15,5.5))
    vlines = [i+0.5 for i in range(len(carbon)-1)]
    mlip_list = carbon['mlip']
    sns.barplot(
        data=carbon, 
        x="mlip", 
        y="MAE", 
        width = 0.6,
        palette=color_dict,   # pick any color scheme
        errorbar=None        # removes those ugly error bars
    )
    plt.xticks(ticks=range(len(mlip_list)), labels=mlip_list, ha='center', rotation=30)
    plt.title(r"Carbon in Fe")
    plt.axhline(0, color="k", linestyle="--", linewidth=1,zorder=0)  # reference line
    for v in vlines:
        plt.axvline(v, color="grey", linestyle="--", linewidth=1,zorder=0)
    plt.ylabel(r"MAE (eV)")
    plt.xlabel("")

    plt.tick_params(axis="x", which="major", pad=5)
    filename = 'fig_mae_carbon.png'
    plt.savefig(filename, transparent=True)

    plt.show()

def sns_df_prep(df, mlip_list):
    sns_df = pd.DataFrame(columns=['index','DFT', 'mlip', 'pred','delta','abs_error'], index=range(len(mlip_list)*len(df)))
    idx = 0
    for i in range(len(df)):
        for mlip in mlip_list:
            sns_df.loc[idx, 'index'] = i
            sns_df.loc[idx, 'DFT'] = df.loc[i, 'DFT']
            sns_df.loc[idx, 'mlip'] =  mlip
            sns_df.loc[idx, 'pred'] = df.loc[i, mlip]
            sns_df.loc[idx, 'delta'] = df.loc[i, mlip+'_delta']
            sns_df.loc[idx, 'abs_error'] = abs(df.loc[i, mlip+'_delta'])
            idx += 1
    return sns_df