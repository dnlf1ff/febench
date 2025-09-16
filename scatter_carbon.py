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

def plot_violin(sns_df, filename):
    plt.figure(figsize=(15,6.5))
    vlines = [i+0.5 for i in range(len(mlip_list)-1)]
    sns.violinplot(
        data=sns_df,
        x="mlip",
        y="delta",          
        inner="quartile",
        cut=2,
        palette=color_dict
    )
    plt.xticks(ticks=range(len(mlip_list)), labels=mlip_list, ha='center', rotation=0)
    ylabel = r"$\Delta E_{\mathrm{bind}}$ (eV)"
    plt.axhline(0, color="k", linestyle="--", linewidth=1)  # reference line
    for v in vlines:
        plt.axvline(v, color="grey", linestyle="--", linewidth=1,zorder=0)
    plt.xlabel("")
    
    # sns.despine()
    
    plt.tick_params(axis="x", which="major", pad=5)
    plt.tight_layout()
    plt.savefig(filename, transparent=True)
    plt.show()
    print(f'Violin plot saved as {filename}')


def plot_box(sns_df,filename):
    plt.figure(figsize=(15,6.5))
    vlines = [i+0.5 for i in range(len(mlip_list)-1)]
    sns.boxplot(
        data=sns_df,
        x="mlip",
        y="delta",
        palette=color_dict,
        whis=(5, 95),
        showfliers=True,
        width=0.6
    )
    plt.xticks(ticks=range(len(mlip_list)), labels=mlip_list, ha='center', rotation=0)
    ylabel = r"$\Delta E_{\mathrm{bind}}$ (eV)"
    plt.axhline(0, color="k", linestyle="--", linewidth=1,zorder=0)  # reference line
    for v in vlines:
        plt.axvline(v, color="grey", linestyle="--", linewidth=1,zorder=0)
    plt.xlabel("")
        
    plt.tick_params(axis="x", which="major", pad=5)
    plt.savefig(filename, transparent=True)
    # plt.close()
    plt.show()
    print(f'Box plot saved as {filename}')

import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica','Arial']
# plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.size'] = 13
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.titlelocation'] = 'right'
plt.rcParams['axes.titlepad'] = 3

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['xtick.labelbottom'] = True    
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.loc'] = 'lower right'
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'in'

plt.rcParams['xtick.major.size'] = 7
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.size'] = 2.5
plt.rcParams['xtick.minor.width'] = 0.3

plt.rcParams['ytick.major.size'] = 7
plt.rcParams['ytick.minor.size'] = 4
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
plt.rcParams['figure.subplot.wspace'] = 0.05
plt.rcParams['figure.subplot.hspace'] = 0.03


color_dict= {
    'ompa-omat24': '#02baae',
    'ompa-mpa': '#72b8b3',
    '7net-omat': '#c1ded8',
    'omni-omat24': '#0073ff',
    'omni-mpa': '#8cbbf5',
    '7net-l3i5': '#0a8a2c',
    '7net-0': '#89c99a',
    'orb_omat': '#d48f3f',
    'orb_mpa': '#d9ba98',
    'esen_omat':'#e8d803',
    'esen_oam': '#f5ea74',
    'dpa-mp': '#d19adf',
    'dpa-omat': '#e829f6',
}
