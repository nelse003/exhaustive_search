import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import to_rgba

if __name__ == "__main__":
    """plot from repeat soak csv"""

    occ_all_df = pd.read_csv( "/dls/science/groups/i04-1/elliot-dev/Work/"
                              "repeat_soak_plots/150919/occ_all.csv")

    # Join FMOPL000662a
    occ_all_df = occ_all_df.replace({'compound':'FMOPL000622a_DSPL'}, 'FMOPL000622a')
    occ_all_df = occ_all_df.replace({'compound': ' FMOPL000622a_DSI_poised'}, 'FMOPL000622a')
    compounds = occ_all_df['compound'].unique()

    ref_types = occ_all_df['refinement'].unique()
    min_b = np.min(occ_all_df['B_factor'].unique())
    max_b = np.max(occ_all_df['B_factor'].unique())

    # Colour map for plots
    cols = {'buster': 'darkgreen',
            'buster_superposed': 'limegreen',
            'exhaustive': 'gold',
            'phenix': 'navy',
            'phenix_superposed': 'cornflowerblue',
            'refmac':'crimson',
            'refmac_superposed':'lightcoral',
            }

    # col_ar = [(0.0, 0.5, 0.0, 1.0),
    # (1.0, 0.0, 0.0, 1.0),
    # (0.0, 0.0, 1.0, 1.0),
    # (0.0, 0.75, 0.75, 1.0),
    # (0.75, 0.0, 0.75, 1.0),
    # (0.75, 0.75, 0, 1.0),
    # (0.0, 0.0, 0.0, 1.0)]
    #
    # cm = LinearSegmentedColormap.from_list("ARGH", col_ar, N=7)

    occ_df_list = []
    keys = []

    # Sort data
    for compound in compounds:
        occ_comp_df = occ_all_df.loc[occ_all_df['compound'] == compound]
        resids = occ_comp_df['resid'].unique()

        for resid in resids:

            fig, ax = plt.subplots(1, 1)

            occ_resid_df = occ_comp_df.loc[occ_comp_df['resid'] == resid]

            for ref_type in ref_types:

                keys.append((ref_type,compound,resid))

                occ_ref_df = occ_resid_df.loc[occ_resid_df['refinement'] == ref_type]

                if len(occ_ref_df['altloc'].unique()) == 1:
                    occ_ref_summary_df = occ_ref_df.groupby(['xtal']).mean()
                else:
                    occ_alt_df = occ_ref_df.groupby(['xtal', 'altloc']).mean()
                    occ_sum = occ_alt_df.groupby(['xtal'])['Occupancy'].sum()
                    occ_ref_summary_df = occ_ref_df.groupby(['xtal']).mean()
                    occ_ref_summary_df['Occupancy'] = occ_sum
                occ_df_list.append(occ_ref_summary_df)

    occ_df = pd.concat(occ_df_list,keys=keys)
    occ_df = occ_df.drop(['resid','Unnamed: 0'],axis='columns')
    occ_df.index.set_names(['program','compound','resid','xtal'], inplace=True)
    occ_df.to_csv("/dls/science/groups/i04-1/elliot-dev/Work/repeat_soak_plots/150919/occ_all_summary.csv")


    # Scatter plots
    for compound in compounds:
        for resid in ['1.0','2.0']:

            if occ_df.query('compound == @compound and resid == @resid').empty:
                continue

            fig, ax = plt.subplots(1, 1)

            refmac_occ_df = occ_df.query('compound == @compound '
                                                'and resid == @resid '
                                                'and program == "refmac"')

            refmac_occ_df =refmac_occ_df.reset_index(level='program',
                                                     drop=True)
            refmac_occ_df =refmac_occ_df.reset_index(level='compound',
                                                     drop=True)
            refmac_occ_df =refmac_occ_df.reset_index(level='resid',
                                                     drop=True)

            for ref_type in cols.keys():

                if ref_type == "refmac":
                    continue

                occ_cmpd_res_ref_df = occ_df.query('compound == @compound '
                                                   'and resid == @resid '
                                                   'and program == @ref_type')
                occ_cmpd_res_ref_df = occ_cmpd_res_ref_df.reset_index(level='program',
                                                          drop=True)
                occ_cmpd_res_ref_df = occ_cmpd_res_ref_df.reset_index(level='compound',
                                                          drop=True)
                occ_cmpd_res_ref_df = occ_cmpd_res_ref_df.reset_index(level='resid',
                                                          drop=True)

                join_df = refmac_occ_df.join(occ_cmpd_res_ref_df,
                                             on='xtal',
                                             how='outer',
                                             lsuffix='_refmac',
                                             rsuffix='_other')



                ax.spines["right"].set_visible(False)
                ax.spines["top"].set_visible(False)


                join_df.plot(kind='scatter',
                             x = 'Occupancy_refmac',
                             y = 'Occupancy_other',
                             label = '{}'.format(ref_type),
                             ax = ax,
                             color = cols[ref_type])

            plt.xlim(0,1.0)
            plt.ylim(0,1.0)
            plt.xlabel("Occupancy Refmac")
            plt.ylabel("Occupancy")
            plt.legend(frameon=False, fontsize=10)
            # plt.legend().set_visible(False)
            plt.savefig("/dls/science/groups/i04-1/elliot-dev/Work/"
            "repeat_soak_plots/occ_scatter_{}_{}.png".format(compound, resid), dpi = 300)
            plt.close()

    mean_occ_df = occ_df.groupby('program').mean()
    std_occ_df = occ_df.groupby('program').std()

    fig = plt.figure(figsize=(7,5.4))
    ax = plt.subplot(111)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    for mean, std in zip(mean_occ_df.iterrows(), std_occ_df.iterrows()):

        ax.errorbar(x=mean[1]['Occupancy'],
                    y=mean[1]['B_factor'],
                    xerr=std[1]['Occupancy'],
                    yerr=std[1]['B_factor'],
                    color=cols[mean[0]],
                    fmt='o',
                    ecolor=to_rgba(cols[mean[0]], 0.4),
                    markersize=8,
                    label=mean[0].replace('_',' ')
                    )

    # Dashed joining lines
    ax.plot([mean_occ_df.loc['buster', 'Occupancy'],
             mean_occ_df.loc['buster_superposed', 'Occupancy']],
            [mean_occ_df.loc['buster', 'B_factor'],
             mean_occ_df.loc['buster_superposed', 'B_factor']],
            linestyle='dashed',
            color=to_rgba('forestgreen', 0.5))

    ax.plot([mean_occ_df.loc['refmac', 'Occupancy'],
             mean_occ_df.loc['refmac_superposed', 'Occupancy']],
            [mean_occ_df.loc['refmac', 'B_factor'],
             mean_occ_df.loc['refmac_superposed', 'B_factor']],
            linestyle='dashed',
            color=to_rgba('crimson', 0.5))

    ax.plot([mean_occ_df.loc['phenix', 'Occupancy'],
             mean_occ_df.loc['phenix_superposed', 'Occupancy']],
            [mean_occ_df.loc['phenix', 'B_factor'],
             mean_occ_df.loc['phenix_superposed', 'B_factor']],
            linestyle='dashed',
            color=to_rgba('navy', 0.5),
            label="Relationship between\nsuperposed and\nnon-superposed\nrefinement methods\n ")

    plt.plot([], [], '-', label="Error bars show\n"
                                "standard deviation\n",
             color=to_rgba('forestgreen',0.4))

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left',
              bbox_to_anchor=(1, 0.5),
              fontsize=8,
              frameon=False)



    plt.xlabel("Mean occupancy across targets")
    plt.ylabel("Mean B factor across targets")
    plt.savefig("/dls/science/groups/i04-1/elliot-dev/"
                "Work/repeat_soak_plots/mean_occ_targets.png",dpi=300)
    plt.close()

    # distribution width plot
    std_df = occ_df.groupby(['program','compound','resid']).std()
    mean_df = occ_df.groupby(['program','compound','resid']).mean()

    var_df = std_df.groupby('program').std()
    width_df = std_df.groupby('program').mean()


    print(width_df.columns.values)
    # Plot each point separately to allow color
    fig = plt.figure(figsize=(7,5.4))
    ax = plt.subplot(111)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    for width, var in zip(width_df.iterrows(), var_df.iterrows()):

        ax.errorbar(x=width[1]['Occupancy'],
                     y=width[1]['B_factor'],
                     xerr=var[1]['Occupancy'],
                     yerr=var[1]['B_factor'],
                     color=cols[var[0]],
                     label=var[0].replace('_', ' '),
                     ecolor=to_rgba(cols[var[0]], 0.4),
                     markersize=8,
                     fmt='o')

    # Dashed joining lines
    ax.plot([width_df.loc['buster','Occupancy'],
                width_df.loc['buster_superposed','Occupancy']],
             [width_df.loc['buster','B_factor'],
                width_df.loc['buster_superposed','B_factor']],
             linestyle='dashed',
             color=to_rgba('forestgreen',0.5))

    ax.plot([width_df.loc['refmac','Occupancy'],
                width_df.loc['refmac_superposed','Occupancy']],
             [width_df.loc['refmac','B_factor'],
                width_df.loc['refmac_superposed','B_factor']],
             linestyle='dashed',
             color=to_rgba('crimson',0.5))

    ax.plot([width_df.loc['phenix','Occupancy'],
                width_df.loc['phenix_superposed','Occupancy']],
             [width_df.loc['phenix','B_factor'],
                width_df.loc['phenix_superposed','B_factor']],
             linestyle='dashed',
             color=to_rgba('navy',0.5),
            label="Relationship between\nsuperposed and\nnon-superposed\nrefinement methods\n ")

    plt.plot([], [], '-', label="Error bars show\n"
                                "varaiability\n"
                                "of distribution width\n"
                                "(standard deviation)\n"
                                "across 8\n"
                                "protein-ligand\n"
                                "repeated soaks\n",
             color=to_rgba('forestgreen',0.4))

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left',
              bbox_to_anchor=(1.0, 0.5),
              fontsize=8,
              frameon=False)
    plt.xlabel("Occupancy distribution width")
    plt.ylabel("B factor distribution width")

    plt.savefig("/dls/science/groups/i04-1/elliot-dev/"
                "Work/repeat_soak_plots/distribution_widths.png",dpi=300)

    # Scatter plots
    for compound in compounds:
        for resid in ['1.0','2.0']:

            if occ_df.query('compound == @compound and resid == @resid').empty:
                continue

            fig, ax = plt.subplots(1, 1)
            plt.xlim(0, 1.05)
            plt.ylim(0, 150)

            for ref_type in cols.keys():

                occ_cmpd_res_ref_df = occ_df.query('compound == @compound '
                                                   'and resid == @resid '
                                                   'and program == @ref_type')
                if occ_cmpd_res_ref_df.empty:
                    continue

                ax.spines["right"].set_visible(False)
                ax.spines["top"].set_visible(False)

                occ_cmpd_res_ref_df.plot(kind='scatter',
                                        x='Occupancy',
                                        y='B_factor',
                                        label='{}'.format(ref_type.replace('_',' ')),
                                        ax=ax,
                                        color=cols[ref_type])
            plt.ylabel("B factor")
            plt.legend(frameon=False, fontsize=12, loc='upper left')
            #plt.legend().set_visible(False)
            plt.savefig("/dls/science/groups/i04-1/elliot-dev/Work/"
                        "repeat_soak_plots/large_leg_scatter_{}_{}.png".format(compound,resid.strip('.0')), dpi=300)
            plt.close()

    # distplots
    for compound in compounds:
        for resid in ['1.0','2.0']:

            if occ_df.query('compound == @compound and resid == @resid').empty:
                continue

            fig, ax = plt.subplots(1, 1)

            for ref_type in cols.keys():
                occ_cmpd_res_ref_df = occ_df.query('compound == @compound '
                                                   'and resid == @resid '
                                                   'and program == @ref_type')

                plt.xlim(0, 1.05)
                plt.ylabel('Frequency density')
                ax.spines["right"].set_visible(False)
                ax.spines["top"].set_visible(False)

                sns.distplot(occ_cmpd_res_ref_df['Occupancy'],
                             hist=False,
                             label="{}".format(ref_type.replace('_',' ')),
                             color=cols[ref_type])

                print(f"{compound} {len(occ_cmpd_res_ref_df)}")

            plt.legend(fontsize=10, loc='best', frameon=False)
            plt.legend().set_visible(False)
            plt.savefig("/dls/science/groups/i04-1/elliot-dev/Work/"
                        "repeat_soak_plots/no_leg_distplot_{}_{}.png".format(compound, resid.strip('.0')), dpi=300)
            plt.close()

