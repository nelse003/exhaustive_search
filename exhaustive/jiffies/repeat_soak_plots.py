import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


if __name__ == "__main__":
    """plot from repeat soak csv"""

    occ_all_df = pd.read_csv( "/dls/science/groups/i04-1/elliot-dev/Work/"
                              "repeat_soak_plots/150919/occ_all.csv")

    compounds = occ_all_df['compound'].unique()
    ref_types = occ_all_df['refinement'].unique()


    # Colour map for plots
    cols = {'refmac':'g',
            'buster':'r',
            'phenix':'b',
            'refmac_superposed':'c',
            'exhaustive':'m',
            'phenix_superposed':'y',
            'buster_superposed':'k',
            }

    for compound in compounds:
        occ_comp_df = occ_all_df.loc[occ_all_df['compound'] == compound]
        resids = occ_comp_df['resid'].unique()

        for resid in resids:

            fig, ax = plt.subplots(1, 1)

            occ_resid_df = occ_comp_df.loc[occ_comp_df['resid'] == resid]

            # Scatter plots
            for ref_type in ref_types:
                occ_ref_df = occ_resid_df.loc[occ_resid_df['refinement'] == ref_type]

                occ_df = occ_ref_df.groupby(['xtal']).mean()

                if occ_df.empty:
                    continue

                occ_df.plot(kind='scatter',
                         x='Occupancy',
                         y='B_factor',
                         label='{} resid {} for n ={}'.format(ref_type, resid, len(occ_df)),
                         ax=ax,
                         color=cols[ref_type])
                plt.legend(fontsize=4, loc='best')

            plt.savefig("/dls/science/groups/i04-1/elliot-dev/Work/"
                        "repeat_soak_plots/scatter_{}_{}.png".format(compound,resid), dpi=300)
            plt.close()

            # Sns distplots
            # for ref_type in ref_types:
            #     occ_ref_df = occ_resid_df.loc[occ_resid_df['refinement'] == ref_type]
            #
            #     occ_df = occ_ref_df.groupby(['xtal']).mean()
            #
            #     if occ_df.empty:
            #         continue

            sns.violinplot(y=occ_resid_df['Occupancy'],
                           x=occ_resid_df['refinement'],
                           label="{} resid {} for n ={}".format(ref_type, resid, len(occ_df)))

            plt.legend(fontsize=4, loc='best')

            plt.savefig("/dls/science/groups/i04-1/elliot-dev/Work/"
                        "repeat_soak_plots/violin_{}_{}.png".format(compound,resid), dpi=300)
            plt.close()