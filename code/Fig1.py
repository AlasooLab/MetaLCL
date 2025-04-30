import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text


def get_annotations_df(df):
    df_unique_eQTL_location_in_var_count_df = df[df['phenotype_id_count'] > 4] # test min
    df_unique_eQTL_location_in_var_count_df = df_unique_eQTL_location_in_var_count_df.sort_values(
        by=["datasets_with_variant", 'phenotype_id_count'], ascending=[False, False])
    df_unique_eQTL_location_in_var_count_df = df_unique_eQTL_location_in_var_count_df[  # # annotatated mannually based on component_group_id
        (~df_unique_eQTL_location_in_var_count_df.duplicated(subset=['component_group_id'], keep='first')) | (
            df_unique_eQTL_location_in_var_count_df['component_group_id'].isnull())]
    df_unique_eQTL_location_in_var_count_df = df_unique_eQTL_location_in_var_count_df.drop_duplicates(
        subset=["eQTL_location_gene_name"], keep="first")
    return df_unique_eQTL_location_in_var_count_df

def make_plot_final(df):
    annotation_labels_df = get_annotations_df(df)
    hue_order = df['colouring'].unique()
    sns.set(style="ticks")
    sns.set_style("whitegrid")
    palette = ["#808080","#e41a1c",
"#377eb8",
"#4daf4a",
"#984ea3",
"#ffff33","#e6ab02"]




    f, (ax_top, ax_top_bottom, ax_bottom) = plt.subplots(ncols=1, nrows=3, sharex=True, gridspec_kw=dict(height_ratios=[0.2,0.5,3.5]), figsize=(18,10)) # 3,0.3 16 10

    sns.scatterplot(data=df[(df["phenotype_id_count"] > 30)], ax=ax_top, y="phenotype_id_count", x="variant_position_cumulative_pos",legend=None)
    sns.scatterplot(data=df[(df["phenotype_id_count"] < 30)], ax=ax_top_bottom, y="phenotype_id_count", x="variant_position_cumulative_pos",legend=None) # , hue_order=hue_order,palette = palette,

    scatterplot = sns.scatterplot(data=df, x="variant_position_cumulative_pos", y="phenotype_pos_cumulative_pos", legend="brief", ax=ax_bottom) #hue="colouring",

    """markers = ['o', 'x']
    legend_labels = ["No cross-mappability", "Cross-mappability"]
    handles = [plt.Line2D([0], [0], marker=markers[i], color='black', markerfacecolor='black', markersize=10) for i in
               range(len(markers))]
    scatterplot_legend = ax_bottom.legend(handles, legend_labels,
                                          loc="lower right", ncol=2)"""


    ax_top.set_ylim(50, 300)
    ax_top_bottom.set_ylim(0, 23)

    ax_top_bottom.set_ylabel('# of trans genes',fontsize=18)
    ax_top.set_ylabel("")
    ax_bottom.set_xlabel('Position of the genetic variant',fontsize=18)
    ax_bottom.set_xticks(df.groupby('variant_chromosome')['variant_position_cumulative_pos'].median())
    ax_bottom.set_xticklabels(sorted(df['variant_chromosome'].unique()), fontsize=10)
    ax_bottom.set_ylabel('Position of the gene',fontsize=18)
    ax_bottom.set_yticks(df.groupby('chromosome')['phenotype_pos_cumulative_pos'].median())
    ax_bottom.set_yticklabels(sorted(df['chromosome'].unique()), fontsize=10)

    ax = ax_top
    d = .003
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)

    ax2 = ax_top_bottom
    kwargs.update(transform=ax2.transAxes)
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    annotations_eQTL_counts_per_eQTL_location_gene = annotation_labels_df.apply(
        lambda p: ax_top_bottom.annotate((p['eQTL_location_gene_name']),
                                    (p[
                                         'variant_position_cumulative_pos'],
                                     p['phenotype_id_count']),fontsize =18),

        axis=1).to_list()
    annotations_eQTL_counts_per_eQTL_location_gene_top = annotation_labels_df.apply(
        lambda p: ax_top.annotate((p['eQTL_location_gene_name']),
                                         (p[
                                              'variant_position_cumulative_pos'],
                                          p['phenotype_id_count']),fontsize=18),

        axis=1).to_list()
    print(annotations_eQTL_counts_per_eQTL_location_gene)
    adjust_text(annotations_eQTL_counts_per_eQTL_location_gene, ax=ax_top_bottom,
                arrowprops={'arrowstyle': '->', 'color': 'red'})
    adjust_text(annotations_eQTL_counts_per_eQTL_location_gene_top, ax=ax_top,
                arrowprops={'arrowstyle': '->', 'color': 'red'})

    plt.subplots_adjust(hspace=0.03)
    plt.savefig(f'FigA.jpeg')


file_to_make_plot = "FigA_data.tsv"
df = pd.read_csv(file_to_make_plot, sep='\t')
make_plot_final(df)

