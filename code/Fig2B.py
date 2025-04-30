import pandas as pd
import matplotlib.pylab as plt
import numpy as np
from adjustText import adjust_text
from statsmodels.stats.multitest import fdrcorrection


fig2b_file = "Fig2B_chr22_18166589_T_C.tsv"
fig2_df = pd.read_csv(fig2b_file, sep='\t')

fig2_df["p_value"] = np.array(10 ** -(fig2_df["log10p"]))
res0 = fdrcorrection(fig2_df["p_value"].to_numpy(), alpha=0.05)
fig2_df["fdr"] = (res0[1])
fdr_genes_path = "chr22_18166589_T_C_meta_fdr_genes.txt"
fdr_genes_df = pd.read_csv(fdr_genes_path, sep='\t')
fdr_genes = fdr_genes_df.phenotype_id.to_list()

# highlight down- or up- regulated genes
down = fig2_df[(fig2_df['pooled_effects']<0)&(fig2_df['fdr']<0.05)]
up = fig2_df[(fig2_df['pooled_effects']>0)&(fig2_df['fdr']<0.05)]

f, (ax, ax2) = plt.subplots(2, 1, sharex=True,gridspec_kw=dict(height_ratios=[0.8,3]))
ax.scatter(x=fig2_df['pooled_effects'], y=fig2_df['fdr'].apply(lambda x:-np.log10(x)), c='grey', s=10)
ax2.scatter(x=fig2_df['pooled_effects'], y=fig2_df['fdr'].apply(lambda x:-np.log10(x)), c='grey', s=10)

ax.scatter(x=down['pooled_effects'],y=down['fdr'].apply(lambda x:-np.log10(x)),s=30,color="red")
ax.scatter(x=up['pooled_effects'],y=up['fdr'].apply(lambda x:-np.log10(x)),s=30,color="red")

ax2.scatter(x=down['pooled_effects'],y=down['fdr'].apply(lambda x:-np.log10(x)),s=30,color="red")
ax2.scatter(x=up['pooled_effects'],y=up['fdr'].apply(lambda x:-np.log10(x)),s=30,color="red")


# zoom in and limit the view to different portions of the data
ax.set_ylim(22,26)  # outliers only
ax2.set_ylim(0, 10)  # most of the data

ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)
ax2.xaxis.tick_bottom()

ax.tick_params(
    axis='x',
    which='both',
    bottom=False,
    top=False,
    labelbottom=False)
d = .015
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)

kwargs.update(transform=ax2.transAxes)
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

texts=[]
for index,r in up.iterrows():
    texts.append(plt.text(x=r['pooled_effects'],y=-np.log10(r['fdr']),s=r['gene_name']))

for index,r in down.iterrows():
    if ( r['phenotype_id'] in fdr_genes and r['pooled_effects'] <= -0.1) or r['gene_name']== "USP41":
        texts.append(plt.text(x=r['pooled_effects'],y=-np.log10(r['fdr']),s=r['gene_name']))

adjust_text(texts,
            arrowprops=dict(arrowstyle="-", color='black', lw=0.4))
fdr_threshold = 0.05
log_fdr_threshold = -np.log10(fdr_threshold)

ax.axhline(y=log_fdr_threshold, color='black', linestyle='--', linewidth=1, label="FDR = 0.05")
ax2.axhline(y=log_fdr_threshold, color='black', linestyle='--', linewidth=1, label="FDR = 0.05")

plt.xlabel("Effect of rs4819670-C",fontsize=16)
plt.ylabel("-logFDR",fontsize=16)


plt.savefig("chr22_18166589_T_C_USP18_Volcano_plot.pdf", dpi=300)
