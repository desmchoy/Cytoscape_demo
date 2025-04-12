#!/usr/bin/env python3

import os
import sys
import re
import time
import argparse
import glob
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
import biomart
import py4cytoscape as p4c
pd.options.mode.chained_assignment = None
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

workdir = os.getcwd()
project = os.path.basename(workdir)

help_text = 'Cytoscape_py4cytoscape.py: Script to run pathway analysis with Cytoscape via py4cytoscape'

parser = argparse.ArgumentParser(description=help_text)
parser.add_argument(
    dest='diffname',
    action='store',
    type=str,
    help='Diff Directory (required)',
    metavar='diffname'
    )
parser.add_argument(
    dest='filedir',
    action='store',
    type=str,
    help='File Directory (required)',
    metavar='filedir'
    )


arguments = parser.parse_args()
diffname = arguments.diffname
filedir = arguments.filedir


outdir = workdir + '/results_Cytoscape/' + diffname
if os.path.exists(outdir):
#    os.makedirs(outdir)
    os.makedirs(outdir + '/tables')
    os.makedirs(outdir + '/images')
    os.makedirs(outdir + '/clusters')
    os.makedirs(outdir + '/images/enrichment')
    os.makedirs(outdir + '/images/networks')


help_text = 'Cytoscape_py4cytoscape.py: Script to run pathway analysis with Cytoscape via py4cytoscape'



### DATA INPUT
files = glob.glob(filedir + '/' + '*.MATS.JCEC.txt')

rMATS_files = [pd.read_csv(f, delimiter = '\t') for f in files]
rMATS_files = [i[['GeneID', 'geneSymbol', 'FDR', 'IncLevelDifference']] for i in rMATS_files]
splice_type = [ os.path.basename(f).replace('.MATS.JCEC.txt', '') for f in files]
for i in range(len(rMATS_files)):
    rMATS_files[i]['Splice'] = splice_type[i]

full_data = pd.concat(rMATS_files, ignore_index=True)



### DATA PRE-PROCESSING
# Include only genes with FDR < 0.05
data = full_data[full_data['FDR'] < 0.05]

## Cytoscape query requires Entrez ID - need biomart conversion
# biomart sorts results by Gene ID - need to sort input to preserve order, followed by deduplication
data = data.sort_values(by = ['GeneID'])
to_query = data['GeneID'].unique().tolist()


# biomart cannot take large queries - need to break the queries into smaller chunks
chunk_size = 100
chunks = math.floor(len(to_query) / chunk_size)
remainder = len(to_query) - chunks
to_iter = []
for i in range(chunks + 1):
    to_iter.append(to_query[i * chunk_size : i * chunk_size + chunk_size])


# Running biomart
server = biomart.BiomartServer('http://www.ensembl.org/biomart')
mart = server.datasets['hsapiens_gene_ensembl']
attributes = ['entrezgene_id', 'ensembl_gene_id', 'external_gene_name']

entrezgene_id = []
ensemble_gene_id = []
external_gene_name = []

for i in to_iter:
    response = mart.search({'attributes': attributes, 'filters': {'ensembl_gene_id': i}})
    decode = response.raw.data.decode('ascii')
    for line in decode.splitlines():
        line = line.split('\t')
        entrezgene_id.append(str(line[0]))
        ensemble_gene_id.append(str(line[1]))
        external_gene_name.append(str(line[2]))


# Need to remove entries with multiple Entrez IDs, or genes without Entrez IDs - also need to tally up the number of times a gene appear in a gene list
conversion = pd.DataFrame({'Entrez': entrezgene_id, 'ID': ensemble_gene_id, 'Name': external_gene_name})
duplicated = conversion[conversion['Name'].isin(conversion[conversion['ID'].duplicated()]['Name'].tolist())]
conversion_dedup = conversion[~conversion['ID'].duplicated()]
conversion_dedup = conversion_dedup[~(conversion_dedup['Entrez'] == '')].reset_index(drop = True)

gene_counts = pd.DataFrame(data.groupby('GeneID').size()).reset_index()
gene_counts.columns = ['ID', 'Count']

conversion_dedup = pd.merge(conversion_dedup, gene_counts, how = 'left', on = 'ID')

#if len(conversion_dedup['Entrez'].tolist()) > 1000:
#    print('Gene list is too long. Not running analysis.')
#    quit()



### CYTOSCAPE
# Running Stringapp to map genes onto network
p4c.cytoscape_version_info()
p4c.session.close_session(False)
#query_str = ','.join(str(v) for v in conversion_dedup['Entrez'])
#string_cmd_list = ['string protein query','query="',query_str,'"', 'species="Homo sapiens"', 'cutoff=0.4']
#string_cmd = " ".join(string_cmd_list)
#p4c.commands.commands_run(string_cmd)

#string_query = conversion_dedup['Entrez'].tolist()
string_cmd = 'string disease query disease="Myelodysplastic syndrome" limit=1000'
p4c.commands.commands_run(string_cmd)



#string_query = conversion_dedup['Entrez'].tolist()[0:1010]
#chunk_size = 100
#chunks = math.floor(len(string_query) / chunk_size)
#remainder = len(string_query) - chunks
#to_add = []
#for i in range(chunks + 1):
#    to_add.append(string_query[i * chunk_size : i * chunk_size + chunk_size])
#
#
#query_str = ','.join(str(v) for v in to_add[0])
#string_cmd_list = ['string protein query','query="',query_str,'"', 'species="Homo sapiens"', 'cutoff=0.4']
#string_cmd = " ".join(string_cmd_list)
#p4c.commands.commands_run(string_cmd)
#
#if len(to_add) > 1:
#    for i in range(1,len(to_add)):
#        query_str = ','.join(str(v) for v in to_add[i])
#        string_cmd_list = ['string add nodes','query="',query_str,'"', 'species="Homo sapiens"', 'cutoff=0.4']
#        string_cmd = " ".join(string_cmd_list)
#        p4c.commands.commands_run(string_cmd)
    

#p4c.map_table_column('display name', 'Human', 'HGNC', 'Entrez Gene')
#p4c.load_table_data(conversion_dedup, data_key_column = 'Entrez', table_key_column = 'Entrez Gene')
p4c.load_table_data(conversion_dedup, data_key_column = 'Name', table_key_column = 'display name')


# Enrichment analysis using Stringapp
p4c.commands.commands_run('string retrieve enrichment allNetSpecies="Homo sapiens"')
p4c.commands.commands_run('string show enrichment')

# Running MCODE Clustering
mcode = p4c.commands.commands_post('mcode cluster degreeCutoff=2 fluff=false fluffNodeDensityCutoff=0.1 haircut=true includeLoops=false kCore=2 maxDepthFromStart=100 network=current nodeScoreCutoff=0.2 scope=NETWORK')



# Extract results tables
#suid = p4c.networks.get_network_suid('STRING network')
suid = p4c.networks.get_network_suid()
nodes_table_t = pd.DataFrame(p4c.commands.cyrest_get(f'networks/{suid}/tables/defaultnode').get('rows'))


# Cluster column is stored as a list - need to convert to text
dummy_df = nodes_table_t[['display name', 'MCODE::Clusters (1)']]
dummy_df['Cluster'] = dummy_df['MCODE::Clusters (1)'].explode().fillna('Unclustered')
dummy_df = dummy_df.drop('MCODE::Clusters (1)', axis = 1)
p4c.load_table_data(dummy_df, data_key_column = 'display name', table_key_column = 'display name')


all_tables = p4c.commands.cyrest_get(f'networks/{suid}/tables')
titles = []
for i in all_tables:
    titles.append(i.get('title'))

nodes_table = pd.DataFrame(p4c.commands.cyrest_get(f'networks/{suid}/tables/defaultnode').get('rows'))
enrichment = pd.DataFrame(all_tables[titles.index('STRING Enrichment: All')].get('rows'))

nodes_table.to_csv(outdir + '/tables/nodes_table.csv')
enrichment.to_csv(outdir + '/tables/enrichment.csv')



# Plot enrichment results
def plot_enrichment(e_file = enrichment, sort_by = 'FDR value', cat = None, xlabel = '-log(FDR)', ylabel = '', xdim = 12, ydim = 5, savedir = os.getcwd()):
    to_plot = e_file.sort_values(sort_by)
    to_plot = to_plot[to_plot['category'] == cat].reset_index(drop = True).iloc[0:20,:]
    fig, ax = plt.subplots(figsize = (xdim, ydim))
    ax.barh(to_plot['term name'] + ': ' + to_plot['description'] + ' (' + to_plot['# genes'].apply(str) + ')', -np.log10(to_plot['FDR value']))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(cat)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.gca().invert_yaxis()
    fig.tight_layout()
    plt.savefig(savedir + '/' + cat.replace(' ', '_').replace(',', '') + '.png')
    plt.close()

for i in enrichment['category'].unique():
    plot_enrichment(cat = i, savedir = outdir + '/images/enrichment')



# Ridge plots for Stringapp scores
tissues = nodes_table.columns[nodes_table.columns.str.startswith('tissue::')]
compartments = nodes_table.columns[nodes_table.columns.str.startswith('compartment::')]

to_plot_tissues = pd.melt(nodes_table[tissues], var_name = 'Tissue', value_vars = tissues, value_name = 'Score').fillna(0)
to_plot_tissues['Tissue'] = to_plot_tissues['Tissue'].str.replace('tissue::', '')

to_plot_compartments = pd.melt(nodes_table[compartments], var_name = 'Compartment', value_vars = compartments, value_name = 'Score').fillna(0)
to_plot_compartments['Compartment'] = to_plot_compartments['Compartment'].str.replace('compartment::', '')


def plot_ridge(dataframe, filename = None, savedir = os.getcwd()):
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
    g = sns.FacetGrid(dataframe, row=dataframe.columns[0], hue=dataframe.columns[0], aspect=15, height=.5, palette=pal)
    g.map(sns.kdeplot, "Score", bw_adjust=.5, clip_on=False, fill=True, alpha=1, linewidth=1.5)
    g.map(sns.kdeplot, "Score", clip_on=False, color="w", lw=2, bw_adjust=.5)
    g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color, ha="left", va="center", transform=ax.transAxes)
    g.map(label, "Score")
    g.figure.subplots_adjust(hspace=-.25)
    g.set_titles("")
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)
    plt.savefig(savedir + '/' + filename + '.png')
    plt.close()


plot_ridge(to_plot_tissues, filename = 'tissues', savedir = outdir + '/images')
plot_ridge(to_plot_compartments, filename = 'compartment', savedir = outdir + '/images')



# Cluster tally
cluster_tally = pd.DataFrame({'Number': pd.DataFrame(nodes_table['Cluster'].apply(str)).groupby('Cluster').size().sort_values(ascending = False)}).reset_index()
cluster_tally = cluster_tally[cluster_tally['Cluster'] != 'nan']
cluster_tally['Percentage'] = 100 * cluster_tally['Number'] / cluster_tally['Number'].sum()

fig, ax = plt.subplots(figsize = (10, 8))
ax.bar(cluster_tally['Cluster'], cluster_tally['Percentage'])
ax.set_xlabel('Cluster')
ax.set_ylabel('Percentage')
ax.set_title('Cluster Size')
ax.tick_params(axis='x', labelrotation=270)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.tight_layout()
plt.savefig(outdir + '/images/clusters.png')



stats = conversion_dedup.merge(nodes_table[['display name', 'Cluster']], left_on = 'Name', right_on = 'display name', how = 'inner')
stats = stats.sort_values(by = 'Count', ascending = False).reset_index(drop = True)
stats50 = stats.iloc[0:50,:]

cluster_count_50 = pd.DataFrame({'Number': stats50.groupby('Cluster').size().sort_values(ascending = False)}).reset_index()
fig, ax = plt.subplots(figsize = (5, 5))
ax.bar(cluster_count_50['Cluster'], cluster_count_50['Number'])
ax.set_xlabel('Cluster')
ax.set_ylabel('Number')
ax.set_title('Clusters of Top 50 Genes (out of ' + str(len(stats)) + ') by Splicing Instances')
ax.tick_params(axis='x', labelrotation=270)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.tight_layout()
plt.savefig(outdir + '/images/top_50_clusters_count.png')


fig, ax = plt.subplots(figsize = (5, 12))
ax.barh(stats50['Name'] + ' [' + stats50['Cluster'].str.replace('Cluster ', 'C') + ']', stats50['Count'])
ax.set_xlabel('Number of Splcing Instances')
ax.set_title('Top 50 Genes')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.gca().invert_yaxis()
fig.tight_layout()
plt.savefig(outdir + '/images/top_50_spliced_genes.png')
plt.close()



# Splicing Volcano
temp_vol = full_data
temp_vol['abs'] = abs(temp_vol['IncLevelDifference'])
temp_vol = temp_vol.loc[temp_vol.groupby('geneSymbol').abs.idxmax()].reset_index(drop = True)


fig, ax = plt.subplots(figsize = (6, 6))
ax.scatter(temp_vol[temp_vol['FDR'] > 0.05]['IncLevelDifference'], -np.log10(temp_vol[temp_vol['FDR'] > 0.05]['FDR']), c = 'lightgrey', s = 1)
ax.scatter(temp_vol[temp_vol['FDR'] < 0.05]['IncLevelDifference'], -np.log10(temp_vol[temp_vol['FDR'] < 0.05]['FDR']), c = 'red', s = 2)

for i in temp_vol.sort_values(by = 'FDR').iloc[0:20,:].index:
    ax.annotate(temp_vol.loc[i,'geneSymbol'], xy = (temp_vol.loc[i,'IncLevelDifference'], -np.log10(temp_vol.loc[i,'FDR'])), fontsize = 8, horizontalalignment = 'center')

ax.set_xlabel('Inclusion Level Difference')
ax.set_ylabel('-log10(FDR)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.tight_layout()
plt.savefig(outdir + '/images/volcano.png')



# Splicing Distribution
splicing_instances = pd.DataFrame(data.groupby('Splice').size().sort_values(ascending = False) * 100 / data.shape[0]).reset_index()
splicing_instances.columns = ['Splice', 'Percentage']

fig, ax = plt.subplots(figsize = (5, 5))
ax.bar(splicing_instances['Splice'], splicing_instances['Percentage'])
ax.set_ylabel('Percentage')
ax.set_title('Alternative Splicings (unique instances, FDR < 0.05)')
ax.tick_params(axis='x', labelrotation=270)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.tight_layout()
plt.savefig(outdir + '/images/splicing_instances.png')





# Visualise network
p4c.layout_network('kamada-kawai')
p4c.set_visual_style('default')
p4c.set_node_shape_default('ELLIPSE', style_name='default')
p4c.lock_node_dimensions(True, style_name='default')
p4c.set_node_size_default(50, style_name='default')
p4c.set_node_color_default('#D3D3D3', style_name='default')
p4c.set_node_border_width_default(2, style_name='default')
p4c.set_node_color_default('#616060', style_name='default')
p4c.set_node_label_mapping('display name', style_name='default')
p4c.set_node_font_size_default(14, style_name='default')

p4c.copy_visual_style('default', 'Targets')
p4c.set_visual_style('Targets')
p4c.set_node_color_mapping(**p4c.gen_node_color_map('Count', p4c.palette_color_brewer_s_YlOrRd(), style_name='Targets'))
p4c.export_image(outdir + '/images/networks/counts', type = 'PNG', height = 1200, width = 1400, overwrite_file = True, force_pre_3_10 = True)

p4c.set_node_color_mapping(**p4c.gen_node_color_map('Cluster', p4c.palette_color_brewer_d_Spectral(), style_name='Targets', mapping_type = 'd'))
p4c.export_image(outdir + '/images/networks/mcode_clusters', type = 'PNG', height = 1200, width = 1400, overwrite_file = True, force_pre_3_10 = True)

scores = nodes_table.columns[nodes_table.columns.str.startswith('tissue::') | nodes_table.columns.str.startswith('compartment::')].tolist()
for i in scores:
    p4c.set_node_color_mapping(**p4c.gen_node_color_map(i, p4c.palette_color_brewer_s_YlOrRd(), style_name='Targets'))
    p4c.export_image(outdir + '/images/networks/' + i.replace('::', '_').replace(' ', '_'), type = 'PNG', height = 1200, width = 1400, overwrite_file = True, force_pre_3_10 = True)



res = p4c.commands.commands_post(f'layout attributes-layout nodeAttribute="Cluster"')
p4c.set_visual_style('Targets')
p4c.set_node_color_mapping(**p4c.gen_node_color_map('Cluster', p4c.palette_color_brewer_d_Spectral(), style_name='Targets', mapping_type = 'd'))
p4c.export_image(outdir + '/images/networks/mcode_clusters_attributes_layout', type = 'PNG', height = 1200, width = 1400, overwrite_file = True, force_pre_3_10 = True)



# Specific visualisation (Not run unless required)
def pivot_enrichment(cat = None):
    temp_subset = enrichment[enrichment['category'] == cat][['genes', 'description']]
    temp_subset = temp_subset.explode('genes')
    temp_subset['Yes'] = 1
    temp_subset = temp_subset[~temp_subset[['genes', 'description']].duplicated()]
    to_return = temp_subset.pivot(index = 'genes', columns = 'description', values = 'Yes')
    to_return.columns = cat + ': ' + to_return.columns
    to_return = to_return.fillna(0)
    return to_return

def add_specific_term(cat = None, term = None):
    temp_df = pivot_enrichment(cat)
    temp_column = temp_df[cat + ': ' + term].to_frame().reset_index()
    p4c.load_table_data(temp_column, data_key_column = 'genes', table_key_column = 'name')

#add_specific_term(cat = 'GO Cellular Component', term = 'Integral component of plasma membrane')
#p4c.set_node_color_mapping(**p4c.gen_node_color_map('GO Cellular Component: Integral component of plasma membrane', p4c.palette_color_brewer_s_YlOrRd(), style_name='Targets'))



# Cluster-level re-analysis
def cluster_reanalysis(clus = None):
    entrez_list = nodes_table[nodes_table['Cluster'] == clus]['display name']
    if len(entrez_list) < 10:
        return 'Cluster is too small.'

    if len(conversion_dedup[conversion_dedup['Name'].isin(entrez_list)]) < 1:
        return 'Irrelevant cluster.'

    os.makedirs(outdir + '/clusters/' + clus.replace(' ', '_'))
    os.makedirs(outdir + '/clusters/' + clus.replace(' ', '_') + '/tables')
    os.makedirs(outdir + '/clusters/' + clus.replace(' ', '_') + '/images')
    os.makedirs(outdir + '/clusters/' + clus.replace(' ', '_') + '/images/enrichment')
    os.makedirs(outdir + '/clusters/' + clus.replace(' ', '_') + '/images/networks')
    p4c.cytoscape_version_info()
    p4c.session.close_session(False)
#    entrez_list = nodes_table[nodes_table['Cluster'] == clus]['display name']
    query_str = ','.join(str(v) for v in entrez_list)
    string_cmd_list = ['string protein query','query="',query_str,'"', 'species="Homo sapiens"', 'cutoff=0.4']
    string_cmd = " ".join(string_cmd_list)
    p4c.commands.commands_run(string_cmd)
    p4c.load_table_data(conversion_dedup[conversion_dedup['Name'].isin(entrez_list)].reset_index(drop = True), data_key_column = 'Name', table_key_column = 'display name')
    p4c.commands.commands_run('string retrieve enrichment allNetSpecies="Homo sapiens"')
    p4c.commands.commands_run('string show enrichment')
#    suid = p4c.networks.get_network_suid('STRING network')
    suid = p4c.networks.get_network_suid()
    all_tables = p4c.commands.cyrest_get(f'networks/{suid}/tables')
    titles = []
    for i in all_tables:
        titles.append(i.get('title'))

    nodes_table_c = pd.DataFrame(p4c.commands.cyrest_get(f'networks/{suid}/tables/defaultnode').get('rows'))
    enrichment_c = pd.DataFrame(all_tables[titles.index('STRING Enrichment: All')].get('rows'))

    nodes_table_c.to_csv(outdir + '/clusters/' + clus.replace(' ', '_') + '/tables/nodes_table.csv')
    enrichment_c.to_csv(outdir + '/clusters/' + clus.replace(' ', '_') + '/tables/enrichment.csv')

    tissues_c = nodes_table_c.columns[nodes_table_c.columns.str.startswith('tissue::')]
    compartments_c = nodes_table_c.columns[nodes_table_c.columns.str.startswith('compartment::')]

    to_plot_tissues_t = pd.melt(nodes_table_c[tissues_c], var_name = 'Tissue', value_vars = tissues_c, value_name = 'Score').fillna(0)
    to_plot_tissues_t['Tissue'] = to_plot_tissues_t['Tissue'].str.replace('tissue::', '')
    to_plot_compartments_t = pd.melt(nodes_table_c[compartments_c], var_name = 'Compartment', value_vars = compartments_c, value_name = 'Score').fillna(0)
    to_plot_compartments_t['Compartment'] = to_plot_compartments_t['Compartment'].str.replace('compartment::', '')
    plot_ridge(to_plot_tissues_t, filename = 'tissues_' + clus.replace(' ', '_')  + '_only', savedir = outdir + '/clusters/' + clus.replace(' ', '_') + '/images')
    plot_ridge(to_plot_compartments_t, filename = 'compartment_' + clus.replace(' ', '_')  + '_only', savedir = outdir + '/clusters/' + clus.replace(' ', '_') + '/images')


    for i in enrichment['category'].unique():
        plot_enrichment(e_file = enrichment_c, cat = i, savedir = outdir + '/clusters/' + clus.replace(' ', '_') + '/images/enrichment')

    p4c.layout_network('kamada-kawai')
    p4c.set_visual_style('default')
    p4c.set_node_shape_default('ELLIPSE', style_name='default')
    p4c.lock_node_dimensions(True, style_name='default')
    p4c.set_node_size_default(50, style_name='default')
    p4c.set_node_color_default('#D3D3D3', style_name='default')
    p4c.set_node_border_width_default(2, style_name='default')
    p4c.set_node_color_default('#616060', style_name='default')
    p4c.set_node_label_mapping('display name', style_name='default')
    p4c.set_node_font_size_default(14, style_name='default')

    p4c.copy_visual_style('default', 'Targets')
    p4c.set_visual_style('Targets')
    p4c.set_node_color_mapping(**p4c.gen_node_color_map('Count', p4c.palette_color_brewer_s_YlOrRd(), style_name='Targets'))
    p4c.export_image(outdir + '/clusters/' + clus.replace(' ', '_') + '/images/networks/counts', type = 'PNG', height = 1200, width = 1400, overwrite_file = True, force_pre_3_10 = True)

    scores = nodes_table_c.columns[nodes_table_c.columns.str.startswith('tissue::') | nodes_table_c.columns.str.startswith('compartment::')].tolist()
    for i in scores:
        p4c.set_node_color_mapping(**p4c.gen_node_color_map(i, p4c.palette_color_brewer_s_YlOrRd(), style_name='Targets'))
        p4c.export_image(outdir + '/clusters/' + clus.replace(' ', '_') + '/images/networks/' + i.replace('::', '_').replace(' ', '_'), type = 'PNG', height = 1200, width = 1400, overwrite_file = True, force_pre_3_10 = True)




for i in nodes_table['Cluster'].unique():
    if str(i) != 'nan':
        cluster_reanalysis(clus = i)

#        to_plot_tissues_c = pd.melt(nodes_table[nodes_table['Cluster'] == i][tissues], var_name = 'Tissue', value_vars = tissues, value_name = 'Score').fillna(0)
#        to_plot_tissues_c['Tissue'] = to_plot_tissues_c['Tissue'].str.replace('tissue::', '')
#        to_plot_compartments_c = pd.melt(nodes_table[nodes_table['Cluster'] == i][compartments], var_name = 'Compartment', value_vars = compartments, value_name = 'Score').fillna(0)
#        to_plot_compartments_c['Compartment'] = to_plot_compartments_c['Compartment'].str.replace('compartment::', '')
#        plot_ridge(to_plot_tissues_c, filename = 'tissues_' + i.replace(' ', '_') + '_global', savedir = outdir + '/clusters/' + i.replace(' ', '_') + '/images')
#        plot_ridge(to_plot_compartments_c, 'compartment_' + i.replace(' ', '_') + '_global', savedir = outdir + '/clusters/' + i.replace(' ', '_') + '/images')


