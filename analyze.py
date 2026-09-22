"""Generate network pharmacology figures from preserved input and database tables."""
from pathlib import Path
import argparse, json, textwrap
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib_venn import venn2

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--data',type=Path,default=ROOT/'数据')
parser.add_argument('--out',type=Path,default=ROOT)
args=parser.parse_args()
D=args.data;OUT=args.out;FIG=OUT/'图件';TAB=OUT/'数据'/'结果'
FIG.mkdir(parents=True,exist_ok=True);TAB.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','DejaVu Serif'],
    'font.size':12,'axes.titlesize':17,'axes.titleweight':'bold','axes.labelsize':13,
    'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,
    'savefig.facecolor':'white','figure.facecolor':'white'})
TEAL='#087F8C';CORAL='#D9795B';PURPLE='#8865A4';INK='#203441';GRAY='#7B8C95'
ORDER=['Luteolin','Apigenin','Kaempferol','Quercetin','Chrysoeriol','Genistein','Isorhamnetin','Biochanin A','Tricetin','Quercitrin']
def csv(df,name): df.to_csv(TAB/(name+'.csv'),index=False,encoding='utf-8-sig')
def save(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=300,bbox_inches='tight',metadata={'Software':''})
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',metadata={'Author':'','Creator':'','Producer':'','CreationDate':None,'ModDate':None})
    plt.close(fig);print(name,flush=True)
def footer(fig,s): fig.text(.5,.025,s,ha='center',va='bottom',fontsize=10,color=INK)
rec=pd.read_csv(D/'compound_target_records.csv',keep_default_na=False)
dis=pd.read_csv(D/'AGA_direct_all_targets.csv')
ann=pd.read_csv(D/'uniprot_annotation.csv')
mapping=pd.read_csv(D/'STRING_mapping.csv')
edges=pd.read_csv(D/'STRING_edges.csv')
en=json.loads((D/'参考'/'STRING_enrichment.json').read_text(encoding='utf-8'))
sen=json.loads((D/'参考'/'STRING_enrichment_target_background.json').read_text(encoding='utf-8'))
assert set(rec.compound)==set(ORDER)
assert len(rec[['source_file','source_row']].drop_duplicates())==1000
assert set(ann['Organism (ID)'])=={9606}
assert set(rec.accession)==set(ann.Entry)
targets=set(rec.gene);disease=set(dis.loc[dis.association_score>0,'gene']);shared=targets&disease
single=set(rec.loc[rec.evidence_scope=='single_protein','gene'])
assert shared==single&disease
csv(pd.DataFrame({'gene':sorted(targets),'has_single_protein_record':[g in single for g in sorted(targets)]}),'成分靶点并集')
csv(dis.loc[dis.gene.isin(shared)].sort_values('association_score',ascending=False),'交集靶点与疾病证据')
pairs=rec.groupby(['compound','gene'],sort=False).agg(
    single_protein_record=('evidence_scope',lambda x:'single_protein' in set(x)),
    complex_member_record=('evidence_scope',lambda x:'complex_member' in set(x)),
    max_record_probability=('probability','max')).reset_index()
csv(pairs,'成分靶点完整关系')
se=pairs[pairs.gene.isin(shared)].copy();csv(se,'成分交集靶点关系')
compound_counts=[]
for c in ORDER:
    x=pairs[pairs.compound==c]
    compound_counts.append({'compound':c,'all_genes':x.gene.nunique(),'intersection_genes':x[x.gene.isin(shared)].gene.nunique(),
       'single_protein_genes':x[x.single_protein_record].gene.nunique()})
csv(pd.DataFrame(compound_counts),'成分靶点计数')

# Area-proportional set overlap. The union includes explicitly annotated complex members.
fig,ax=plt.subplots(figsize=(9,7))
v=venn2([targets,disease],set_labels=(f'10 flavonoids\n{len(targets)} mapped genes',f'Androgenetic alopecia\n{len(disease)} associated genes'),set_colors=(TEAL,CORAL),alpha=.65,ax=ax)
for t in v.subset_labels:
    if t is not None:t.set_fontsize(20);t.set_fontweight('bold')
for t in v.set_labels:t.set_fontsize(13)
ax.set_title('Shared genes between flavonoids and AGA',pad=24)
footer(fig,'Open Targets MONDO_0005339; direct associations, score > 0; Homo sapiens.\nMapped input includes complex members; excluding complex-only genes retains the same 36 shared genes.')
fig.subplots_adjust(bottom=.18,top=.86);save(fig,'01_靶点交集_Venn')

# Tripartite provenance network with separately styled complex membership edges.
genes=sorted(shared)
pos={g:(5*np.cos(t),5*np.sin(t)) for g,t in zip(genes,np.linspace(np.pi/2,np.pi/2-2*np.pi,len(genes),endpoint=False))}
pos.update({c:(2.45*np.cos(t),2.45*np.sin(t)) for c,t in zip(ORDER,np.linspace(np.pi/2,np.pi/2-2*np.pi,len(ORDER),endpoint=False))})
source='Platycladi Cacumen\nPCDENs';pos[source]=(0,0)
fig,ax=plt.subplots(figsize=(14,14))
for r in se.itertuples():
    x,y=pos[r.compound],pos[r.gene]
    ax.plot([x[0],y[0]],[x[1],y[1]],color=TEAL if r.single_protein_record else PURPLE,
       alpha=.22 if r.single_protein_record else .8,lw=.7 if r.single_protein_record else 1.2,
       ls='-' if r.single_protein_record else '--',zorder=1)
for c in ORDER:
    ax.plot([0,pos[c][0]],[0,pos[c][1]],color=CORAL,alpha=.65,lw=1.5,zorder=1)
ax.scatter([pos[g][0] for g in genes],[pos[g][1] for g in genes],s=150,color=CORAL,edgecolor='white',zorder=3)
for g in genes:
    x,y=pos[g];ax.text(x*1.09,y*1.09,g,fontsize=11.5,ha='center',va='center',fontweight='bold',color=INK)
for c in ORDER:
    x,y=pos[c];ax.text(x,y,c,ha='center',va='center',fontsize=11.5,fontweight='bold',color=TEAL,
        bbox={'boxstyle':'round,pad=.45','fc':'white','ec':TEAL,'lw':1.5},zorder=4)
ax.scatter([0],[0],s=6000,marker='D',color='#D9EBD9',edgecolor='#57845B',zorder=3)
ax.text(0,0,source,ha='center',va='center',fontsize=12,fontweight='bold',color=INK,zorder=4)
ax.set(xlim=(-6,6),ylim=(-6,6),aspect='equal');ax.axis('off')
ax.set_title('Herb–flavonoid–AGA shared-gene network\n'+f'1 source | 10 flavonoids | 36 genes | {len(se)} compound–gene links',pad=18)
ax.legend(handles=[Line2D([0],[0],color=CORAL,label='Source–compound'),Line2D([0],[0],color=TEAL,label='Single-protein prediction'),
 Line2D([0],[0],color=PURPLE,ls='--',label='Complex-member association only')],loc='lower center',bbox_to_anchor=(.5,-.055),ncol=3,frameon=False,fontsize=11)
footer(fig,'Compound identities follow the supplied files and Figure 2.8C labels.\nComplex-member links do not establish independent binding to each subunit; no abundance or activity is inferred.')
fig.subplots_adjust(bottom=.13,top=.90);save(fig,'02_药材成分交集靶点网络')

# Membership matrix gives an unambiguous readable view of every network edge.
mat=np.zeros((len(genes),len(ORDER)),dtype=int)
for r in se.itertuples():mat[genes.index(r.gene),ORDER.index(r.compound)]=2 if r.single_protein_record else 1
fig,ax=plt.subplots(figsize=(10,13))
ax.imshow(mat,cmap=ListedColormap(['#F2F4F5','#C7AEDB',TEAL]),norm=BoundaryNorm([-.5,.5,1.5,2.5],3),aspect='auto')
ax.set_xticks(range(10),ORDER,rotation=45,ha='right');ax.set_yticks(range(len(genes)),genes,fontsize=11)
ax.set_xticks(np.arange(-.5,10,1),minor=True);ax.set_yticks(np.arange(-.5,len(genes),1),minor=True)
ax.grid(which='minor',color='white',linewidth=1);ax.tick_params(which='minor',bottom=False,left=False)
ax.set_title('Compound–shared-gene membership',pad=22)
ax.legend(handles=[Patch(color='#F2F4F5',label='No supplied record'),Patch(color='#C7AEDB',label='Complex-member only'),Patch(color=TEAL,label='Single-protein record')],
 loc='upper center',bbox_to_anchor=(.5,-.16),ncol=1,frameon=False,fontsize=11)
footer(fig,'Cells encode supplied prediction records, not expression, abundance, or experimental activity.')
fig.subplots_adjust(left=.14,bottom=.23,top=.93,right=.97);save(fig,'03_成分交集靶点关系热图')

# No additional interactors: network nodes are the actual STRING edge endpoints.
assert edges.score.min()>=.7
g=nx.from_pandas_edgelist(edges,'preferredName_A','preferredName_B',edge_attr='score')
assert set(g)<=shared
assert g.number_of_edges()==len(edges) and nx.number_of_selfloops(g)==0
isolated=sorted(shared-set(g));csv(pd.DataFrame({'gene':isolated,'reason':'No edge at STRING score >= 0.700'}),'无高置信度连边靶点')
bet=nx.betweenness_centrality(g,normalized=True,weight=None)
deg=dict(g.degree());close=nx.closeness_centrality(g)
central=pd.DataFrame([{'gene':x,'degree':deg.get(x,0),'betweenness':bet.get(x,0),'closeness':close.get(x,0),'in_connected_network':x in g} for x in sorted(shared)])
central=central.sort_values(['degree','gene'],ascending=[False,True])
cutoff=int(central.iloc[9].degree);central['core_candidate']=central.degree>=cutoff
csv(central,'PPI中心性与核心候选');core=central.loc[central.core_candidate,'gene'].tolist()
csv(edges,'PPI关联边');csv(central[central.core_candidate],'核心候选靶点')
nx.write_graphml(g,TAB/'PPI网络.graphml')
def network_plot(graph,title,stem,note):
    fig,ax=plt.subplots(figsize=(12,11))
    layout=nx.kamada_kawai_layout(graph,weight=None)
    # Separate nearby node centres without changing the graph or its statistics.
    node_list=list(graph)
    for _ in range(100):
        moved=False
        for i,a in enumerate(node_list):
            for b in node_list[i+1:]:
                delta=layout[a]-layout[b];dist=np.linalg.norm(delta)
                if dist<.255:
                    shift=delta/max(dist,1e-9)*(.255-dist)*.51
                    layout[a]+=shift;layout[b]-=shift;moved=True
        if not moved:break
    ds=np.array([deg[n] for n in graph])
    nx.draw_networkx_edges(graph,layout,ax=ax,edge_color=GRAY,alpha=.34,width=[1.1+2*(d['score']-.7) for a,b,d in graph.edges(data=True)])
    nodes=nx.draw_networkx_nodes(graph,layout,ax=ax,node_size=650+ds*70,node_color=ds,cmap=plt.cm.YlOrRd,vmin=0,vmax=max(deg.values()),edgecolors='white',linewidths=1.5)
    for n,(x,y) in layout.items():
        ax.text(x,y,n,ha='center',va='center',fontsize=10.5,fontweight='bold',
            color='white' if deg[n]>=12 else INK,zorder=4)
    ax.set_title(title,pad=18);ax.margins(.18);ax.axis('off')
    cb=fig.colorbar(nodes,ax=ax,shrink=.55,pad=.02);cb.set_label('Degree in the full PPI network')
    footer(fig,note);fig.subplots_adjust(top=.90,bottom=.13,left=.03,right=.9);save(fig,stem)
network_plot(g,f'High-confidence protein association network\n{len(g)} connected proteins | {g.number_of_edges()} edges','04_PPI高置信度网络',
 f'STRING v12.0; Homo sapiens; functional associations; score >= 0.700; no added proteins.\n{len(isolated)} shared genes have no qualifying edge and are retained in the accompanying table. Edges are not exclusively physical binding.')
fig,ax=plt.subplots(figsize=(9,7))
c=central[central.core_candidate].iloc[::-1]
ax.barh(c.gene,c.degree,color=TEAL,height=.68)
for i,v in enumerate(c.degree):ax.text(v+.25,i,str(v),va='center',fontsize=12)
ax.set_xlim(0,max(deg.values())+3);ax.set_xlabel('Unweighted degree');ax.set_title(f'Core candidates ranked by network degree\nTop 10 positions including ties: {len(core)} genes',pad=18)
ax.xaxis.get_major_locator().set_params(integer=True);ax.grid(axis='x',alpha=.18);ax.set_axisbelow(True)
footer(fig,f'Degree cutoff >= {cutoff}; ties at the 10th position retained. Centrality is a prioritisation criterion, not experimental validation.')
fig.subplots_adjust(left=.16,bottom=.16,top=.84,right=.95);save(fig,'05_核心候选靶点排名')
cg=g.subgraph(core).copy()
network_plot(cg,f'Core-candidate induced subnetwork\n{len(cg)} proteins | {cg.number_of_edges()} edges','06_核心候选靶点子网络',
 f'Candidates selected by degree >= {cutoff} in the full network, including all cutoff ties.\nNode size and colour use full-network degree; edges retain the original STRING score threshold.')

# Retain the API's reported enrichment FDR; do not re-adjust a pre-filtered result list.
enrich=pd.DataFrame(en)
for col in ['inputGenes','preferredNames']:
    enrich[col]=enrich[col].map(lambda x:';'.join(x))
enrich['query_gene_count']=len(shared);enrich['gene_ratio']=enrich.number_of_genes/len(shared)
enrich['background']='STRING v12.0 human proteome (default)'
csv(enrich,'STRING返回富集条目')
sens=pd.DataFrame(sen)
for col in ['inputGenes','preferredNames']:
    if col in sens:sens[col]=sens[col].map(lambda x:';'.join(x))
sens['background']='214 input-derived proteins';csv(sens,'候选靶点背景敏感性富集')
catmap=[('Process','GO biological process','07_GO生物过程',20),('Component','GO cellular component','08_GO细胞组分',15),
 ('Function','GO molecular function','09_GO分子功能',15),('KEGG','KEGG pathways','10_KEGG通路',20)]
for cat,title,stem,limit in catmap:
    table=enrich[enrich.category==cat].sort_values(['fdr','p_value','term'])
    csv(table,stem+'_全部返回条目')
    top=table[table.fdr<.05].head(limit).iloc[::-1].copy();csv(top.iloc[::-1],stem+'_展示条目')
    fig,ax=plt.subplots(figsize=(13,10 if len(top)>15 else 8.5))
    if top.empty:
        ax.text(.5,.5,'No terms pass FDR < 0.05',ha='center',va='center',transform=ax.transAxes);ax.axis('off')
    else:
        xx=top.gene_ratio.to_numpy();yy=np.arange(len(top));color=-np.log10(top.fdr.to_numpy())
        dots=ax.scatter(xx,yy,s=top.number_of_genes.to_numpy()*17,c=color,cmap='YlOrRd',edgecolors='#805038',linewidths=.4)
        labels=[textwrap.fill(s,48) for s in top.description]
        ax.set_yticks(yy,labels,fontsize=11);ax.set_ylim(-.8,len(top)-.2)
        ax.set_xlabel('Gene ratio (overlap / 36 input genes)');ax.grid(axis='x',alpha=.20);ax.set_axisbelow(True)
        ax.set_xlim(max(0,xx.min()-.055),min(1,xx.max()+.065))
        cb=fig.colorbar(dots,ax=ax,shrink=.60,pad=.04);cb.set_label(r'$-\log_{10}$(FDR)')
        ns=sorted(set([int(top.number_of_genes.min()),int(np.median(top.number_of_genes)),int(top.number_of_genes.max())]))
        handles=[ax.scatter([],[],s=n*17,color='#B1B9BD',edgecolors='#805038',linewidths=.4,label=str(n)) for n in ns]
        ax.legend(handles=handles,title='Gene count',loc='upper left',bbox_to_anchor=(1.22,.18),frameon=False,labelspacing=1.3,fontsize=10,title_fontsize=11)
    ax.set_title(title+f' enrichment\nTop {len(top)} terms by FDR',pad=18)
    footer(fig,'36 shared genes; STRING v12.0; default human-proteome background; reported FDR < 0.05.\nExploratory over-representation only; candidate-protein background sensitivity results are supplied separately.')
    fig.subplots_adjust(left=.42,right=.86,bottom=.14,top=.88);save(fig,stem)

notes=[
 ('Disease','Open Targets MONDO_0005339, androgenetic alopecia; direct associations only; enableIndirect=false; score > 0.'),
 ('Disease scope','904 genes retained, including 899 protein-coding genes and 5 lncRNAs. None of the 36 shared genes is a lncRNA.'),
 ('Website count','The website default includes the descendant MONDO_0007184 and displays 964 targets; the direct API query yields 904. These are different scopes.'),
 ('Input','10 supplied SwissTargetPrediction tables; 100 rows each. All supplied rows retained; no compound re-screening or new prediction.'),
 ('Compound identity','Names follow supplied filenames and thesis Figure 2.8C. insorehamnetin is mapped to Isorhamnetin by the supplied figure. Original prediction structures were not supplied or independently verified.'),
 ('Complexes','32 source rows describe complexes. Member links are marked explicitly and do not establish independent subunit binding. The 214-gene union includes 12 complex-only genes; the 36-gene intersection is unchanged on excluding them.'),
 ('Probability','Supplied probabilities are retained as prediction metadata, not used as efficacy, abundance, expression or calibrated between-compound weights.'),
 ('PPI','STRING v12.0 functional network; human 9606; combined score >= 0.700; no added interactors; 25 connected nodes and 83 undirected edges. Eleven genes without qualifying edges remain in the tables.'),
 ('Core candidates',f'Unweighted full-network degree; top 10 positions with all ties included: degree >= {cutoff}, {len(core)} candidates. Betweenness and closeness are supporting metrics.'),
 ('Enrichment','All 36 shared genes, including PPI isolates, entered STRING enrichment. Default human-proteome background. Reported STRING FDR used without re-adjusting the returned subset; FDR < 0.05 shown.'),
 ('Returned terms','Tables contain terms returned by STRING, not all tested terms. Unreturned terms must not be treated as zero p-value or FDR=1.'),
 ('Sensitivity','Using all 214 supplied candidate proteins as background, STRING returned only Cell development (GO:0048468), FDR 0.0078. No KEGG term was returned. Whole-proteome enrichment is sensitive to target ascertainment and is not evidence of AGA-specific therapeutic effects.'),
 ('Interpretation','Networks are computational associations. No activation/inhibition direction, efficacy, animal validation, novel mechanism or independent binding is established by these analyses.'),
 ('Reference','Thesis Figure 2.8C identifies the ten compounds only. No thesis figure or experimental value was reused to generate these figures.')]
csv(pd.DataFrame(notes,columns=['item','detail']),'分析参数与解释边界')
counts={'compounds':10,'source_rows':1000,'complex_rows':32,'mapped_genes':len(targets),'disease_genes':len(disease),'intersection':len(shared),
    'compound_shared_links':len(se),'ppi_nodes':len(g),'ppi_edges':g.number_of_edges(),'ppi_isolates':len(isolated),'core_candidates':len(core),'core_degree_cutoff':cutoff,'figure_count':10}
csv(pd.DataFrame(counts.items(),columns=['metric','value']),'分析计数')
print(json.dumps(counts),flush=True)
