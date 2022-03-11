#下载最新版的KEGG信息，并且解析好
import pandas as pd
def KEGG_info(filename):
    with open(filename) as f:
        dic = {}
        for line in f:
            if line.startswith('C'):
                pathway_id = line.split(' ')[4]
                dic[pathway_id]=[]
            elif line.startswith('D'):
                gene_id = line.split(' ')[6]
                dic[pathway_id].append(gene_id)
    with open('/Users/caita/Desktop/第七题/KEGG2gene.txt','wt') as f:
        for pathway_id,gene_id in dic.items():
            for gene in gene_id:
                f.write(pathway_id+'\t'+gene+'\n')
filepath ='/Users/caita/Desktop/第七题/hsa00001.keg'
KEGG_info(filepath)