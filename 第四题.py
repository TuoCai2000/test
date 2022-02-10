###hg38每条染色体的基因、转录本分布
#统计每条染色体基因数、转录本数、内含子数、外显子数
def gene_count(filename):
    dic = {} #空字典
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                continue
            gtf = f.readlines()
            for i in range(1,23):
                dic[i] = {}  # 空字典里嵌套空字典
                dic[i]['gene']=0
                dic[i]['transcript']=0
                dic[i]['exon']=0
                dic[i]['CDS']=0
                dic[i]['protein_coding_gene']=0
            for k in ['X','Y','MT']:
                dic[k] = {}
                dic[k]['gene']=0
                dic[k]['transcript']=0
                dic[k]['exon']=0
                dic[k]['CDS']=0
                dic[k]['protein_coding_gene']=0
            for j in range(len(gtf)):
                gtf[j]=gtf[j].split('\t')
                name=gtf[j][0]#染色体号
                if name=='Y' or name=='X' or name=='MT':
                    type=gtf[j][2]
                    if type=='gene':
                        dic[name]['gene']+=1
                        attributes = gtf[j][8].split(';')[4]
                        if 'protein_coding' in attributes:
                            dic[name]['protein_coding_gene'] += 1
                    if type=='transcript':
                        dic[name]['transcript']+=1
                    if type=='exon':
                        dic[name]['exon']+=1
                    if type=='CDS':
                        dic[name]['exon']+=1
                else:
                    name=int(name)
                    type=gtf[j][2]
                    if type=='gene':
                        dic[name]['gene']+=1
                        attributes = gtf[j][8].split(';')[4]
                        if 'protein_coding' in attributes:
                            dic[name]['protein_coding_gene']+=1
                    if type=='transcript':
                        dic[name]['transcript']+=1
                    if type=='exon':
                        dic[name]['exon']+=1
                    if type=='CDS':
                        dic[name]['CDS']+=1
            for m,n in dic.items():
                print(m,n)
file = '/Users/caita/Desktop/第四题/Homo_sapiens.GRCh38.87.chr.gtf'
gene_count(file)


