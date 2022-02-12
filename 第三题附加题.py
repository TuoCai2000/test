#蛋白编码区域的GC含量会比基因组其它区域高吗
#列出蛋白编码序列的文件
list = []
import re
cds_pos = open('/Users/caita/Desktop/第三题/cds_pos.txt','w')
with open('/Users/caita/Desktop/第三题/CCDS.current.txt','r')as f1,open('/Users/caita/Desktop/第三题/chr1.fa','r') as f2:
    # 对chr1.fa进行处理
    lines = f2.readlines()
    a = ''
    for line in lines:
        if line.startswith('>'):
            continue
        else:
            a += line.strip('\n')  # 去掉每行末尾的换行符,进行叠加
    # a则为chr1的完整序列
    #对CCDS文件进行处理
    for line in f1: #逐行读取CCDS文件
        if line.startswith('#'):
            continue #去除第一行表头信息
        line = line.rstrip() #把空格去掉
        lst = line.split('\t') #以制表符分隔
        if lst[-2]=='-':
            continue #如果是-，则无外显子，跳过该行
        lst[-2] = re.sub('\[|\]','',lst[-2]) #将倒数第二列的[]改成空格
        exons = lst[9].split(',') #把每行（该基因）的多个外显子
        if lst[0]=='1':
            list.append(lst[9])
    # 取出chr1 所有CDS位置进行遍历
    for i in range(len(list)):
        pos = list[i].split(',')
        for j in range(len(pos)):
            pos[j] = pos[j].split('-')
            start = pos[j][0].replace('"','').replace(' ','')
            end = pos[j][1].replace('"','').replace(' ','')
            #cds_seq = a[start:end]
            #print(cds_seq)
            cds_pos.write(start+'\t'+end+'\n') #输出所有的外显子坐标结果
cds_pos.close()
cds_seq = ''
with open('/Users/caita/Desktop/第三题/cds_pos.txt','r') as f3:
    lines = f3.readlines()
    for line in lines:
        line = line.split('\t')
        start = int(line[0])-1
        end = int(line[1])-1 #注意碱基位置和索引的不同
        cds_seq += a[start:end]
    cds_gcpercent = (cds_seq.count('G')+cds_seq.count('g')+cds_seq.count('C')+cds_seq.count('c'))\
                     /len(cds_seq)
    cds_gcpercent =round(cds_gcpercent,2)
    print(cds_gcpercent)
