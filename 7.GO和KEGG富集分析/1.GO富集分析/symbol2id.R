BiocManager::install("AnnotationDbi")
install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

library("BiocManager")
library("org.Hs.eg.db")          #引用包
inputFile="Drug_Disease.txt"     #输入文件名称
outFile="id.txt"                 #输出文件名称
setwd("C:\\Users\\chenq\\My Work\\Network Pharmacology\\7.GO和KEGG富集分析\\1.GO富集分析")          #设置工作目录
                                        
rt=read.table(inputFile,sep="\t",check.names=F,header=F)       #读取文件
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    #找出基因对应id
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
colnames(out)=c("symbol","entrezID")
write.table(out,file=outFile,sep="\t",quote=F,row.names=F)     #输出结果


