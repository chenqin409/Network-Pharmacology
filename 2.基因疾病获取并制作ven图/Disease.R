#install.packages("venn")

library(venn)              #引用包
outFile="Disease.txt"      #输出文件
setwd("C:\\Users\\chenq\\My Work\\Network Pharmacology\\2.基因疾病获取并制作ven图")    #设置工作目录
files=dir()                        #获取所有文件
files=grep("txt$",files,value=T)   #提取.txt的文件
geneList=list()

#读取txt文件中的基因信息，保存到geneList
for(i in 1:length(files)){
    inputFile=files[i]
	if(inputFile==outFile){next}
    rt=read.table(inputFile,header=F)        #读取输入文件
    geneNames=as.vector(rt[,1])              #读取基因名称
    geneNames=gsub("^ | $","",geneNames)     #去掉基因收尾空格
    uniqGene=unique(geneNames)               #基因取unique
    header=unlist(strsplit(inputFile,"\\.|\\-"))
    geneList[[header[1]]]=uniqGene
    uniqLength=length(uniqGene)
    print(paste(header[1],uniqLength,sep=" "))
}

#绘制venn图
mycol=c("#029149","#E0367A","#5D90BA","#431A3D","#91612D","#FFD121","#D8D155","#223D6C","#D20A13","#088247","#11AA4D","#7A142C","#5D90BA","#64495D","#7CC767")
pdf(file="venn.pdf",width=5,height=5)
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F,ilabels = "counts")
dev.off()

#保存并集基因
unionGenes=Reduce(union,geneList)
write.table(file=outFile,unionGenes,sep="\t",quote=F,col.names=F,row.names=F)

