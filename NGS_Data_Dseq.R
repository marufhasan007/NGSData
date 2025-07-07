#library(GMD)
library(genefilter)
library(arrayQualityMetrics)
library(affy)
library(qvalue)
setwd("Z:/FBNProjectGroups/I3-NGS-DexaPhos/Maruf_VitaminD/RNA-seq/DESeq2")

tissue=c("Liv")     #J K
#sex=c("M")
effect=c("Group")
outlier=c("None")                   ##None
anno<-read.table("annotated.ensembl.gene.name.txt",header=FALSE,sep="\t")
colnames(anno)<-c("SSC_ID","gene")
pheno<-read.table("Liver_AFBI-UVB_pheno_subset.csv",header=TRUE,sep=",")

f.path<-"Z:/FBNProjectGroups/I3-NGS-DexaPhos/Maruf_VitaminD/RNA-seq/htseq.gene.count"
file<-as.data.frame(dir(f.path))
f<-as.data.frame(file[grep("GCount.txt",file[,1]),])
for (u in 1:dim(f)[1]) {
data1 <- read.table(paste0(f.path,"/",f[u,1]),header=FALSE,sep="\t")
colnames(data1)<-c("V1",gsub(".GCount.txt","",f[u,1]))
if (u==1) {
  x<-as.data.frame(data1)
}else {
x<-(x,data1)
}
}
y<-x[grep("^[^__]",x[,1]),]
row.names(y)<-y[,1];
y<-y[,-1]

write.csv(y,"raw_counts_all.csv")
#h<-subset(y,select=c(colnames(y)[grep(paste("^",paste(tissue,collapse="|^"),sep=""),colnames(y))]))
#h<-subset(h,select=c(colnames(h)[!grepl(paste(paste(paste0(outlier,"$"),collapse="|"),sep=""),colnames(h))]))
#dim(h)
#colnames(y)<-gsub("^.*_","",colnames(y))
h<-y

p<-data.frame(ID=colnames(h))
#p$sample<-gsub("[A-Z]+","",p[,1])

i<-merge(p,pheno,by.x = "ID",by.y="RNAseq", sort=FALSE,all.y=FALSE)

i$Group<-as.factor(i$Group)
i$Birth.wt<-as.factor(i$Birth.wt)
i$avdg.14.fin<-as.factor(i$avdg.14.fin)

###DESeq2####
library('DESeq2')
dds1<-DESeqDataSetFromMatrix(countData=h,colData=i,design=as.formula(paste("~",paste(effect,collapse="+"))))
#dds1<-DESeqDataSetFromMatrix(countData=h,colData=i,design=~Group + FCR14.fin)

#m1<-model.matrix(~ group + Sex + Di?t.Sau:Sau.n,colData(dds1))
#all.zero <- apply(m1, 2, function(x) all(x==0))
#idx <- which(all.zero)
#m1 <- m1[,-idx]

keep <- rowSums(counts(dds1) >= 50) >= 4  #number of counts, number of observations
dds1 <- dds1[keep,]
dim(dds1)

dds<-DESeq(dds1,test="Wald",fitType="local", minReplicatesForReplace=3)
#dds<-DESeq(dds1,test="LRT",fitType="local",reduced= ~1, minReplicatesForReplace=3)

plotDispEsts(dds)
#e = ExpressionSet(assay(vst(dds)), AnnotatedDataFrame(as.data.frame(colData(dds))))
#arrayQualityMetrics(e, outdir=paste("Quality_",paste(tissue,collapse="+"),"DESeq",sep="_"),intgroup=c("UVB"), force=T)

res<-results(dds)
#res<-results(dds,contrast=c("group","HH","HM"),lfcThreshold = 0)#,addMLE = T)
res<-res[order(res$padj),]
head(res)
#mcols(res,use.names=TRUE)

mcols(dds)$maxCooks <- apply(assays(dds)[["cooks"]], 1, max)
plot(mcols(dds)$baseMean, mcols(dds)$maxCooks)
# this requires you not filter or order the 'res' object
stopifnot(all.equal(rownames(dds), rownames(res)))
plot(res$log2FoldChange, mcols(dds)$maxCooks)

cooksCutoff<-quantile(mcols(dds)$maxCooks,0.9999,na.rm=T) 
abline(cooksCutoff,0)

res$pvalue[mcols(dds)$maxCooks > cooksCutoff] <- NA
# optionally, also mean filtering:
#meanFilter<-quantile(res$baseMean,0.0005)
#res$pvalue[res$baseMean < meanFilter] <- NA
res$padj <- p.adjust(res$pvalue, method="BH")

#p<-as.data.frame(res$pvalue)
#qVals_hilfs <- c(); pqVals <- matrix(nrow=dim(p)[1], ncol=0);
#qVals_hilfs <- qvalue(as.vector(p[,1]),pi0.method="bootstrap")
#pqVals  <- cbind(pqVals, qVals_hilfs[[4]], qVals_hilfs[[3]])
#res<-cbind(as.data.frame(res),"qvalue"=pqVals[,2])
res<-merge(as.data.frame(res),anno,by.x=0,by.y="SSC_ID",sort=FALSE)
colnames(res)[1]<-"SSC_ID"

res<-res[order(res$pvalue),]
head(res,50)

cat(dim(subset(res,res$padj<0.05))[1]); "/n"
write.table(res, sep = ",", file=paste(paste(tissue,collapse="+"),paste(effect,collapse="+"),paste0("-",outlier,collapse=""),"anno_DESeq2.csv",sep="_"),row.names=F)
