library(pheatmap)
library(data.table)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

rank <- args[1]

print(paste0("working on rank ",rank))


mat <- fread(paste0("rank",rank,".png.order.matrix"))
dat <- as.data.frame(mat)
dim(dat)
colnames(dat)<-colnames(mat)
si = read.table("GDC_identifiers_all.tsv", header=T, sep = "\t")
si_distinct<-si %>% distinct(associated_entities.entity_submitter_id, .keep_all = TRUE)
n<-colnames(mat)

group<-data.frame()
f<-data.frame()
for (i in 1:length(n)) {
  de<-si_distinct[si_distinct$File.Name==paste0(n[i],"_atacseq_gdc_realn.bam"),]$Project
  #print(paste0(de,n[i]))
  if (length(de)==0) {
    group<-rbind(group,"ALL")
    f<-rbind(f,"ALL")
  }
  else{
    group<-rbind(group,de)
    f<-rbind(f,"TCGA")
  }
}
flist<-cbind(group,f)
names(flist)<-c("cancertypes","group")
rownames(flist)<-n
dim(flist)
print(head(flist))

###BRCA subtypes
bca<-si_distinct[si_distinct$Project=="TCGA-BRCA",]
dim(bca)
bcasubtypes<-read.table("tcga-bcasubtypes", header=F, sep = "\t")
dim(bcasubtypes)
fn<-data.frame()
for (i in 1:nrow(bcasubtypes)) {
  de<-si_distinct[si_distinct$cases.case_id==bcasubtypes$V1[i],]$File.Name
  fn <- rbind(fn, de)
}
bcasubtypes<-cbind(bcasubtypes,fn)
colnames(bcasubtypes)<-c("caseid","subtype","filename")
extra<-data.frame(cbind("b0f8d698-a30e-4d8d-b0a2-a5a01fac8406","LumB","00a595c5-0565-4289-878a-80ff4ba89d97_atacseq_gdc_realn.bam"))
colnames(extra)<-c("caseid","subtype","filename")
bcasubtypes<-rbind(bcasubtypes,extra)
dim(bcasubtypes)

stype<-data.frame()
ni<-n
for (i in 1:length(ni)) {
  de<-bcasubtypes[bcasubtypes$filename==paste0(n[i],"_atacseq_gdc_realn.bam"),]$subtype
  print(de)
  if (length(de)==0) {
    stype <- rbind(stype, NA)
  }
  else{
    stype <- rbind(stype, de)
  }
}
names(stype)<-"bca_subtype"

print(stype)






Categories <- data.frame(cancertypes = factor(flist$cancertypes),bcasubtypes=factor(stype$bca_subtype))
rownames(Categories)<-n
cl<-c("TCGA-ACC"="#ff7d7d","TCGA-BRCA"="#fe0000","TCGA-SKCM"="#a30d0e","TCGA-CESC"="#612e01","TCGA-GBM"="#a37856","TCGA-COAD"="#fe6801","TCGA-PCPG"="#ffb55c","TCGA-BLCA"="#fbe016","TCGA-HNSC"="#c1b52b","TCGA-MESO"="#9bc48a","TCGA-KIRP"="#82fd57", "TCGA-ESCA"="#026b02","TCGA-LUSC"="#1fbf81", "TCGA-TGCT"="#83ffff","TCGA-LIHC"="#43b9f9", "TCGA-LUAD"="#002bfe","TCGA-PRAD"="#262b6b","TCGA-LGG"="#a983f2","TCGA-KIRC"="#7e2ad8","TCGA-THCA"="#aa059f","TCGA-STAD"="#e26ddf","TCGA-UCEC"="#fbc0ff","TCGA-CHOL"="#757575")
bcasub = c("Basal"="indianred4","Her2"="royalblue2","LumA"="indianred1","LumB"="olivedrab4","NA"="grey","Normal"="yellow")

CategoriesColor <- list(cancertypes = cl,bcasubtypes=bcasub)

png(paste0("bca-rank",rank,"ordered.png" ))
pheatmap(dat,cluster_rows = F,cluster_cols = F, main=paste0("NMF rank",rank), color=colorRampPalette(c( "white", "navy"))(50),show_rownames=F,show_colnames=F,annotation_col = Categories,annotation_colors = CategoriesColor,border_color =NA,treeheight_row=0)

#png(paste0("rank",rank,"ordered.png" ))
dev.off()

