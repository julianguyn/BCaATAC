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
Categories <- data.frame(cancertypes = factor(flist$cancertypes))
rownames(Categories)<-n
cl<-c("TCGA-ACC"="#ff7d7d","TCGA-BRCA"="#fe0000","TCGA-SKCM"="#a30d0e","TCGA-CESC"="#612e01","TCGA-GBM"="#a37856","TCGA-COAD"="#fe6801","TCGA-PCPG"="#ffb55c","TCGA-BLCA"="#fbe016","TCGA-HNSC"="#c1b52b","TCGA-MESO"="#9bc48a","TCGA-KIRP"="#82fd57", "TCGA-ESCA"="#026b02","TCGA-LUSC"="#1fbf81", "TCGA-TGCT"="#83ffff","TCGA-LIHC"="#43b9f9", "TCGA-LUAD"="#002bfe","TCGA-PRAD"="#262b6b","TCGA-LGG"="#a983f2","TCGA-KIRC"="#7e2ad8","TCGA-THCA"="#aa059f","TCGA-STAD"="#e26ddf","TCGA-UCEC"="#fbc0ff","TCGA-CHOL"="#757575")
CategoriesColor <- list(cancertypes = cl)

png(paste0("rank",rank,"ordered.png" ))
pheatmap(dat,cluster_rows = F,cluster_cols = F, main=paste0("NMF rank",rank), color=colorRampPalette(c( "white", "navy"))(50),show_rownames=F,show_colnames=F,annotation_col = Categories,annotation_colors = CategoriesColor,border_color =NA,treeheight_row=0)

#png(paste0("rank",rank,"ordered.png" ))
dev.off()

