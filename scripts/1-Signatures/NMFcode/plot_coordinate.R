library(data.table)
f<-fread("rank21.ordered.coordinate.matrix")
m<-as.matrix(f)
rownames(m)<-c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10",
                "s11","s12","s13","s14","s15","s16","s17","s18","s19","s20","s21")

library("reshape2")
tmp <- melt(m)
head(tmp)

library("ggplot2")

##DHS colors
cl<-c("s1"="#fee423","s2"="#f8a059","s3"="#ed1f23","s4"="#2cb149","s5"="#719748","s6"="#5a833c","s7"="#2fc2d7","s8"="#436bb4","s9"="#099487","s10"="#9b4f9d","s11"="#6b55a3","s12"="#536e7e","s13"="#1e2959","s14"="#be502a","s15"="#713119","s16"="#c2c2c2","s17"="#55FF55","s18"="#84A98C","s19"="#FFCCD5","s20"="#fbc0ff","s21"="#a2a832")

##DHS colors+mix
#cl<-c("s1"="#c2c2c2","s2"="#f8a059","s3"="#ed1f23","s4"="#2cb149","s5"="#e26ddf","s6"="#9f66fa","s7"="#2fc2d7","s8"="#436bb4","s9"="#099487","s10"="#c2c2c2","s11"="#6b55a3","s12"="#536e7e","s13"="#1e2959","s14"="#c2c2c2","s15"="#713119","s16"="#fee423","s17"="#9b4f9d","s18"="#be502a","s19"="#FFCCD5","s20"="#00BBF9")

##TCGA cancer color
# cl<-c("s1"="#c2c2c2","s2"="#82fd57","s3"="#fe6801","s4"="#a983f2","s5"="#fe0000","s6"="#aa059f","s7"="#262b6b","s8"="#43b9f9","s9"="#1fbf81","s10"="#9b4f9d","s11"="#ffb55c","s12"="#a30d0e","s13"="#83ffff","s14"="#be502a","s15"="#ff7d7d","s16"="#fe0000","s17"="#fbc0ff","s18"="#9bc48a","s19"="#002bfe","s20"="#e26ddf")

##TCGA cancer color + mix same color
# cl<-c("s1"="#c2c2c2","s2"="#82fd57","s3"="#fe6801","s4"="#a983f2","s5"="#fe0000","s6"="#aa059f","s7"="#262b6b","s8"="#43b9f9","s9"="#1fbf81","s10"="#c2c2c2","s11"="#ffb55c","s12"="#a30d0e","s13"="#83ffff","s14"="#c2c2c2","s15"="#ff7d7d","s16"="#fe0000","s17"="#fbc0ff","s18"="#9bc48a","s19"="#002bfe","s20"="#e26ddf")


gg<-ggplot(tmp, aes(x=Var2, y=value, fill=Var1))+geom_bar(stat = "identity",width=1)+ scale_fill_manual(values = cl)+ theme_minimal()+ theme(legend.position = "bottom",axis.text.x = element_blank())+ xlab("Chromatin Variants") + ylab("Value")+guides(fill=guide_legend(title="Signatures"))

ggsave(file.path("new.rank21.coordinates.plot.png"),gg,width = 30,height = 12,units = "cm")
#ggsave(file.path("new.rank21.coordinates.plot.pdf"),gg,device = "pdf",width = 30,height = 12,units = "cm")
