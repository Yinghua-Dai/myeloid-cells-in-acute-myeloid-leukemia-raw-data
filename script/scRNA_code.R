rm(list = ls())
suppressMessages(library(clusterProfiler))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(magrittr))
suppressMessages(library(gtools))
suppressMessages(library(stringr))
suppressMessages(library(Matrix))
suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))
suppressMessages(library(data.table))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggpubr))
suppressMessages(library(scales))
suppressMessages(library(ggsci))
suppressMessages(library(sctransform))
suppressMessages(library(harmony))
suppressMessages(library(tidydr))
suppressMessages(library(celldex))
suppressMessages(library(pheatmap))
suppressMessages(library(clustree))
suppressMessages(library(xlsx))
suppressMessages(library(gridExtra))
suppressMessages(library(ggthemes))
suppressMessages(library(ggnewscale))
suppressMessages(library(CellChat))
suppressMessages(library(ggpubr))
suppressMessages(library(patchwork))
suppressMessages(library(monocle))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(ggrepel))
suppressMessages(library(DoubletFinder))
suppressMessages(library(future))
scRNA.fis <- list.files('00_origin_datas/GEO/GSE116256_RAW/')
scRNA.fis <- scRNA.fis[endsWith(scRNA.fis, 'dem.txt.gz')]
AMLs <- scRNA.fis[endsWith(scRNA.fis, 'D0.dem.txt.gz')]
BMs <- c("GSM3587996_BM1.dem.txt.gz",
         "GSM3587997_BM2.dem.txt.gz",
         "GSM3587998_BM3.dem.txt.gz",
         "GSM3588000_BM4.dem.txt.gz")
scRNA.fis <- c(AMLs, BMs)



GSE116256_cli <- data.frame()
datalist <- list()
tmp=NULL
for (i in scRNA.fis) {
  print(i)
  
  infos <- unlist(strsplit(i, '_'))
  infos1 <- infos[1]
  if (startsWith(infos[2], 'AM')) {
    infos2 <- 'AML'
  } else {
    infos2 <- 'BM'
  }
  
  tmp <- read.delim(paste0('00_origin_datas/GEO/GSE116256_RAW/', i),
                    check.names = F, stringsAsFactors = F, row.names = 1)
  colnames(tmp) <- paste0(infos1, '_', str_split_fixed(colnames(tmp), '_', 2)[, 2])
  
  
  infos3 <- colnames(tmp)
  infos.dat <- data.frame(Samples = infos1,
                          Type = infos2,
                          Barcode = infos3)
  GSE116256_cli <- rbind(GSE116256_cli, infos.dat)
  tmp <- CreateSeuratObject(counts=tmp,
                            project = infos1,
                            min.cells = 3, 
                            min.features = 200)
  tmp$orig.ident <- infos1
  datalist[[i]] <- tmp
}

names(datalist)
names(datalist) <- str_split_fixed(names(datalist), '_', 3)[, 1]
names(datalist)

rownames(GSE116256_cli) <- GSE116256_cli$Barcode
colors=c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","black")


for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  #sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  datalist[[i]] <- sce
  rm(sce)
}
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
rm(datalist)
unique(sce@meta.data$orig.ident)
head(sce@meta.data)

sum(table(sce@meta.data$orig.ident))  20362
table(sce$orig.ident)
  
scRNA = subset(sce, subset=nFeature_RNA>300&nCount_RNA <= 100000 &percent.mt<10)
sum(table(scRNA@meta.data$orig.ident ))  20323

scRNA = SCTransform(scRNA, vars.to.regress="percent.mt", verbose=FALSE)
scRNA = RunPCA(scRNA, verbose=FALSE)
scRNA = RunHarmony(scRNA, group.by.vars="orig.ident", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
##############################################################################################################
saveRDS(scRNA,file = "02_scRNA/scRNA.RDS")
#scRNA=readRDS("02_scRNA/scRNA.RDS")
ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA, dims=1:30, reduction="harmony")
scRNA_umap <- RunUMAP(scRNA, dims=1:30, reduction="harmony")
#scRNA_tsne <- RunTSNE(scRNA, dims=1:30, reduction="harmony")

######################################  
p1 = DimPlot(scRNA_umap , reduction="umap", group.by="orig.ident", pt.size=1)+
  theme(legend.position="right", plot.title=element_blank())+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed"))
#+scale_color_manual(values=colors)+theme(text=element_text(family="Times"))

p2= VlnPlot(scRNA_umap, features=c("nFeature_RNA"), #cols=colors
            pt.size=0)+NoLegend()+theme(text=element_text(size=12,face = "plain"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3= VlnPlot(scRNA_umap, features=c("nCount_RNA"), #cols=colors
            pt.size=0)+NoLegend()+theme(text=element_text(size=12,face = "plain"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p4= VlnPlot(scRNA_umap, features=c("percent.mt"), #cols=colors
            pt.size=0)+NoLegend()+theme(text=element_text(size=12,face = "plain"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


resolution=0.1
mydata <- FindClusters(scRNA_umap , resolution=0.3)

u_plot_clusters=DimPlot(mydata, reduction="umap",label = T,label.size = 8,group.by="seurat_clusters", pt.size=1)+
  theme(legend.position="right", plot.title=element_blank())+
  theme(text=element_text(family="Times"))+scale_fill_lancet()+
  #scale_color_manual(values=colors)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed"))
u_plot_clusters
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "02_scRNA/All_celltype_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
markers <- read.table("02_scRNA/All_celltype_DEG.txt",header = T,check.names = F,fill=T,sep = "\t")
markers["pct.diff"]=markers$pct.1-markers$pct.2
#################################################################################

Top5 <- markers %>% group_by(cluster) %>% slice_max(n =50, order_by =avg_logFC )
Top51 <- markers %>% group_by(cluster) %>% slice_max(n =50, order_by =pct.diff )
Cellmarker <- read.table(header = T,check.names = F,fill=T,sep = "\t")
Cellmarker <- Cellmarker[,c(9,7)]
colnames(Top51)[7]="marker"
colnames(Top5)[7]="marker"
Cellmarker2 <- merge(Cellmarker,Top5,by="marker")
write.table(Cellmarker2, "02_scRNA/Cellmarker2.txt", col.names=T, row.names=F, quote=F, sep="\t")



u_plot_clusters
VlnPlot(mydata, features=c( "CD69","JUN","CD44"), pt.size=0)+NoLegend()+
  theme(axis.title.x=element_blank())

#0：Myeloid cell："FCN1","LYZ","S100A8"
#1:LSPC( leukemia stem/progenitor cells):"SPINK2","CD34"
#2：T cell："CCL5","CD2","CD3D"
#3：Myeloid cell："FCN1","LYZ","S100A8"
#4：MKI67+ progenitor cell："MKI67","ASPM","BIRC5"
#5:Erythroid cell:"HBB","HBA1","HBA2"
#6:Neutrophil:"AZU1","CTSG","ELANE"
#7：T cell："CCL5","CD2","CD3D"
#8:Plasma cell:"MZB1","JCHAIN"
#9:LSPC( leukemia stem/progenitor cells):"SPINK2","CD34"
#10:B cell:"CD79A","CD37"
#11:NK/T cell:"CD69","JUN"
#12：Myeloid cell："FCN1","LYZ","S100A8"
#13:NK/T cell:"CD69","JUN"






cell_label = c(
  "Myeloid cell",
  "LSPC",
  "T cell",
 "Myeloid cell",
  "MKI67+ progenitor cell",
  "Erythroid cell",
  "Neutrophil",
  "T cell",
 "Plasma cell",
  "LSPC",
  "B cell",
  "NK/T cell",
 "Myeloid cell",
"NK/T cell"
)
#################################################################################

names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)


u_plot=DimPlot(mydata, reduction="umap",label = T,label.size = 3.5,group.by="cell_type", pt.size=1)+
  theme(legend.position="right", plot.title=element_blank())+
  theme(text=element_text(family="Arial"))+NoLegend()+
  scale_color_manual(values=colors)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed"))


genes = unique(c(
  "FCN1","LYZ","S100A8"
 ,"SPINK2","CD34"
 ,"CCL5","CD2","CD3D"
  ,"FCN1","LYZ","S100A8"
  ,"MKI67","ASPM","BIRC5"
 ,"HBB","HBA1","HBA2"
 ,"AZU1","CTSG","ELANE"
  ,"CCL5","CD2","CD3D"
 ,"MZB1","JCHAIN"
  ,"SPINK2","CD34"
 ,"CD79A","CD37"
 ,"CD69","JUN"
  ,"FCN1","LYZ","S100A8"
  ,"CD69","JUN"
  
))

diff.marker.dotplot1= DotPlot(object = mydata, features = unique(genes),
                              dot.scale =6,
                              dot.min = 0,
                              scale =T)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('')+ylab('')+coord_flip()+scale_color_gradientn(colors=c( "grey","snow", "orange"))

cell.prop<-as.data.frame(prop.table(table(Idents(mydata), as.vector(mydata$orig.ident))))
colnames(cell.prop)<-c("Cell_type","Samples","proportion")
cell.prop$"group"=""
cell.prop$"group"[which(cell.prop$Samples%in%substr(AMLs,0,10))]="AML"
cell.prop$"group"[which(cell.prop$Samples%in%substr(BMs,0,10))]="BM"

bili=ggplot(cell.prop,aes(Samples,proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=colors)+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0,'cm'),axis.text.x = element_text(angle = 45,hjust =1))+
  guides(fill=guide_legend(title=NULL))

############AUcell#########
library(AUCell)
library(GSEABase)
countexp<-mydata@assays$RNA@data

cells_rankings <- AUCell_buildRankings(as.matrix(countexp)) #rank
pcdgenesets=read.table("01_DEGs//diffSig_mRNA_log2FC1.5_P.Vale0.05.txt",header = T, fill = TRUE)
geneSets <- list('Ultrasound'=unique(rownames(pcdgenesets))) #signature read
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) #calc

HP_AUCell <- as.numeric(getAUC(cells_AUC)['Ultrasound', ])
mydata$HP_AUCell  <- HP_AUCell
saveRDS(mydata,"02_scRNA/mydata.rds")
#mydata <- readRDS("02_scRNA/mydata.rds")

sc_df=data.frame(mydata@meta.data, mydata@reductions$umap@cell.embeddings)
head(sc_df)


hp.boxplot=ggplot(sc_df) +
  geom_boxplot(aes(x = cell_type, y = HP_AUCell,fill=cell_type))+
  #scale_fill_manual(values =colors)+
  xlab('')+ylab('Metabolic AUCell')+ theme_classic()+
  theme(text = element_text(family = 'Times',size = 12),legend.position = 'none',
        axis.text.x = element_text(color = "black", size = 12,angle = 30,hjust = 1))
hp.boxplot
sc_df$"group"=""
sc_df$"group"[which(sc_df$orig.ident%in%substr(AMLs,0,10))]="AML"
sc_df$"group"[which(sc_df$orig.ident%in%substr(BMs,0,10))]="BM"

hp.boxplot1=ggplot(sc_df,aes(x = cell_type , y = HP_AUCell,fill=group)) +
  geom_boxplot()+
  scale_fill_manual(values =c("#E64B35","#4DBBD5"))+
  xlab('')+ylab('AUCell Score')+ theme_classic()+
  theme(text = element_text(family = 'Times',size = 12),
        axis.text.x = element_text(color = "black", size = 12,angle = 30,hjust = 1))+stat_compare_means(label = "p.signif", method = "t.test")+theme (legend.title=element_blank ())

ggsave('02_scRNA/hp.boxplot.pdf',hp.boxplot,height = 6,width = 6)

FigureS1=mg_merge_plot(p2,p3,p4,u_plot_clusters,ncol=2,nrow=2,labels = LETTERS[1:4])

ggsave('PDFs/Figure S1.pdf',FigureS1,height = 10,width = 12)
ggsave('PDFs/Figure S1.jpg',FigureS1,height = 10,width = 12)

Figure1=mg_merge_plot(u_plot,diff.marker.dotplot1,bili,hp.boxplot1,nrow=2,ncol=2,labels = c('A','B',"C",'D'))
ggsave('PDFs/Fig2.pdf',Figure1,height = 15,width = 16)
ggsave('PDFs/Fig2.jpg',Figure1,height = 15,width = 16)

save.image("project_scRNA.Rdata")




table(mydata$cell_type)
immuce_cells=subset(mydata,cell_type%in%c('Myeloid cell'))
meta_sample=immuce_cells@meta.data


immuce_cells=SCTransform(immuce_cells, vars.to.regress="percent.mt", verbose=FALSE)

#immuce_cells <- NormalizeData(immuce_cells)
#immuce_cells <- FindVariableFeatures(immuce_cells, selection.method = "vst", nfeatures = 2000)
#immuce_cells <- ScaleData(immuce_cells)
immuce_cells <- RunPCA(immuce_cells, verbose=FALSE, assay="SCT")
#immuce_cells = RunHarmony(immuce_cells, group.by.vars="Samples", max.iter.harmony=50, lambda=0.5, assay.use="SCT")

ElbowPlot(immuce_cells,ndims = 50)+theme(text = element_text(family = 'Times',size = 12))
dev.off()
immuce_cells <- FindNeighbors(immuce_cells, dims=1:30, reduction="harmony")
#immuce_cells <- RunUMAP(immuce_cells, dims=1:30, reduction="harmony")
immuce_cells <- RunTSNE(immuce_cells, dims=1:30, reduction="harmony")

#DimPlot(immuce_cells, reduction="umap", group.by="orig.ident", pt.size=1)+theme(legend.position="right", plot.title=element_blank())
DimPlot(immuce_cells, reduction="tsne", group.by="orig.ident", pt.size=1)+theme(legend.position="right", plot.title=element_blank())


immuce_cells <- FindClusters(immuce_cells, resolution=0.1)
#Umap_Plot=UMAPPlot(immuce_cells, pt.size=1, label=T)
tsne_state=TSNEPlot(immuce_cells, pt.size=1, label=T)
markers <- FindAllMarkers(immuce_cells, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.5)
write.table(markers, "03_Myeloid/Myeloid_All_celltype_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
markers <- read.table("03_Myeloid/Myeloid_All_celltype_DEG.txt",header = T,check.names = F,fill=T,sep = "\t")
markers["pct.diff"]=markers$pct.1-markers$pct.2
#################################################################################

Top51 <- markers %>% group_by(cluster) %>% slice_max(n =5, order_by =pct.diff )
Cellmarker <- read.table(header = T,check.names = F,fill=T,sep = "\t")
Cellmarker <- Cellmarker[,c(9,7)]
colnames(Top51)[7]="marker"
Cellmarker2 <- merge(Cellmarker,Top51,by="marker")
write.table(Cellmarker2, "03_Myeloid/Myeloid_Cellmarker2.txt", col.names=T, row.names=F, quote=F, sep="\t")



VlnPlot(immuce_cells, features=c("ACP5","GPNMB","SPP1"), pt.size=0)+NoLegend()+theme(axis.title.x=element_blank())


#0:M2 macrophage:"CD209","FOLR2","SIGLEC1"
#1:Monocyte:"AQP9","S100A8","VCAN"
#2:Macrophage:"ACP5","GPNMB","SPP1"
#3:Dendritic cell:"CD1C","CD1E","CLEC10A"

Idents(immuce_cells)="HP_AUCell"
TSNE_Plot=TSNEPlot(immuce_cells, pt.size=1, label=F, label.size=5)
immuce_cells$AUCell_group=ifelse( immuce_cells$HP_AUCell > median(immuce_cells$HP_AUCell),'High','Low')

LRG_group_umap=DimPlot(immuce_cells,group.by='AUCell_group',reduction="tsne",label = F,pt.size = 0.2)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12))
LRG_group_umap
ggsave('results/02.scRNA/LRG_group_umap.pdf',LRG_group_umap,height = 5,width = 5.5)


gs=c(
  "CD209","FOLR2","SIGLEC1"
  ,"AQP9","S100A8","VCAN"
  ,"ACP5","GPNMB","SPP1"
  ,"CD1C","CD1E","CLEC10A"
)
diff.marker.dotplot_subcell= DotPlot(object = immuce_cells, features = unique(as.vector(Top51$marker)),
                                     dot.scale =6,
                                     #dot.min = 0,
                                     scale =T)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('')+ylab('')+coord_flip()+scale_color_gradientn(colors=c( "grey","snow", "#0077c8"))

markers2 <- FindAllMarkers(immuce_cells, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.5)
write.table(markers2, "02_Myeloid/Myeloid_All_celltype_DEG2.txt", col.names=T, row.names=T, quote=F, sep="\t")
#VlnPlot_Macrophage=VlnPlot(immuce_cells, features=gs, pt.size=0, cols=c("#DF8F44FF","#DF8F44FF","#DF8F44FF"))+NoLegend()+theme(axis.title.x=element_blank())
saveRDS(immuce_cells,file = '02_Myeloid/immuce_cells.rds')
