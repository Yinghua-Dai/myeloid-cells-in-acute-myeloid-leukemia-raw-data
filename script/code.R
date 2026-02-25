rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(devtools))
suppressMessages(library(customLayout))
suppressMessages(library(stringr))
suppressMessages(library(ConsensusClusterPlus))
suppressMessages(library(tidydr))
suppressMessages(library(openxlsx))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))
suppressMessages(library(pheatmap))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(GSVA))
suppressMessages(library(GSEABase))
suppressMessages(library(fgsea))
suppressMessages(library(corrplot))
suppressMessages(library(colorspace))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(maftools))
suppressMessages(library(vegan))
suppressMessages(library(forcats))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(rstatix))
suppressMessages(library(ggstatsplot))
suppressMessages(library(ggcor))
suppressMessages(library(ggstance))
suppressMessages(library(tidyverse))
suppressMessages(library(GOplot))
suppressMessages(library(caret))
suppressMessages(library(writexl))
suppressMessages(library(rcartocolor))
suppressMessages(library(ggcorrplot))
suppressMessages(library(psych))
suppressMessages(library(clusterProfiler))
suppressMessages(library(dplyr))
suppressMessages(library(cols4all))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(scales))
suppressMessages(library(oncoPredict))
suppressMessages(library(gghalves))
suppressMessages(library(cowplot))
suppressMessages(library(IOBR))
suppressMessages(library(estimate))
suppressMessages(library(UpSetR))
suppressMessages(library(ggbiplot))
suppressMessages(library(ggsci))
suppressMessages(library(WGCNA))
suppressMessages(library(circlize))
suppressMessages(library(rJava))
suppressMessages(library(xlsxjars))
suppressMessages(library(xlsx))
suppressMessages(library(glmnet))
suppressMessages(library(tidyr))
suppressMessages(library(pROC))
suppressMessages(library(ROCR))
suppressMessages(library(patchwork))


##############DEG##########
GSE10212_exp=readMatrix('00_origin_datas/Preprocessed/GSE10212_exp.txt')
GSE10212_cli=readMatrix('00_origin_datas/Preprocessed/GSE10212_cli.txt')

table(GSE10212_cli$group)

identical(rownames(GSE10212_cli),colnames(GSE10212_exp))
geo.limma=mg_limma_DEG(exp = GSE10212_exp,
                       group = GSE10212_cli$group,
                       ulab = 'Case',dlab = 'Control')
geo.limma$Summary

geo.degs=geo.limma$DEG[which(geo.limma$DEG$P.Value<0.05 & abs(geo.limma$DEG$logFC) > log2(1.5)),]
head(geo.degs)
geo.degs$group=ifelse(geo.degs$logFC>0,'up','down')
geo.degs$gene=rownames(geo.degs)
table(geo.degs$group)
dim(geo.degs)
write.table(geo.degs,file="01_DEGs//diffSig_mRNA_log2FC1.5_P.Vale0.05.txt",sep="\t",quote=F)


fig2a=my_volcano (geo.limma,p_cutoff = 0.05,fc_cutoff = log2(1.5),col = c("#FDB462","#8BACD1",'grey'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        text = element_text(color = "black",family = 'Times',size = 14),
        axis.text = element_text(color = "black",family = 'Times',size = 14),
        legend.position = 'top')
fig2a


GEO_expression=GSE10212_exp[rownames(geo.degs),]
GEO_expression1 <- apply(GEO_expression, 1, scale)
rownames(GEO_expression1) <- colnames(GEO_expression)
GEO_expression2 <- as.data.frame(t(GEO_expression1))
GEO_expression2 <-GEO_expression2[,rownames(GSE10212_cli)]
identical(rownames(GSE10212_cli), colnames(GEO_expression2))

Subtype.color=c("#FDB462","#8BACD1")
names(Subtype.color)=c('Case','Control')


GSE10212_cli1=as.matrix(GSE10212_cli[,2])
rownames(GSE10212_cli1)=rownames(GSE10212_cli)
colnames(GSE10212_cli1)=colnames(GSE10212_cli)[2]
GSE10212_cli1=as.data.frame(GSE10212_cli1)
column_ha=HeatmapAnnotation(df = GSE10212_cli1
                            , na_col = "grey"
                            , annotation_height = unit(0.01, "mm")
                            , gap = unit(1, 'mm')
                            ,col = list(group=Subtype.color)
                                        
                            )


identical(rownames(GSE10212_cli1), colnames(GEO_expression2))

heatmap_plot=Heatmap(as.matrix(GEO_expression2),
                     col = colorRamp2(c(-3, -1.5, 0, 1.5, 3),c("#226ED1", "#7AA8E2", "white", "#F09090","#E01010")), 
                     #col = circlize::colorRamp2(c(-3, 0, 3), c('navy', 'white', 'red')),
                     name = "Expression",top_annotation = column_ha,
                     show_row_names = F,show_column_names = F,
                     clustering_method_rows = "complete",
                     row_names_gp = gpar(fontsize = 10),
                     column_names_gp = gpar(fontsize = 10),
                     column_names_rot = 45,
                     cluster_columns = F,
                     cluster_rows = T,
                     show_row_dend = F,
                     #  width = ncol(exp)*unit(cell_size, "mm"),
                     use_raster = F,
                     ##`use_raster` is automatically set to TRUE for a matrix with more than 2000 rows.
                     #  row_names_max_width = row_name_width,
                     heatmap_legend_param = list(direction = "vertical",
                                                 legend_width = unit(3.05,"cm"),
                                                 legend_height = unit(2.8, "cm"),
                                                 title_position = "lefttop-rot"))
pdf(file="01_DEGs/heatmap_plot.pdf",width=7,height=6)
heatmap_plot
dev.off()

enrichment=mg_clusterProfiler(as.vector(geo.degs$gene))
kegg_dot=enrichplot::dotplot(enrichment$KEGG)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))+ scale_x_continuous(limits = c(0.05, 0.11))

bp_dot=enrichplot::dotplot(enrichment$GO_BP)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))+ scale_x_continuous(limits = c(0.04, 0.07))
Fig2=mg_merge_plot(fig2a, grid.grabExpr(draw(heatmap_plot)),kegg_dot,bp_dot,nrow = 2,ncol = 2,labels=LETTERS[c(1:4)])
ggsave('PDFs/Fig1.pdf',Fig2,height = 11,width = 12.5)
ggsave('PDFs/Fig1.jpg',Fig2,height = 11,width = 12.5)


s.deg <- readRDS("02_scRNA/deg.RDS")
markers.deg <- readRDS("02_scRNA/markers.deg.RDS")
hub=intersect(rownames(markers.deg),s.deg$gene)
length(hub)
#write.table(hub,file="03_DEGs//hub.txt",sep="\t",quote=F,row.names = F,col.names = F)
############lasso##############

tcga.exp=readMatrix('00_origin_datas/Preprocessed/tcga.exp.txt')
tcga.cli=readMatrix('00_origin_datas/Preprocessed/tcga.cli.txt')
identical(as.vector(tcga.cli$sampleID),as.vector(colnames(tcga.exp)))

tcga.t.exp_use=as.data.frame(tcga.exp)
tcga.subtype.cli=tcga.cli
identical(colnames(tcga.t.exp_use),as.vector(tcga.subtype.cli$sampleID))
colnames(tcga.subtype.cli)[1]="Samples"
rownames(tcga.subtype.cli)=tcga.subtype.cli$Samples


tcga.cox=cox_batch(t(scale(t(tcga.t.exp_use[hub,])))
                   ,time =  tcga.subtype.cli$OS.time/365
                   ,event =tcga.subtype.cli$OS)
dim(tcga.cox)

table(tcga.cox$p.value<0.05)
table(tcga.cox$p.value<0.01)
table(tcga.cox$p.value<0.001)
writeMatrix(tcga.cox,outpath = '03_Lasso/tcga.cox.txt')



p.cutoff=0.01
tcga.cox_use=tcga.cox
tcga.cox_use$coef=log(tcga.cox_use$HR)
tcga.cox_use$Gene=rownames(tcga.cox_use)
tcga.cox_use$type=rep('None',nrow(tcga.cox_use))
tcga.cox_use$type[which(tcga.cox_use$p.value<p.cutoff & tcga.cox_use$coef>0)]='Risk'
tcga.cox_use$type[which(tcga.cox_use$p.value<p.cutoff & tcga.cox_use$coef<0)]='Protective'
table(tcga.cox_use$type)

######### lasso
tcga.gene.sig=rownames(tcga.cox)[which(tcga.cox$p.value<p.cutoff)]
length(tcga.gene.sig)

table(tcga.cox_use$type)


tcga.cox_forVis=tcga.cox_use
tcga.cox_forVis=tcga.cox_forVis[which(tcga.cox_forVis$type %in% c('Risk','Protective')),]
tcga.cox_forVis$p.value=-log10(tcga.cox_forVis$p.value)
range(tcga.cox_forVis$p.value)


#################### LASSO
table(tcga.cox$p.value<0.05)
table(tcga.cox$p.value<0.01)
table(tcga.cox$p.value<0.001)

tcga.exp.sig=tcga.t.exp_use[tcga.gene.sig,]
tcga.exp.sig=t(tcga.exp.sig)
dim(tcga.exp.sig)


dim(tcga.exp.sig)
options(ggrepel.max.hnscerlaps = Inf)
tcga.subtype.cli$Samples=as.vector(tcga.subtype.cli$Samples)
identical(rownames(tcga.exp.sig),tcga.subtype.cli$Samples)
tcga.lasso.res=mg_lasso_cox_use(tcga.exp.sig
                                , time = tcga.subtype.cli$OS.time/365
                                , event = tcga.subtype.cli$OS
                                , nfolds = 5
                                , lambda.min = T
                                , figLabels=c('B','C'))
tcga.lasso.res$Genes

tcga.lasso.res$lambda

tcga.lasso.res$plot

tcga.exp.for.cox=tcga.t.exp_use[match(tcga.lasso.res$Genes,row.names(tcga.t.exp_use)),]
dim(tcga.exp.for.cox)
identical(colnames(tcga.exp.for.cox),tcga.subtype.cli$Samples)

lst.modl=createCoxModel_use((t(tcga.exp.for.cox))
                            , time = tcga.subtype.cli$OS.time/365
                            , event = tcga.subtype.cli$OS
                            , isStep =T)
lst.modl$Cox
lst.modl$Genes
lst.modl$fmla

lst.modl.Coef=lst.modl$Coef
names(lst.modl.Coef)=lst.modl$Genes
lst.modl.Coef

tcga.risk.score=lst.modl$Score
#tcga.risk.score=scale(tcga.risk.score)[,1]
tcga.risk.score=mosaic::zscore(tcga.risk.score)

range(tcga.risk.score)

lst.modl$Coef

gene.coef=data.frame(Gene=lst.modl$Genes,Coef=lst.modl$Coef)
gene.coef$Type=ifelse(lst.modl$Coef>0,'Risk','Protective')
gene.coef$Type=factor(gene.coef$Type,levels=c('Risk','Protective'))
table(gene.coef$Type)

fig3c=gene.coef %>% 
  ggplot(aes(reorder(Gene, Coef), Coef)) +
  geom_col(aes(fill = Type)) +
  geom_text(aes(label=round(Coef,digits = 3)),color="black",hjust = "left")+
  
  coord_flip() +
  scale_fill_manual(values=pal_nejm(alpha = 0.9)(8)[c(1,2)]) +
  coord_flip() +
  labs(x = "") +
  labs(y = "Lasso Cox coefficient") +
  theme_classic()+theme(legend.position = c(0,1))
# theme(axis.text.y = element_text(angle = 0, hjust = 1),legend.position="top")


fig3AB=mg_plot_lasso_use(fit = tcga.lasso.res$Mode1
                         , cv_fit = tcga.lasso.res$Model2
                         , show_text = F
                         , figLabels = c('A', 'B'))
fig3AB
fig3abc=mg_merge_plot(fig3AB,fig3c,nrow = 1,ncol = 2,widths = c(2,1))
#savePDF('PDFs/Fig7AB.pdf',fig7A,height = 4,width = 9)
#savePDF('PDFs/fig3abc.pdf',fig3abc,height = 5,width = 15)


tcga.exp.forCox<- cbind(time=tcga.subtype.cli$OS.time/365,
                        status=tcga.subtype.cli$OS,
                        t(tcga.t.exp_use)[rownames(tcga.subtype.cli), lst.modl$Genes])



dim(tcga.exp.forCox)

fmla <- as.formula(paste0("Surv(time, status) ~",paste0(lst.modl$Genes,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga.exp.forCox))
fig3d=survminer::ggforest(cox,data=tcga.exp.forCox,noDigits = 3)
fig3abd=mg_merge_plot(fig3AB,fig3d,nrow = 1,ncol =2,widths = c(2,1))

#savePDF('PDFs/fig3abc.pdf',fig3abcd,height = 5,width = 16)


############### TCGA
#tcga.cutoff <- survminer::surv_cutpoint(data.frame(time=tcga.subtype.cli$OS.time/365, event=tcga.subtype.cli$OS, risk=tcga.risk.score),time = "time", event = "event",variables = c("risk"))
#tcga.cutoff=tcga.cutoff$cutpoint$cutpoint
#tcga.cutoff=median(tcga.risk.score)
tcga.cutoff=0
identical(colnames(tcga.exp.for.cox),tcga.subtype.cli$Samples)
risk.group.color=c("#d15034","#264a5f")
names(risk.group.color)=c('High','Low')
tcga.roc=plotCoxModel_Batch_use(riskScore = tcga.risk.score
                                ,dat = t(tcga.exp.for.cox[match(lst.modl$Genes, row.names(tcga.exp.for.cox)), ])
                                , time = tcga.subtype.cli$OS.time/365
                                , event = tcga.subtype.cli$OS
                                , cutoff = tcga.cutoff
                                , labs = c('High','Low')
                                , title = 'RiskType'
                                , hetColor = c('#BF1D2D', 'white', '#293890')
                                
                                , pal = risk.group.color
                                , mks = c(1:5))
tcga.roc1=tcga.roc[[2]]
#pdf('PDFs/fig3e2.pdf',height = 6,width = 6)
tcga.roc1

dev.off()

tcga.group=ifelse(tcga.risk.score>tcga.cutoff,'High','Low')
tcga.group=data.frame(tcga.group)
colnames(tcga.group)='group'
table(tcga.group$group)
write.table(cbind(tcga.risk.score,tcga.group),file = '03_Lasso/tcga.group.txt',sep='\t',quote = F)


#GSE10358
GSE10358.t.exp=readMatrix('00_origin_datas/Preprocessed/GSE10358.t.exp.txt')
GSE10358.t.cli=readMatrix('00_origin_datas/Preprocessed/GSE10358.t.cli.txt')
GSE10358.t.cli$"geo_accession"=rownames(GSE10358.t.cli)
identical(colnames(GSE10358.t.exp),as.vector(GSE10358.t.cli$geo_accession))

match(lst.modl$Genes,row.names(GSE10358.t.exp))
length(lst.modl$Genes)
GSE10358.t.cli.os=GSE10358.t.cli
GSE10358.t.cli.os$Samples =as.vector(GSE10358.t.cli.os$Samples )
identical(GSE10358.t.cli.os$geo_accession , colnames(GSE10358.t.exp)) 
GSE10358.model.dat=GSE10358.t.exp[match(lst.modl$Genes,row.names(GSE10358.t.exp)),]

GSE10358.risk.score=predictScoreByCoxModel(coxModel = lst.modl
                                           ,(t(GSE10358.model.dat)))
#GSE10358.risk.score=scale(GSE10358.risk.score)
GSE10358.risk.score=mosaic::zscore(GSE10358.risk.score)

lst.modl$fmla
#lst.vd.mod1$fmla
#GSE10358.cutoff <- survminer::surv_cutpoint(data.frame(time=as.numeric(GSE10358.t.cli.os$OS.time),
#                                                       event=GSE10358.t.cli.os$OS, risk=GSE10358.risk.score),                                                        time =  "time", event = "event", variables = c("risk"))
#GSE10358.cutoff=GSE10358.cutoff$cutpoint$cutpoint
GSE10358.cutoff=0


GSE10358.roc=plotCoxModel_Batch_use(riskScore = GSE10358.risk.score
                                    , dat = t(GSE10358.t.exp[intersect(lst.modl$Genes, row.names(GSE10358.t.exp)),])
                                    , time = as.numeric(GSE10358.t.cli.os$OS.time/365) 
                                    , event = as.numeric(GSE10358.t.cli.os$OS)
                                    , cutoff = GSE10358.cutoff
                                    , labs = c('High','Low')
                                    , title = 'RiskType'
                                    , hetColor = c('#BF1D2D', 'white', '#293890')
                                    , pal = risk.group.color
                                    , mks = c(1:5))
GSE10358.roc1=GSE10358.roc[[2]]
#pdf('PDFs/fig3g2.pdf',height = 6,width = 6)

GSE10358.roc1
dev.off()
GSE10358.group=ifelse(GSE10358.risk.score>GSE10358.cutoff,'High','Low')
GSE10358.group=data.frame(GSE10358.group)
colnames(GSE10358.group)='group'
table(GSE10358.group)

write.table(cbind(GSE10358.risk.score,GSE10358.group),file = '03_Lasso//GSE10358.group.txt',sep='\t',quote = F)






#GSE37642
GSE37642.t.exp=readMatrix('00_origin_datas/Preprocessed/GSE37642.t.exp.txt')
GSE37642.t.cli=readMatrix('00_origin_datas/Preprocessed/GSE37642.t.cli.txt')
GSE37642.t.cli$"geo_accession"=rownames(GSE37642.t.cli)
identical(colnames(GSE37642.t.exp),as.vector(GSE37642.t.cli$geo_accession))

match(lst.modl$Genes,row.names(GSE37642.t.exp))
length(lst.modl$Genes)
GSE37642.t.cli.os=GSE37642.t.cli
GSE37642.t.cli.os$Samples =as.vector(GSE37642.t.cli.os$Samples )
identical(GSE37642.t.cli.os$geo_accession , colnames(GSE37642.t.exp))
GSE37642.model.dat=GSE37642.t.exp[match(lst.modl$Genes,row.names(GSE37642.t.exp)),]

GSE37642.risk.score=predictScoreByCoxModel(coxModel = lst.modl
                                           ,(t(GSE37642.model.dat)))
#GSE37642.risk.score=scale(GSE37642.risk.score)
GSE37642.risk.score=mosaic::zscore(GSE37642.risk.score)

lst.modl$fmla
#lst.vd.mod1$fmla
#GSE37642.cutoff <- survminer::surv_cutpoint(data.frame(time=as.numeric(GSE37642.t.cli.os$OS.time),
#                                                       event=GSE37642.t.cli.os$OS, risk=GSE37642.risk.score),                                                        time =  "time", event = "event", variables = c("risk"))
#GSE37642.cutoff=GSE37642.cutoff$cutpoint$cutpoint
GSE37642.cutoff=0


GSE37642.roc=plotCoxModel_Batch_use(riskScore = GSE37642.risk.score
                                    , dat = t(GSE37642.t.exp[intersect(lst.modl$Genes, row.names(GSE37642.t.exp)),])
                                    , time = as.numeric(GSE37642.t.cli.os$OS.time/365) 
                                    , event = as.numeric(GSE37642.t.cli.os$OS)
                                    , cutoff = GSE37642.cutoff
                                    , labs = c('High','Low')
                                    , title = 'RiskType'
                                    , hetColor = c('#BF1D2D', 'white', '#293890')
                                    , pal = risk.group.color
                                    , mks = c(1:5))
GSE37642.roc1=GSE37642.roc[[2]]
#pdf('PDFs/fig3g2.pdf',height = 6,width = 6)

GSE37642.roc1
dev.off()
GSE37642.group=ifelse(GSE37642.risk.score>GSE37642.cutoff,'High','Low')
GSE37642.group=data.frame(GSE37642.group)
colnames(GSE37642.group)='group'
table(GSE37642.group)

write.table(cbind(GSE37642.risk.score,GSE37642.group),file = '03_Lasso//GSE37642.group.txt',sep='\t',quote = F)





figROC=mg_merge_plot(tcga.roc1,GSE10358.roc1,GSE37642.roc1,nrow = 1,ncol = 3)

Fig3=mg_merge_plot(fig3abd,figROC,nrow = 2,ncol = 1,heights = c(1,2))
ggsave("PDFs/Fig4.pdf",Fig3, width = 12, height = 10)
identical(rownames(as.data.frame(tcga.risk.score)),rownames(tcga.subtype.cli))
tcga.risktype.cli=data.frame(tcga.subtype.cli,Riskscore=tcga.risk.score)
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>tcga.cutoff,'High','Low')
table(tcga.risktype.cli$Risktype)
write.table(tcga.risktype.cli,file = '03_Lasso/tcga.risktype.cli.txt',sep='\t',quote = F)





tcga.risktype.cli=read.table('03_Lasso/tcga.risktype.cli.txt',header = T,check.names = F,fill=T,sep = "\t")
tcga.exp=readMatrix('00_origin_datas/Preprocessed/tcga.exp.txt')
tcga.cli=readMatrix('00_origin_datas/Preprocessed/tcga.cli.txt')
identical(as.vector(tcga.risktype.cli$Samples),colnames(tcga.exp))

risk.group.color1=c("#d15034","#264a5f")
names(risk.group.color1)=c('High','Low')

library(estimate)
#### ESTIMATE
#tcga.exp.estimate<-deconvo_estimate(eset=tcga.exp)
#save(tcga.exp.estimate,file='05_imm/tcga.exp.estimate.RData')
load('05_imm/tcga.exp.estimate.RData')
tcga.exp.estimate=get.IOBR.immu.format(tcga.exp.estimate)

############ MCP-counter 
#tcga.exp.mcp<-deconvo_mcpcounter(eset=as.matrix(tcga.exp))
#save(tcga.exp.mcp,file='05_imm/tcga.exp.mcp.RData')
load('05_imm/tcga.exp.mcp.RData')
tcga.exp.mcp=get.IOBR.immu.format(tcga.exp.mcp)

### CIBERSORT
#tcga.exp.cibersort<-deconvo_cibersort(eset=tcga.exp,arrays=F)
#save(tcga.exp.cibersort,file='05_imm/tcga.exp.cibersort.RData')
load('05_imm/tcga.exp.cibersort.RData')
tcga.exp.cibersort=get.IOBR.immu.format(tcga.exp.cibersort)

#######sssGSEA#######
#geo.immu.ssgsea=immu_ssgsea(exp = tcga.exp)
##save(geo.immu.ssgsea,file='05_imm/geo.immu.ssgsea.RData')
load('05_imm/geo.immu.ssgsea.RData')



#
tcga.t.estimate=tcga.exp.estimate[rownames(tcga.risktype.cli),1:3]
tcga.t.mcp=tcga.exp.mcp[rownames(tcga.risktype.cli),]
tcga.t.cibersort=tcga.exp.cibersort[rownames(tcga.risktype.cli),1:22]
tcga.t.ssGSEA28=as.data.frame(geo.immu.ssgsea[rownames(tcga.risktype.cli),])

fig5a=get_PlotMutiBoxplot(tcga.t.estimate,tcga.risktype.cli
                          ,group_cols = risk.group.color1
                          ,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')

fig5a=groupViolin(tcga.t.estimate,
                  tcga.risktype.cli$Risktype,
                  ylab = 'Score',
                  group_col=risk.group.color1)

fig5a





fig5d=get_PlotMutiBoxplot(tcga.exp.mcp,tcga.risktype.cli
                          ,group_cols = risk.group.color1
                          ,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)

fig5d=groupViolin(tcga.exp.mcp,
                  tcga.risktype.cli$Risktype,
                  ylab = 'Score',
                  group_col=risk.group.color1)
fig5d





fig5e=get_PlotMutiBoxplot(tcga.t.cibersort,tcga.risktype.cli
                          ,group_cols = risk.group.color
                          ,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')
fig5e


fig5f=get_PlotMutiBoxplot(tcga.t.ssGSEA28,tcga.risktype.cli   ,group_cols = risk.group.color1,ylab = 'Score' ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')

fig5f=groupViolin(tcga.t.ssGSEA28,
                  tcga.risktype.cli$Risktype,
                  ylab = 'Score',
                  group_col=risk.group.color1)

fig5f

fig5ab=mg_merge_plot(fig5a,fig5d,nrow = 1,ncol = 2
                     ,labels = LETTERS[2:3])



fig5=mg_merge_plot(fig5f,fig5ab,nrow = 2,ncol = 1 ,labels =c("A","") )


ggsave('PDFs/Fig5a.pdf',fig5,height = 12,width = 15)
ggsave('PDFs/Fig5a.jpg',fig5,height = 12,width = 15)



library('oncoPredict')
tcga.risktype.cli=read.table('03_Lasso/tcga.risktype.cli.txt',header = T,check.names = F,fill=T,sep = "\t")
tcga.exp=readMatrix('00_origin_datas/Preprocessed/tcga.exp.txt')
tcga.cli=readMatrix('00_origin_datas/Preprocessed/tcga.cli.txt')
identical(as.vector(tcga.risktype.cli$Samples),colnames(tcga.exp))
##ctrl + shift + C
drug_exp=as.matrix(tcga.exp)
# GDSC2_Expr = readRDS(file=file.path(dir,'Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
# GDSC2_Res = readRDS(file = file.path(dir,"Training Data/GDSC2_Res.rds"))
# GDSC2_Res <- exp(GDSC2_Res)
# calcPhenotype(trainingExprData = as.matrix(GDSC2_Expr),
#              trainingPtype = as.matrix(GDSC2_Res),
#              testExprData = as.matrix(drug_exp),
#              batchCorrect = 'eb',  #   "eb" for ComBat
#              powerTransformPhenotype = TRUE,
#              removeLowVaryingGenes = 0.2,
#              minNumSamples = 10,
#             printOutput = TRUE,
#             removeLowVaringGenesFrom = 'rawData' )

tcga.drug.ic50=read.csv('05_imm/calcPhenotype_Output/DrugPredictions.csv',row.names = 1)
dim(tcga.drug.ic50)

head(tcga.drug.ic50)
dim(tcga.drug.ic50)
all(rownames(tcga.drug.ic50)==rownames(tcga.risktype.cli))

cr=psych::corr.test(y=data.frame(tcga.drug.ic50[rownames(tcga.risktype.cli),],check.names = F),
                    x=data.frame(t(tcga.exp)[rownames(tcga.risktype.cli),lst.modl$Genes],RiskScore=tcga.risktype.cli$Riskscore))

df_cor=cr$r
df_pval=cr$p
#df_cor=round(df_cor,2)


inds=which(abs(df_cor[nrow(df_cor),])>0.3 & df_pval[nrow(df_cor),]<0.05)
length(inds)
df_cor=df_cor[,inds]
df_pval=df_pval[,inds]
### Pivot data from wide to long
library(tidyverse)
g = pivot_longer(data=rownames_to_column(as.data.frame(df_cor),var = "from"),
                 cols = 2:(ncol(df_cor)+1), ## Columns to pivot into longer format
                 names_to = "to",
                 values_to = "cor")
gp = pivot_longer(data=rownames_to_column(as.data.frame(df_pval)),
                  cols = 2:(ncol(df_pval)+1),
                  names_to = "gene",
                  values_to = "p")
all(g$from==gp$rowname & g$to==gp$gene)
g$p.adj = gp$p

###################
df=g
df <- df %>%
  mutate(col = cut(cor, breaks = c(-1, 0, 1),
                   labels = c("negative", "positive")),
         p.signif = cut(p.adj, breaks = c(0,0.0001, 0.001, 0.01, 0.05,1),
                        labels = c("****", "**", "**","*",""),
                        right = FALSE, include.lowest = TRUE))
df=data.frame(df)
head(df)
mode(df$p.adj)="numeric"
mode(df$cor)="numeric"

# df=df[which(df$p.adj<0.05 & abs(df$cor)>0.4),]
length(unique(df$to))
head(df)
writeMatrix(df,outpath = '05_imm/tcga.drug.model.cor.txt')

corr.mat=pivot_wider(df[,c(1,2,3)],names_from  ="from",values_from ='cor')
p.mat=pivot_wider(df[,c(1,2,4)],names_from  ="from",values_from ='p.adj')
corr.mat=data.frame(corr.mat)
p.mat=data.frame(p.mat)
corr.mat[1:4,1:5]
rownames(corr.mat)=corr.mat$to
corr.mat=corr.mat[-1]
rownames(p.mat)=p.mat$to
p.mat=p.mat[-1]
dim(corr.mat)

corr.mat[1:4,1:5]
####################
library(corrplot)
pdf('05_imm/drug_plot.pdf',height = 10,width = 10,onefile = F)
drug_plot=corrplot(corr = as.matrix(corr.mat),
         p.mat = as.matrix(p.mat),
         mar = c(0,0,0,0),
         col=colorRampPalette(c('#3969AC', 'white','#E73F74'))(50),
         tl.srt = 90,tl.cex = 1,tl.col = 'black',tl.offset = 0.5,
         cl.pos = c("b","r","n")[1],cl.align.text = 'l',cl.length = 5,
         cl.ratio = 0.1,cl.cex = 0.8,
         addgrid.col = 'white',
         method = "pie",
         insig = 'label_sig',
         sig.level=c(0.001,0.01,0.05),
         pch.cex=1,is.corr=T,xpd=T)
dev.off()
rownames(corr.mat)


tcga_tmb <- mg_getTCGATMBByCode('LAML')
tcga_tmb <- as.data.frame(tcga_tmb)
head(tcga_tmb)
table(substr(tcga.risktype.cli$Samples,14,15))
tcga_tmb$Sample <- paste0(tcga_tmb$Sample, '-03')
tcga.risktype.cli$Samples=substr(tcga.risktype.cli$Samples,0,15)
#colnames(tcga.risktype.cli)[1]="Sample"
tcga_tmb1 <- merge(tcga_tmb, tcga.risktype.cli[, c("Samples", "Risktype",'Riskscore',"OS.time","OS")],
                   by.x = 'Sample', by.y = 'Samples')
head(tcga_tmb1)
tcga_tmb1$logTMB <- log2(tcga_tmb1$TMB + 1)


# tcga_tmb_plot <- wb_beeswarm_plot(tcga_tmb1[, c("Risktype", "logTMB")],
#                                   ylab = 'log2(Tumor mutation burden)',
#                                   show_compare = T,
#                                   xlab = "Risktype",
#                                   title = 'TCGA',
#                                   col = risk.group.color)

# tcga_tmb_plot <- mg_violin(tcga_tmb1[, c("Risktype", "logTMB")]
#                            ,melt = T
#                            ,xlab = 'Risktype'
#                            ,legend.pos = 'tl'
#                            ,ylab = 'log2(Tumor mutation burden)')
tcga_tmb_plot <-ggboxplot(tcga_tmb1, x = "Risktype", y = "logTMB",
                fill = "Risktype", palette = c("#d15034", "#264a5f")) +
                stat_compare_means(method = "t.test") +
                labs(x = "Risktype", y = "logTMB")+geom_point(position = position_jitter(width = 0.2))
tcga_tmb_plot
ggsave(plot = tcga_tmb_plot,
       filename = 'PDFs/tcga_tmb_plot.pdf',
       width = 5, height = 5)
#######TMB##########
rownames(tcga_tmb1)=tcga_tmb1$Sample
head(tcga_tmb1)


tcga.tmb.cli=tcga_tmb1
library(survival)
library(survminer)
head(tcga.tmb.cli)
tmb.cutoff<-surv_cutpoint(tcga.tmb.cli,
                          time="OS.time",
                          event="OS",
                          variables=c("TMB"))
summary(tmb.cutoff)
tcga.tmb.cli$type <- ifelse(tcga.tmb.cli$TMB > tmb.cutoff$cutpoint$cutpoint, 'TMB-High', 'TMB-Low')
tcga.tmb.cli$ri_tmb=rep('none',nrow(tcga.tmb.cli))
tcga.tmb.cli$ri_tmb[which(tcga.tmb.cli$Risktype=='High' & tcga.tmb.cli$type=='TMB-High')]='H-Risk & H-TMB'
tcga.tmb.cli$ri_tmb[which(tcga.tmb.cli$Risktype=='Low' & tcga.tmb.cli$type=='TMB-Low')]='L-Risk & L-TMB'
tcga.tmb.cli$ri_tmb[which(tcga.tmb.cli$Risktype=='Low' & tcga.tmb.cli$type=='TMB-High')]='L-Risk & H-TMB'
tcga.tmb.cli$ri_tmb[which(tcga.tmb.cli$Risktype=='High' & tcga.tmb.cli$type=='TMB-Low')]='H-Risk & L-TMB'
tcga.tmb.cli$ri_tmb[which(tcga.tmb.cli$ri_tmb=='none')]=NA
table(tcga.tmb.cli$ri_tmb)
fig6b=ggsurvplot(fit=survfit(Surv(OS.time, OS) ~ type,
                             data = data.frame(OS.time = tcga.tmb.cli$OS.time/365
                                               , OS = tcga.tmb.cli$OS
                                               , type=tcga.tmb.cli$ri_tmb)),
                 data=data.frame(OS.time = tcga.tmb.cli$OS.time/365
                                 , OS = tcga.tmb.cli$OS
                                 , type=tcga.tmb.cli$ri_tmb),
                 conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,surv.median.line = 'hv',
                 title='TMB & Risktype',ggtheme=theme_classic(),
                 linetype = c("solid", "dashed","strata")[1],
                 palette =c("#F09496","#5F97C6","#F4DFDD","#A3CDEA"),
                 legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                 #legend = c(1,0.75), 
                 legend.title = "")#,legend.labs =  c('TMB-High','TMB-Low')

Fig6b=fig6b$plot
Fig6ab=mg_merge_plot(tcga_tmb_plot,Fig6b,nrow=1,ncol = 2)
ggsave('PDFs/Fig6ab.pdf',Fig6ab,height = 5,width = 10)


tcga.maf=getTCGAMAFByCode('LAML')
tcga.risktype.use=tcga.risktype.cli[,c('Samples','Risktype')]
table(tcga.risktype.use$Risktype)
colnames(tcga.risktype.use)[1]='Tumor_Sample_Barcode'
tcga.risktype.use$Tumor_Sample_Barcode=substr(tcga.risktype.use$Tumor_Sample_Barcode,1,12)
tcga.risktype.use.high=tcga.risktype.use[which(tcga.risktype.use$Risktype=='High'),]
tcga.risktype.use.low=tcga.risktype.use[which(tcga.risktype.use$Risktype=='Low'),]

write.table(tcga.risktype.use.high,file='06_risktype.mut/tcga.risktype.use.high.txt',row.names = F)
write.table(tcga.risktype.use.low,file='06_risktype.mut/tcga.risktype.use.low.txt',row.names = F)

tcga.maf.high=subsetMaf(tcga.maf,tsb=intersect(as.vector(tcga.maf@data$Tumor_Sample_Barcode),tcga.risktype.use.high$Tumor_Sample_Barcode))
tcga.maf.high<-read.maf(tcga.maf.high@data,isTCGA=T,clinicalData = '06_risktype.mut/tcga.risktype.use.high.txt')
tcga.maf.high@clinical.data

tcga.maf.low=subsetMaf(tcga.maf,tsb=intersect(as.vector(tcga.maf@data$Tumor_Sample_Barcode),tcga.risktype.use.low$Tumor_Sample_Barcode))
tcga.maf.low<-read.maf(tcga.maf.low@data,isTCGA=T,clinicalData = '06_risktype.mut/tcga.risktype.use.low.txt')
tcga.maf.low@clinical.data
dev.off()

pdf('06_risktype.mut/tcga.maf.high.pdf',height = 6,width = 7,onefile = F)
oncoplot(maf = tcga.maf.high,top = 20,sortByAnnotation = T)
dev.off()
pdf('06_risktype.mut/tcga.maf.low.pdf',height = 6,width = 7,onefile = F)
oncoplot(maf = tcga.maf.low,top =20,sortByAnnotation = T)
dev.off()

save.image(file = 'project.RData')





